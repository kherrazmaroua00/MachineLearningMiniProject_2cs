import time
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import pandas as pd
from pathlib import Path

# 1. Load datasets
df = pd.read_csv("molecules_expanded_raw.csv")
sources = pd.read_csv("molecules_expanded_raw_sources.csv")
df = df.merge(sources, on="molecule_chembl_id", how="left")

target_keys = (
    df.loc[df["source"] == "zenodo_new", "inchi_key"]
    .dropna()
    .unique()
    .tolist()
)

CACHE_FILE = Path("chembl_inchikey_lookup_cache.csv")
if CACHE_FILE.exists():
    cache_df = pd.read_csv(CACHE_FILE)
    # Keep only rows that already finished (or re-verify unfound ones if needed)
    done_keys = set(cache_df["inchi_key"])
else:
    cache_df = pd.DataFrame(
        columns=["inchi_key", "real_chembl_id", "pref_name", "max_phase", "first_approval"]
    )
    done_keys = set()

remaining_keys = [k for k in target_keys if k not in done_keys]
print(f"Total unique InChIKeys: {len(target_keys)} | Remaining: {len(remaining_keys)}")

session = requests.Session()
retries = Retry(total=5, backoff_factor=1.5, status_forcelist=[429, 500, 502, 503, 504])
session.mount("https://", HTTPAdapter(max_retries=retries))

# STEP 1: Map InChIKeys -> ChEMBL IDs + PubChem Titles via PubChem PUG REST POST (200 keys/req)
PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/property/Title/JSON"
CHEMBL_PK_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"

def resolve_batch_via_pubchem_and_chembl(batch_keys):
    # 1a. Get PubChem Titles (real chemical/drug names) in bulk via POST
    pc_titles = {}
    try:
        r = session.post(PUBCHEM_URL, data={"inchikey": ",".join(batch_keys)}, timeout=30)
        if r.status_code == 200:
            # Map back via individual or xref if needed; PubChem returns in order or with CID
            props = r.json().get("PropertyTable", {}).get("Properties", [])
            # Note: PubChem POST doesn't echo input InChIKey unless requested via InChIKey property
    except Exception:
        pass

    # Request InChIKey + Title together so mapping is 100% exact
    pc_url_exact = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/property/InChIKey,Title/JSON"
    try:
        r = session.post(pc_url_exact, data={"inchikey": ",".join(batch_keys)}, timeout=30)
        if r.status_code == 200:
            for item in r.json().get("PropertyTable", {}).get("Properties", []):
                ik = item.get("InChIKey")
                title = item.get("Title")
                if ik and title and not title.startswith("SCHEMBL"):
                    pc_titles[ik] = title
    except Exception as e:
        print(f"PubChem batch warning: {e}")

    # 1b. Query ChEMBL in small safe batches of 10 with fallback to single-key on 500
    chembl_hits = {}
    for j in range(0, len(batch_keys), 10):
        sub = batch_keys[j : j + 10]
        params = {
            "molecule_structures__standard_inchi_key__in": ",".join(sub),
            "only": "molecule_chembl_id,pref_name,max_phase,first_approval,molecule_structures",
            "limit": 10,
        }
        try:
            resp = session.get(CHEMBL_PK_URL, params=params, timeout=(10, 30))
            if resp.status_code == 200:
                for mol in resp.json().get("molecules", []):
                    ik = (mol.get("molecule_structures") or {}).get("standard_inchi_key")
                    if ik:
                        chembl_hits[ik] = mol
            else:
                raise RuntimeError(f"HTTP {resp.status_code}")
        except Exception:
            # Fallback: query one-by-one for this 10-item slice so a single bad key never kills the run
            for single_k in sub:
                try:
                    r1 = session.get(
                        CHEMBL_PK_URL,
                        params={
                            "molecule_structures__standard_inchi_key": single_k,
                            "only": "molecule_chembl_id,pref_name,max_phase,first_approval",
                        },
                        timeout=(10, 20),
                    )
                    if r1.status_code == 200:
                        mols = r1.json().get("molecules", [])
                        if mols:
                            chembl_hits[single_k] = mols[0]
                except Exception:
                    pass
                time.sleep(0.05)

    rows = []
    for k in batch_keys:
        mol = chembl_hits.get(k, {})
        real_name = mol.get("pref_name") or pc_titles.get(k)
        rows.append({
            "inchi_key": k,
            "real_chembl_id": mol.get("molecule_chembl_id"),
            "pref_name": real_name,
            "max_phase": mol.get("max_phase"),
            "first_approval": mol.get("first_approval"),
        })
    return rows

BATCH_SIZE = 100
new_records = []

for i in range(0, len(remaining_keys), BATCH_SIZE):
    batch = remaining_keys[i : i + BATCH_SIZE]
    rows = resolve_batch_via_pubchem_and_chembl(batch)
    new_records.extend(rows)

    # Save checkpoint every 100 keys so progress is never lost
    cache_df = pd.concat([cache_df, pd.DataFrame(new_records)], ignore_index=True)
    cache_df.drop_duplicates(subset=["inchi_key"], keep="last", inplace=True)
    cache_df.to_csv(CACHE_FILE, index=False)
    new_records = []

    done_total = len(done_keys) + i + len(batch)
    print(f"Saved {min(done_total, len(target_keys))}/{len(target_keys)} keys...")
    time.sleep(0.2)

# 2. Merge real values cleanly with explicit float64 casting (fixes FutureWarning)
lookup = cache_df.set_index("inchi_key")
mask_new = df["source"] == "zenodo_new"

df.loc[mask_new, "name"] = df.loc[mask_new, "inchi_key"].map(lookup["pref_name"])
df.loc[mask_new, "max_phase"] = pd.to_numeric(
    df.loc[mask_new, "inchi_key"].map(lookup["max_phase"]), errors="coerce"
)
df.loc[mask_new, "first_approval"] = pd.to_numeric(
    df.loc[mask_new, "inchi_key"].map(lookup["first_approval"]), errors="coerce"
)
df["resolved_chembl_id"] = df["inchi_key"].map(lookup["real_chembl_id"]).fillna(df["molecule_chembl_id"])

df.to_csv("molecules_expanded_enriched.csv", index=False)
print("Done! Saved to molecules_expanded_enriched.csv")
print(df.loc[mask_new, ["name", "max_phase", "first_approval"]].notna().sum())