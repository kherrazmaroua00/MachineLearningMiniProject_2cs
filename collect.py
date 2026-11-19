"""
DrugScope AI - Data Collection (collection only, no cleaning)
Robust version: automatic retries + resume from saved files.

Setup:   pip install chembl_webresource_client pandas tqdm
Run:     python collect.py
If it stops, just run it again: finished steps are skipped.
"""
import time
import datetime
from pathlib import Path

import pandas as pd
from tqdm import tqdm
from chembl_webresource_client.new_client import new_client

try:
    from chembl_webresource_client.settings import Settings
    Settings.Instance().TIMEOUT = 60
except Exception:
    pass

RAW = Path("data/raw")
RAW.mkdir(parents=True, exist_ok=True)

# ---- Selection criteria ----
MIN_PHASE = 1                       # 1-3 = clinical, 4 = approved
ACTIVITY_TYPES = ["IC50", "Ki", "Kd"]
ACTION_TYPES = ["INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR"]


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def fetch(make_query, retries=6):
    """Run a ChEMBL query; on network errors wait and retry."""
    for attempt in range(retries):
        try:
            return list(make_query())
        except Exception as e:
            wait = min(5 * 2 ** attempt, 60)
            print(f"\n   ! {type(e).__name__}; retry {attempt + 1}/{retries} in {wait}s")
            time.sleep(wait)
    raise RuntimeError("ChEMBL kept failing. Wait a few minutes and run again.")


def load_or_none(name):
    p = RAW / name
    return pd.read_csv(p) if p.exists() and p.stat().st_size > 0 else None


def main():
    counts = {}

    # ChEMBL version (for documentation)
    chembl_version = "unknown - write it from https://www.ebi.ac.uk/chembl/"
    try:
        chembl_version = new_client.status.get()["chembl_db_version"]
    except Exception:
        try:
            chembl_version = new_client.status.status.get()["chembl_db_version"]
        except Exception:
            pass
    print("ChEMBL version:", chembl_version)

    # 1. Kinase targets
    t_df = load_or_none("targets_raw.csv")
    if t_df is None:
        print("[1/5] Fetching kinase targets...")
        t_df = pd.DataFrame(fetch(lambda: new_client.target.filter(
            target_type="SINGLE PROTEIN",
            organism="Homo sapiens",
            pref_name__icontains="kinase",
        ).only(["target_chembl_id", "pref_name", "organism", "target_type"])))
        t_df.to_csv(RAW / "targets_raw.csv", index=False)
    else:
        print("[1/5] targets_raw.csv found, skipping")
    kinase_ids = t_df["target_chembl_id"].unique().tolist()
    counts["targets_raw.csv"] = len(t_df)
    print(f"      {len(kinase_ids)} kinase targets")

    # 2. Mechanisms
    m_df = load_or_none("mechanisms_raw.csv")
    if m_df is None:
        print("[2/5] Fetching mechanisms of action...")
        mech = []
        for ch in tqdm(list(chunks(kinase_ids, 20))):
            mech += fetch(lambda ch=ch: new_client.mechanism.filter(
                target_chembl_id__in=ch, action_type__in=ACTION_TYPES
            ).only(["molecule_chembl_id", "target_chembl_id", "action_type",
                    "mechanism_of_action", "max_phase"]))
            time.sleep(1)
        m_df = pd.DataFrame(mech)
        m_df.to_csv(RAW / "mechanisms_raw.csv", index=False)
    else:
        print("[2/5] mechanisms_raw.csv found, skipping")
    mol_ids = m_df["molecule_chembl_id"].dropna().unique().tolist()
    counts["mechanisms_raw.csv"] = len(m_df)
    print(f"      {len(mol_ids)} candidate molecules")

    # 3. Molecule records
    mol_df = load_or_none("molecules_raw.csv")
    if mol_df is None:
        print("[3/5] Fetching molecule records...")
        rows = []
        for ch in tqdm(list(chunks(mol_ids, 20))):
            for m in fetch(lambda ch=ch: new_client.molecule.filter(molecule_chembl_id__in=ch)):
                props = m.get("molecule_properties") or {}
                structs = m.get("molecule_structures") or {}
                rows.append({
                    "molecule_chembl_id": m.get("molecule_chembl_id"),
                    "name": m.get("pref_name"),
                    "max_phase": m.get("max_phase"),
                    "first_approval": m.get("first_approval"),
                    "molecule_type": m.get("molecule_type"),
                    "smiles": structs.get("canonical_smiles"),
                    "inchi_key": structs.get("standard_inchi_key"),
                    "mol_weight": props.get("full_mwt"),
                    "alogp": props.get("alogp"),
                    "hbd": props.get("hbd"),
                    "hba": props.get("hba"),
                    "psa": props.get("psa"),
                    "rtb": props.get("rtb"),
                    "aromatic_rings": props.get("aromatic_rings"),
                    "ro5_violations": props.get("num_ro5_violations"),
                })
            time.sleep(1)
        mol_df = pd.DataFrame(rows)
        mol_df.to_csv(RAW / "molecules_raw.csv", index=False)
    else:
        print("[3/5] molecules_raw.csv found, skipping")
    counts["molecules_raw.csv"] = len(mol_df)
    phase = pd.to_numeric(mol_df["max_phase"], errors="coerce")
    selected_ids = mol_df.loc[phase >= MIN_PHASE, "molecule_chembl_id"].tolist()
    n_approved = int((phase == 4).sum())
    print(f"      {len(mol_df)} molecules, {len(selected_ids)} with max_phase >= {MIN_PHASE}, "
          f"{n_approved} approved")

    # 4. Bioactivities
    a_df = load_or_none("activities_raw.csv")
    if a_df is None:
        print("[4/5] Fetching bioactivities (slowest step)...")
        acts = []
        for ch in tqdm(list(chunks(selected_ids, 10))):
            acts += fetch(lambda ch=ch: new_client.activity.filter(
                molecule_chembl_id__in=ch,
                target_chembl_id__in=kinase_ids,
                standard_type__in=ACTIVITY_TYPES,
                standard_relation="=",
                pchembl_value__isnull=False,
            ).only(["molecule_chembl_id", "target_chembl_id", "target_pref_name",
                    "standard_type", "standard_value", "standard_units",
                    "pchembl_value", "assay_chembl_id"]))
            time.sleep(1)
        a_df = pd.DataFrame(acts)
        a_df.to_csv(RAW / "activities_raw.csv", index=False)
    else:
        print("[4/5] activities_raw.csv found, skipping")
    counts["activities_raw.csv"] = len(a_df)

    # 5. Indications
    i_df = load_or_none("indications_raw.csv")
    if i_df is None:
        print("[5/5] Fetching therapeutic indications...")
        inds = []
        for ch in tqdm(list(chunks(selected_ids, 20))):
            inds += fetch(lambda ch=ch: new_client.drug_indication.filter(
                molecule_chembl_id__in=ch
            ).only(["molecule_chembl_id", "mesh_heading", "efo_term",
                    "max_phase_for_ind"]))
            time.sleep(1)
        i_df = pd.DataFrame(inds)
        i_df.to_csv(RAW / "indications_raw.csv", index=False)
    else:
        print("[5/5] indications_raw.csv found, skipping")
    counts["indications_raw.csv"] = len(i_df)

    # Documentation
    doc = [
        "# DrugScope AI - Data Sources Documentation",
        "",
        f"- Collection date: {datetime.date.today()}",
        "- Source: ChEMBL (https://www.ebi.ac.uk/chembl/)",
        f"- ChEMBL version: {chembl_version}",
        "- Access method: chembl_webresource_client (official ChEMBL Python client)",
        "- PubChem: not used (add a note here if used later)",
        "",
        "## Selection criteria",
        "- Targets: human, SINGLE PROTEIN, name contains 'kinase'",
        f"- Molecules: mechanism action type in {ACTION_TYPES} against those targets",
        f"- Activities/indications collected for molecules with max_phase >= {MIN_PHASE} "
        "(4 = approved, 1-3 = clinical/exploratory)",
        f"- Activity types: {ACTIVITY_TYPES}, relation '=', pChEMBL value present",
        "",
        "## Raw files and record counts",
    ] + [f"- {name}: {n} rows" for name, n in counts.items()] + [
        "",
        f"- Distinct candidate molecules: {len(mol_df)}",
        f"- Selected molecules (max_phase >= {MIN_PHASE}): {len(selected_ids)}",
        f"- Approved drugs (max_phase = 4): {n_approved}",
        "",
        "## Known limitations",
        "- Kinase targets selected by name keyword, not by protein classification.",
        "- Not every molecule has activity or indication data.",
        "- Data is raw: duplicates, invalid SMILES and missing values are NOT handled "
        "(cleaning is a separate task).",
    ]
    Path("data/DATA_SOURCES.md").write_text("\n".join(doc), encoding="utf-8")
    print("\nDone. Raw files in data/raw/, documentation in data/DATA_SOURCES.md")


if __name__ == "__main__":
    main()