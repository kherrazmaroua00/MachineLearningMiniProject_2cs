import pandas as pd
from chembl_webresource_client.new_client import new_client

# 1. Kinase targets (same rule as before)
targets = new_client.target.filter(
    organism="Homo sapiens", target_type="SINGLE PROTEIN",
    pref_name__icontains="kinase"
).only(["target_chembl_id", "pref_name"])
target_ids = [t["target_chembl_id"] for t in targets]

# 2. Activities: no phase filter, so preclinical compounds are included
mol_ids = set()
for i in range(0, len(target_ids), 20):
    acts = new_client.activity.filter(
        target_chembl_id__in=target_ids[i:i+20],
        standard_type__in=["IC50", "Ki", "Kd"],
        standard_relation="=",
        pchembl_value__isnull=False,
    ).only(["molecule_chembl_id"])
    mol_ids.update(a["molecule_chembl_id"] for a in acts)
print(len(mol_ids), "molecules found")

# 3. Fetch the same columns as molecules_raw.csv
ids = list(mol_ids)
rows = []
for i in range(0, len(ids), 50):
    for m in new_client.molecule.filter(molecule_chembl_id__in=ids[i:i+50]):
        p = m.get("molecule_properties") or {}
        s = m.get("molecule_structures") or {}
        rows.append({
            "molecule_chembl_id": m["molecule_chembl_id"],
            "name": m.get("pref_name"),
            "max_phase": m.get("max_phase"),
            "first_approval": m.get("first_approval"),
            "molecule_type": m.get("molecule_type"),
            "smiles": s.get("canonical_smiles"),
            "inchi_key": s.get("standard_inchi_key"),
            "mol_weight": p.get("full_mwt"),
            "alogp": p.get("alogp"),
            "hbd": p.get("hbd"), "hba": p.get("hba"),
            "psa": p.get("psa"), "rtb": p.get("rtb"),
            "aromatic_rings": p.get("aromatic_rings"),
            "ro5_violations": p.get("num_ro5_violations"),
        })

df = pd.DataFrame(rows).drop_duplicates("inchi_key").dropna(subset=["smiles"])
df.to_csv("molecules_expanded_raw.csv", index=False)