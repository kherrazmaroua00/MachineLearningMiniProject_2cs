import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog("rdApp.*")

P = r"C:\Users\ibrah\Downloads\ML\MINIPROJET"
orig = pd.read_csv(rf"{P}\data\raw\molecules_raw.csv")
kept = set(pd.read_csv(rf"{P}\molecules_expanded_raw.csv")["molecule_chembl_id"])
chooser = rdMolStandardize.LargestFragmentChooser()

def reason(smi):
    if not isinstance(smi, str) or not smi.strip():
        return "no SMILES"
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return "invalid SMILES"
    try:
        Chem.MolToInchiKey(chooser.choose(m))
    except Exception:
        return "InChIKey failed"
    return "duplicate of another molecule"

lost = orig[~orig["molecule_chembl_id"].isin(kept)].copy()
lost["reason"] = lost["smiles"].apply(reason)
print(lost["reason"].value_counts())
print("\nPhase of lost molecules:\n", lost["max_phase"].value_counts())
lost[["molecule_chembl_id", "name", "max_phase", "molecule_type", "reason"]] \
    .to_csv(rf"{P}\lost_molecules.csv", index=False)
print(lost[["molecule_chembl_id", "name", "max_phase", "reason"]].head(20))