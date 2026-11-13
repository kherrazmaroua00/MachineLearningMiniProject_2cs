# eda_utils.py — shared loading helper for the three EDA notebooks
from pathlib import Path
import pandas as pd

_MODULE_DIR = Path(__file__).resolve().parent
_CANDIDATES = [
    _MODULE_DIR / "data/raw",
    Path("data/raw"),
    Path("../data/raw"),
    Path("../../data/raw"),
]

def _find_raw():
    for c in _CANDIDATES:
        if (Path(c) / "molecules_raw.csv").exists():
            return Path(c)
    raise FileNotFoundError(
        "molecules_raw.csv not found. Looked in: " + ", ".join(_CANDIDATES)
        + f". Current folder: {Path.cwd()}."
    )

RAW = _find_raw()

def load_raw():
    """Return the five raw tables: mol, act, ind, mech, tgt."""
    names = ["molecules", "activities", "indications", "mechanisms", "targets"]
    return tuple(pd.read_csv(RAW / f"{n}_raw.csv") for n in names)


# ---------------------------------------------------------------------------
# Parent structures (in memory only — nothing is saved). Same logic as 2.4 of notebook 01.
# ---------------------------------------------------------------------------
DESC = ["mw", "logp", "hbd", "hba", "tpsa", "rotb", "arom_rings"]

def parent_table(mol):
    """Small molecules -> largest fragment, neutralized -> one row per parent structure.
    Descriptors are recomputed with RDKit on the parent. Returns a new DataFrame (in memory)."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize
    RDLogger.DisableLog("rdApp.*")
    lfc, unc = rdMolStandardize.LargestFragmentChooser(), rdMolStandardize.Uncharger()

    small = mol[(mol.molecule_type == "Small molecule") & mol.smiles.notna()].copy()
    small["parent"] = small.smiles.apply(lambda s: unc.uncharge(lfc.choose(Chem.MolFromSmiles(s))))
    small["parent_ik"] = small.parent.apply(Chem.MolToInchiKey)
    small["n_fragments"] = small.smiles.str.count(r"\.") + 1

    g = small.groupby("parent_ik")
    small["n_merged"] = g.molecule_chembl_id.transform("count")
    small["max_phase_group"] = g.max_phase.transform("max")
    small["first_approval_group"] = g.first_approval.transform("min")

    kept = small.sort_values(["parent_ik", "n_fragments"]).drop_duplicates("parent_ik").copy()
    kept["max_phase"] = kept.max_phase_group
    kept["first_approval"] = kept.first_approval_group
    P = kept.parent
    kept["mw"] = P.apply(Descriptors.MolWt)
    kept["logp"] = P.apply(Descriptors.MolLogP)
    kept["hbd"] = P.apply(rdMolDescriptors.CalcNumHBD)
    kept["hba"] = P.apply(rdMolDescriptors.CalcNumHBA)
    kept["tpsa"] = P.apply(rdMolDescriptors.CalcTPSA)
    kept["rotb"] = P.apply(rdMolDescriptors.CalcNumRotatableBonds)
    kept["arom_rings"] = P.apply(rdMolDescriptors.CalcNumAromaticRings)
    kept["approved"] = kept.max_phase == 4
    keep_cols = ["molecule_chembl_id", "name", "max_phase", "first_approval", "approved", "n_merged",
                 "parent", "parent_ik"] + DESC
    return kept[keep_cols].reset_index(drop=True)