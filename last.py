import pandas as pd

df = pd.read_csv("molecules_expanded_final.csv", low_memory=False)

# Preserve the full IUPAC / chemical title in a dedicated column
df["iupac_name"] = df["name"]

# Keep 'name' concise (<= 40 chars, matching the original 407 ChEMBL rows)
mask_too_long = df["name"].astype(str).str.len() > 40
df.loc[mask_too_long, "name"] = df.loc[mask_too_long, "resolved_chembl_id"]

df.to_csv("molecules_expanded_clean_names.csv", index=False)