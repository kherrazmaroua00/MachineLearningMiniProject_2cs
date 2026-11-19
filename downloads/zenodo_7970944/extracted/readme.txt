----------------------------
Date: 2023_04 of 04-Apr-2023
----------------------------


Description of the files
------------------------

1. human_PKI_active.tsv - human PKIs with -log potency of at least 5 
2. mouse_PKI_active.tsv - mouse PKIs with -log potency of ar least 5 
3. human_PKI_inactive.tsv - designated human PKIs with log potency < 5 
4. mouse_PKI_inactive.tsv - designated mouse PKIs with log potency < 5

Description of the file contents
--------------------------------

Compound_new_ID : internal ID assigned to each compound
nonstereo_aromatic_smile: standardized compound SMILES
Uniprot_ID: UniProtKB/Swiss-Prot human and mouse kinase IDs
pref_name: preferred/full name of the kinase family 
activity_id: reference to the binding data retrieved by ChEMBL and/or BindingDB
mean_log: mean potency (negative decadic logarithmic values of IC50, Ki or Kd). NaN if the measure is qualitative
selected_stvalue: standard value used
ORGANISM: organism of the kinase

Only for human_PKI_active.tsv:
CPKI: True if the compound contains one of the 14 investigated warheads



