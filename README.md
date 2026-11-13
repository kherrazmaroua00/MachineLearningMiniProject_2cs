# Machine Learning Mini Project

Exploratory data analysis of kinase-inhibitor data from ChEMBL. The project
examines molecular structures alongside their biological activities, targets,
mechanisms, and disease indications.

## Project contents

```text
data/
├── raw/          Original CSV datasets
├── dictionary/   Column-level data dictionaries
└── reports/      Generated data-quality and statistical reports
notebooks/
├── 01_overview_and_data_dictionary.ipynb
├── 02_data_quality.ipynb
├── 03_statistical_analysis.ipynb
├── 04_dataset_specific_analysis.ipynb
└── eda_utils.py
```

## Analysis workflow

1. **Overview and data dictionary** — inspects the five datasets and documents
   their columns, types, missing values, and intended roles.
2. **Data quality** — evaluates missing values, duplicates, invalid values, and
   inconsistencies without modifying the raw data.
3. **Statistical analysis** — explores distributions, relationships, and
   outliers, including molecular descriptors.
4. **Dataset-specific analysis** — investigates the characteristics of each
   biological dataset and prepares context for later modelling.

The raw datasets contain:

| Dataset | Rows | Description |
|---|---:|---|
| `molecules` | 407 | ChEMBL compounds and molecular descriptors |
| `activities` | 15,655 | Molecule–kinase activity measurements |
| `indications` | 2,864 | Molecule–disease associations |
| `mechanisms` | 610 | Molecule–target mechanism records |
| `targets` | 548 | Kinase protein targets |

## Requirements

Python 3.9+ with:

- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `rdkit`
- Jupyter Notebook or JupyterLab

Install the commonly available packages with:

```bash
pip install pandas numpy matplotlib seaborn jupyter
```

Install RDKit using the recommended package manager for your platform (for
example, `conda install -c conda-forge rdkit`).

## Running the notebooks

Open the repository in Jupyter and run the notebooks in numerical order:

```bash
jupyter notebook
```

Run each notebook from its own file in `notebooks/`. The shared
`notebooks/eda_utils.py` helper locates `data/raw/` and loads the five raw
tables. The analysis is exploratory: raw CSV files are not cleaned or
overwritten. Reports are written to `data/reports/`.

## Data provenance

The datasets are derived from ChEMBL kinase-inhibitor data for educational
analysis. Consult the original ChEMBL documentation and data-release
information before using the data for scientific or clinical conclusions.

## Scope and limitations

This repository contains exploratory analysis, not a production machine
learning pipeline. Results depend on the source data, selected descriptors,
missing-value patterns, and the limitations documented in the notebooks.
