# MachineLearningMiniProject_2cs

## About the project

This project, also referred to as **DrugScope AI**, is an exploratory machine
learning project built around kinase-inhibitor data from ChEMBL. It studies
small-molecule compounds and connects their chemical properties with biological
activity, protein targets, mechanisms of action, and therapeutic indications.

The objective is to understand the structure and quality of the data before
building a downstream clustering or predictive model. The analysis is
deliberately exploratory: it documents patterns, limitations, and potential
preprocessing decisions without changing the original datasets.

## Objectives

- Understand the structure and meaning of each dataset.
- Build reusable data dictionaries for all source tables.
- Identify missing values, duplicates, invalid values, and inconsistencies.
- Describe molecular properties and their statistical distributions.
- Explore relationships between molecular descriptors and biological data.
- Investigate each dataset independently before combining results.
- Produce documented, reproducible evidence for future modelling work.

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

### Dataset relationships

The tables are linked through ChEMBL identifiers:

- `molecules` contains compounds and chemical descriptors.
- `activities` records measured molecule–target activity values.
- `targets` describes the biological targets associated with activities.
- `mechanisms` records known molecule–target mechanisms of action.
- `indications` records molecule–disease associations.

`molecules` is the primary source for future molecular clustering. The other
tables provide biological context that can be used to interpret clusters.

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

The row counts describe the current project snapshot and may change if the
source data is refreshed.

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

For a reproducible environment, create and activate a dedicated environment
before installing packages:

```bash
conda create -n drugscope python=3.10
conda activate drugscope
conda install -c conda-forge rdkit
pip install pandas numpy matplotlib seaborn jupyter
```

## Running the notebooks

Open the repository in Jupyter and run the notebooks in numerical order:

```bash
jupyter notebook
```

Run each notebook from its own file in `notebooks/`. The shared
`notebooks/eda_utils.py` helper locates `data/raw/` and loads the five raw
tables. The analysis is exploratory: raw CSV files are not cleaned or
overwritten. Reports are written to `data/reports/`.

## Outputs

The notebooks create or document the following outputs:

- One data dictionary per source dataset.
- Missing-value, duplicate, and data-quality reports.
- Univariate numerical statistics.
- A molecular descriptor correlation matrix.
- An outlier analysis report.
- Dataset-specific visualizations and observations inside the notebooks.

Existing CSV outputs are included in `data/dictionary/` and `data/reports/` so
the results can be reviewed without rerunning every notebook.

## Data provenance

The datasets are derived from ChEMBL kinase-inhibitor data for educational
analysis. ChEMBL is a manually curated database of bioactive molecules and
their properties. Consult the [ChEMBL documentation](https://www.ebi.ac.uk/chembl/)
and the relevant data-release information for complete provenance and
licensing details.

## Scope and limitations

This repository contains exploratory analysis, not a production machine
learning pipeline. Results depend on the source data, selected descriptors,
missing-value patterns, and the limitations documented in the notebooks.
Activity measurements may come from different experiments and conditions, so
they should not be treated as directly interchangeable without additional
normalization and domain review. Nothing in this repository is intended for
clinical decision-making.

## Contributing

Suggestions, corrections, and improvements are welcome. When extending the
analysis, preserve the raw data, document new assumptions, and place generated
reports in the appropriate directory.
