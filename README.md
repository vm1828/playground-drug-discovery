# playground-drug-discovery

A computational drug-discovery pipeline for **acetylcholinesterase (AChE)**
inhibitors. It pulls IC50 bioactivity data from
[ChEMBL](https://www.ebi.ac.uk/chembl/), engineers Lipinski descriptors and
881-bit PaDEL PubChem fingerprints, explores the chemical space, and trains a
QSAR regression model that predicts pIC50 from molecular structure.

> Based on **Part 1** of the *Data Professor* "Computational Drug Discovery"
> tutorial series ([YouTube](https://www.youtube.com/watch?v=plVLRashaA8)).

The reusable logic lives in the importable `cdd` package (`src/cdd/`); the
notebook (`notebooks/bioactivity_qsar.ipynb`) is a thin narrative on top of it.
The same notebook runs unchanged in **VS Code / local Jupyter** and in
**Google Colab**.

> **Background.** This is a QSAR (Quantitative Structure–Activity Relationship)
> workflow: represent each molecule as a fixed-length numeric vector
> (fingerprints + physicochemical descriptors), then regress a potency value
> (pIC50) against it. It follows the standard ChEMBL → Lipinski → PaDEL →
> scikit-learn recipe, refactored into a tested, de-vendored project.

---

## Table of contents

- [Quick start](#quick-start)
- [Project layout](#project-layout)
- [The pipeline, step by step](#the-pipeline-step-by-step)
- [The `cdd` package API](#the-cdd-package-api)
- [Dataset](#dataset)
- [The `padel/` folder](#the-padel-folder)
- [Regenerating artifacts](#regenerating-artifacts)
- [Libraries](#libraries)
- [Design decisions](#design-decisions)
- [Troubleshooting](#troubleshooting)
- [Glossary](docs/GLOSSARY.md) — concise definitions of every term & concept

---

## Quick start

### Locally (VS Code / Jupyter)

```bash
./init.sh                    # creates .venv, installs deps, registers a kernel
```

Then open `notebooks/bioactivity_qsar.ipynb` and select the
**Python (drug-discovery)** kernel. The cached CSVs in `data/` let the whole
notebook run **offline** — no ChEMBL calls, no PaDEL recompute, no Java needed.

### Google Colab

1. Copy this project folder to your Google Drive.
2. Open `notebooks/bioactivity_qsar.ipynb` in Colab.
3. In the first (Setup) cell, set `PROJECT_ROOT` to the project's path on your
   Drive (default: `/gdrive/MyDrive/playground-drug-discovery`).

The Setup cell auto-detects Colab, mounts Drive, and `pip install`s
`requirements.txt`. Java is already available in Colab.

---

## Project layout

```text
playground-drug-discovery/
├── README.md                     this file
├── requirements.txt              Python dependencies
├── init.sh                       local environment bootstrap (.venv + kernel)
├── .gitignore
│
├── src/cdd/                      importable pipeline package
│   ├── __init__.py               exposes config, data, chem, analysis, features
│   ├── config.py                 paths & constants (single source of truth)
│   ├── data.py                   ChEMBL retrieval + preprocessing
│   ├── chem.py                   Lipinski descriptors + pIC50
│   ├── analysis.py               EDA plots + Mann-Whitney U tests
│   └── features.py               PaDEL fingerprints + aligned (X, y)
│
├── notebooks/
│   └── bioactivity_qsar.ipynb    end-to-end narrative (drives the cdd package)
│
├── data/                         cached datasets (checked in, run offline)
│   ├── data_raw.csv              raw ChEMBL IC50 records for CHEMBL220
│   ├── data_preprocessed.csv     cleaned + Lipinski + pIC50 table
│   └── PaDEL_descriptors_output.csv   881-bit PubChem fingerprints per molecule
│
├── padel/
│   └── PubchemFingerprinter.xml  PaDEL descriptor-type preset (PubChem only)
│
└── plots/                        generated figures (PDF, git-ignored)
```

Why paths resolve everywhere: `config.PROJECT_ROOT` is derived from
`config.py`'s own location (`Path(__file__).resolve().parents[2]`), so every
data/plot path is correct regardless of the current working directory — whether
you run from the repo root, from `notebooks/`, or from a Colab Drive mount.

---

## The pipeline, step by step

| # | Stage | What happens | Module |
| --- | ------- | -------------- | -------- |
| 1 | **Data** | Query ChEMBL for IC50 records of Human AChE (`CHEMBL220`). Cached in `data/data_raw.csv`. | `data.search_targets`, `data.fetch_bioactivity` |
| 2 | **Preprocess** | Drop rows missing IC50/SMILES, coerce IC50 to numeric, replace zeros with `1e-10`, label each compound, drop the ambiguous "intermediate" class. | `data.preprocess` |
| 3 | **Descriptors** | Add RDKit Lipinski descriptors (`MolWt`, `MolLogP`, `NumHDonors`, `NumHAcceptors`) and convert IC50 (nM) → pIC50. | `chem.build_features_table` |
| 4 | **EDA** | Class balance, MW-vs-LogP chemical space, and a Mann-Whitney U test per descriptor (active vs inactive). Figures saved to `plots/`. | `analysis.*` |
| 5 | **Fingerprints** | Compute 881-bit PubChem fingerprints with PaDEL. Cached in `data/PaDEL_descriptors_output.csv`. | `features.compute_padel_fingerprints` |
| 6 | **Assemble** | Align fingerprints (`X`) to labels (`y = pIC50`) **by molecule id**, drop low-variance features. | `features.prepare_xy` |
| 7 | **Model** | LazyPredict screen of many regressors, then train a Random Forest QSAR model (~0.73 test R²). | notebook |

**Bioactivity labelling** (IC50 in nM):

- `active`   — IC50 ≤ 1000 nM
- `inactive` — IC50 ≥ 10000 nM
- `intermediate` — everything in between (**dropped**)

**pIC50** = −log₁₀(IC50 in molar) = −log₁₀(IC50_nM × 10⁻⁹). Higher
pIC50 = more potent. Regressing pIC50 (rather than raw IC50) gives a
well-scaled, roughly linear target.

On the cached dataset: **9415 raw records → 6220 modelled compounds**
(3793 active / 2427 inactive), 137 fingerprint features after variance
filtering, RandomForest **test R² ≈ 0.73**.

---

## The `cdd` package API

Import once (the notebook's Setup cell adds `src/` to `sys.path`):

```python
from cdd import config, data, chem, analysis, features
```

### `config` — paths & constants

All tunables in one place: `PROJECT_ROOT`, `DATA_DIR`, `PLOTS_DIR`,
`RAW_CSV`, `PREPROCESSED_CSV`, `PADEL_OUTPUT_CSV`, `SMI_FILE`,
`PUBCHEM_FINGERPRINTER_XML`, `TARGET_CHEMBL_ID` (`"CHEMBL220"`),
`LIPINSKI_DESCRIPTORS`, `ACTIVE_THRESHOLD_NM` (1000),
`INACTIVE_THRESHOLD_NM` (10000).

### `data` — retrieval & preprocessing

| Function | Purpose |
| ---------- | --------- |
| `search_targets(query="acetylcholinesterase")` | Search ChEMBL targets (network). |
| `fetch_bioactivity(target_chembl_id="CHEMBL220")` | Download IC50 records (network). |
| `classify_bioactivity(ic50_nm)` | Label a Series active/inactive/intermediate. |
| `preprocess(df)` | Clean raw records → tidy `active`/`inactive` frame. |

### `chem` — molecular features

| Function | Purpose |
| ---------- | --------- |
| `lipinski_descriptors(smiles)` | RDKit descriptors as a DataFrame. |
| `to_pic50(ic50_nm)` | IC50 (nM) → pIC50. |
| `build_features_table(df)` | Attach Lipinski + pIC50, drop raw IC50. |

### `analysis` — EDA (figures saved to `plots/` by default, `save=False` to skip)

| Function | Purpose |
| ---------- | --------- |
| `plot_bioactivity_counts(df)` | Active vs inactive count plot. |
| `plot_mw_vs_logp(df)` | MW-vs-LogP scatter, sized by pIC50. |
| `mannwhitney(df, descriptor)` | Returns `(U, p, reject_H0)`. |
| `plot_descriptor(df, descriptor)` | Boxplot annotated with the U-test result. |

### `features` — fingerprints & model matrix

| Function | Purpose |
| ---------- | --------- |
| `write_smi(df)` | Write a `.smi` (SMILES, id) file for PaDEL. |
| `compute_padel_fingerprints(smi_path, output_csv)` | Run PaDEL via padelpy (needs Java). |
| `load_fingerprints(padel_csv)` | Load fingerprints indexed by molecule id (deduplicated). |
| `prepare_xy(features_table, fingerprints=None)` | **Id-aligned** `(X, y)`, low-variance features dropped. |

---

## Dataset

### Where the data comes from

Everything originates from **[ChEMBL](https://www.ebi.ac.uk/chembl/)**, the
EMBL-EBI open database of bioactive molecules. It is pulled over ChEMBL's public
REST API through the official `chembl_webresource_client` — no manual downloads,
no local database. The exact lineage, as coded in `cdd.data`:

1. **Find the target.** Search ChEMBL targets for `"acetylcholinesterase"` and
   take **Human Acetylcholinesterase**, ChEMBL id **`CHEMBL220`** (organism
   *Homo sapiens*, tax id 9606).

   ```python
   new_client.target.search("acetylcholinesterase")   # → pick CHEMBL220
   ```

2. **Pull its bioactivities.** Fetch every activity record for that target where
   the assay endpoint is **IC50** (half-maximal inhibitory concentration — the
   concentration of compound that inhibits the enzyme by 50%; lower = more
   potent).

   ```python
   new_client.activity.filter(target_chembl_id="CHEMBL220").filter(standard_type="IC50")
   ```

The result — one row per measured compound/assay, 46 ChEMBL columns — is cached
verbatim as `data/data_raw.csv`. Each record's key fields are
`molecule_chembl_id`, `canonical_smiles` (the structure), and `standard_value` /
`standard_units` (the IC50, in nM).

> **Provenance & licence.** ChEMBL data is released under
> [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/). IC50 values
> here come from many primary literature assays (see the `document_journal` /
> `document_year` columns), so absolute potencies carry inter-assay variability —
> a known limitation of pooled ChEMBL QSAR datasets.

### Cached files in `data/`

All three CSVs are checked in so the notebook runs **offline** (no network, no
Java). They are cache/output and can be regenerated any time (see
[Regenerating artifacts](#regenerating-artifacts)).

| File | Rows × cols | What it is | Produced by |
| ------ | ------------- | ------------ | ------------- |
| `data_raw.csv` | 9415 × 46 | Raw ChEMBL IC50 records for CHEMBL220, unmodified. | `data.fetch_bioactivity` |
| `data_preprocessed.csv` | 6220 × 8 | Cleaned, labelled compounds + Lipinski descriptors + pIC50. | `data.preprocess` → `chem.build_features_table` |
| `PaDEL_descriptors_output.csv` | 6220 × 882 | `Name` (ChEMBL id) + 881 PubChem fingerprint bits (`PubchemFP0…880`). | `features.compute_padel_fingerprints` |

### From 9415 raw records to 6220 modelled compounds

`data.preprocess` applies, in order:

1. Drop rows missing `standard_value` (IC50) or `canonical_smiles`.
2. Coerce IC50 to numeric; drop what won't parse.
3. Replace IC50 `0` with `1e-10` (keeps `pIC50` finite).
4. **Label** each compound from its IC50 in nM:
   `active` ≤ 1000, `inactive` ≥ 10000, else `intermediate`.
5. **Drop** the `intermediate` class (the two ends give a cleaner signal).

Yielding **6220 compounds — 3793 active / 2427 inactive**. Note 6220 rows span
only 5018 *unique* molecules: the same compound can have several IC50
measurements from different assays. This is why fingerprints are joined to labels
**by molecule id** rather than row position (see
[Design decisions](#design-decisions)).

### Columns in `data_preprocessed.csv`

| Column | Meaning |
| -------- | --------- |
| `molecule_chembl_id` | ChEMBL compound id. |
| `canonical_smiles` | Structure as a SMILES string. |
| `bioactivity_class` | `active` / `inactive`. |
| `MolWt` | Molecular weight (Da). |
| `MolLogP` | Calculated octanol–water partition coefficient (lipophilicity). |
| `NumHDonors` | Hydrogen-bond donor count. |
| `NumHAcceptors` | Hydrogen-bond acceptor count. |
| `pIC50` | −log₁₀(IC50 in M); the regression **target**. |

---

## The `padel/` folder

Contains a single 4.6 KB file, `PubchemFingerprinter.xml` — a PaDEL-Descriptor
**"descriptor types" preset**. It is an on/off checklist telling PaDEL *which*
descriptors to compute; every entry is `false` except the PubChem fingerprinter:

```xml
<Descriptor name="PubchemFingerprinter" value="true"/>
```

So this preset means *"compute only the 881-bit PubChem fingerprints"* — the
columns you see in `data/PaDEL_descriptors_output.csv`.

It lives in the repo (rather than the original 54 MB `PaDEL/` jar tree) because
`padelpy` bundles the PaDEL-Descriptor **jar** but not these XML presets.
`features.compute_padel_fingerprints()` passes it through:

```python
padeldescriptor(..., descriptortypes=str(config.PUBCHEM_FINGERPRINTER_XML),
                fingerprints=True, ...)
```

You only touch it when **recomputing** fingerprints; the cached CSV makes it
optional for normal runs. The original project shipped ~15 such presets (MACCS,
Substructure, Klekota-Roth, …); only the PubChem one is kept because that is all
this pipeline uses.

---

## Regenerating artifacts

Everything cached can be rebuilt. All snippets assume the notebook's imports.

**Re-download from ChEMBL** (needs network):

```python
raw = data.fetch_bioactivity("CHEMBL220")
raw.to_csv(config.RAW_CSV, index=False)
```

**Recompute PaDEL fingerprints** (needs Java):

```python
df  = chem.build_features_table(data.preprocess(pd.read_csv(config.RAW_CSV)))
smi = features.write_smi(df)                 # writes data/molecule.smi
features.compute_padel_fingerprints(smi)     # writes PaDEL_descriptors_output.csv
```

**Regenerate all plots**: just re-run the EDA cells; each `analysis.plot_*`
function writes its PDF to `plots/` (git-ignored).

---

## Libraries

Installed by `init.sh` / `requirements.txt`. Two are version-pinned to match the
original project; the rest float on their latest compatible release.

| Library | Role in this project | Version |
| --------- | ---------------------- | --------- |
| [`chembl_webresource_client`](https://github.com/chembl/chembl_webresource_client) | Talks to the ChEMBL REST API — target search + activity retrieval (`data.py`). | `==0.10.9` |
| [RDKit](https://www.rdkit.org/) | Cheminformatics toolkit: parse SMILES, compute Lipinski descriptors (`chem.py`). | `==2024.9.5` |
| [`padelpy`](https://github.com/ecrl/padelpy) | Python wrapper that bundles the **PaDEL-Descriptor** jar to compute PubChem fingerprints (`features.py`). Needs Java. | latest |
| [scikit-learn](https://scikit-learn.org/) | `VarianceThreshold` feature selection, `train_test_split`, `RandomForestRegressor`. | latest |
| [LazyPredict](https://github.com/shankarpandala/lazypredict) | One-line screen of ~40 regressors to pick a model family. | latest |
| [pandas](https://pandas.pydata.org/) / [NumPy](https://numpy.org/) | Data frames and numerics throughout. | latest |
| [SciPy](https://scipy.org/) | `mannwhitneyu` U-test in `analysis.py`. | latest |
| [seaborn](https://seaborn.pydata.org/) / [Matplotlib](https://matplotlib.org/) | All EDA and result plots. | latest |
| ipykernel | Registers the local Jupyter kernel (`init.sh`). | latest |

**External runtime — Java.** `padelpy` shells out to a Java jar, so a **JRE** is
required *only to recompute fingerprints* (`sudo apt-get install default-jre`;
Colab already has it). Because `data/PaDEL_descriptors_output.csv` is cached, the
full notebook runs end-to-end with no Java installed.

---

## Methods

The modelling choices, and where each lives:

- **Bioactivity binning** (`data.classify_bioactivity`) — IC50 thresholds of
  1000 / 10000 nM split compounds into `active` / `inactive`, dropping the
  ambiguous middle. A standard way to turn continuous potency into clean classes
  for EDA.
- **Lipinski descriptors** (`chem.lipinski_descriptors`) — `MolWt`, `MolLogP`,
  `NumHDonors`, `NumHAcceptors`: the four properties of Lipinski's *Rule of Five*,
  a rule-of-thumb for oral drug-likeness. Used here as interpretable features and
  EDA axes.
- **pIC50 target** (`chem.to_pic50`) — regressing −log₁₀(IC50 in M) instead of
  raw IC50 gives an evenly-scaled, roughly linear target (e.g. 100 nM → 7.0).
- **Mann-Whitney U test** (`analysis.mannwhitney`) — a non-parametric test of
  whether each descriptor differs between active and inactive compounds (no
  normality assumption). On this data `MolWt`, `MolLogP`, `NumHDonors` differ
  significantly; `NumHAcceptors` does not.
- **PubChem fingerprints** (`features.compute_padel_fingerprints`) — each
  molecule is encoded as an **881-bit** binary substructure fingerprint via PaDEL.
  These bits are the model's input features.
- **Variance-threshold selection** (`features.prepare_xy`) — drop near-constant
  bits (same value in >80% of molecules), trimming 881 → **137** informative
  features.
- **Model selection** — LazyPredict screens many regressors, then a
  **Random Forest** (`n_estimators=100`) is trained on an 80/20 split. Test
  **R² ≈ 0.73** on pIC50 — a solid QSAR baseline.

---

## Design decisions

- **De-vendored.** The original copy carried a 54 MB `PaDEL/` folder — the same
  jars duplicated 8×, a `__MACOSX/` junk tree, and a `padel.zip`. Replaced by the
  `padelpy` pip dependency (which supplies the jar) plus one 4.6 KB XML preset.
  Repo shrank **122 MB → 19 MB** (the remainder is the cached data).
- **Logic in a package, not the notebook.** Reusable functions live in
  `src/cdd/` so they are importable, diffable, and testable; the notebook stays a
  readable narrative instead of Colab `#@title` form cells.
- **Portable Setup cell.** Detects Colab via `try: import google.colab`; locally
  it walks up from the notebook to find `src/cdd`. One notebook, both worlds.
- **Id-based X/y alignment (bug fix).** The original paired the PaDEL matrix and
  the labels **by row position**. PaDEL does **not** preserve input order, so the
  labels were misaligned to the fingerprints. `features.prepare_xy` aligns them
  **by molecule id** (fingerprints depend only on structure, so duplicate
  measurements of the same molecule map to the same, correct fingerprint).
- **Cached data checked in.** `data/*.csv` is committed so the notebook is
  reproducible offline; `plots/*.pdf` is git-ignored as regenerated output.

---

## Troubleshooting

- **`ModuleNotFoundError: cdd`** — the Setup cell didn't find `src/cdd`. Run the
  notebook from within the repo (locally), or set `PROJECT_ROOT` correctly in the
  Setup cell (Colab).
- **PaDEL / padelpy errors** — you're recomputing fingerprints without Java.
  Install a JRE, or skip that cell and use the cached
  `data/PaDEL_descriptors_output.csv`.
- **`X` and `y` length mismatch / all-NaN fingerprints** — the PaDEL CSV is stale
  relative to `data_raw.csv`. Recompute fingerprints (see above).
- **ChEMBL request errors** — network/rate limiting. Use the cached
  `data/data_raw.csv` instead of the live-query cell.
