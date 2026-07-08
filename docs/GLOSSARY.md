# Glossary

Concise definitions of every domain term, metric, and tool used in this
project. Grouped by theme; skim the header you need.

See also: [README](../README.md) for how these fit together in the pipeline.

---

## Drug discovery & biology

| Term | Meaning |
| --- | --- |
| **Drug discovery** | The process of identifying chemical compounds that produce a desired biological effect (here, inhibiting an enzyme). |
| **Target** | The biological molecule (usually a protein) a drug acts on. Ours is an enzyme. |
| **Acetylcholinesterase (AChE)** | The enzyme this project targets. It breaks down the neurotransmitter acetylcholine; inhibiting it is a strategy for Alzheimer's and other conditions. |
| **Inhibitor** | A compound that reduces a target's activity. We look for AChE inhibitors. |
| **Enzyme** | A protein that catalyses a biochemical reaction. |
| **Ligand / compound / molecule** | A small chemical that binds the target. Used interchangeably here for the things we measure and model. |
| **Assay** | A lab experiment that measures how strongly a compound affects the target. |
| **Bioactivity** | A compound's measured effect on the target — the raw signal we model. |
| **Potency** | How little compound is needed for an effect. High potency = low IC50 = high pIC50. |
| **Organism / *Homo sapiens* / tax id 9606** | The species the target comes from. We use human AChE. |

## Potency measures

| Term | Meaning |
| --- | --- |
| **IC50** | Half-maximal inhibitory concentration: the compound concentration that inhibits the target by 50%. **Lower = more potent.** Measured here in nanomolar. |
| **nM (nanomolar)** | Concentration unit, 10⁻⁹ mol/L. The units IC50 is reported in. |
| **M (molar)** | Concentration unit, mol/L. pIC50 is defined on IC50 expressed in molar. |
| **pIC50** | −log₁₀(IC50 in molar) = −log₁₀(IC50_nM × 10⁻⁹). **Higher = more potent** (e.g. 100 nM → 7.0). The regression **target**; better-scaled and more linear than raw IC50. |
| **Bioactivity class** | Discrete potency label derived from IC50: `active` (≤ 1000 nM), `inactive` (≥ 10000 nM), or `intermediate` (in between, **dropped**). |
| **Active / inactive / intermediate** | The three bins above. Only the two extremes are kept — they give a cleaner signal for EDA and modelling. |

## Cheminformatics

| Term | Meaning |
| --- | --- |
| **Cheminformatics** | Using computation to represent and analyse chemical structures. |
| **SMILES** | Simplified Molecular-Input Line-Entry System: a molecule's structure encoded as a text string (column `canonical_smiles`). |
| **Canonical SMILES** | A single standardised SMILES per molecule, so identical structures get identical strings. |
| **Molecular descriptor** | Any numeric property computed from a structure (e.g. molecular weight). Used as model features / EDA axes. |
| **Fingerprint** | A fixed-length bit vector where each bit flags the presence/absence of a substructure. Turns a molecule into a numeric feature vector. |
| **PubChem fingerprint** | The specific **881-bit** substructure fingerprint standard used here (columns `PubchemFP0…880`). Computed by PaDEL. |
| **Lipinski descriptors** | The four properties in Lipinski's *Rule of Five*: `MolWt`, `MolLogP`, `NumHDonors`, `NumHAcceptors`. Used as interpretable features. |
| **Rule of Five** | Lipinski's rule-of-thumb for oral drug-likeness (limits on weight, lipophilicity, and H-bond donors/acceptors). |
| **MolWt** | Molecular weight, in daltons (Da). |
| **MolLogP** | Calculated octanol–water partition coefficient — a measure of **lipophilicity** (fat- vs water-solubility). |
| **NumHDonors / NumHAcceptors** | Counts of hydrogen-bond donor / acceptor atoms in the molecule. |
| **Salt removal / nitro standardisation / aromaticity detection** | Structure clean-up steps PaDEL applies before computing fingerprints (`removesalt`, `standardizenitro`, `detectaromaticity`). |

## Modelling, statistics & ML

| Term | Meaning |
| --- | --- |
| **QSAR** | Quantitative Structure–Activity Relationship: model that predicts a potency value from a molecule's numeric structure representation. The core method of this project. |
| **Feature / X** | The model inputs — the (variance-filtered) fingerprint bits. |
| **Target / label / y** | The value being predicted — `pIC50`. |
| **Regression** | Predicting a continuous number (pIC50), as opposed to classification. |
| **EDA** | Exploratory Data Analysis: plots and tests to understand the data before modelling. |
| **Chemical space** | The spread of molecules across descriptor axes; visualised here as MolWt vs MolLogP. |
| **Class balance** | The ratio of active to inactive compounds (3793 / 2427 here). |
| **Mann-Whitney U test** | Non-parametric test of whether a descriptor's distribution differs between active and inactive compounds (no normality assumption). |
| **U statistic** | The Mann-Whitney test's rank-based statistic. |
| **p-value** | Probability of seeing a difference this large by chance; small p ⇒ a real difference. |
| **α (alpha)** | Significance cutoff (0.05). If p ≤ α we reject the null hypothesis. |
| **H0 (null hypothesis)** | "No difference between active and inactive." "Reject H0" = the descriptor differs significantly. |
| **VarianceThreshold** | Feature selection that drops near-constant bits (same value in > 80% of molecules), trimming 881 → 137 features. |
| **train_test_split** | Splitting data into training and held-out test sets (80/20 here). |
| **Random Forest** | The chosen QSAR model: an ensemble of decision trees (`n_estimators=100`). Test **R² ≈ 0.73**. |
| **R² (coefficient of determination)** | Regression quality score; 1.0 is perfect, 0 is no better than the mean. |
| **LazyPredict** | Utility that trains ~40 regressors at once to quickly compare model families before committing. |
| **Id-based alignment** | Joining fingerprints (X) to labels (y) by `molecule_chembl_id`, not row position — because PaDEL doesn't preserve input order. |

## Data & project artifacts

| Term | Meaning |
| --- | --- |
| **ChEMBL** | EMBL-EBI's open database of bioactive molecules; the source of all raw data, via its REST API. |
| **CHEMBL220** | The ChEMBL id of the specific target used: Human Acetylcholinesterase. |
| **molecule_chembl_id** | ChEMBL's unique id per compound; the key used to align fingerprints to labels. |
| **standard_value / standard_units / standard_type** | ChEMBL activity columns: the numeric measurement, its units (nM), and its type (`IC50`). |
| **data_raw.csv** | Unmodified ChEMBL IC50 records for CHEMBL220 (9415 × 46). |
| **data_preprocessed.csv** | Cleaned, labelled compounds + Lipinski descriptors + pIC50 (6220 × 8). |
| **PaDEL_descriptors_output.csv** | `Name` (ChEMBL id) + 881 PubChem fingerprint bits per molecule. |
| **.smi file** | Tab-separated (SMILES, id) input file that PaDEL reads. |
| **Descriptor-type preset** | PaDEL XML checklist selecting which descriptors to compute; ours (`PubchemFingerprinter.xml`) enables only PubChem fingerprints. |
| **cdd package** | The importable pipeline code in `src/cdd/` (`config`, `data`, `chem`, `analysis`, `features`). |
| **De-vendored** | Removing bundled third-party binaries (the 54 MB PaDEL jars) in favour of a pip dependency. |

## Tooling

| Term | Meaning |
| --- | --- |
| **PaDEL-Descriptor** | Java program that computes molecular descriptors and fingerprints. |
| **padelpy** | Python wrapper that bundles the PaDEL jar; used to compute fingerprints. Requires Java. |
| **RDKit** | Cheminformatics toolkit: parses SMILES and computes the Lipinski descriptors. |
| **chembl_webresource_client** | Official Python client for the ChEMBL REST API. |
| **scikit-learn** | ML library: `VarianceThreshold`, `train_test_split`, `RandomForestRegressor`. |
| **SciPy** | Provides the `mannwhitneyu` test. |
| **seaborn / Matplotlib** | Plotting libraries for all EDA and result figures. |
| **JRE (Java Runtime)** | Needed only to *recompute* fingerprints (PaDEL is a Java program); the cached CSV avoids it for normal runs. |
