"""Project paths and constants.

Paths are derived from this file's location, so they resolve correctly no
matter the working directory (VS Code, a script, or Colab).
"""
from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data"
PLOTS_DIR = PROJECT_ROOT / "plots"
PADEL_DIR = PROJECT_ROOT / "padel"

RAW_CSV = DATA_DIR / "data_raw.csv"
PREPROCESSED_CSV = DATA_DIR / "data_preprocessed.csv"
PADEL_OUTPUT_CSV = DATA_DIR / "PaDEL_descriptors_output.csv"
SMI_FILE = DATA_DIR / "molecule.smi"

# PaDEL PubChem fingerprint descriptor-type preset (881-bit fingerprints).
PUBCHEM_FINGERPRINTER_XML = PADEL_DIR / "PubchemFingerprinter.xml"

# ChEMBL target: Human Acetylcholinesterase.
TARGET_QUERY = "acetylcholinesterase"
TARGET_CHEMBL_ID = "CHEMBL220"

# Lipinski Rule-of-Five descriptors.
LIPINSKI_DESCRIPTORS = ["MolWt", "MolLogP", "NumHDonors", "NumHAcceptors"]

# IC50 (nM) thresholds used to label compounds by bioactivity.
ACTIVE_THRESHOLD_NM = 1_000
INACTIVE_THRESHOLD_NM = 10_000
