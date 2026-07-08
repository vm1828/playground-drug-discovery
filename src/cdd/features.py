"""Molecular fingerprints (PaDEL) and feature-matrix preparation."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from sklearn.feature_selection import VarianceThreshold

from . import config

# Default VarianceThreshold cutoff: drop binary features that are the same value
# in more than 80% of molecules.
DEFAULT_VARIANCE_THRESHOLD = 0.8 * (1 - 0.8)


def write_smi(df: pd.DataFrame, path: Path | None = None) -> Path:
    """Write a tab-separated ``.smi`` file (SMILES, id) for PaDEL."""
    path = path or config.SMI_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    df[["canonical_smiles", "molecule_chembl_id"]].to_csv(
        path, sep="\t", index=False, header=False
    )
    return path


def compute_padel_fingerprints(
    smi_path: Path | None = None, output_csv: Path | None = None
) -> Path:
    """Compute 881-bit PubChem fingerprints via padelpy (requires Java).

    padelpy bundles the PaDEL-Descriptor jar, so nothing is vendored here beyond
    the small XML descriptor-type preset in ``padel/``.
    """
    from padelpy import padeldescriptor

    smi_path = smi_path or config.SMI_FILE
    output_csv = output_csv or config.PADEL_OUTPUT_CSV
    padeldescriptor(
        mol_dir=str(smi_path),
        d_file=str(output_csv),
        descriptortypes=str(config.PUBCHEM_FINGERPRINTER_XML),
        fingerprints=True,
        removesalt=True,
        standardizenitro=True,
        detectaromaticity=True,
        threads=-1,
    )
    return output_csv


def load_fingerprints(padel_csv: Path | None = None) -> pd.DataFrame:
    """Load PaDEL fingerprints indexed by molecule id.

    PaDEL does not preserve input order, and a molecule may appear multiple
    times (repeated IC50 measurements). Fingerprints depend only on structure,
    so duplicates are identical; we keep one row per id to allow alignment.
    """
    padel_csv = padel_csv or config.PADEL_OUTPUT_CSV
    fingerprints = pd.read_csv(padel_csv).drop_duplicates("Name").set_index("Name")
    fingerprints.index.name = "molecule_chembl_id"
    return fingerprints


def prepare_xy(
    features_table: pd.DataFrame,
    fingerprints: pd.DataFrame | None = None,
    target: str = "pIC50",
    variance_threshold: float = DEFAULT_VARIANCE_THRESHOLD,
) -> tuple[pd.DataFrame, pd.Series]:
    """Build aligned (X, y) for modelling.

    X (fingerprints) is aligned to y (``target``) row-by-row **by molecule id**,
    then low-variance features are dropped. This avoids the positional-alignment
    pitfall of pairing PaDEL output with labels by row number.
    """
    if fingerprints is None:
        fingerprints = load_fingerprints()

    aligned = fingerprints.reindex(features_table["molecule_chembl_id"].to_numpy())
    if aligned.isna().any().any():
        missing = int(aligned.isna().any(axis=1).sum())
        raise ValueError(f"{missing} molecules have no fingerprint; recompute PaDEL.")

    y = features_table[target].reset_index(drop=True)
    selector = VarianceThreshold(threshold=variance_threshold)
    kept = selector.fit_transform(aligned)
    X = pd.DataFrame(kept, columns=selector.get_feature_names_out())
    return X, y
