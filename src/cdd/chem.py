"""Cheminformatics features: Lipinski descriptors and pIC50."""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from rdkit.Chem import Descriptors, MolFromSmiles

from . import config


def lipinski_descriptors(
    smiles: pd.Series, descriptors: Sequence[str] | None = None
) -> pd.DataFrame:
    """Compute RDKit Lipinski descriptors for a series of SMILES strings."""
    descriptors = list(descriptors or config.LIPINSKI_DESCRIPTORS)
    molecules = smiles.map(MolFromSmiles)
    return pd.DataFrame(
        {name: molecules.map(getattr(Descriptors, name)) for name in descriptors},
        index=smiles.index,
    )


def to_pic50(ic50_nm: pd.Series) -> pd.Series:
    """Convert IC50 in nanomolar to pIC50 = -log10(IC50 in molar)."""
    return ic50_nm.map(lambda value: -np.log10(value * 1e-9))


def build_features_table(df: pd.DataFrame) -> pd.DataFrame:
    """Attach Lipinski descriptors and pIC50, dropping the raw IC50 column.

    Expects the preprocessed frame from :func:`cdd.data.preprocess`.
    """
    out = df.assign(**lipinski_descriptors(df["canonical_smiles"]))
    out["pIC50"] = to_pic50(out["standard_value"])
    return out.drop(columns=["standard_value"])
