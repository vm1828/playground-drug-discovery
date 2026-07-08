"""Retrieve and preprocess ChEMBL bioactivity data."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config


def search_targets(query: str = config.TARGET_QUERY) -> pd.DataFrame:
    """Search ChEMBL for targets matching ``query`` (requires network access)."""
    from chembl_webresource_client.new_client import new_client

    return pd.DataFrame.from_dict(new_client.target.search(query))


def fetch_bioactivity(target_chembl_id: str = config.TARGET_CHEMBL_ID) -> pd.DataFrame:
    """Download IC50 bioactivity records for a target (requires network access)."""
    from chembl_webresource_client.new_client import new_client

    activity = new_client.activity
    records = activity.filter(target_chembl_id=target_chembl_id).filter(
        standard_type="IC50"
    )
    return pd.DataFrame.from_dict(records)


def classify_bioactivity(ic50_nm: pd.Series) -> pd.Series:
    """Label IC50 values (nM) as active / inactive / intermediate."""
    conditions = [
        ic50_nm <= config.ACTIVE_THRESHOLD_NM,
        ic50_nm >= config.INACTIVE_THRESHOLD_NM,
    ]
    labels = np.select(conditions, ["active", "inactive"], default="intermediate")
    return pd.Series(labels, index=ic50_nm.index)


def preprocess(df: pd.DataFrame) -> pd.DataFrame:
    """Clean raw ChEMBL records and keep only active/inactive compounds.

    Drops rows missing an IC50 or SMILES, coerces IC50 to numeric, replaces
    zeros with a tiny value (so pIC50 stays finite), labels each compound, and
    discards the ambiguous "intermediate" class.
    """
    df = df.dropna(subset=["standard_value", "canonical_smiles"]).copy()
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value"])
    df["standard_value"] = df["standard_value"].replace(0, 1e-10)

    df["bioactivity_class"] = classify_bioactivity(df["standard_value"])
    df = df[df["bioactivity_class"] != "intermediate"]

    columns = [
        "molecule_chembl_id",
        "canonical_smiles",
        "bioactivity_class",
        "standard_value",
    ]
    return df[columns].reset_index(drop=True)
