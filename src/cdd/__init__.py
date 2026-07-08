"""Computational drug discovery pipeline for acetylcholinesterase inhibitors.

Modules:
    config    - project paths and constants
    data      - ChEMBL retrieval and preprocessing
    chem      - Lipinski descriptors and pIC50
    analysis  - exploratory plots and Mann-Whitney tests
    features  - PaDEL fingerprints and feature-matrix prep
"""
from __future__ import annotations

from . import analysis, chem, config, data, features

__all__ = ["analysis", "chem", "config", "data", "features"]
