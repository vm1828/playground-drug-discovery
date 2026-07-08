"""Exploratory plots and the Mann-Whitney U test for descriptor differences."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

from . import config


def _save(figure: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=300, bbox_inches="tight")


def plot_bioactivity_counts(df: pd.DataFrame, save: bool = True) -> plt.Axes:
    """Frequency of the active vs inactive classes."""
    ax = sns.countplot(x="bioactivity_class", data=df)
    if save:
        _save(ax.figure, config.PLOTS_DIR / "count_bioactivity_class.pdf")
    return ax


def plot_mw_vs_logp(df: pd.DataFrame, save: bool = True) -> plt.Axes:
    """Chemical space: molecular weight vs LogP, sized by potency."""
    ax = sns.scatterplot(
        x="MolWt", y="MolLogP", data=df,
        hue="bioactivity_class", size="pIC50", alpha=0.7,
    )
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    if save:
        _save(ax.figure, config.PLOTS_DIR / "scatter_MW_LogP.pdf")
    return ax


def mannwhitney(
    df: pd.DataFrame, descriptor: str, alpha: float = 0.05
) -> tuple[float, float, bool]:
    """Compare a descriptor across active/inactive classes.

    Returns ``(U statistic, p-value, reject_null)``.
    """
    is_active = df["bioactivity_class"] == "active"
    stat, p = mannwhitneyu(df[is_active][descriptor], df[~is_active][descriptor])
    return stat, p, p <= alpha


def plot_descriptor(
    df: pd.DataFrame, descriptor: str, alpha: float = 0.05, save: bool = True
) -> plt.Axes:
    """Boxplot of a descriptor by class, annotated with the Mann-Whitney result."""
    fig, ax = plt.subplots(figsize=(3, 3))
    sns.boxplot(x="bioactivity_class", y=descriptor, data=df, ax=ax)

    stat, p, reject = mannwhitney(df, descriptor, alpha)
    annotation = (
        f"U = {stat:.2f}, p = {p:.3e}, α = {alpha}\n"
        f"{'Reject H0' if reject else 'Fail to reject H0'}"
    )
    ax.text(
        0.5, -0.6, annotation, transform=ax.transAxes, ha="center", fontsize=10,
        bbox=dict(facecolor="white", alpha=0.6),
    )
    if save:
        _save(fig, config.PLOTS_DIR / f"box_mannwhitney_{descriptor}.pdf")
    return ax
