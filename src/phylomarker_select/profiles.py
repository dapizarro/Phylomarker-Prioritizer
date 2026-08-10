"""Perfiles de panel y ranking por perfil."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .scoring import percentile_score


PROFILE_WEIGHTS = {
    "core_complete": {
        "coverage_percentile": 0.35,
        "alignment_percentile": 0.25,
        "trimming_percentile": 0.15,
        "information_percentile": 0.15,
        "low_bias_percentile": 0.10,
    },
    "backbone_balanced": {
        "coverage_percentile": 0.30,
        "clade_balance_percentile": 0.25,
        "alignment_percentile": 0.20,
        "trimming_percentile": 0.15,
        "information_percentile": 0.10,
    },
    "deep_robust": {
        "alignment_percentile": 0.25,
        "coverage_percentile": 0.20,
        "trimming_percentile": 0.20,
        "rate_percentile": 0.20,
        "low_bias_percentile": 0.15,
    },
    "low_bias": {
        "low_bias_percentile": 0.35,
        "alignment_percentile": 0.25,
        "trimming_percentile": 0.20,
        "coverage_percentile": 0.10,
        "information_percentile": 0.10,
    },
    "diverse_rate": {
        "coverage_percentile": 0.25,
        "alignment_percentile": 0.25,
        "trimming_percentile": 0.20,
        "information_percentile": 0.20,
        "low_bias_percentile": 0.10,
    },
    "occupancy_only": {
        "coverage_percentile": 1.00,
    },
    "information_only": {
        "information_percentile": 1.00,
    },
}


PROFILE_EXCLUDED_TRIMMING_CLASSES = {
    "core_complete": {"extreme"},
    "backbone_balanced": {"extreme"},
    "deep_robust": {"extreme", "signal_sensitive"},
    "low_bias": {"extreme"},
    "diverse_rate": {"extreme"},
    "occupancy_only": {"extreme"},
    "information_only": set(),
    "random_matched": set(),
}


def calculate_profile_ranking(
    scored: pd.DataFrame,
    profile: str,
) -> pd.DataFrame:
    if profile not in PROFILE_WEIGHTS:
        raise ValueError(
            f"Unknown panel profile: {profile}"
        )

    ranked = scored.copy()

    percentile_sources = {
        "coverage_percentile": ("coverage_score", True),
        "clade_balance_percentile": ("clade_balance_score", True),
        "alignment_percentile": ("alignment_score", True),
        "information_percentile": ("information_score", True),
        "rate_percentile": ("rate_score", True),
        "low_bias_percentile": ("bias_penalty", False),
        "trimming_percentile": ("trimming_stability_score", True),
    }

    for percentile_column, (source_column, higher_is_better) in (
        percentile_sources.items()
    ):
        if (
            percentile_column not in ranked.columns
            and source_column in ranked.columns
        ):
            ranked[percentile_column] = percentile_score(
                ranked[source_column],
                higher_is_better=higher_is_better,
            )

    weights = PROFILE_WEIGHTS[profile]

    total_score = np.zeros(
        len(ranked),
        dtype=float,
    )

    for feature, weight in weights.items():
        total_score += (
            weight
            * pd.to_numeric(
                ranked[feature],
                errors="coerce",
            )
            .fillna(0)
            .to_numpy()
        )

    ranked["profile"] = profile
    ranked["profile_score"] = total_score
    ranked["profile_eligible"] = ranked["eligible"].astype(bool)

    excluded_classes = PROFILE_EXCLUDED_TRIMMING_CLASSES.get(
        profile,
        set(),
    )

    if excluded_classes and "trimming_class" in ranked.columns:
        ranked["profile_eligible"] &= ~ranked[
            "trimming_class"
        ].isin(excluded_classes)

    ranked.loc[
        ~ranked["profile_eligible"],
        "profile_score",
    ] = -np.inf

    ranked = ranked.sort_values(
        [
            "profile_score",
            "gene_id",
        ],
        ascending=[
            False,
            True,
        ],
    ).reset_index(drop=True)

    ranked["rank"] = np.arange(
        1,
        len(ranked) + 1,
    )

    return ranked
