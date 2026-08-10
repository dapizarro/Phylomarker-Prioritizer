"""Estandarizacion robusta, percentiles y puntuaciones biologicas."""
from __future__ import annotations

import numpy as np
import pandas as pd


def robust_standardize(
    series: pd.Series,
    higher_is_better: bool = True,
) -> pd.Series:
    numeric = pd.to_numeric(
        series,
        errors="coerce",
    )

    median = numeric.median()
    mad = (numeric - median).abs().median()

    if pd.isna(mad) or mad == 0:
        standardized = pd.Series(
            np.zeros(len(numeric)),
            index=numeric.index,
        )
    else:
        standardized = (
            numeric - median
        ) / (1.4826 * mad)

    standardized = standardized.clip(-4, 4)

    return (
        standardized
        if higher_is_better
        else -standardized
    )


def percentile_score(
    series: pd.Series,
    higher_is_better: bool = True,
) -> pd.Series:
    """Transform a numeric series into comparable scores in [0, 1]."""
    numeric = pd.to_numeric(series, errors="coerce")
    observed = numeric.notna().sum()

    if observed <= 1 or numeric.nunique(dropna=True) <= 1:
        return pd.Series(0.5, index=series.index, dtype=float)

    scores = numeric.rank(
        method="average",
        pct=True,
        na_option="keep",
    )

    if not higher_is_better:
        scores = 1.0 - scores

    return scores.fillna(0.5).clip(0.0, 1.0)


def add_biological_scores(
    metrics: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    scored = metrics.copy()

    scored["coverage_score"] = (
        0.45
        * robust_standardize(
            scored["taxon_occupancy"]
        )
        + 0.35
        * robust_standardize(
            scored["sequence_completeness"]
        )
        + 0.20
        * robust_standardize(
            scored["cell_occupancy"]
        )
    )

    if {
        "min_group_sequence_completeness",
        "sd_group_sequence_completeness",
        "worst_group_gap_fraction",
    }.issubset(scored.columns):
        scored["clade_balance_score"] = (
            0.50
            * robust_standardize(
                scored["min_group_sequence_completeness"]
            )
            + 0.30
            * robust_standardize(
                scored["sd_group_sequence_completeness"],
                higher_is_better=False,
            )
            + 0.20
            * robust_standardize(
                scored["worst_group_gap_fraction"],
                higher_is_better=False,
            )
        )
    else:
        scored["clade_balance_score"] = (
            0.60
            * robust_standardize(
                scored["mean_group_occupancy"]
            )
            + 0.40
            * robust_standardize(
                scored[
                    "min_replicated_group_occupancy"
                ]
            )
        )

    scored["alignment_score"] = (
        0.45
        * robust_standardize(
            scored["alignment_length"]
        )
        + 0.30
        * robust_standardize(
            scored["gap_fraction"],
            higher_is_better=False,
        )
        + 0.25
        * robust_standardize(
            scored["ambiguous_fraction"],
            higher_is_better=False,
        )
    )

    scored["information_score"] = (
        0.45
        * robust_standardize(
            scored["pis_per_length"]
        )
        + 0.25
        * robust_standardize(
            scored["variable_sites_per_length"]
        )
        + 0.20
        * robust_standardize(
            scored["mean_entropy"]
        )
        + 0.10
        * robust_standardize(
            scored["parsimony_informative_sites"]
        )
    )

    pairwise_distances = pd.to_numeric(
        scored["mean_pairwise_distance"],
        errors="coerce",
    )

    median_rate = pairwise_distances.median()
    rate_deviation = (
        pairwise_distances - median_rate
    ).abs()

    scored["rate_score"] = robust_standardize(
        rate_deviation,
        higher_is_better=False,
    )

    scored["bias_penalty"] = (
        0.55
        * robust_standardize(
            scored["composition_variability"]
        )
        + 0.30
        * robust_standardize(
            scored["gap_fraction"]
        )
        + 0.15
        * robust_standardize(
            scored["ambiguous_fraction"]
        )
    )

    scored["coverage_percentile"] = percentile_score(
        scored["coverage_score"]
    )
    scored["clade_balance_percentile"] = percentile_score(
        scored["clade_balance_score"]
    )
    scored["alignment_percentile"] = percentile_score(
        scored["alignment_score"]
    )
    scored["information_percentile"] = percentile_score(
        scored["information_score"]
    )
    scored["rate_percentile"] = percentile_score(
        scored["rate_score"]
    )
    scored["low_bias_percentile"] = percentile_score(
        scored["bias_penalty"],
        higher_is_better=False,
    )
    scored["trimming_percentile"] = percentile_score(
        scored["trimming_stability_score"]
        if "trimming_stability_score" in scored.columns
        else pd.Series(1.0, index=scored.index)
    )

    eligibility = config.get(
        "eligibility",
        {},
    )

    minimum_taxa = int(
        eligibility.get(
            "min_taxa",
            4,
        )
    )

    minimum_occupancy = float(
        eligibility.get(
            "min_taxon_occupancy",
            0.5,
        )
    )

    minimum_length = int(
        eligibility.get(
            "min_alignment_length",
            80,
        )
    )

    minimum_pis = int(
        eligibility.get(
            "min_informative_sites",
            2,
        )
    )

    maximum_gaps = float(
        eligibility.get(
            "max_gap_fraction",
            0.6,
        )
    )

    maximum_ambiguous = float(
        eligibility.get(
            "max_ambiguous_fraction",
            0.1,
        )
    )

    scored["eligible"] = (
        (scored["n_taxa"] >= minimum_taxa)
        & (
            scored["taxon_occupancy"]
            >= minimum_occupancy
        )
        & (
            scored["alignment_length"]
            >= minimum_length
        )
        & (
            scored["parsimony_informative_sites"]
            >= minimum_pis
        )
        & (
            scored["gap_fraction"]
            <= maximum_gaps
        )
        & (
            scored["ambiguous_fraction"]
            <= maximum_ambiguous
        )
    )

    exclusion_reasons: list[str] = []

    for row in scored.itertuples(index=False):
        failures: list[str] = []

        if row.n_taxa < minimum_taxa:
            failures.append("min_taxa")

        if row.taxon_occupancy < minimum_occupancy:
            failures.append("taxon_occupancy")

        if row.alignment_length < minimum_length:
            failures.append("alignment_length")

        if row.parsimony_informative_sites < minimum_pis:
            failures.append("informative_sites")

        if row.gap_fraction > maximum_gaps:
            failures.append("gap_fraction")

        if row.ambiguous_fraction > maximum_ambiguous:
            failures.append("ambiguous_fraction")

        exclusion_reasons.append(
            ";".join(failures)
        )

    scored["exclusion_reason"] = exclusion_reasons

    return scored
