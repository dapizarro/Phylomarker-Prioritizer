"""Optimizadores greedy de panel con penalizacion de redundancia."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import RobustScaler


def feature_vectors(
    frame: pd.DataFrame,
) -> np.ndarray:
    features = [
        "taxon_occupancy",
        "min_replicated_group_occupancy",
        "alignment_length",
        "pis_per_length",
        "mean_pairwise_distance",
        "composition_variability",
        "retained_length_fraction",
    ]

    features = [
        feature
        for feature in features
        if feature in frame.columns
        and pd.to_numeric(
            frame[feature],
            errors="coerce",
        ).nunique(dropna=True) > 1
    ]

    if not features:
        return np.zeros((len(frame), 1), dtype=float)

    matrix = frame[features].apply(
        pd.to_numeric,
        errors="coerce",
    )

    imputed = SimpleImputer(
        strategy="median"
    ).fit_transform(matrix)

    return RobustScaler().fit_transform(
        imputed
    )


def optimize_panel_greedily(
    ranked: pd.DataFrame,
    panel_size: int,
    redundancy_penalty: float,
    candidate_top_fraction: float = 0.10,
    candidate_minimum_pool_size: int = 50,
    maximum_score_drop: float = 0.20,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligibility_column = (
        "profile_eligible"
        if "profile_eligible" in ranked.columns
        else "eligible"
    )
    candidates = ranked[
        ranked[eligibility_column]
    ].copy().reset_index(drop=True)

    if candidates.empty:
        return candidates, pd.DataFrame()

    candidates = candidates.sort_values(
        ["profile_score", "gene_id"],
        ascending=[False, True],
    ).reset_index(drop=True)

    pool_size = max(
        int(candidate_minimum_pool_size),
        int(math.ceil(len(candidates) * candidate_top_fraction)),
        int(panel_size * 5),
    )
    pool_size = min(pool_size, len(candidates))
    candidates = candidates.head(pool_size).copy()

    best_score = float(candidates["profile_score"].max())
    candidates = candidates[
        candidates["profile_score"]
        >= best_score - maximum_score_drop
    ].copy().reset_index(drop=True)

    if len(candidates) < panel_size:
        candidates = ranked[
            ranked[eligibility_column]
        ].sort_values(
            ["profile_score", "gene_id"],
            ascending=[False, True],
        ).head(max(panel_size, pool_size)).copy().reset_index(drop=True)

    requested_size = min(panel_size, len(candidates))

    vectors = feature_vectors(candidates)

    norms = np.linalg.norm(
        vectors,
        axis=1,
    )

    norms[norms == 0] = 1.0

    normalized = vectors / norms[:, None]
    similarity_matrix = normalized @ normalized.T

    selected_indices: list[int] = []
    trace_rows: list[dict] = []

    for step in range(requested_size):
        best_index: int | None = None
        best_gain = -np.inf
        best_redundancy = 0.0

        for index, row in candidates.iterrows():
            if index in selected_indices:
                continue

            redundancy = (
                max(
                    0.0,
                    max(
                        float(
                            similarity_matrix[
                                index,
                                selected_index,
                            ]
                        )
                        for selected_index in selected_indices
                    ),
                )
                if selected_indices
                else 0.0
            )

            marginal_gain = (
                float(row["profile_score"])
                - redundancy_penalty
                * redundancy
            )

            if marginal_gain > best_gain:
                best_gain = marginal_gain
                best_index = index
                best_redundancy = redundancy

        if best_index is None:
            break

        selected_indices.append(best_index)

        selected_row = candidates.iloc[
            best_index
        ]

        trace_rows.append(
            {
                "step": step + 1,
                "gene_id": selected_row["gene_id"],
                "profile_score": selected_row[
                    "profile_score"
                ],
                "maximum_redundancy": best_redundancy,
                "redundancy_penalty": redundancy_penalty,
                "marginal_gain": best_gain,
                "trimming_class": selected_row.get(
                    "trimming_class",
                    "unknown",
                ),
                "trimming_stability_score": selected_row.get(
                    "trimming_stability_score",
                    np.nan,
                ),
                "candidate_pool_size": len(candidates),
            }
        )

    panel = candidates.iloc[
        selected_indices
    ].copy()

    panel["panel_order"] = np.arange(
        1,
        len(panel) + 1,
    )

    return panel, pd.DataFrame(trace_rows)



def optimize_diverse_rate_panel(
    ranked: pd.DataFrame,
    panel_size: int,
    redundancy_penalty: float,
    candidate_top_fraction: float = 0.10,
    candidate_minimum_pool_size: int = 50,
    maximum_score_drop: float = 0.20,
    number_rate_bins: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Select a quality-controlled panel with balanced coverage of rate bins.

    Genes are divided into quantile bins using mean_pairwise_distance.
    Selection proceeds round-robin across bins while applying the same
    redundancy penalty used by the standard greedy optimizer.
    """
    eligibility_column = (
        "profile_eligible"
        if "profile_eligible" in ranked.columns
        else "eligible"
    )
    candidates = ranked[
        ranked[eligibility_column]
    ].copy().reset_index(drop=True)

    if candidates.empty:
        return candidates, pd.DataFrame()

    candidates = candidates.sort_values(
        ["profile_score", "gene_id"],
        ascending=[False, True],
    ).reset_index(drop=True)

    pool_size = max(
        int(candidate_minimum_pool_size),
        int(math.ceil(len(candidates) * candidate_top_fraction)),
        int(panel_size * 5),
    )
    pool_size = min(pool_size, len(candidates))
    quality_pool = candidates.head(pool_size).copy()

    best_score = float(quality_pool["profile_score"].max())
    quality_pool = quality_pool[
        quality_pool["profile_score"]
        >= best_score - maximum_score_drop
    ].copy()

    if len(quality_pool) < panel_size:
        quality_pool = candidates.head(
            max(panel_size, pool_size)
        ).copy()

    rate_values = pd.to_numeric(
        quality_pool["mean_pairwise_distance"],
        errors="coerce",
    )
    valid_rate_count = int(rate_values.notna().sum())
    q = min(
        int(number_rate_bins),
        max(1, valid_rate_count),
    )

    if q <= 1 or rate_values.nunique(dropna=True) <= 1:
        panel, trace = optimize_panel_greedily(
            ranked=ranked,
            panel_size=panel_size,
            redundancy_penalty=redundancy_penalty,
            candidate_top_fraction=candidate_top_fraction,
            candidate_minimum_pool_size=candidate_minimum_pool_size,
            maximum_score_drop=maximum_score_drop,
        )
        panel["rate_bin"] = "single_bin"
        if not trace.empty:
            trace["rate_bin"] = "single_bin"
        return panel, trace

    try:
        rate_bins = pd.qcut(
            rate_values,
            q=q,
            labels=False,
            duplicates="drop",
        )
    except ValueError:
        rate_bins = pd.Series(
            0,
            index=quality_pool.index,
            dtype="Int64",
        )

    quality_pool = quality_pool.reset_index(drop=True)
    quality_pool["rate_bin"] = pd.Series(
        rate_bins,
        index=quality_pool.index,
    ).astype("Int64")

    if quality_pool["rate_bin"].isna().any():
        median_bin = int(
            quality_pool["rate_bin"].dropna().median()
        )
        quality_pool["rate_bin"] = (
            quality_pool["rate_bin"].fillna(median_bin)
        )

    available_bins = sorted(
        int(value)
        for value in quality_pool["rate_bin"].unique()
    )
    requested_size = min(panel_size, len(quality_pool))
    base_quota = requested_size // len(available_bins)
    remainder = requested_size % len(available_bins)

    quotas = {
        rate_bin: base_quota
        + (1 if position < remainder else 0)
        for position, rate_bin in enumerate(available_bins)
    }

    vectors = feature_vectors(quality_pool)
    norms = np.linalg.norm(vectors, axis=1)
    norms[norms == 0] = 1.0
    normalized = vectors / norms[:, None]
    similarity_matrix = normalized @ normalized.T

    selected_indices: list[int] = []
    selected_by_bin = {rate_bin: 0 for rate_bin in available_bins}
    trace_rows: list[dict] = []

    def choose_best(allowed_indices: list[int]) -> tuple[int | None, float, float]:
        best_index: int | None = None
        best_gain = -np.inf
        best_redundancy = 0.0

        for index in allowed_indices:
            if index in selected_indices:
                continue

            row = quality_pool.iloc[index]
            redundancy = (
                max(
                    0.0,
                    max(
                        float(similarity_matrix[index, selected])
                        for selected in selected_indices
                    ),
                )
                if selected_indices
                else 0.0
            )
            gain = (
                float(row["profile_score"])
                - redundancy_penalty * redundancy
            )

            if (
                gain > best_gain
                or (
                    math.isclose(gain, best_gain)
                    and best_index is not None
                    and str(row["gene_id"])
                    < str(quality_pool.iloc[best_index]["gene_id"])
                )
            ):
                best_index = index
                best_gain = gain
                best_redundancy = redundancy

        return best_index, best_gain, best_redundancy

    while len(selected_indices) < requested_size:
        progress = False

        for rate_bin in available_bins:
            if len(selected_indices) >= requested_size:
                break
            if selected_by_bin[rate_bin] >= quotas[rate_bin]:
                continue

            allowed = quality_pool.index[
                quality_pool["rate_bin"] == rate_bin
            ].tolist()
            best_index, best_gain, best_redundancy = choose_best(allowed)

            if best_index is None:
                continue

            selected_indices.append(best_index)
            selected_by_bin[rate_bin] += 1
            selected_row = quality_pool.iloc[best_index]

            trace_rows.append(
                {
                    "step": len(selected_indices),
                    "gene_id": selected_row["gene_id"],
                    "profile_score": selected_row["profile_score"],
                    "maximum_redundancy": best_redundancy,
                    "redundancy_penalty": redundancy_penalty,
                    "marginal_gain": best_gain,
                    "trimming_class": selected_row.get(
                        "trimming_class",
                        "unknown",
                    ),
                    "trimming_stability_score": selected_row.get(
                        "trimming_stability_score",
                        np.nan,
                    ),
                    "candidate_pool_size": len(quality_pool),
                    "rate_bin": int(rate_bin),
                    "rate_value": selected_row.get(
                        "mean_pairwise_distance",
                        np.nan,
                    ),
                }
            )
            progress = True

        if not progress:
            remaining = [
                index
                for index in quality_pool.index
                if index not in selected_indices
            ]
            best_index, best_gain, best_redundancy = choose_best(remaining)
            if best_index is None:
                break

            selected_indices.append(best_index)
            selected_row = quality_pool.iloc[best_index]
            rate_bin = int(selected_row["rate_bin"])
            selected_by_bin[rate_bin] += 1

            trace_rows.append(
                {
                    "step": len(selected_indices),
                    "gene_id": selected_row["gene_id"],
                    "profile_score": selected_row["profile_score"],
                    "maximum_redundancy": best_redundancy,
                    "redundancy_penalty": redundancy_penalty,
                    "marginal_gain": best_gain,
                    "trimming_class": selected_row.get(
                        "trimming_class",
                        "unknown",
                    ),
                    "trimming_stability_score": selected_row.get(
                        "trimming_stability_score",
                        np.nan,
                    ),
                    "candidate_pool_size": len(quality_pool),
                    "rate_bin": rate_bin,
                    "rate_value": selected_row.get(
                        "mean_pairwise_distance",
                        np.nan,
                    ),
                }
            )

    panel = quality_pool.iloc[selected_indices].copy()
    panel["panel_order"] = np.arange(1, len(panel) + 1)

    return panel, pd.DataFrame(trace_rows)
