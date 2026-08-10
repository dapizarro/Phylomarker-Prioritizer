"""PCA exploratorio. No alimenta el ranking."""
from __future__ import annotations

import logging
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import RobustScaler

from .layout import OutputLayout

LOGGER = logging.getLogger("phylomarker-select")


def run_exploratory_pca(
    scored: pd.DataFrame,
    layout: OutputLayout,
) -> None:
    features = [
        "taxon_occupancy",
        "sequence_completeness",
        "gap_fraction",
        "ambiguous_fraction",
        "alignment_length",
        "pis_per_length",
        "mean_entropy",
        "mean_pairwise_distance",
        "composition_variability",
        "mean_group_occupancy",
    ]

    available_features = [
        feature
        for feature in features
        if feature in scored.columns
        and pd.to_numeric(
            scored[feature],
            errors="coerce",
        ).nunique() > 1
    ]

    if (
        len(available_features) < 2
        or len(scored) < 3
    ):
        LOGGER.warning(
            "PCA skipped because there is insufficient variation."
        )
        return

    matrix = scored[
        available_features
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    imputed = SimpleImputer(
        strategy="median"
    ).fit_transform(matrix)

    scaled = RobustScaler().fit_transform(
        imputed
    )

    number_components = min(
        5,
        scaled.shape[0],
        scaled.shape[1],
    )

    pca = PCA(
        n_components=number_components,
        random_state=0,
    )

    scores = pca.fit_transform(scaled)
    loadings = pca.components_.T.copy()

    for component_index in range(
        loadings.shape[1]
    ):
        pivot = int(
            np.argmax(
                np.abs(
                    loadings[:, component_index]
                )
            )
        )

        if loadings[
            pivot,
            component_index,
        ] < 0:
            loadings[
                :,
                component_index,
            ] *= -1

            scores[
                :,
                component_index,
            ] *= -1

    layout.pca_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    score_frame = pd.DataFrame(
        scores,
        columns=[
            f"PC{index + 1}"
            for index in range(
                number_components
            )
        ],
    )

    score_frame.insert(
        0,
        "gene_id",
        scored["gene_id"].values,
    )

    score_frame.to_csv(
        layout.pca_scores_table,
        sep="\t",
        index=False,
    )

    loading_frame = pd.DataFrame(
        loadings,
        index=available_features,
        columns=[
            f"PC{index + 1}"
            for index in range(
                number_components
            )
        ],
    )

    loading_frame.index.name = "metric"

    loading_frame.to_csv(
        layout.pca_loadings_table,
        sep="\t",
    )

    variance_frame = pd.DataFrame(
        {
            "component": [
                f"PC{index + 1}"
                for index in range(
                    number_components
                )
            ],
            "explained_variance_ratio": (
                pca.explained_variance_ratio_
            ),
            "cumulative_explained_variance": (
                np.cumsum(
                    pca.explained_variance_ratio_
                )
            ),
        }
    )

    variance_frame.to_csv(
        layout.pca_explained_variance_table,
        sep="\t",
        index=False,
    )
