"""Construccion, resumen y exportacion de paneles."""
from __future__ import annotations

import logging
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .layout import OutputLayout
from .optimize import optimize_diverse_rate_panel, optimize_panel_greedily
from .profiles import calculate_profile_ranking

LOGGER = logging.getLogger("phylomarker-select")


def summarize_panel(panel: pd.DataFrame) -> pd.DataFrame:
    """Return one-row panel-level diagnostics."""
    if panel.empty:
        return pd.DataFrame(
            [{
                "n_genes": 0,
                "mean_profile_score": np.nan,
                "minimum_profile_score": np.nan,
                "total_alignment_length": 0,
                "total_pis": 0,
                "minimum_taxon_occupancy": np.nan,
                "mean_trimming_stability_score": np.nan,
                "maximum_composition_variability": np.nan,
            }]
        )

    return pd.DataFrame(
        [{
            "n_genes": int(len(panel)),
            "mean_profile_score": float(
                pd.to_numeric(
                    panel.get("profile_score"),
                    errors="coerce",
                ).mean()
            ),
            "minimum_profile_score": float(
                pd.to_numeric(
                    panel.get("profile_score"),
                    errors="coerce",
                ).min()
            ),
            "total_alignment_length": int(
                pd.to_numeric(
                    panel.get("alignment_length"),
                    errors="coerce",
                ).fillna(0).sum()
            ),
            "total_pis": int(
                pd.to_numeric(
                    panel.get("parsimony_informative_sites"),
                    errors="coerce",
                ).fillna(0).sum()
            ),
            "minimum_taxon_occupancy": float(
                pd.to_numeric(
                    panel.get("taxon_occupancy"),
                    errors="coerce",
                ).min()
            ),
            "mean_trimming_stability_score": float(
                pd.to_numeric(
                    panel.get("trimming_stability_score"),
                    errors="coerce",
                ).mean()
            ),
            "maximum_composition_variability": float(
                pd.to_numeric(
                    panel.get("composition_variability"),
                    errors="coerce",
                ).max()
            ),
            "minimum_group_sequence_completeness": float(
                pd.to_numeric(
                    panel.get("min_group_sequence_completeness"),
                    errors="coerce",
                ).min()
            )
            if "min_group_sequence_completeness" in panel.columns
            else np.nan,
            "maximum_group_completeness_sd": float(
                pd.to_numeric(
                    panel.get("sd_group_sequence_completeness"),
                    errors="coerce",
                ).max()
            )
            if "sd_group_sequence_completeness" in panel.columns
            else np.nan,
            "rate_min": float(
                pd.to_numeric(
                    panel.get("mean_pairwise_distance"),
                    errors="coerce",
                ).min()
            ),
            "rate_median": float(
                pd.to_numeric(
                    panel.get("mean_pairwise_distance"),
                    errors="coerce",
                ).median()
            ),
            "rate_max": float(
                pd.to_numeric(
                    panel.get("mean_pairwise_distance"),
                    errors="coerce",
                ).max()
            ),
            "rate_sd": float(
                pd.to_numeric(
                    panel.get("mean_pairwise_distance"),
                    errors="coerce",
                ).std(ddof=1)
            ),
            "n_rate_bins": int(
                panel["rate_bin"].nunique(dropna=True)
            )
            if "rate_bin" in panel.columns
            else 0,
        }]
    )


def random_panel(
    scored: pd.DataFrame,
    panel_size: int,
    seed: int,
) -> pd.DataFrame:
    candidates = scored[
        scored["eligible"]
    ].copy()

    if candidates.empty:
        return candidates

    random_generator = np.random.default_rng(
        seed
    )

    selected_indices = random_generator.choice(
        candidates.index.to_numpy(),
        size=min(
            panel_size,
            len(candidates),
        ),
        replace=False,
    )

    panel = candidates.loc[
        selected_indices
    ].copy()

    panel = panel.sort_values(
        "gene_id"
    ).reset_index(drop=True)

    panel["panel_order"] = np.arange(
        1,
        len(panel) + 1,
    )

    panel["profile"] = "random_matched"
    panel["profile_score"] = np.nan

    return panel


def export_panel_alignments(
    panel: pd.DataFrame,
    alignment_directory: Path,
    destination: Path,
    sequence_type: str,
) -> None:
    extension = ".faa" if sequence_type == "protein" else ".fna"

    destination.mkdir(
        parents=True,
        exist_ok=True,
    )

    for gene_id in panel["gene_id"]:
        candidates = [
            alignment_directory
            / f"{gene_id}.trimmed{extension}",
            alignment_directory
            / f"{gene_id}.aln{extension}",
        ]

        source = next(
            (
                candidate
                for candidate in candidates
                if candidate.is_file()
            ),
            None,
        )

        if source is None:
            LOGGER.warning(
                "Alignment not found for marker %s",
                gene_id,
            )
            continue

        shutil.copy2(
            source,
            destination / source.name,
        )


def create_panels(
    scored: pd.DataFrame,
    layout: OutputLayout,
    config: dict,
    trimming_enabled: bool,
) -> None:
    sequence_type = layout.sequence_type

    alignment_directory = layout.analysis_alignment_directory(
        trimming_enabled
    )

    panel_config = config.get(
        "panels",
        {},
    )

    profiles = panel_config.get(
        "profiles",
        [
            "core_complete",
            "backbone_balanced",
            "deep_robust",
            "low_bias",
            "diverse_rate",
            "occupancy_only",
            "information_only",
            "random_matched",
        ],
    )

    sizes = [
        int(value)
        for value in panel_config.get(
            "sizes",
            [25, 50, 100],
        )
    ]

    redundancy_penalty = float(
        panel_config.get(
            "redundancy_penalty",
            0.10,
        )
    )

    candidate_top_fraction = float(
        panel_config.get("candidate_top_fraction", 0.10)
    )
    candidate_minimum_pool_size = int(
        panel_config.get("candidate_minimum_pool_size", 50)
    )
    maximum_score_drop = float(
        panel_config.get("maximum_score_drop", 0.20)
    )

    random_seed = int(
        config.get(
            "project",
            {},
        ).get(
            "random_seed",
            20260729,
        )
    )

    layout.panels_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    layout.rankings_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    for profile in profiles:
        if profile == "random_matched":
            for size in sizes:
                panel = random_panel(
                    scored,
                    size,
                    random_seed + size,
                )

                layout.panel_directory(
                    profile,
                    size,
                ).mkdir(
                    parents=True,
                    exist_ok=True,
                )

                panel.to_csv(
                    layout.panel_genes_table(profile, size),
                    sep="\t",
                    index=False,
                )

                summarize_panel(panel).to_csv(
                    layout.panel_summary_table(profile, size),
                    sep="\t",
                    index=False,
                )

                panel["gene_id"].to_csv(
                    layout.panel_genes_file(profile, size),
                    index=False,
                    header=False,
                )

                export_panel_alignments(
                    panel,
                    alignment_directory,
                    layout.panel_alignments_directory(
                        profile,
                        size,
                    ),
                    sequence_type,
                )

            continue

        ranked = calculate_profile_ranking(
            scored,
            profile,
        )

        ranked.to_csv(
            layout.profile_ranking_table(profile),
            sep="\t",
            index=False,
        )

        for size in sizes:
            if profile == "diverse_rate":
                panel, trace = optimize_diverse_rate_panel(
                    ranked,
                    size,
                    redundancy_penalty,
                    candidate_top_fraction,
                    candidate_minimum_pool_size,
                    maximum_score_drop,
                    number_rate_bins=int(
                        panel_config.get("diverse_rate_bins", 5)
                    ),
                )
            else:
                panel, trace = optimize_panel_greedily(
                    ranked,
                    size,
                    redundancy_penalty,
                    candidate_top_fraction,
                    candidate_minimum_pool_size,
                    maximum_score_drop,
                )

            layout.panel_directory(
                profile,
                size,
            ).mkdir(
                parents=True,
                exist_ok=True,
            )

            panel.to_csv(
                layout.panel_genes_table(profile, size),
                sep="\t",
                index=False,
            )

            summarize_panel(panel).to_csv(
                layout.panel_summary_table(profile, size),
                sep="\t",
                index=False,
            )

            trace.to_csv(
                layout.selection_trace_table(profile, size),
                sep="\t",
                index=False,
            )

            panel["gene_id"].to_csv(
                layout.panel_genes_file(profile, size),
                index=False,
                header=False,
            )

            export_panel_alignments(
                panel,
                alignment_directory,
                layout.panel_alignments_directory(
                    profile,
                    size,
                ),
                sequence_type,
            )

            manifest = {
                "schema_version": "0.1",
                "profile": profile,
                "requested_size": size,
                "actual_size": int(
                    len(panel)
                ),
                "sequence_type": sequence_type,
                "redundancy_penalty": (
                    redundancy_penalty
                ),
                "candidate_top_fraction": candidate_top_fraction,
                "candidate_minimum_pool_size": candidate_minimum_pool_size,
                "maximum_score_drop": maximum_score_drop,
                "ranking_uses_pca": False,
            }

            with layout.panel_manifest(
                profile,
                size,
            ).open(
                "w",
                encoding="utf-8",
            ) as handle:
                yaml.safe_dump(
                    manifest,
                    handle,
                    sort_keys=False,
                )
