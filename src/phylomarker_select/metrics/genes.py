"""Tabla de metricas por gen."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..constants import MISSING_CHARACTERS
from .distance import composition_variability, mean_pairwise_p_distance
from .sites import alignment_matrix, variable_sites_and_pis
from .trimming import calculate_trimming_metrics


def calculate_metrics(
    alignment_directory: Path,
    output_directory: Path,
    metadata: pd.DataFrame,
    sample_id_column: str,
    sequence_type: str,
    balance_level: str,
    raw_alignment_directory: Path | None = None,
    trimming_enabled: bool = False,
) -> pd.DataFrame:
    extension = ".faa" if sequence_type == "protein" else ".fna"

    paths = sorted(
        list(
            alignment_directory.glob(
                f"*.trimmed{extension}"
            )
        )
        + list(
            alignment_directory.glob(
                f"*.aln{extension}"
            )
        )
    )

    if not paths:
        raise ValueError(
            f"No alignments found in {alignment_directory}"
        )

    metadata_index = metadata.set_index(
        sample_id_column,
        drop=False,
    )

    total_taxa = len(metadata)
    metric_rows: list[dict] = []
    group_rows: list[dict] = []

    for path in paths:
        gene_id = path.name
        gene_id = gene_id.removesuffix(
            f".trimmed{extension}"
        )
        gene_id = gene_id.removesuffix(
            f".aln{extension}"
        )

        taxa, matrix = alignment_matrix(path)

        number_taxa = len(taxa)
        alignment_length = matrix.shape[1]

        missing_mask = np.isin(
            matrix,
            list(MISSING_CHARACTERS),
        )

        gap_mask = matrix == "-"

        taxon_occupancy = number_taxa / total_taxa
        sequence_completeness = float(
            1.0 - missing_mask.mean()
        )

        cell_occupancy = (
            taxon_occupancy
            * sequence_completeness
        )

        gap_fraction = float(gap_mask.mean())

        ambiguous_fraction = float(
            np.isin(
                matrix,
                ["?", "X", "N"],
            ).mean()
        )

        (
            variable_sites,
            informative_sites,
            entropy,
        ) = variable_sites_and_pis(matrix)

        pairwise_distance = mean_pairwise_p_distance(
            matrix
        )

        composition_variation = composition_variability(
            matrix
        )

        current_group_rows: list[dict] = []

        if balance_level in metadata_index.columns:
            for group_name, group in metadata.groupby(
                balance_level,
                dropna=False,
            ):
                members = set(
                    group[sample_id_column]
                )

                represented = sum(
                    taxon in members
                    for taxon in taxa
                )

                occupancy = (
                    represented / len(members)
                    if members
                    else float("nan")
                )

                taxon_to_index = {
                    taxon: index
                    for index, taxon in enumerate(taxa)
                }

                member_completeness: list[float] = []
                member_gap_fractions: list[float] = []

                for member in members:
                    if member not in taxon_to_index:
                        member_completeness.append(0.0)
                        member_gap_fractions.append(1.0)
                        continue

                    member_index = taxon_to_index[member]
                    member_missing = missing_mask[member_index]
                    member_gaps = gap_mask[member_index]

                    member_completeness.append(
                        float(1.0 - member_missing.mean())
                    )
                    member_gap_fractions.append(
                        float(member_gaps.mean())
                    )

                group_sequence_completeness = (
                    float(np.mean(member_completeness))
                    if member_completeness
                    else float("nan")
                )
                group_gap_fraction = (
                    float(np.mean(member_gap_fractions))
                    if member_gap_fractions
                    else float("nan")
                )

                group_row = {
                    "gene_id": gene_id,
                    "group_level": balance_level,
                    "group": group_name or "unknown",
                    "n_group_taxa": len(members),
                    "n_present": represented,
                    "group_occupancy": occupancy,
                    "group_sequence_completeness":
                        group_sequence_completeness,
                    "group_gap_fraction": group_gap_fraction,
                    "singleton_group": len(members) == 1,
                }

                group_rows.append(group_row)
                current_group_rows.append(group_row)

        replicated_occupancies = [
            row["group_occupancy"]
            for row in current_group_rows
            if not row["singleton_group"]
        ]

        all_group_occupancies = [
            row["group_occupancy"]
            for row in current_group_rows
        ]

        all_group_completeness = [
            row["group_sequence_completeness"]
            for row in current_group_rows
            if np.isfinite(row["group_sequence_completeness"])
        ]

        replicated_group_completeness = [
            row["group_sequence_completeness"]
            for row in current_group_rows
            if (
                not row["singleton_group"]
                and np.isfinite(
                    row["group_sequence_completeness"]
                )
            )
        ]

        all_group_gap_fractions = [
            row["group_gap_fraction"]
            for row in current_group_rows
            if np.isfinite(row["group_gap_fraction"])
        ]

        trimming_metrics: dict[str, object]

        if trimming_enabled:
            if raw_alignment_directory is None:
                raise ValueError(
                    "raw_alignment_directory is required when trimming "
                    "is enabled"
                )

            raw_alignment_path = (
                raw_alignment_directory
                / f"{gene_id}.aln{extension}"
            )

            if not raw_alignment_path.is_file():
                raise FileNotFoundError(
                    f"Raw alignment not found: {raw_alignment_path}"
                )

            trimming_metrics = calculate_trimming_metrics(
                raw_alignment_path=raw_alignment_path,
                trimmed_alignment_path=path,
            )
        else:
            trimming_metrics = {
                "raw_alignment_length": int(alignment_length),
                "trimmed_alignment_length": int(alignment_length),
                "retained_length_fraction": 1.0,
                "raw_variable_sites": int(variable_sites),
                "trimmed_variable_sites": int(variable_sites),
                "variable_sites_retained_fraction": (
                    1.0 if variable_sites > 0 else None
                ),
                "raw_pis": int(informative_sites),
                "trimmed_pis": int(informative_sites),
                "pis_retained_fraction": (
                    1.0 if informative_sites > 0 else None
                ),
                "raw_mean_entropy": float(entropy),
                "trimmed_mean_entropy": float(entropy),
                "trimming_class": "not_trimmed",
                "trimming_stability_score": 1.0,
            }

        metric_rows.append(
            {
                "gene_id": gene_id,
                **trimming_metrics,
                "alignment_file": str(path),
                "n_taxa": number_taxa,
                "n_taxa_total": total_taxa,
                "alignment_length": alignment_length,
                "taxon_occupancy": taxon_occupancy,
                "sequence_completeness": sequence_completeness,
                "cell_occupancy": cell_occupancy,
                "gap_fraction": gap_fraction,
                "ambiguous_fraction": ambiguous_fraction,
                "variable_sites": variable_sites,
                "parsimony_informative_sites": informative_sites,
                "variable_sites_per_length": (
                    variable_sites / alignment_length
                ),
                "pis_per_length": (
                    informative_sites / alignment_length
                ),
                "mean_entropy": entropy,
                "mean_pairwise_distance": pairwise_distance,
                "composition_variability": composition_variation,
                "mean_group_occupancy": (
                    float(
                        np.mean(
                            all_group_occupancies
                        )
                    )
                    if all_group_occupancies
                    else float("nan")
                ),
                "min_replicated_group_occupancy": (
                    float(
                        np.min(
                            replicated_occupancies
                        )
                    )
                    if replicated_occupancies
                    else float("nan")
                ),
                "n_groups_present": sum(
                    occupancy > 0
                    for occupancy in all_group_occupancies
                ),
                "mean_group_sequence_completeness": (
                    float(np.mean(all_group_completeness))
                    if all_group_completeness
                    else float("nan")
                ),
                "min_group_sequence_completeness": (
                    float(np.min(all_group_completeness))
                    if all_group_completeness
                    else float("nan")
                ),
                "min_replicated_group_sequence_completeness": (
                    float(np.min(replicated_group_completeness))
                    if replicated_group_completeness
                    else float("nan")
                ),
                "sd_group_sequence_completeness": (
                    float(np.std(all_group_completeness, ddof=0))
                    if all_group_completeness
                    else float("nan")
                ),
                "worst_group_gap_fraction": (
                    float(np.max(all_group_gap_fractions))
                    if all_group_gap_fractions
                    else float("nan")
                ),
            }
        )

    metrics = pd.DataFrame(metric_rows)
    metrics = metrics.drop_duplicates(
        "gene_id",
        keep="first",
    )

    metrics_directory = output_directory / "metrics"
    metrics_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    metrics.to_csv(
        metrics_directory / "per_gene_metrics.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(group_rows).to_csv(
        metrics_directory
        / f"per_gene_{balance_level}_occupancy.tsv",
        sep="\t",
        index=False,
    )

    return metrics
