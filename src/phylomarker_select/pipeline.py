"""Orquestacion de las etapas del flujo."""
from __future__ import annotations

import logging
from pathlib import Path

import yaml

from .align import align_markers, trim_alignments
from .busco import discover_busco_runs, extract_markers, validate_runs
from .layout import OutputLayout
from .metadata import load_metadata
from .metrics.genes import calculate_metrics
from .panels import create_panels
from .pca import run_exploratory_pca
from .provenance import write_provenance
from .report import create_html_report
from .scoring import add_biological_scores

LOGGER = logging.getLogger("phylomarker-select")


def load_yaml(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    with path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    if not isinstance(config, dict):
        raise ValueError("The YAML configuration must contain a mapping.")

    return config


def run_pipeline(
    configuration_path: Path,
) -> None:
    config = load_yaml(
        configuration_path
    )

    project_config = config.get(
        "project",
        {},
    )

    input_config = config.get(
        "inputs",
        {},
    )

    alignment_config = config.get(
        "alignment",
        {},
    )

    trimming_config = config.get(
        "trimming",
        {},
    )

    taxonomy_config = config.get(
        "taxonomy",
        {},
    )

    busco_directory = Path(
        input_config["busco_directory"]
    )

    metadata_path = Path(
        input_config["metadata"]
    )

    output_directory = Path(
        input_config.get(
            "output",
            "results",
        )
    )

    sample_id_column = input_config.get(
        "sample_id_column",
        "sample_ID",
    )

    sequence_type = project_config.get(
        "sequence_type",
        "protein",
    )

    if sequence_type not in {
        "protein",
        "nucleotide",
    }:
        raise ValueError(
            "sequence_type must be 'protein' or 'nucleotide'"
        )

    balance_level = taxonomy_config.get(
        "balance_level",
        "order",
    )

    layout = OutputLayout(
        root=output_directory,
        sequence_type=sequence_type,
        balance_level=balance_level,
    )

    trimming_enabled = bool(
        trimming_config.get(
            "enabled",
            True,
        )
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    write_provenance(
        layout,
        config,
    )

    metadata = load_metadata(
        metadata_path,
        sample_id_column,
    )

    runs = discover_busco_runs(
        busco_directory
    )

    validated_runs, warnings = validate_runs(
        runs,
        metadata,
        sample_id_column,
    )

    layout.validation_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    runs.to_csv(
        layout.discovered_runs_table,
        sep="\t",
        index=False,
    )

    validated_runs.to_csv(
        layout.validated_runs_table,
        sep="\t",
        index=False,
    )

    warnings.to_csv(
        layout.warnings_table,
        sep="\t",
        index=False,
    )

    if (
        not warnings.empty
        and (
            warnings["severity"] == "error"
        ).any()
    ):
        raise RuntimeError(
            "Input validation produced errors. "
            "Inspect validation/warnings.tsv."
        )

    extract_markers(
        validated_runs,
        layout,
    )

    align_markers(
        layout=layout,
        mafft_executable=alignment_config.get(
            "mafft_executable",
            "mafft",
        ),
        threads=int(
            alignment_config.get(
                "threads_per_gene",
                2,
            )
        ),
        strategy=alignment_config.get(
            "strategy",
            "auto",
        ),
    )

    if trimming_enabled:
        trim_alignments(
            layout=layout,
            trimal_executable=trimming_config.get(
                "trimal_executable",
                "trimal",
            ),
            mode=trimming_config.get(
                "mode",
                "automated1",
            ),
        )

    metrics = calculate_metrics(
        layout=layout,
        metadata=metadata,
        sample_id_column=sample_id_column,
        trimming_enabled=trimming_enabled,
    )

    scored = add_biological_scores(
        metrics,
        config,
    )

    layout.rankings_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    scored.to_csv(
        layout.all_gene_scores_table,
        sep="\t",
        index=False,
    )

    if config.get(
        "pca",
        {},
    ).get(
        "enabled",
        True,
    ):
        run_exploratory_pca(
            scored,
            layout,
        )

    create_panels(
        scored=scored,
        layout=layout,
        config=config,
        trimming_enabled=trimming_enabled,
    )

    create_html_report(
        layout=layout,
        runs=validated_runs,
        warnings=warnings,
        metrics=metrics,
        scored=scored,
        config=config,
    )

    LOGGER.info(
        "Analysis complete: %s",
        output_directory,
    )

    LOGGER.info(
        "HTML report: %s",
        layout.report_index,
    )
