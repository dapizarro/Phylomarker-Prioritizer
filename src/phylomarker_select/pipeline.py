"""Orquestacion de las etapas del flujo."""
from __future__ import annotations

import logging
from pathlib import Path

from .align import align_markers, trim_alignments
from .busco import discover_busco_runs, extract_markers, validate_runs
from .config import SelectConfig
from .layout import OutputLayout
from .metadata import load_metadata
from .metrics.genes import calculate_metrics
from .panels import create_panels
from .pca import run_exploratory_pca
from .provenance import write_provenance
from .report import create_html_report
from .scoring import add_biological_scores

LOGGER = logging.getLogger("phylomarker-select")


def run_pipeline(
    configuration_path: Path,
) -> None:
    config = SelectConfig.load(
        configuration_path
    )

    busco_directory = config.inputs.busco_directory
    metadata_path = config.inputs.metadata
    output_directory = config.inputs.output
    sample_id_column = config.inputs.sample_id_column

    layout = OutputLayout(
        root=output_directory,
        sequence_type=config.project.sequence_type,
        balance_level=config.taxonomy.balance_level,
    )

    trimming_enabled = config.trimming.enabled

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    write_provenance(
        layout,
        config.raw,
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
        mafft_executable=config.alignment.mafft_executable,
        threads=config.alignment.threads_per_gene,
        strategy=config.alignment.strategy,
    )

    if trimming_enabled:
        trim_alignments(
            layout=layout,
            trimal_executable=config.trimming.trimal_executable,
            mode=config.trimming.mode,
        )

    metrics = calculate_metrics(
        layout=layout,
        metadata=metadata,
        sample_id_column=sample_id_column,
        trimming_enabled=trimming_enabled,
    )

    scored = add_biological_scores(
        metrics,
        config.eligibility,
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

    if config.pca.enabled:
        run_exploratory_pca(
            scored,
            layout,
        )

    create_panels(
        scored=scored,
        layout=layout,
        panel_config=config.panels,
        random_seed=config.project.random_seed,
        trimming_enabled=trimming_enabled,
    )

    create_html_report(
        layout=layout,
        runs=validated_runs,
        warnings=warnings,
        metrics=metrics,
        scored=scored,
        config=config.raw,
    )

    LOGGER.info(
        "Analysis complete: %s",
        output_directory,
    )

    LOGGER.info(
        "HTML report: %s",
        layout.report_index,
    )
