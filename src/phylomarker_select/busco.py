"""Descubrimiento, validacion y extraccion de las corridas BUSCO."""
from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

import pandas as pd

from .constants import VALID_NUCLEOTIDE, VALID_PROTEIN
from .fasta import FastaRecord, first_fasta_record, write_fasta
from .layout import OutputLayout
from .utils import sha256_file

LOGGER = logging.getLogger("phylomarker-select")


def find_single_copy_directory(run_directory: Path) -> Path | None:
    candidates = [
        run_directory / "single_copy_busco_sequences",
        run_directory / "busco_sequences" / "single_copy_busco_sequences",
    ]

    candidates.extend(
        path
        for path in run_directory.rglob("single_copy_busco_sequences")
        if path.is_dir()
    )

    seen: set[Path] = set()

    for candidate in candidates:
        if not candidate.is_dir():
            continue

        resolved = candidate.resolve()

        if resolved in seen:
            continue

        seen.add(resolved)
        return candidate

    return None


def discover_busco_runs(busco_root: Path) -> pd.DataFrame:
    if not busco_root.is_dir():
        raise NotADirectoryError(
            f"BUSCO root directory not found: {busco_root}"
        )

    rows: list[dict] = []

    for run_directory in sorted(busco_root.rglob("run_*")):
        if not run_directory.is_dir():
            continue

        sample_id = run_directory.name.removeprefix("run_")
        sequence_directory = find_single_copy_directory(run_directory)

        protein_files = (
            sorted(sequence_directory.glob("*.faa"))
            if sequence_directory
            else []
        )

        nucleotide_files = (
            sorted(sequence_directory.glob("*.fna"))
            if sequence_directory
            else []
        )

        summaries = sorted(run_directory.glob("short_summary*.txt"))
        full_tables = sorted(run_directory.glob("full_table*.tsv"))

        rows.append(
            {
                "sample_ID": sample_id,
                "run_directory": str(run_directory),
                "single_copy_directory": (
                    str(sequence_directory)
                    if sequence_directory is not None
                    else ""
                ),
                "n_faa": len(protein_files),
                "n_fna": len(nucleotide_files),
                "summary_file": (
                    str(summaries[0]) if summaries else ""
                ),
                "full_table_file": (
                    str(full_tables[0]) if full_tables else ""
                ),
            }
        )

    if not rows:
        raise ValueError(
            f"No BUSCO run_* directories were found under {busco_root}"
        )

    runs = pd.DataFrame(rows)

    duplicates = runs.loc[
        runs["sample_ID"].duplicated(keep=False),
        "sample_ID",
    ].tolist()

    if duplicates:
        raise ValueError(
            "Duplicated BUSCO run identifiers: "
            + ", ".join(sorted(set(duplicates)))
        )

    return runs


def validate_runs(
    runs: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_id_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    run_ids = set(runs["sample_ID"])
    metadata_ids = set(metadata[sample_id_column])

    warnings: list[dict] = []

    for sample_id in sorted(run_ids - metadata_ids):
        warnings.append(
            {
                "severity": "warning",
                "sample_ID": sample_id,
                "code": "RUN_WITHOUT_METADATA",
                "message": (
                    "BUSCO run exists but no matching metadata row was found."
                ),
            }
        )

    for sample_id in sorted(metadata_ids - run_ids):
        warnings.append(
            {
                "severity": "warning",
                "sample_ID": sample_id,
                "code": "METADATA_WITHOUT_RUN",
                "message": (
                    "Metadata row exists but no matching BUSCO run was found."
                ),
            }
        )

    for row in runs.itertuples(index=False):
        if not row.single_copy_directory:
            warnings.append(
                {
                    "severity": "error",
                    "sample_ID": row.sample_ID,
                    "code": "MISSING_SINGLE_COPY_DIRECTORY",
                    "message": (
                        "single_copy_busco_sequences directory was not found."
                    ),
                }
            )

        if row.n_faa == 0:
            warnings.append(
                {
                    "severity": "warning",
                    "sample_ID": row.sample_ID,
                    "code": "NO_PROTEIN_FASTA",
                    "message": "No .faa files were detected.",
                }
            )

        if row.n_fna == 0:
            warnings.append(
                {
                    "severity": "warning",
                    "sample_ID": row.sample_ID,
                    "code": "NO_NUCLEOTIDE_FASTA",
                    "message": "No .fna files were detected.",
                }
            )

    common_ids = sorted(run_ids & metadata_ids)

    validated = runs[
        runs["sample_ID"].isin(common_ids)
    ].copy()

    validated = validated.merge(
        metadata,
        left_on="sample_ID",
        right_on=sample_id_column,
        how="left",
        suffixes=("", "_metadata"),
    )

    warning_frame = pd.DataFrame(
        warnings,
        columns=[
            "severity",
            "sample_ID",
            "code",
            "message",
        ],
    )

    return validated, warning_frame


def validate_sequence(
    sequence: str,
    sequence_type: str,
) -> tuple[int, int]:
    valid_characters = (
        VALID_PROTEIN
        if sequence_type == "protein"
        else VALID_NUCLEOTIDE
    )

    invalid_characters = sum(
        character not in valid_characters
        for character in sequence
    )

    internal_stops = (
        sequence[:-1].count("*")
        if sequence_type == "protein" and sequence
        else 0
    )

    return invalid_characters, internal_stops


def extract_markers(
    validated_runs: pd.DataFrame,
    layout: OutputLayout,
) -> pd.DataFrame:
    sequence_type = layout.sequence_type
    extension = layout.extension

    sequence_output = layout.per_gene_sequences_directory

    sequence_output.mkdir(parents=True, exist_ok=True)

    records_by_gene: dict[str, list[FastaRecord]] = defaultdict(list)
    provenance_rows: list[dict] = []

    for row in validated_runs.itertuples(index=False):
        sample_id = row.sample_ID
        source_directory = Path(row.single_copy_directory)

        for source_path in sorted(
            source_directory.glob(f"*{extension}")
        ):
            gene_id = source_path.stem
            record = first_fasta_record(source_path)

            invalid_count, internal_stops = validate_sequence(
                record.sequence,
                sequence_type,
            )

            selected = invalid_count == 0 and internal_stops == 0

            provenance_rows.append(
                {
                    "gene_id": gene_id,
                    "sample_ID": sample_id,
                    "sequence_type": sequence_type,
                    "source_file": str(source_path),
                    "original_header": record.identifier,
                    "sequence_length": len(record.sequence),
                    "invalid_character_count": invalid_count,
                    "internal_stop_count": internal_stops,
                    "selected": selected,
                    "selection_reason": (
                        "single_copy_busco_sequence"
                        if selected
                        else "invalid_sequence"
                    ),
                    "sha256": sha256_file(source_path),
                }
            )

            if selected:
                records_by_gene[gene_id].append(
                    FastaRecord(
                        identifier=sample_id,
                        sequence=record.sequence,
                    )
                )

    for gene_id, records in records_by_gene.items():
        write_fasta(
            records,
            sequence_output / f"{gene_id}{extension}",
        )

    provenance = pd.DataFrame(provenance_rows)

    layout.validation_directory.mkdir(parents=True, exist_ok=True)

    provenance.to_csv(
        layout.sequence_provenance_table,
        sep="\t",
        index=False,
    )

    return provenance
