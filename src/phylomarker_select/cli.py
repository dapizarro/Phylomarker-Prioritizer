from __future__ import annotations

import argparse
import hashlib
import html
import logging
import math
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
import pandas as pd
import yaml
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import RobustScaler


LOGGER = logging.getLogger("phylomarker-select")

MISSING_CHARACTERS = {"-", "?", "X", "N"}
VALID_PROTEIN = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*-?")
VALID_NUCLEOTIDE = set("ACGTUNRYKMSWBDHVX*-?")


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


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    sequence: str


def configure_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def load_yaml(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    with path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    if not isinstance(config, dict):
        raise ValueError("The YAML configuration must contain a mapping.")

    return config


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def read_fasta(path: Path) -> Iterator[FastaRecord]:
    identifier: str | None = None
    sequence_parts: list[str] = []

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if identifier is not None:
                    yield FastaRecord(
                        identifier=identifier,
                        sequence="".join(sequence_parts).upper(),
                    )

                identifier = line[1:].split()[0]
                sequence_parts = []
            else:
                if identifier is None:
                    raise ValueError(
                        f"Sequence found before FASTA header in {path}"
                    )

                sequence_parts.append(line.replace(" ", ""))

    if identifier is not None:
        yield FastaRecord(
            identifier=identifier,
            sequence="".join(sequence_parts).upper(),
        )


def first_fasta_record(path: Path) -> FastaRecord:
    records = list(read_fasta(path))

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    if len(records) > 1:
        LOGGER.warning(
            "Multiple FASTA records found in %s; using the first one.",
            path,
        )

    return records[0]


def write_fasta(
    records: Iterable[FastaRecord],
    path: Path,
    width: int = 80,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(path.suffix + ".tmp")

    with temporary_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.identifier}\n")

            for start in range(0, len(record.sequence), width):
                handle.write(record.sequence[start : start + width] + "\n")

    temporary_path.replace(path)


def load_metadata(
    path: Path,
    sample_id_column: str,
) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Metadata file not found: {path}")

    metadata = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    if sample_id_column not in metadata.columns:
        raise ValueError(
            f"Required metadata column '{sample_id_column}' was not found. "
            f"Available columns: {', '.join(metadata.columns)}"
        )

    metadata[sample_id_column] = (
        metadata[sample_id_column].astype(str).str.strip()
    )

    if (metadata[sample_id_column] == "").any():
        raise ValueError("Empty sample identifiers were found.")

    duplicates = metadata.loc[
        metadata[sample_id_column].duplicated(keep=False),
        sample_id_column,
    ].tolist()

    if duplicates:
        raise ValueError(
            "Duplicated sample identifiers: "
            + ", ".join(sorted(set(duplicates)))
        )

    return metadata


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
    output_directory: Path,
    sequence_type: str,
) -> pd.DataFrame:
    extension = ".faa" if sequence_type == "protein" else ".fna"

    sequence_output = (
        output_directory
        / "sequences"
        / sequence_type
        / "per_gene"
    )

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

    validation_directory = output_directory / "validation"
    validation_directory.mkdir(parents=True, exist_ok=True)

    provenance.to_csv(
        validation_directory
        / f"{sequence_type}_sequence_provenance.tsv",
        sep="\t",
        index=False,
    )

    return provenance


def require_executable(name: str) -> str:
    executable = shutil.which(name)

    if executable is None:
        raise RuntimeError(
            f"Required executable '{name}' was not found in PATH."
        )

    return executable


def run_to_file(
    command: list[str],
    output_path: Path,
) -> None:
    LOGGER.debug("Running command: %s", " ".join(command))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(
        output_path.suffix + ".tmp"
    )

    with temporary_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            command,
            check=True,
            stdout=handle,
        )

    temporary_path.replace(output_path)


def align_markers(
    input_directory: Path,
    output_directory: Path,
    sequence_type: str,
    mafft_executable: str,
    threads: int,
    strategy: str,
) -> None:
    extension = ".faa" if sequence_type == "protein" else ".fna"

    aligned_directory = (
        output_directory
        / "alignments"
        / sequence_type
        / "untrimmed"
    )

    aligned_directory.mkdir(parents=True, exist_ok=True)

    mafft = require_executable(mafft_executable)

    for input_path in sorted(input_directory.glob(f"*{extension}")):
        output_path = (
            aligned_directory
            / f"{input_path.stem}.aln{extension}"
        )

        if output_path.is_file() and output_path.stat().st_size > 0:
            LOGGER.info(
                "Existing alignment retained: %s",
                input_path.stem,
            )
            continue

        command = [
            mafft,
            "--thread",
            str(threads),
        ]

        if strategy == "auto":
            command.append("--auto")
        elif strategy == "retree2":
            command.extend(
                [
                    "--retree",
                    "2",
                    "--maxiterate",
                    "0",
                ]
            )
        else:
            raise ValueError(
                f"Unsupported MAFFT strategy: {strategy}"
            )

        command.append(str(input_path))
        run_to_file(command, output_path)


def trim_alignments(
    input_directory: Path,
    output_directory: Path,
    sequence_type: str,
    trimal_executable: str,
    mode: str,
) -> None:
    extension = ".faa" if sequence_type == "protein" else ".fna"

    trimmed_directory = (
        output_directory
        / "alignments"
        / sequence_type
        / "trimmed"
    )

    trimmed_directory.mkdir(parents=True, exist_ok=True)

    trimal = require_executable(trimal_executable)

    for input_path in sorted(
        input_directory.glob(f"*.aln{extension}")
    ):
        gene_id = input_path.name.removesuffix(
            f".aln{extension}"
        )

        output_path = (
            trimmed_directory
            / f"{gene_id}.trimmed{extension}"
        )

        if output_path.is_file() and output_path.stat().st_size > 0:
            continue

        command = [
            trimal,
            "-in",
            str(input_path),
            "-out",
            str(output_path),
        ]

        if mode == "automated1":
            command.append("-automated1")
        elif mode == "gappyout":
            command.append("-gappyout")
        else:
            raise ValueError(
                f"Unsupported trimAl mode: {mode}"
            )

        subprocess.run(command, check=True)


def alignment_matrix(
    path: Path,
) -> tuple[list[str], np.ndarray]:
    records = list(read_fasta(path))

    if not records:
        raise ValueError(f"Empty alignment: {path}")

    lengths = {
        len(record.sequence)
        for record in records
    }

    if len(lengths) != 1:
        raise ValueError(
            f"Sequences with unequal lengths found in {path}"
        )

    taxa = [
        record.identifier
        for record in records
    ]

    matrix = np.asarray(
        [
            list(record.sequence)
            for record in records
        ],
        dtype="U1",
    )

    return taxa, matrix


def site_entropy(column: np.ndarray) -> float:
    usable = [
        character
        for character in column.tolist()
        if character not in MISSING_CHARACTERS
    ]

    if not usable:
        return 0.0

    counts = Counter(usable)
    total = len(usable)

    return float(
        -sum(
            (count / total)
            * math.log(count / total, 2)
            for count in counts.values()
        )
    )


def variable_sites_and_pis(
    matrix: np.ndarray,
) -> tuple[int, int, float]:
    variable_sites = 0
    informative_sites = 0
    entropies: list[float] = []

    for column in matrix.T:
        usable = [
            character
            for character in column.tolist()
            if character not in MISSING_CHARACTERS
        ]

        counts = Counter(usable)

        if len(counts) >= 2:
            variable_sites += 1

        if sum(
            count >= 2
            for count in counts.values()
        ) >= 2:
            informative_sites += 1

        entropies.append(site_entropy(column))

    mean_entropy = (
        float(np.mean(entropies))
        if entropies
        else 0.0
    )

    return (
        variable_sites,
        informative_sites,
        mean_entropy,
    )



def safe_fraction(
    numerator: int | float,
    denominator: int | float,
) -> float | None:
    """Return numerator / denominator, or None for a zero denominator."""
    if denominator <= 0:
        return None
    return float(numerator / denominator)


def classify_trimming(
    retained_length_fraction: float,
    pis_retained_fraction: float | None,
) -> str:
    """Classify gene sensitivity to alignment trimming."""
    if retained_length_fraction < 0.25:
        return "extreme"

    if (
        pis_retained_fraction is not None
        and pis_retained_fraction < 0.50
    ):
        return "extreme"

    if (
        retained_length_fraction >= 0.75
        and (
            pis_retained_fraction is None
            or pis_retained_fraction >= 0.80
        )
    ):
        return "stable"

    if (
        pis_retained_fraction is None
        or pis_retained_fraction >= 0.80
    ):
        return "signal_preserved"

    return "signal_sensitive"


def calculate_trimming_stability_score(
    retained_length_fraction: float,
    pis_retained_fraction: float | None,
) -> float:
    """Return a bounded score in which preservation of PIS has more weight."""
    length_score = min(
        1.0,
        max(0.0, retained_length_fraction / 0.75),
    )

    if pis_retained_fraction is None:
        score = length_score
    else:
        information_score = min(
            1.0,
            max(0.0, pis_retained_fraction / 0.80),
        )
        score = 0.40 * length_score + 0.60 * information_score

    if retained_length_fraction < 0.25:
        score *= 0.25

    if (
        pis_retained_fraction is not None
        and pis_retained_fraction < 0.50
    ):
        score *= 0.50

    return float(min(1.0, max(0.0, score)))


def calculate_trimming_metrics(
    raw_alignment_path: Path,
    trimmed_alignment_path: Path,
) -> dict[str, object]:
    """Compare raw and trimmed alignments generated in the same run."""
    raw_taxa, raw_matrix = alignment_matrix(raw_alignment_path)
    trimmed_taxa, trimmed_matrix = alignment_matrix(
        trimmed_alignment_path
    )

    if set(raw_taxa) != set(trimmed_taxa):
        raise ValueError(
            "Raw and trimmed alignments contain different taxa for "
            f"{raw_alignment_path.name}"
        )

    raw_variable, raw_pis, raw_entropy = variable_sites_and_pis(
        raw_matrix
    )
    trimmed_variable, trimmed_pis, trimmed_entropy = (
        variable_sites_and_pis(trimmed_matrix)
    )

    raw_length = int(raw_matrix.shape[1])
    trimmed_length = int(trimmed_matrix.shape[1])
    retained_length_fraction = safe_fraction(
        trimmed_length,
        raw_length,
    )

    if retained_length_fraction is None:
        raise ValueError(
            f"Raw alignment has zero length: {raw_alignment_path}"
        )

    if trimmed_pis > raw_pis:
        raise ValueError(
            "Trimmed alignment contains more PIS than its raw "
            f"alignment: {raw_alignment_path.name}; "
            f"raw={raw_pis}, trimmed={trimmed_pis}"
        )

    variable_retained = safe_fraction(
        trimmed_variable,
        raw_variable,
    )
    pis_retained = safe_fraction(
        trimmed_pis,
        raw_pis,
    )
    trimming_class = classify_trimming(
        retained_length_fraction,
        pis_retained,
    )

    return {
        "raw_alignment_length": raw_length,
        "trimmed_alignment_length": trimmed_length,
        "retained_length_fraction": retained_length_fraction,
        "raw_variable_sites": int(raw_variable),
        "trimmed_variable_sites": int(trimmed_variable),
        "variable_sites_retained_fraction": variable_retained,
        "raw_pis": int(raw_pis),
        "trimmed_pis": int(trimmed_pis),
        "pis_retained_fraction": pis_retained,
        "raw_mean_entropy": float(raw_entropy),
        "trimmed_mean_entropy": float(trimmed_entropy),
        "trimming_class": trimming_class,
        "trimming_stability_score":
            calculate_trimming_stability_score(
                retained_length_fraction,
                pis_retained,
            ),
    }

def mean_pairwise_p_distance(
    matrix: np.ndarray,
) -> float:
    distances: list[float] = []

    for left_index in range(matrix.shape[0]):
        for right_index in range(
            left_index + 1,
            matrix.shape[0],
        ):
            left = matrix[left_index]
            right = matrix[right_index]

            valid = np.asarray(
                [
                    left_character not in MISSING_CHARACTERS
                    and right_character not in MISSING_CHARACTERS
                    for left_character, right_character in zip(
                        left,
                        right,
                    )
                ],
                dtype=bool,
            )

            if not valid.any():
                continue

            distances.append(
                float(
                    np.mean(
                        left[valid] != right[valid]
                    )
                )
            )

    return (
        float(np.mean(distances))
        if distances
        else 0.0
    )


def composition_variability(
    matrix: np.ndarray,
) -> float:
    alphabet = sorted(
        set(matrix.flatten().tolist())
        - MISSING_CHARACTERS
        - {"*"}
    )

    if not alphabet:
        return 0.0

    frequencies: list[list[float]] = []

    for row in matrix:
        usable = [
            character
            for character in row.tolist()
            if character not in MISSING_CHARACTERS
            and character != "*"
        ]

        if not usable:
            continue

        counts = Counter(usable)

        frequencies.append(
            [
                counts[character] / len(usable)
                for character in alphabet
            ]
        )

    if len(frequencies) < 2:
        return 0.0

    frequency_array = np.asarray(
        frequencies,
        dtype=float,
    )

    return float(
        np.mean(
            np.std(
                frequency_array,
                axis=0,
            )
        )
    )


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


def run_exploratory_pca(
    scored: pd.DataFrame,
    output_directory: Path,
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

    pca_directory = output_directory / "pca"
    pca_directory.mkdir(
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
        pca_directory / "pca_scores.tsv",
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
        pca_directory / "pca_loadings.tsv",
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
        pca_directory
        / "pca_explained_variance.tsv",
        sep="\t",
        index=False,
    )


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
    alignment_directory: Path,
    output_directory: Path,
    config: dict,
    sequence_type: str,
) -> None:
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

    panels_directory = (
        output_directory / "panels"
    )

    rankings_directory = (
        output_directory / "rankings"
    )

    panels_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    rankings_directory.mkdir(
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

                panel_directory = (
                    panels_directory
                    / profile
                    / f"n{size}"
                )

                panel_directory.mkdir(
                    parents=True,
                    exist_ok=True,
                )

                panel.to_csv(
                    panel_directory
                    / "panel_genes.tsv",
                    sep="\t",
                    index=False,
                )

                summarize_panel(panel).to_csv(
                    panel_directory / "panel_summary.tsv",
                    sep="\t",
                    index=False,
                )

                panel["gene_id"].to_csv(
                    panel_directory
                    / "genes.txt",
                    index=False,
                    header=False,
                )

                export_panel_alignments(
                    panel,
                    alignment_directory,
                    panel_directory
                    / "alignments",
                    sequence_type,
                )

            continue

        ranked = calculate_profile_ranking(
            scored,
            profile,
        )

        ranked.to_csv(
            rankings_directory
            / f"{profile}.tsv",
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

            panel_directory = (
                panels_directory
                / profile
                / f"n{size}"
            )

            panel_directory.mkdir(
                parents=True,
                exist_ok=True,
            )

            panel.to_csv(
                panel_directory
                / "panel_genes.tsv",
                sep="\t",
                index=False,
            )

            summarize_panel(panel).to_csv(
                panel_directory / "panel_summary.tsv",
                sep="\t",
                index=False,
            )

            trace.to_csv(
                panel_directory
                / "selection_trace.tsv",
                sep="\t",
                index=False,
            )

            panel["gene_id"].to_csv(
                panel_directory
                / "genes.txt",
                index=False,
                header=False,
            )

            export_panel_alignments(
                panel,
                alignment_directory,
                panel_directory
                / "alignments",
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

            with (
                panel_directory / "panel.yaml"
            ).open(
                "w",
                encoding="utf-8",
            ) as handle:
                yaml.safe_dump(
                    manifest,
                    handle,
                    sort_keys=False,
                )


def create_html_report(
    output_directory: Path,
    runs: pd.DataFrame,
    warnings: pd.DataFrame,
    metrics: pd.DataFrame,
    scored: pd.DataFrame,
    config: dict,
) -> None:
    report_directory = (
        output_directory / "report"
    )

    report_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    eligible_count = int(
        scored["eligible"].sum()
    )

    excluded_count = int(
        (~scored["eligible"]).sum()
    )

    top_complete = scored.sort_values(
        "cell_occupancy",
        ascending=False,
    ).head(15)

    table_rows = "\n".join(
        (
            "<tr>"
            f"<td>{html.escape(str(row.gene_id))}</td>"
            f"<td>{row.cell_occupancy:.3f}</td>"
            f"<td>{row.parsimony_informative_sites}</td>"
            f"<td>{row.gap_fraction:.3f}</td>"
            f"<td>{row.composition_variability:.4f}</td>"
            "</tr>"
        )
        for row in top_complete.itertuples(
            index=False
        )
    )

    warning_count = len(warnings)

    report = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>PhyloMarker Select report</title>
<style>
body {{
    max-width: 1100px;
    margin: 40px auto;
    padding: 0 24px;
    font-family: system-ui, sans-serif;
    color: #17202a;
    line-height: 1.5;
}}
h1, h2 {{
    color: #17324d;
}}
.cards {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 16px;
}}
.card {{
    border: 1px solid #d9e2ec;
    border-radius: 10px;
    padding: 16px;
    background: #f8fafc;
}}
.value {{
    font-size: 2rem;
    font-weight: 700;
}}
table {{
    border-collapse: collapse;
    width: 100%;
}}
th, td {{
    border-bottom: 1px solid #d9e2ec;
    padding: 8px;
    text-align: left;
}}
.notice {{
    border-left: 4px solid #d97706;
    padding: 12px;
    background: #fff7ed;
}}
pre {{
    white-space: pre-wrap;
    background: #f1f5f9;
    padding: 16px;
}}
</style>
</head>
<body>
<h1>PhyloMarker Select</h1>

<p>
Evolution-aware marker characterization and panel optimization.
PCA is exploratory and is not used as a biological quality ranking.
</p>

<div class="cards">
<div class="card">
<div class="value">{len(runs)}</div>
Validated BUSCO runs
</div>

<div class="card">
<div class="value">{len(metrics)}</div>
Aligned markers
</div>

<div class="card">
<div class="value">{eligible_count}</div>
Eligible markers
</div>

<div class="card">
<div class="value">{excluded_count}</div>
Excluded markers
</div>

<div class="card">
<div class="value">{warning_count}</div>
Validation warnings
</div>
</div>

<h2>Scientific warning</h2>

<div class="notice">
Marker suitability depends on taxonomic sampling and evolutionary objective.
High occupancy does not guarantee phylogenetic signal, and many informative
sites do not guarantee an unbiased or correct gene tree. Final panel quality
must be evaluated using independent phylogenetic analyses.
</div>

<h2>Most complete markers</h2>

<table>
<thead>
<tr>
<th>Gene</th>
<th>Cell occupancy</th>
<th>PIS</th>
<th>Gap fraction</th>
<th>Composition variability</th>
</tr>
</thead>
<tbody>
{table_rows}
</tbody>
</table>

<h2>Interpretation</h2>

<ul>
<li>PCA describes multivariate structure; it does not rank biological quality.</li>
<li>Singleton taxonomic groups are reported separately from replicated groups.</li>
<li>Eligibility thresholds are hard constraints and cannot be compensated by a high score.</li>
<li>Panel optimization penalizes redundant genes.</li>
<li>Random and single-criterion panels are generated as controls.</li>
</ul>

<h2>Resolved configuration</h2>

<pre>{html.escape(yaml.safe_dump(config, sort_keys=False))}</pre>
</body>
</html>
"""

    (
        report_directory / "index.html"
    ).write_text(
        report,
        encoding="utf-8",
    )


def write_provenance(
    output_directory: Path,
    config: dict,
) -> None:
    provenance_directory = (
        output_directory / "provenance"
    )

    provenance_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    versions = [
        {
            "software": "python",
            "version": sys.version.replace(
                "\n",
                " ",
            ),
        },
        {
            "software": "numpy",
            "version": np.__version__,
        },
        {
            "software": "pandas",
            "version": pd.__version__,
        },
    ]

    pd.DataFrame(versions).to_csv(
        provenance_directory
        / "software_versions.tsv",
        sep="\t",
        index=False,
    )

    with (
        provenance_directory
        / "resolved_config.yaml"
    ).open(
        "w",
        encoding="utf-8",
    ) as handle:
        yaml.safe_dump(
            config,
            handle,
            sort_keys=False,
        )


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

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    write_provenance(
        output_directory,
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

    validation_directory = (
        output_directory / "validation"
    )

    validation_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    runs.to_csv(
        validation_directory
        / "discovered_runs.tsv",
        sep="\t",
        index=False,
    )

    validated_runs.to_csv(
        validation_directory
        / "validated_runs.tsv",
        sep="\t",
        index=False,
    )

    warnings.to_csv(
        validation_directory
        / "warnings.tsv",
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
        output_directory,
        sequence_type,
    )

    extension = (
        ".faa"
        if sequence_type == "protein"
        else ".fna"
    )

    per_gene_directory = (
        output_directory
        / "sequences"
        / sequence_type
        / "per_gene"
    )

    align_markers(
        input_directory=per_gene_directory,
        output_directory=output_directory,
        sequence_type=sequence_type,
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

    untrimmed_directory = (
        output_directory
        / "alignments"
        / sequence_type
        / "untrimmed"
    )

    analysis_alignment_directory = (
        untrimmed_directory
    )

    if trimming_config.get(
        "enabled",
        True,
    ):
        trim_alignments(
            input_directory=untrimmed_directory,
            output_directory=output_directory,
            sequence_type=sequence_type,
            trimal_executable=trimming_config.get(
                "trimal_executable",
                "trimal",
            ),
            mode=trimming_config.get(
                "mode",
                "automated1",
            ),
        )

        analysis_alignment_directory = (
            output_directory
            / "alignments"
            / sequence_type
            / "trimmed"
        )

    metrics = calculate_metrics(
        alignment_directory=analysis_alignment_directory,
        output_directory=output_directory,
        metadata=metadata,
        sample_id_column=sample_id_column,
        sequence_type=sequence_type,
        balance_level=balance_level,
        raw_alignment_directory=untrimmed_directory,
        trimming_enabled=bool(
            trimming_config.get("enabled", True)
        ),
    )

    scored = add_biological_scores(
        metrics,
        config,
    )

    ranking_directory = (
        output_directory / "rankings"
    )

    ranking_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    scored.to_csv(
        ranking_directory
        / "all_gene_scores.tsv",
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
            output_directory,
        )

    create_panels(
        scored=scored,
        alignment_directory=analysis_alignment_directory,
        output_directory=output_directory,
        config=config,
        sequence_type=sequence_type,
    )

    create_html_report(
        output_directory=output_directory,
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
        output_directory
        / "report"
        / "index.html",
    )


def validate_command(
    args: argparse.Namespace,
) -> None:
    metadata = load_metadata(
        Path(args.metadata),
        args.sample_id_column,
    )

    runs = discover_busco_runs(
        Path(args.busco_directory)
    )

    validated, warnings = validate_runs(
        runs,
        metadata,
        args.sample_id_column,
    )

    output_directory = Path(
        args.output
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    runs.to_csv(
        output_directory
        / "discovered_runs.tsv",
        sep="\t",
        index=False,
    )

    validated.to_csv(
        output_directory
        / "validated_runs.tsv",
        sep="\t",
        index=False,
    )

    warnings.to_csv(
        output_directory
        / "warnings.tsv",
        sep="\t",
        index=False,
    )

    print(
        f"Detected BUSCO runs: {len(runs)}"
    )

    print(
        f"Validated BUSCO runs: {len(validated)}"
    )

    print(
        f"Warnings or errors: {len(warnings)}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="phylomarker-select",
        description=(
            "Evolution-aware selection of "
            "phylogenomic marker panels from BUSCO outputs."
        ),
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable detailed logging.",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    run_parser = subparsers.add_parser(
        "run",
        help="Run the complete workflow.",
    )

    run_parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="YAML configuration file.",
    )

    validate_parser = subparsers.add_parser(
        "validate",
        help="Validate BUSCO runs and metadata.",
    )

    validate_parser.add_argument(
        "--busco-directory",
        required=True,
    )

    validate_parser.add_argument(
        "--metadata",
        required=True,
    )

    validate_parser.add_argument(
        "--sample-id-column",
        default="sample_ID",
    )

    validate_parser.add_argument(
        "--output",
        default="validation_results",
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    configure_logging(
        args.verbose
    )

    try:
        if args.command == "run":
            run_pipeline(
                args.config
            )
        elif args.command == "validate":
            validate_command(
                args
            )
        else:
            parser.error(
                f"Unsupported command: {args.command}"
            )
    except Exception as error:
        LOGGER.error("%s", error)

        if args.verbose:
            raise

        sys.exit(1)


if __name__ == "__main__":
    main()
