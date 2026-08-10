"""Alineamiento con MAFFT y recorte con trimAl."""
from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

LOGGER = logging.getLogger("phylomarker-select")


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
