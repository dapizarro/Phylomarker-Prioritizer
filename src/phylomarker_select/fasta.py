"""Lectura y escritura de FASTA."""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

LOGGER = logging.getLogger("phylomarker-select")


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    sequence: str


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
