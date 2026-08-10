"""Metricas por columna del alineamiento."""
from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

import numpy as np

from ..constants import MISSING_CHARACTERS
from ..fasta import read_fasta


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
