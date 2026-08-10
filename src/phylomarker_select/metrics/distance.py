"""Distancia p media y variabilidad composicional."""
from __future__ import annotations

from collections import Counter

import numpy as np

from ..constants import MISSING_CHARACTERS


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
