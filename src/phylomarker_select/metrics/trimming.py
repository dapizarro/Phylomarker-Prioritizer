"""Sensibilidad de cada gen al recorte del alineamiento."""
from __future__ import annotations

from pathlib import Path

from .sites import alignment_matrix, safe_fraction, variable_sites_and_pis


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
