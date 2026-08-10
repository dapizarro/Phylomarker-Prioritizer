"""Carga de la tabla de metadatos de muestras."""
from __future__ import annotations

from pathlib import Path

import pandas as pd


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
