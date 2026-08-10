"""Registro de configuracion resuelta y versiones de software."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


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
