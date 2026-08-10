"""Registro de configuracion resuelta y versiones de software."""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import yaml

from .layout import OutputLayout
from .utils import executable_version


def write_provenance(
    layout: OutputLayout,
    config: dict,
    mafft_executable: str = "mafft",
    trimal_executable: str = "trimal",
) -> None:
    layout.provenance_directory.mkdir(
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
        #  MAFFT y trimAl determinan los alineamientos, y por tanto todo lo
        #  que se calcula despues. Sin su version la procedencia no permite
        #  interpretar una corrida.
        {
            "software": "mafft",
            "version": executable_version(
                mafft_executable,
            ),
        },
        {
            "software": "trimal",
            "version": executable_version(
                trimal_executable,
            ),
        },
    ]

    pd.DataFrame(versions).to_csv(
        layout.software_versions_table,
        sep="\t",
        index=False,
    )

    with layout.resolved_config_file.open(
        "w",
        encoding="utf-8",
    ) as handle:
        yaml.safe_dump(
            config,
            handle,
            sort_keys=False,
        )
