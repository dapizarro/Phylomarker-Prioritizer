"""Utilidades compartidas sin dependencias del dominio."""
from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path


def executable_version(
    executable: str,
    argument: str = "--version",
) -> str:
    """Primera linea util de `<executable> --version`.

    Nunca lanza: la procedencia se escribe antes de comprobar que las
    herramientas externas existan, y un fallo aqui no debe adelantarse al
    mensaje claro de `require_executable`.

    MAFFT escribe su version en stderr, asi que se capturan ambos flujos.
    """
    try:
        completed = subprocess.run(
            [executable, argument],
            capture_output=True,
            text=True,
            timeout=20,
        )
    except (OSError, subprocess.SubprocessError) as error:
        return f"unavailable: {error}"

    lines = [
        line.strip()
        for line in (
            (completed.stdout or "")
            + "\n"
            + (completed.stderr or "")
        ).splitlines()
        if line.strip()
    ]

    return (
        lines[0]
        if lines
        else f"unavailable: returncode={completed.returncode}"
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()
