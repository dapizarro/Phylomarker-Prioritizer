"""PhyloMarker Select — etapa 1: selección de paneles de marcadores.

Superficie pública estable. El resto de módulos son detalle de implementación
y pueden reorganizarse; ver `.claude/docs/select.md`.
"""

from .config import SelectConfig
from .layout import OutputLayout
from .pipeline import run_pipeline
from .profiles import PROFILES, Profile

__version__ = "0.4.0"

__all__ = [
    "OutputLayout",
    "PROFILES",
    "Profile",
    "SelectConfig",
    "__version__",
    "run_pipeline",
]
