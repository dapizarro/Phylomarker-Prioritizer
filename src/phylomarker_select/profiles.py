"""Registro de perfiles de panel y ranking por perfil.

Anadir un perfil es anadir una entrada a PROFILES. Antes habia que tocar tres
sitios: PROFILE_WEIGHTS, PROFILE_EXCLUDED_TRIMMING_CLASSES y las ramas por
nombre literal de `create_panels`.

CUIDADO al editar los pesos: `calculate_profile_ranking` acumula la puntuacion
iterando `weights.items()`, asi que **el orden de insercion decide el orden de
la suma en coma flotante**. Reordenar las claves de un perfil puede mover los
ultimos bits de `profile_score` y cambiar desempates del ranking.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

import numpy as np
import pandas as pd

from .scoring import percentile_score

# Optimizadores admitidos, tal como los despacha `create_panels`.
GREEDY = "greedy"
DIVERSE_RATE = "diverse_rate"
RANDOM = "random"


@dataclass(frozen=True)
class Profile:
    """Un objetivo evolutivo explicito y como se construye su panel."""

    name: str
    weights: Mapping[str, float] = field(default_factory=dict)
    excluded_trimming_classes: frozenset[str] = frozenset()
    optimizer: str = GREEDY
    #  `random_matched` no pasa por el ranking, asi que no hay traza que
    #  escribir ni manifiesto de parametros de optimizacion.
    writes_trace: bool = True
    is_negative_control: bool = False

    @property
    def is_rankable(self) -> bool:
        """Si tiene pesos, `calculate_profile_ranking` puede puntuarlo."""
        return bool(self.weights)


PROFILES: dict[str, Profile] = {
    profile.name: profile
    for profile in (
        Profile(
            name="core_complete",
            weights={
                "coverage_percentile": 0.35,
                "alignment_percentile": 0.25,
                "trimming_percentile": 0.15,
                "information_percentile": 0.15,
                "low_bias_percentile": 0.10,
            },
            excluded_trimming_classes=frozenset({"extreme"}),
        ),
        Profile(
            name="backbone_balanced",
            weights={
                "coverage_percentile": 0.30,
                "clade_balance_percentile": 0.25,
                "alignment_percentile": 0.20,
                "trimming_percentile": 0.15,
                "information_percentile": 0.10,
            },
            excluded_trimming_classes=frozenset({"extreme"}),
        ),
        Profile(
            name="deep_robust",
            weights={
                "alignment_percentile": 0.25,
                "coverage_percentile": 0.20,
                "trimming_percentile": 0.20,
                "rate_percentile": 0.20,
                "low_bias_percentile": 0.15,
            },
            excluded_trimming_classes=frozenset(
                {"extreme", "signal_sensitive"}
            ),
        ),
        Profile(
            name="low_bias",
            weights={
                "low_bias_percentile": 0.35,
                "alignment_percentile": 0.25,
                "trimming_percentile": 0.20,
                "coverage_percentile": 0.10,
                "information_percentile": 0.10,
            },
            excluded_trimming_classes=frozenset({"extreme"}),
        ),
        Profile(
            name="diverse_rate",
            weights={
                "coverage_percentile": 0.25,
                "alignment_percentile": 0.25,
                "trimming_percentile": 0.20,
                "information_percentile": 0.20,
                "low_bias_percentile": 0.10,
            },
            excluded_trimming_classes=frozenset({"extreme"}),
            optimizer=DIVERSE_RATE,
        ),
        # --- controles negativos: no son perfiles a "mejorar" ---
        Profile(
            name="occupancy_only",
            weights={
                "coverage_percentile": 1.00,
            },
            excluded_trimming_classes=frozenset({"extreme"}),
            is_negative_control=True,
        ),
        Profile(
            name="information_only",
            weights={
                "information_percentile": 1.00,
            },
            excluded_trimming_classes=frozenset(),
            is_negative_control=True,
        ),
        Profile(
            name="random_matched",
            excluded_trimming_classes=frozenset(),
            optimizer=RANDOM,
            writes_trace=False,
            is_negative_control=True,
        ),
    )
}


def get_profile(name: str) -> Profile:
    """Perfil del registro. Lanza ValueError si no existe."""
    if name not in PROFILES:
        raise ValueError(
            f"Unknown panel profile: {name}"
        )

    return PROFILES[name]


def calculate_profile_ranking(
    scored: pd.DataFrame,
    profile: str,
) -> pd.DataFrame:
    #  `random_matched` no tiene pesos y nunca ha pasado por aqui: antes no
    #  estaba en PROFILE_WEIGHTS y lanzaba este mismo error.
    selected = PROFILES.get(profile)

    if selected is None or not selected.is_rankable:
        raise ValueError(
            f"Unknown panel profile: {profile}"
        )

    ranked = scored.copy()

    percentile_sources = {
        "coverage_percentile": ("coverage_score", True),
        "clade_balance_percentile": ("clade_balance_score", True),
        "alignment_percentile": ("alignment_score", True),
        "information_percentile": ("information_score", True),
        "rate_percentile": ("rate_score", True),
        "low_bias_percentile": ("bias_penalty", False),
        "trimming_percentile": ("trimming_stability_score", True),
    }

    for percentile_column, (source_column, higher_is_better) in (
        percentile_sources.items()
    ):
        if (
            percentile_column not in ranked.columns
            and source_column in ranked.columns
        ):
            ranked[percentile_column] = percentile_score(
                ranked[source_column],
                higher_is_better=higher_is_better,
            )

    weights = selected.weights

    total_score = np.zeros(
        len(ranked),
        dtype=float,
    )

    for feature, weight in weights.items():
        total_score += (
            weight
            * pd.to_numeric(
                ranked[feature],
                errors="coerce",
            )
            .fillna(0)
            .to_numpy()
        )

    ranked["profile"] = profile
    ranked["profile_score"] = total_score
    ranked["profile_eligible"] = ranked["eligible"].astype(bool)

    excluded_classes = selected.excluded_trimming_classes

    if excluded_classes and "trimming_class" in ranked.columns:
        ranked["profile_eligible"] &= ~ranked[
            "trimming_class"
        ].isin(excluded_classes)

    ranked.loc[
        ~ranked["profile_eligible"],
        "profile_score",
    ] = -np.inf

    ranked = ranked.sort_values(
        [
            "profile_score",
            "gene_id",
        ],
        ascending=[
            False,
            True,
        ],
    ).reset_index(drop=True)

    ranked["rank"] = np.arange(
        1,
        len(ranked) + 1,
    )

    return ranked
