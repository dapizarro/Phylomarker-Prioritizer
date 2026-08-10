"""Esquema tipado de la configuracion YAML.

Punto unico donde viven los valores por defecto. Antes estaban repartidos como
literales en las llamadas a `config.get(clave, default)` de `run_pipeline`,
`create_panels` y `add_biological_scores`.

Las claves desconocidas se ignoran, igual que antes. No se valida el YAML de
forma estricta a proposito: `configs/select.dikarya.yaml` contiene claves que
el codigo nunca ha leido (ver `TrimmingConfig.retain_untrimmed` y
`PcaConfig.use_for_ranking`) y rechazarlas romperia configuraciones existentes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

DEFAULT_PROFILES: tuple[str, ...] = (
    "core_complete",
    "backbone_balanced",
    "deep_robust",
    "low_bias",
    "diverse_rate",
    "occupancy_only",
    "information_only",
    "random_matched",
)

VALID_SEQUENCE_TYPES = frozenset({"protein", "nucleotide"})


def load_yaml(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(
            f"Configuration file not found: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        data = yaml.safe_load(handle)

    if not isinstance(data, dict):
        raise ValueError(
            "Configuration root must be a mapping"
        )

    return data


def _section(config: dict, name: str) -> dict:
    """Subseccion del YAML, o vacia si no esta."""
    value = config.get(name, {})
    return value if isinstance(value, dict) else {}


@dataclass(frozen=True)
class ProjectConfig:
    name: str | None = None
    sequence_type: str = "protein"
    random_seed: int = 20260729

    @classmethod
    def from_mapping(cls, section: dict) -> ProjectConfig:
        sequence_type = str(
            section.get("sequence_type", "protein")
        )

        if sequence_type not in VALID_SEQUENCE_TYPES:
            raise ValueError(
                "sequence_type must be 'protein' or 'nucleotide'"
            )

        return cls(
            name=section.get("name"),
            sequence_type=sequence_type,
            random_seed=int(
                section.get("random_seed", 20260729)
            ),
        )


@dataclass(frozen=True)
class InputsConfig:
    busco_directory: Path
    metadata: Path
    sample_id_column: str = "sample_ID"
    output: Path = Path("results")

    @classmethod
    def from_mapping(cls, section: dict) -> InputsConfig:
        # Requeridas: sin valor por defecto, igual que antes.
        return cls(
            busco_directory=Path(
                section["busco_directory"]
            ),
            metadata=Path(section["metadata"]),
            sample_id_column=str(
                section.get("sample_id_column", "sample_ID")
            ),
            output=Path(
                section.get("output", "results")
            ),
        )


@dataclass(frozen=True)
class TaxonomyConfig:
    balance_level: str = "order"

    @classmethod
    def from_mapping(cls, section: dict) -> TaxonomyConfig:
        return cls(
            balance_level=str(
                section.get("balance_level", "order")
            ),
        )


@dataclass(frozen=True)
class AlignmentConfig:
    mafft_executable: str = "mafft"
    strategy: str = "auto"
    threads_per_gene: int = 2

    @classmethod
    def from_mapping(cls, section: dict) -> AlignmentConfig:
        return cls(
            mafft_executable=str(
                section.get("mafft_executable", "mafft")
            ),
            strategy=str(
                section.get("strategy", "auto")
            ),
            threads_per_gene=int(
                section.get("threads_per_gene", 2)
            ),
        )


@dataclass(frozen=True)
class TrimmingConfig:
    enabled: bool = True
    trimal_executable: str = "trimal"
    mode: str = "automated1"
    # Descriptiva: los alineamientos sin recortar se conservan siempre, porque
    # `calculate_metrics` los necesita para la estabilidad al recorte. El
    # codigo nunca ha leido esta clave.
    retain_untrimmed: bool = True

    @classmethod
    def from_mapping(cls, section: dict) -> TrimmingConfig:
        return cls(
            enabled=bool(section.get("enabled", True)),
            trimal_executable=str(
                section.get("trimal_executable", "trimal")
            ),
            mode=str(
                section.get("mode", "automated1")
            ),
            retain_untrimmed=bool(
                section.get("retain_untrimmed", True)
            ),
        )


@dataclass(frozen=True)
class EligibilityConfig:
    """Puerta dura aplicada antes de cualquier ranking."""

    min_taxa: int = 4
    min_taxon_occupancy: float = 0.5
    min_alignment_length: int = 80
    min_informative_sites: int = 2
    max_gap_fraction: float = 0.6
    max_ambiguous_fraction: float = 0.1

    @classmethod
    def from_mapping(cls, section: dict) -> EligibilityConfig:
        return cls(
            min_taxa=int(section.get("min_taxa", 4)),
            min_taxon_occupancy=float(
                section.get("min_taxon_occupancy", 0.5)
            ),
            min_alignment_length=int(
                section.get("min_alignment_length", 80)
            ),
            min_informative_sites=int(
                section.get("min_informative_sites", 2)
            ),
            max_gap_fraction=float(
                section.get("max_gap_fraction", 0.6)
            ),
            max_ambiguous_fraction=float(
                section.get("max_ambiguous_fraction", 0.1)
            ),
        )


@dataclass(frozen=True)
class PcaConfig:
    enabled: bool = True
    # Descriptiva: el PCA es exploratorio y PC1 nunca entra en una puntuacion
    # de calidad. Lo garantiza la estructura del codigo, no esta clave, que
    # nunca se ha leido.
    use_for_ranking: bool = False

    @classmethod
    def from_mapping(cls, section: dict) -> PcaConfig:
        return cls(
            enabled=bool(section.get("enabled", True)),
            use_for_ranking=bool(
                section.get("use_for_ranking", False)
            ),
        )


@dataclass(frozen=True)
class PanelsConfig:
    profiles: tuple[str, ...] = DEFAULT_PROFILES
    sizes: tuple[int, ...] = (25, 50, 100)
    redundancy_penalty: float = 0.10
    candidate_top_fraction: float = 0.10
    candidate_minimum_pool_size: int = 50
    maximum_score_drop: float = 0.20
    diverse_rate_bins: int = 5

    @classmethod
    def from_mapping(cls, section: dict) -> PanelsConfig:
        return cls(
            profiles=tuple(
                section.get("profiles", DEFAULT_PROFILES)
            ),
            sizes=tuple(
                int(value)
                for value in section.get(
                    "sizes",
                    (25, 50, 100),
                )
            ),
            redundancy_penalty=float(
                section.get("redundancy_penalty", 0.10)
            ),
            candidate_top_fraction=float(
                section.get("candidate_top_fraction", 0.10)
            ),
            candidate_minimum_pool_size=int(
                section.get("candidate_minimum_pool_size", 50)
            ),
            maximum_score_drop=float(
                section.get("maximum_score_drop", 0.20)
            ),
            diverse_rate_bins=int(
                section.get("diverse_rate_bins", 5)
            ),
        )


@dataclass(frozen=True)
class SelectConfig:
    inputs: InputsConfig
    project: ProjectConfig = field(
        default_factory=ProjectConfig
    )
    taxonomy: TaxonomyConfig = field(
        default_factory=TaxonomyConfig
    )
    alignment: AlignmentConfig = field(
        default_factory=AlignmentConfig
    )
    trimming: TrimmingConfig = field(
        default_factory=TrimmingConfig
    )
    eligibility: EligibilityConfig = field(
        default_factory=EligibilityConfig
    )
    pca: PcaConfig = field(default_factory=PcaConfig)
    panels: PanelsConfig = field(
        default_factory=PanelsConfig
    )
    # El YAML tal cual se leyo. `create_html_report` y `write_provenance` lo
    # vuelcan entero con `yaml.safe_dump`, asi que debe conservarse literal.
    raw: dict = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, config: dict) -> SelectConfig:
        return cls(
            project=ProjectConfig.from_mapping(
                _section(config, "project")
            ),
            inputs=InputsConfig.from_mapping(
                _section(config, "inputs")
            ),
            taxonomy=TaxonomyConfig.from_mapping(
                _section(config, "taxonomy")
            ),
            alignment=AlignmentConfig.from_mapping(
                _section(config, "alignment")
            ),
            trimming=TrimmingConfig.from_mapping(
                _section(config, "trimming")
            ),
            eligibility=EligibilityConfig.from_mapping(
                _section(config, "eligibility")
            ),
            pca=PcaConfig.from_mapping(
                _section(config, "pca")
            ),
            panels=PanelsConfig.from_mapping(
                _section(config, "panels")
            ),
            raw=config,
        )

    @classmethod
    def load(cls, path: Path) -> SelectConfig:
        return cls.from_mapping(load_yaml(path))
