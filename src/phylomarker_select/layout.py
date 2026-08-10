"""Rutas del directorio de salida.

Punto unico donde se decide la estructura de `inputs.output`. Antes cada
fragmento (`validation`, `alignments/<tipo>/trimmed`, `rankings`, ...) estaba
escrito a mano en la funcion que lo usaba y recalculado otra vez en el
orquestador.

Dos de estas rutas son contrato con la etapa 2 y estan marcadas como tal:
`panel_genes_file()` y `trimmed_alignment()`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class OutputLayout:
    """Estructura del arbol de resultados de una corrida de seleccion."""

    root: Path
    sequence_type: str
    balance_level: str

    @property
    def extension(self) -> str:
        """Extension de fichero de secuencia, con punto inicial."""
        return (
            ".faa"
            if self.sequence_type == "protein"
            else ".fna"
        )

    # --- directorios ---------------------------------------------------

    @property
    def validation_directory(self) -> Path:
        return self.root / "validation"

    @property
    def per_gene_sequences_directory(self) -> Path:
        return (
            self.root
            / "sequences"
            / self.sequence_type
            / "per_gene"
        )

    @property
    def untrimmed_directory(self) -> Path:
        return (
            self.root
            / "alignments"
            / self.sequence_type
            / "untrimmed"
        )

    @property
    def trimmed_directory(self) -> Path:
        return (
            self.root
            / "alignments"
            / self.sequence_type
            / "trimmed"
        )

    @property
    def metrics_directory(self) -> Path:
        return self.root / "metrics"

    @property
    def rankings_directory(self) -> Path:
        return self.root / "rankings"

    @property
    def pca_directory(self) -> Path:
        return self.root / "pca"

    @property
    def panels_directory(self) -> Path:
        return self.root / "panels"

    @property
    def report_directory(self) -> Path:
        return self.root / "report"

    @property
    def provenance_directory(self) -> Path:
        return self.root / "provenance"

    # --- validacion ----------------------------------------------------

    @property
    def discovered_runs_table(self) -> Path:
        return self.validation_directory / "discovered_runs.tsv"

    @property
    def validated_runs_table(self) -> Path:
        return self.validation_directory / "validated_runs.tsv"

    @property
    def warnings_table(self) -> Path:
        return self.validation_directory / "warnings.tsv"

    @property
    def sequence_provenance_table(self) -> Path:
        return (
            self.validation_directory
            / f"{self.sequence_type}_sequence_provenance.tsv"
        )

    # --- alineamientos -------------------------------------------------

    def untrimmed_alignment(self, gene_id: str) -> Path:
        return (
            self.untrimmed_directory
            / f"{gene_id}.aln{self.extension}"
        )

    def trimmed_alignment(self, gene_id: str) -> Path:
        """Alineamiento recortado de un gen.

        CONTRATO con la etapa 2: `phylomarker_phylogeny.core.alignment_path()`
        reconstruye esta misma ruta. Cambiarla rompe la etapa 2 en silencio.
        """
        return (
            self.trimmed_directory
            / f"{gene_id}.trimmed{self.extension}"
        )

    def analysis_alignment_directory(
        self,
        trimming_enabled: bool,
    ) -> Path:
        """Alineamientos sobre los que se calculan metricas y paneles."""
        return (
            self.trimmed_directory
            if trimming_enabled
            else self.untrimmed_directory
        )

    # --- metricas y rankings -------------------------------------------

    @property
    def per_gene_metrics_table(self) -> Path:
        return self.metrics_directory / "per_gene_metrics.tsv"

    @property
    def group_occupancy_table(self) -> Path:
        return (
            self.metrics_directory
            / f"per_gene_{self.balance_level}_occupancy.tsv"
        )

    @property
    def all_gene_scores_table(self) -> Path:
        return self.rankings_directory / "all_gene_scores.tsv"

    def profile_ranking_table(self, profile: str) -> Path:
        return self.rankings_directory / f"{profile}.tsv"

    # --- PCA -----------------------------------------------------------

    @property
    def pca_scores_table(self) -> Path:
        return self.pca_directory / "pca_scores.tsv"

    @property
    def pca_loadings_table(self) -> Path:
        return self.pca_directory / "pca_loadings.tsv"

    @property
    def pca_explained_variance_table(self) -> Path:
        return self.pca_directory / "pca_explained_variance.tsv"

    # --- paneles -------------------------------------------------------

    def panel_directory(self, profile: str, size: int) -> Path:
        return self.panels_directory / profile / f"n{size}"

    def panel_genes_file(self, profile: str, size: int) -> Path:
        """Lista plana de genes del panel, una linea por gen, sin cabecera.

        CONTRATO con la etapa 2: `phylomarker_phylogeny.core.panel_genes()`
        lee este fichero y exige exactamente `size` genes sin duplicados.
        """
        return self.panel_directory(profile, size) / "genes.txt"

    def panel_genes_table(self, profile: str, size: int) -> Path:
        return self.panel_directory(profile, size) / "panel_genes.tsv"

    def panel_summary_table(self, profile: str, size: int) -> Path:
        return self.panel_directory(profile, size) / "panel_summary.tsv"

    def selection_trace_table(self, profile: str, size: int) -> Path:
        return self.panel_directory(profile, size) / "selection_trace.tsv"

    def panel_manifest(self, profile: str, size: int) -> Path:
        return self.panel_directory(profile, size) / "panel.yaml"

    def panel_alignments_directory(
        self,
        profile: str,
        size: int,
    ) -> Path:
        return self.panel_directory(profile, size) / "alignments"

    # --- informe y procedencia -----------------------------------------

    @property
    def report_index(self) -> Path:
        return self.report_directory / "index.html"

    @property
    def software_versions_table(self) -> Path:
        return self.provenance_directory / "software_versions.tsv"

    @property
    def resolved_config_file(self) -> Path:
        return self.provenance_directory / "resolved_config.yaml"
