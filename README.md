# PhyloMarker Prioritizer

Selección de paneles de marcadores filogenómicos a partir de salidas BUSCO, con
validación filogenética independiente.

El programa tiene dos etapas encadenadas:

| Etapa | Comando | Qué hace |
|---|---|---|
| 1. Selección | `phylomarker-select` | Caracteriza cada marcador BUSCO y construye paneles optimizados para objetivos evolutivos explícitos |
| 2. Validación | `phylomarker-phylogeny` | Evalúa esos paneles con árboles génicos, concatenación particionada, ASTRAL, gCF/sCF, jackknife y comparación topológica |

La etapa 2 consume directamente el directorio de resultados de la etapa 1: lee
`panels/<perfil>/n<tamaño>/genes.txt` y los alineamientos asociados.

## Instalación

```bash
mamba env create -f environment.yml
mamba activate phylomarker
python -m pip install -e ".[test]"
pytest
```

Herramientas externas requeridas en `PATH`:

- etapa 1: `mafft`, `trimal`
- etapa 2: `iqtree2` (o `iqtree3`), `java` y un `.jar` de ASTRAL

## Uso

### Etapa 1 — selección

```bash
phylomarker-select validate \
    --busco-directory /ruta/dicarya_BUSCO \
    --metadata data/samples.tsv \
    --output validation_dicarya

phylomarker-select --verbose run --config configs/select.dikarya.yaml
```

Se ejecuta de principio a fin; no es reanudable. Produce:

```text
<output>/
├── provenance/          resolved_config.yaml, software_versions.tsv
├── validation/
├── sequences/protein/per_gene/
├── alignments/protein/{untrimmed,trimmed}/
├── metrics/             per_gene_metrics.tsv, per_gene_order_occupancy.tsv
├── pca/                 scores, loadings, varianza explicada
├── rankings/            all_gene_scores.tsv, <perfil>.tsv
├── panels/<perfil>/n<N>/  genes.txt, panel_genes.tsv, panel_summary.tsv,
│                          selection_trace.tsv, panel.yaml, alignments/
└── report/index.html
```

### Etapa 2 — validación filogenética

```bash
phylomarker-phylogeny local-run --config configs/phylogeny.local.n10.yaml
```

`local-run` encadena los 10 pasos y es **reanudable**: omite cualquier
`.treefile`, `.cf.tree` o árbol ASTRAL ya existente y no vacío.

Para el cluster, cada paso se invoca por separado. Los pasos caros
(`gene-trees`, `concatenated`, `astral`, `concordance`, `jackknife`) **sólo
escriben manifiestos y arrays SLURM** salvo que añadas `--execute`:

```bash
phylomarker-phylogeny validate    --config configs/phylogeny.cluster.yaml
phylomarker-phylogeny prepare     --config configs/phylogeny.cluster.yaml
phylomarker-phylogeny gene-trees  --config configs/phylogeny.cluster.yaml
# → sbatch <output>/slurm/gene_trees.sbatch
```

`scripts/run_cluster_pipeline.sh` automatiza el encadenamiento con dependencias
SLURM. Detalles en [`docs/cluster.md`](docs/cluster.md).

## Principio científico

Un marcador no es universalmente bueno. Su utilidad depende de la escala
evolutiva, el muestreo taxonómico, la representación de clados, el occupancy, la
completitud, la calidad del alineamiento, la tasa evolutiva, la información
filogenética, los sesgos composicionales y la redundancia con el resto del panel.

Decisiones de diseño que el código sostiene y los tests protegen:

- **Los paneles no son los N primeros del ranking.** La selección es una ganancia
  marginal voraz `Δ_g(P) = S_g − λ·max_{h∈P} similitud(g,h)`, restringida a un
  pool de candidatos de alta calidad.
- **El PCA es exploratorio.** PC1 no es una puntuación de calidad y no ordena los
  genes (`pca.use_for_ranking: false`).
- **La elegibilidad es una barrera dura previa al ranking.** Una puntuación alta
  de información nunca rescata un gen inelegible.
- **`diverse_rate`** reparte el panel entre cinco cuantiles de
  `mean_pairwise_distance` en lugar de concentrarlo en un régimen de tasas.
- **`occupancy_only`, `information_only` y `random_matched` son controles
  negativos deliberados** para comprobar si el enfoque multicriterio supera a
  estrategias simples o aleatorias.
- **Reutilización de modelos.** Con `concat_partition_strategy: gene_models`,
  ModelFinder se ejecuta una sola vez durante `gene-trees`; la concatenación
  reaprovecha el mejor modelo BIC por gen en lugar de repetir la selección.

Documentación detallada de cada etapa en [`docs/select.md`](docs/select.md) y
[`docs/phylogeny.md`](docs/phylogeny.md).

## Estado del análisis Dikarya

Conjunto piloto de 14 genomas de hongos liquenizados (`data/samples.tsv`).

- Etapa 1 completa. Corrida canónica: `results_dicarya_trimmed_v4`.
- Etapa 2 parcial: pasos 1–8 terminados (46 genes únicos en 6 paneles de 10);
  el jackknife se interrumpió en 31 de 120 réplicas. `local-run` reanuda.
- RF concatenación–ASTRAL: 0 en `core_complete`, `diverse_rate`,
  `information_only` y `random_matched`; 0,091 en `deep_robust` y
  `occupancy_only`.

Los datos pesados no están versionados. Ver [`data/README.md`](data/README.md).

## Limitaciones

Este software no demuestra que un panel produzca el árbol de especies correcto,
ni que esté libre de paralogía oculta o introgresión, ni distingue ILS de error
de estimación. RF, gCF y sCF miden aspectos diferentes: una RF pequeña no
descarta ILS, y un gCF bajo puede deberse a error de árbol génico, señal
insuficiente, paralogía, composición u otros procesos. Las métricas son
descriptivas del alineamiento observado, no garantías de calidad.

## Licencia

MIT. Ver [`LICENSE`](LICENSE).
