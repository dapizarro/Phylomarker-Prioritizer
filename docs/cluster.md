# PhyloMarker Phylogeny 0.2

Pipeline reproducible para evaluar paneles de **PhyloMarker Select** mediante
árboles génicos, concatenación particionada, ASTRAL, factores de concordancia,
jackknife y comparación topológica.

## Correcciones e integración de v0.2

- `prepare` genera un NEXUS válido con `charset`; ya no escribe `AA, gen = ...`.
- `gene-trees` ejecuta ModelFinder por gen cuando `gene_model: MFP`.
- `gene-models` extrae el modelo BIC de los `.iqtree` ya generados.
- `concatenated` admite tres estrategias:
  - `gene_models`: modelo específico por gen, recomendada;
  - `model_finder_merge`: `MFP+MERGE`;
  - `fixed`: modelo global fijo.
- `concordance` usa automáticamente el árbol y las particiones de la estrategia elegida.
- `compare-trees` calcula RF no enraizada con DendroPy, evitando el problema observado con IQ-TREE 2.0.7.
- Los directorios de salida se crean antes de ejecutar los comandos.
- Los resultados originales no se sobrescriben cuando se usa `gene_models`.

## Instalación

```bash
mamba env create -f environment.yml
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phylomarker-phylogeny
python -m pip install -e ".[test]"
pytest -v
```

## Flujo local o SLURM

```bash
phylomarker-phylogeny validate --config phylogeny.cluster.yaml
phylomarker-phylogeny prepare --config phylogeny.cluster.yaml
phylomarker-phylogeny gene-trees --config phylogeny.cluster.yaml
```

Tras terminar los árboles génicos:

```bash
phylomarker-phylogeny gene-models --config phylogeny.cluster.yaml
phylomarker-phylogeny concatenated --config phylogeny.cluster.yaml
phylomarker-phylogeny astral --config phylogeny.cluster.yaml
```

Tras terminar concatenación y ASTRAL:

```bash
phylomarker-phylogeny concordance --config phylogeny.cluster.yaml
phylomarker-phylogeny compare-trees --config phylogeny.cluster.yaml
phylomarker-phylogeny summarize --config phylogeny.cluster.yaml
```

Los comandos costosos generan manifiestos y arrays SLURM. Añada `--execute` solo para
una prueba local.

## Configuración recomendada

```yaml
iqtree:
  executable: iqtree2
  gene_model: MFP
  concat_partition_strategy: gene_models
  missing_model_policy: error
  fallback_model: LG+G
```

Con `gene_models`, el pipeline reutiliza el ModelFinder de los árboles génicos:
no repite la selección del modelo.

## Salidas principales

```text
summaries/gene_models.tsv
panels/<perfil>/<perfil>.gene_models.nex
concatenated_gene_models/<perfil>/<perfil>.treefile
concordance_gene_models/<perfil>/<perfil>.cf.tree
summaries/tree_comparisons.tsv
```

## Precaución científica

RF, gCF y sCF miden aspectos diferentes. Una RF pequeña no demuestra ausencia de
ILS, y gCF bajo puede deberse a error de árbol génico, señal insuficiente,
paralogía, composición, introgressión u otros procesos.
