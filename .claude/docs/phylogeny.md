# Etapa 2 — `phylomarker-phylogeny` (referencia de implementación)

Referencia interna para **trabajar sobre el código**. La prosa de usuario está en
`docs/phylogeny.md` y `docs/cluster.md`; este documento no la repite.

Dos módulos, deliberadamente separados:

- `src/phylomarker_phylogeny/cli.py` (375 líneas) — una función delgada por
  subcomando, más el `argparse`.
- `src/phylomarker_phylogeny/core.py` (220 líneas) — FASTA, concatenación,
  particiones NEXUS, construcción de comandos, manifiestos, SLURM, RF, jackknife.

`cli.py` hace `from .core import *`, así que todo lo de `core` está en su espacio
de nombres sin prefijo.

---

## 1. El principio rector: *manifest-first*

Ningún subcomando ejecuta nada por defecto. `gene-trees`, `concatenated`,
`astral`, `concordance` y `jackknife` escriben **un manifiesto de comandos** y
**un array de SLURM**, y solo ejecutan si se pasa `--execute`. `local-run` lo
fuerza a `True` (cli.py:329).

Esto es lo que permite que **un mismo config sirva en portátil y en cluster**. Al
añadir un subcomando nuevo, respetar el patrón: construir comandos → `write_manifest`
→ `emit_slurm` → ejecutar solo bajo `--execute`.

```
comandos  ──→ manifests/<paso>.commands.txt
          └─→ slurm/<paso>.sbatch        (#SBATCH --array=1-N, sed -n "${TASK_ID}p")
```

`write_slurm_array` (core.py:187) genera un array donde cada tarea extrae su línea
del manifiesto con `sed`. Por eso **un comando por línea** es un invariante del
formato del manifiesto.

## 2. Reanudabilidad por existencia de salida

No hay base de datos de estado. `run_command` (core.py:150) llama a
`command_expected_output` (core.py:139), que deduce el fichero de terminación a
partir del propio comando:

| Patrón en el comando | Salida esperada |
|---|---|
| `--prefix P` con `--gcf` o `--scf` | `P.cf.tree` |
| `--prefix P` (resto) | `P.treefile` |
| `-o F` (ASTRAL) | `F` |

Si ese fichero existe y **no está vacío**, se salta el comando y se registra
`SKIPPED completed output:`. Consecuencias prácticas:

- Un `local-run` interrumpido se reanuda solo. Es exactamente el mecanismo que
  recupera el jackknife parado a mitad.
- Un fichero de salida truncado pero no vacío se toma por completo. Ante sospecha
  de corrupción, borrar el fichero a mano.
- Un comando sin `--prefix` ni `-o` devuelve `None` y **siempre** se re-ejecuta.

## 3. Acoplamiento con la etapa 1

`context()` (cli.py:7) resuelve el config y devuelve la tupla que usan todos los
subcomandos: `(cfg, base, out, sel, sequence_type, profiles, size, taxa)`.

`inputs.selection_results` (`sel`) apunta al `inputs.output` de la etapa 1. Dos
funciones de `core` leen ese árbol directamente:

- `panel_genes(sel, profile, size)` (core.py:48) →
  `<sel>/panels/<profile>/n<size>/genes.txt`.
  Valida que haya **exactamente** `size` genes y sin duplicados.
- `alignment_path(sel, gene, seq)` (core.py:56) →
  `<sel>/alignments/<seq>/trimmed/<gene>.trimmed.<faa|fna>`.

**Cambiar el layout de salida de la etapa 1 rompe la etapa 2 en silencio.** Las
rutas están codificadas aquí, no derivadas de ningún manifiesto compartido.

`inputs.metadata` es opcional; si está, fija el **orden de taxones** de todas las
concatenaciones y hace que un taxón desconocido en un alineamiento aborte
(cli.py:30 y core.py:70).

## 4. Subcomandos

Registrados en `parser()` (cli.py:356). Solo los cinco marcados aceptan `--execute`.

| Subcomando | Función | Qué hace | `--execute` |
|---|---|---|---|
| `validate` | 21 | Lee cada alineamiento del panel, comprueba taxones contra metadata → `validation.tsv` | — |
| `prepare` | 38 | Enlaza/copia alineamientos, concatena, escribe `partitions.nex` y `concatenation_stats.json` | — |
| `gene-trees` | 58 | Un IQ-TREE por gen **único** entre todos los perfiles | ✓ |
| `gene-models` | 78 | Raspa el mejor modelo BIC de cada `.iqtree` → `gene_models.nex` | — |
| `concatenated` | 105 | Árbol concatenado particionado | ✓ |
| `astral` | 147 | Recoge árboles génicos y lanza el jar de ASTRAL | ✓ |
| `concordance` | 164 | gCF/sCF de IQ-TREE sobre el árbol concatenado | ✓ |
| `compare-trees` | 191 | RF concatenado vs ASTRAL → `tree_comparisons.tsv` | — |
| `jackknife` | 210 | Réplicas deterministas con submuestras de genes | ✓ |
| `jackknife-summary` | 236 | RF de cada réplica contra el árbol completo | — |
| `summarize` | 295 | Estado del pipeline + `provenance.json` | — |
| `local-run` / `all` | 327 / 353 | Los 10 pasos en orden, con `--execute` forzado | (implícito) |

`gene_trees` deduplica genes entre perfiles (`unique_genes`, cli.py:61): un gen
compartido por varios paneles se calcula **una sola vez**, y los perfiles
comparten el resultado. Es lo que hace barata la comparación entre perfiles.

## 5. Estrategias de partición y reutilización de modelos

`iqtree.concat_partition_strategy` gobierna `concatenated` y `concordance`:

| Estrategia | Partición usada | Modelo pasado a IQ-TREE | Raíz de salida |
|---|---|---|---|
| `gene_models` *(recomendada)* | `<profile>.gene_models.nex` | ninguno (va dentro del NEXUS) | `concatenated_gene_models/`, `concordance_gene_models/` |
| `model_finder_merge` | `<profile>.partitions.nex` | `MFP+MERGE` | `concatenated/`, `concordance/` |
| `fixed` | `<profile>.partitions.nex` | `iqtree.concat_fixed_model` | `concatenated/`, `concordance/` |

Con `gene_models`, **ModelFinder corre una sola vez**, en `gene-trees` (con
`gene_model: MFP`). Después `gene-models` lee cada `<gene>.iqtree` con
`read_iqtree_best_model` (core.py:95), que prueba tres regex en orden
(`Best-fit model according to BIC:`, `Best-fit substitution model:`,
`Model of substitution:`) y falla con un mensaje que recuerda usar `-m MFP`.

`missing_model_policy`:

- `error` *(defecto y recomendado)* — si algún gen no tiene modelo, aborta y
  remite a `summaries/gene_models.tsv`. **Falla ruidosamente en vez de usar en
  silencio un modelo por defecto.**
- `fallback` — sustituye por `iqtree.fallback_model` y lo marca en la columna
  `status`.

Las raíces separadas (`concatenated_gene_models/` frente a `concatenated/`)
existen para que **cambiar de estrategia nunca sobrescriba una corrida anterior**.
Ambas pueden coexistir.

`write_gene_model_nexus` (core.py:103) emite `charset` + `charpartition
gene_models = MODELO:gen, …;` y aborta si falta el modelo de algún gen.

## 6. Concatenación y particiones

`concatenate` (core.py:62):

- `read_fasta` (core.py:22) exige alineamiento real: **todas las secuencias con la
  misma longitud**, sin IDs duplicados, no vacío. Secuencia en mayúsculas.
- Un taxón ausente de un gen se rellena con `-` de la longitud de ese gen (74).
- Un taxón presente en un alineamiento pero ausente de la metadata es error (70).
- Devuelve `{n_taxa, n_genes, alignment_length, partitions}`, serializado en
  `concatenation_stats.json`. `gene-models` y `summarize` releen ese JSON en vez
  de recalcular.

`write_charset_nexus` (core.py:82) emite un NEXUS **solo de charsets** — válido
para `iqtree -p`, sin bloque de datos.

## 7. Distancia RF: DendroPy, no IQ-TREE

`compare_unrooted_trees` (core.py:116) usa `dendropy.calculate.treecompare`.
**No sustituir por el RF de IQ-TREE**: existe un bug en IQ-TREE 2.0.7 que motivó
este cambio.

Los árboles se leen con `rooting="force-unrooted"`, `preserve_underscores=True`,
`suppress_internal_node_taxa=True` y un `TaxonNamespace` compartido — los tres
últimos son necesarios para que las biparticiones sean comparables. Si los
conjuntos de hojas difieren, se lanza un error que lista las diferencias en vez
de devolver un número sin sentido.

Devuelve `rf`, `splits_only_first`, `splits_only_second`, `rf_max_binary`
(`2(n−3)`) y `rf_normalized`.

## 8. Jackknife

`deterministic_subsets` (core.py:201): `k = max(2, min(n, round(n · keep_fraction)))`,
muestreo con `random.Random(seed)`, **subconjuntos únicos** (se descartan
repetidos), hasta `n_replicates · 100` intentos. La semilla efectiva es
`seed + índice_del_perfil` (cli.py:216), de modo que cada perfil recibe
submuestras distintas pero reproducibles.

Cada réplica concatena su propio subconjunto y corre IQ-TREE con
`jackknife.model` (por defecto `MFP+MERGE`) y semilla `seed + r`.

`jackknife-summary` compara cada `tree.treefile` con el árbol concatenado
completo. Las réplicas que faltan se registran con `status: "missing"` y no
contaminan las medias — por eso el resumen es útil sobre una corrida a medias.
Escribe `jackknife_replicates.tsv` (por réplica) y `jackknife_summary.tsv`
(agregado: `identical_fraction`, `rf_normalized_{mean,median,max}`).

> `MFP+MERGE` repite la selección **y fusión** de modelos en cada réplica: es el
> cuello de botella en portátil. Para probar en local: `model: LG+G` y
> `replicates: 5`. Ver el comentario en `configs/phylogeny.local.n10.yaml:57`.

## 9. Layout de salida

```
<inputs.output>/
├── validation.tsv, validation_warnings.tsv
├── panels/<profile>/
│   ├── genes/                      symlinks (o copias si prepare.copy_alignments)
│   ├── gene_manifest.tsv
│   ├── <profile>.concat.<faa|fna>
│   ├── <profile>.partitions.nex
│   ├── <profile>.gene_models.nex   (tras gene-models)
│   └── concatenation_stats.json
├── gene_trees/by_gene/<gene>/<gene>.{treefile,iqtree,log,…}
├── concatenated[_gene_models]/<profile>/<profile>.treefile
├── astral/<profile>/{gene_trees.tre,species_tree.tre}
├── concordance[_gene_models]/<profile>/<profile>.cf.tree
├── jackknife/<profile>/replicate_NNN/{*.faa,partitions.nex,tree.treefile}
├── jackknife/jackknife_manifest.tsv
├── manifests/*.commands.txt
├── slurm/*.sbatch
├── logs/*.log
└── summaries/
    ├── gene_models.tsv
    ├── tree_comparisons.tsv
    ├── jackknife_replicates.tsv, jackknife_summary.tsv
    ├── pipeline_status.tsv
    └── provenance.json
```

## 10. Estilo y trampas

- Estilo **compacto y terso**: punto y coma, funciones de una línea, imports
  agrupados. No arrastrar aquí el estilo vertical de la etapa 1.
- `run_command` ejecuta con `bash -lc`, así que el comando pasa por el shell de
  login — el entorno del `.bashrc` está disponible, pero también sus efectos.
- `main()` (cli.py:369) traga la excepción y sale con código 1 salvo que se pase
  `--debug` en la línea de órdenes, que la vuelve a lanzar con traza completa.
  `--debug` no está declarado en el `argparse`: se detecta mirando `sys.argv`.
- `astral.jar` no es un paquete de conda: la ruta va en el config
  (`configs/phylogeny.local.n10.yaml:54`) y se comprueba antes de emitir comandos.
- El ejecutable de IQ-TREE se toma de `iqtree.executable` (`iqtree2` en el config
  local, `iqtree3` como defecto en código). Conviene verificar cuál resuelve el
  PATH: puede haber un `iqtree2` del sistema fuera del entorno conda.
- Tests: `tests/phylogeny/test_core.py` (69 líneas) y `test_v03.py` (21 líneas).
