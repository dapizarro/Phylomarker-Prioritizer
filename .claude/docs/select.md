# Etapa 1 — `phylomarker-select` (referencia de implementación)

Referencia interna para **trabajar sobre el código**. La prosa orientada a usuario
(instalación, interpretación biológica, límites científicos) está en `docs/select.md`;
este documento no la repite.

Hasta agosto de 2026 todo vivía en un único `cli.py` de ~3.500 líneas. Hoy son 21
módulos. El grafo de llamadas es un DAG sin ciclos y no hay estado mutable a nivel
de módulo, así que se puede leer de abajo arriba.

Sigue sin haber **ejecución parcial ni reanudable**: `run_pipeline()` recorre las
etapas de principio a fin en cada invocación.

---

## 1. Mapa de módulos

```
phylomarker_select/
├── cli.py            argparse, main, configure_logging, validate_command
├── pipeline.py       run_pipeline — orquestación pura, sin lógica propia
├── config.py         esquema tipado del YAML  ← todos los valores por defecto
├── layout.py         OutputLayout             ← todas las rutas de salida
├── profiles.py       registro PROFILES        ← definición de cada perfil
├── constants.py      alfabetos y caracteres ausentes
├── utils.py          sha256_file
├── fasta.py          FastaRecord, read/write/first_fasta_record
├── metadata.py       load_metadata
├── busco.py          discover_busco_runs, validate_runs, extract_markers
├── align.py          MAFFT y trimAl
├── metrics/
│   ├── sites.py      alignment_matrix, site_entropy, variable_sites_and_pis
│   ├── distance.py   mean_pairwise_p_distance, composition_variability
│   ├── trimming.py   classify_trimming y la estabilidad al recorte
│   └── genes.py      calculate_metrics — la tabla por gen
├── scoring.py        robust_standardize, percentile_score, add_biological_scores
├── pca.py            run_exploratory_pca
├── optimize.py       feature_vectors y los dos greedy
├── panels.py         create_panels, summarize_panel, random_panel, export
├── report.py         create_html_report
└── provenance.py     write_provenance
```

**Los tres módulos de arriba son el sitio donde tocar para extender.** Un valor por
defecto nuevo va en `config.py`; una ruta de salida nueva, en `layout.py`; un perfil
nuevo, en `profiles.py`. Si te encuentras escribiendo un literal de ruta o un
`get(clave, default)` fuera de esos tres, es señal de que el cambio va en el sitio
equivocado.

---

## 2. Cadena de ejecución

`run_pipeline()` (`pipeline.py`) está dirigido por el YAML. Orden real:

| # | Llamada | Módulo | Escribe en `<out>/` |
|---|---|---|---|
| 0 | `SelectConfig.load` | `config` | — |
| 1 | `write_provenance` | `provenance` | `provenance/resolved_config.yaml`, `provenance/software_versions.tsv` |
| 2 | `load_metadata` | `metadata` | — |
| 3 | `discover_busco_runs` | `busco` | `validation/discovered_runs.tsv` |
| 4 | `validate_runs` | `busco` | `validation/validated_runs.tsv`, `validation/warnings.tsv` |
| 5 | `extract_markers` | `busco` | `sequences/<tipo>/per_gene/`, `validation/<tipo>_sequence_provenance.tsv` |
| 6 | `align_markers` | `align` | `alignments/<tipo>/untrimmed/` |
| 7 | `trim_alignments` | `align` | `alignments/<tipo>/trimmed/` |
| 8 | `calculate_metrics` | `metrics.genes` | `metrics/per_gene_metrics.tsv`, `metrics/per_gene_<nivel>_occupancy.tsv` |
| 9 | `add_biological_scores` | `scoring` | `rankings/all_gene_scores.tsv` |
| 10 | `run_exploratory_pca` | `pca` | `pca/pca_*.tsv` |
| 11 | `create_panels` | `panels` | `rankings/<perfil>.tsv`, `panels/<perfil>/n<tamaño>/` |
| 12 | `create_html_report` | `report` | `report/index.html` |

**Punto de parada duro**: tras el paso 4, si `warnings` contiene alguna fila con
`severity == "error"`, se lanza `RuntimeError` y no se continúa. Los avisos
no-error se registran y se sigue.

**Qué alineamientos se analizan.** `layout.analysis_alignment_directory(trimming_enabled)`
devuelve `trimmed/` o `untrimmed/`. Los sin recortar **siempre** se conservan,
porque `calculate_metrics` los necesita para la estabilidad al recorte.

---

## 3. Configuración (`config.py`)

Ocho dataclasses congeladas por sección más `SelectConfig.load(path)`. Cada valor
por defecto aparece **una sola vez**, en el campo del dataclass.

`busco_directory` y `metadata` son requeridos (acceso con `[]`, lanzan `KeyError`);
el resto tiene default. `SelectConfig.raw` guarda el dict crudo porque
`create_html_report` y `write_provenance` lo vuelcan entero con `yaml.safe_dump`.

`sequence_type` se valida **al cargar**, no a mitad del flujo: un valor inválido
falla antes de crear el directorio de salida.

**Claves que el YAML declara y el código nunca ha leído** — están modeladas como
campos documentados, no eliminadas, porque `configs/select.dikarya.yaml` las trae
y rechazarlas rompería configuraciones existentes:

| Clave | Por qué no hace nada |
|---|---|
| `trimming.retain_untrimmed` | Los sin recortar se conservan siempre; nunca se borran. |
| `pca.use_for_ranking` | Lo garantiza la estructura del código, no la clave. |
| `project.name` | Puramente documental. |

No hay validación estricta de claves desconocidas, a propósito.

---

## 4. Rutas de salida (`layout.py`)

`OutputLayout(root, sequence_type, balance_level)`, congelado, con una propiedad
por directorio y por fichero. Antes cada fragmento estaba escrito a mano en la
función que lo usaba **y recalculado otra vez** en el orquestador.

**Dos accesores son contrato con la etapa 2** y llevan docstring que lo dice:

- `panel_genes_file(profile, size)` → lo lee `phylomarker_phylogeny.core.panel_genes()`
- `trimmed_alignment(gene)` → lo lee `phylomarker_phylogeny.core.alignment_path()`

Convención de nombres: sin recortar `<gen>.aln.faa`, recortado `<gen>.trimmed.faa`.

**Excepción deliberada**: `cli.validate_command` escribe sus tres TSV en **plano**
bajo `--output`, no bajo `validation/`. Son dos convenciones distintas con los
mismos nombres de fichero. No unificar sin decidir cuál gana.

---

## 5. Métricas por gen (`metrics/`)

`calculate_metrics` (`metrics/genes.py`) escribe una fila por gen. Familias:

- **Cobertura** — `n_taxa`, `taxon_occupancy`, `sequence_completeness`, `cell_occupancy`.
- **Alineamiento** — `alignment_length`, `gap_fraction`, `ambiguous_fraction`.
- **Información** — `variable_sites`, `parsimony_informative_sites`,
  `variable_sites_per_length`, `pis_per_length`, `mean_entropy` (`metrics/sites.py`).
- **Tasa y sesgo** — `mean_pairwise_distance`, `composition_variability`
  (`metrics/distance.py`).
- **Balance de clados** — agregados por `taxonomy.balance_level` (por defecto
  `order`). El desglose va a `metrics/per_gene_<nivel>_occupancy.tsv`.
- **Estabilidad al recorte** — `raw_*` vs `trimmed_*`, `retained_length_fraction`,
  `pis_retained_fraction`, `trimming_class`, `trimming_stability_score`.

### Clasificación del recorte (`metrics/trimming.py`)

Cuatro clases, en este orden:

| Clase | Condición |
|---|---|
| `extreme` | `retained_length_fraction < 0.25` **o** `pis_retained_fraction < 0.50` |
| `stable` | `retained_length_fraction ≥ 0.75` y `pis_retained_fraction ≥ 0.80` |
| `signal_preserved` | `pis_retained_fraction ≥ 0.80` (longitud entre 0,25 y 0,75) |
| `signal_sensitive` | resto |

Sin recorte se fuerza `trimming_class = "not_trimmed"` y score `1.0`.

`calculate_trimming_stability_score` da **más peso a conservar PIS que longitud**:
`0.40 · length_score + 0.60 · information_score`, con penalizaciones multiplicativas
(`×0.25` si la longitud cae bajo 0,25; `×0.50` si los PIS caen bajo 0,50).

---

## 6. Puntuación biológica (`scoring.py`)

Dos transformaciones con papeles distintos — es la confusión más fácil de cometer:

- **`robust_standardize`** — `(x − mediana) / (1.4826 · MAD)`, recortado a `[−4, 4]`.
  Combina **métricas crudas dentro** de una dimensión.
- **`percentile_score`** — rango percentil en `[0, 1]`. Combina **dimensiones entre
  sí** en el ranking por perfil. Con ≤1 valor observado devuelve 0,5 constante.

### Las seis dimensiones

| Dimensión | Composición |
|---|---|
| `coverage_score` | 0,45·`taxon_occupancy` + 0,35·`sequence_completeness` + 0,20·`cell_occupancy` |
| `clade_balance_score` | 0,50·`min_group_sequence_completeness` + 0,30·(−`sd_group_sequence_completeness`) + 0,20·(−`worst_group_gap_fraction`) |
| `alignment_score` | 0,45·`alignment_length` + 0,30·(−`gap_fraction`) + 0,25·(−`ambiguous_fraction`) |
| `information_score` | 0,45·`pis_per_length` + 0,25·`variable_sites_per_length` + 0,20·`mean_entropy` + 0,10·`parsimony_informative_sites` |
| `rate_score` | −\|`mean_pairwise_distance` − mediana\| estandarizado: premia la **cercanía a la tasa mediana**, no la tasa alta ni la baja |
| `bias_penalty` | 0,55·`composition_variability` + 0,30·`gap_fraction` + 0,15·`ambiguous_fraction` (mayor = peor) |

`clade_balance_score` tiene **fallback**: sin las tres columnas de grupo finas usa
0,60·`mean_group_occupancy` + 0,40·`min_replicated_group_occupancy`.

De ahí salen siete percentiles, incluido `low_bias_percentile` (= percentil
**invertido** de `bias_penalty`).

### Elegibilidad — puerta dura previa al ranking

Se calcula en `add_biological_scores`, que recibe un `EligibilityConfig`. Seis
condiciones en AND:

| Campo | Defecto | En `configs/select.dikarya.yaml` |
|---|---|---|
| `min_taxa` | 4 | 7 |
| `min_taxon_occupancy` | 0.5 | 0.50 |
| `min_alignment_length` | 80 | 100 |
| `min_informative_sites` | 2 | 3 |
| `max_gap_fraction` | 0.6 | 0.60 |
| `max_ambiguous_fraction` | 0.1 | 0.10 |

La columna `exclusion_reason` acumula los filtros incumplidos separados por `;`:
es el sitio al que mirar cuando un gen esperado no aparece.

---

## 7. Perfiles (`profiles.py`)

`PROFILES: dict[str, Profile]`. Cada `Profile` declara:

```python
name, weights, excluded_trimming_classes, optimizer, writes_trace, is_negative_control
```

`profile_score = Σ peso · percentil`. Los ocho perfiles:

| Perfil | Pesos | Optimizador |
|---|---|---|
| `core_complete` | coverage 0,35 · alignment 0,25 · trimming 0,15 · information 0,15 · low_bias 0,10 | greedy |
| `backbone_balanced` | coverage 0,30 · clade_balance 0,25 · alignment 0,20 · trimming 0,15 · information 0,10 | greedy |
| `deep_robust` | alignment 0,25 · coverage 0,20 · trimming 0,20 · rate 0,20 · low_bias 0,15 | greedy |
| `low_bias` | low_bias 0,35 · alignment 0,25 · trimming 0,20 · coverage 0,10 · information 0,10 | greedy |
| `diverse_rate` | coverage 0,25 · alignment 0,25 · trimming 0,20 · information 0,20 · low_bias 0,10 | diverse_rate |
| `occupancy_only` | coverage 1,00 — **control negativo** | greedy |
| `information_only` | information 1,00 — **control negativo** | greedy |
| `random_matched` | sin pesos — **control negativo** | random |

> ⚠️ **Al editar los pesos, no reordenes las claves.** `calculate_profile_ranking`
> acumula `total_score += peso · percentil` iterando `weights.items()`, así que el
> orden de inserción decide el orden de la suma en coma flotante. Reordenar puede
> mover los últimos bits de `profile_score` y cambiar desempates del ranking.

**Segundo filtro, por perfil**: `excluded_trimming_classes`. Todos los perfiles
reales excluyen `extreme`; `deep_robust` excluye además `signal_sensitive`.
`information_only` y `random_matched` **no excluyen nada** — es deliberado.

`random_matched` no tiene pesos, así que `Profile.is_rankable` es falso y
`calculate_profile_ranking` lanza `ValueError` si se le pasa, igual que antes.
Tampoco escribe `selection_trace.tsv` ni `panel.yaml` (`writes_trace=False`).

Los genes no elegibles reciben `profile_score = −inf` en lugar de eliminarse, para
que el TSV conserve todos los genes. Orden determinista por
`(profile_score desc, gene_id asc)`.

---

## 8. Optimización de paneles (`optimize.py`)

Los paneles **no son top-N**. Esto es lo que hay que preservar.

### `optimize_panel_greedily`

1. **Filtrar** por `profile_eligible`.
2. **Pool**: `max(candidate_minimum_pool_size, ceil(n · candidate_top_fraction), panel_size · 5)`.
3. **Corte por calidad**: descartar los que caen más de `maximum_score_drop` bajo el
   mejor `profile_score`. Si eso deja menos genes que `panel_size`, **se revierte al
   pool sin recortar** — la red que evita paneles incompletos.
4. **Greedy**, maximizando

   ```
   Δ_g(P) = profile_score(g) − redundancy_penalty · max_{h∈P} similarity(g, h)
   ```

   `similarity` es coseno sobre `feature_vectors`: `taxon_occupancy`,
   `min_replicated_group_occupancy`, `alignment_length`, `pis_per_length`,
   `mean_pairwise_distance`, `composition_variability`, `retained_length_fraction`
   — imputadas con mediana y escaladas con `RobustScaler`. Se descartan las
   features con un solo valor único.

Devuelve `(panel, trace)`; la traza es `selection_trace.tsv`.

### `optimize_diverse_rate_panel`

Mismo pool y penalizador, pero reparte plazas entre **bins de tasa**:

- `pd.qcut` sobre `mean_pairwise_distance` en `panels.diverse_rate_bins` cuantiles
  (5 por defecto), `duplicates="drop"`.
- Cuotas: `base = size // n_bins`, el resto a los primeros bins.
- Round-robin; dentro de cada bin, la misma ganancia marginal. Si una pasada no
  progresa, rellena desde el resto sin restricción de bin.
- **Degradación**: con ≤1 bin efectivo o tasa constante, delega en el greedy
  estándar y marca `rate_bin = "single_bin"`.

Diferencia menor: `diverse_rate` desempata explícitamente por `gene_id`; el greedy
estándar depende del orden de iteración, que ya es determinista por la ordenación
previa.

---

## 9. Salidas — contrato con la etapa 2

```
<inputs.output>/
├── validation/     discovered_runs.tsv, validated_runs.tsv, warnings.tsv
├── sequences/<tipo>/per_gene/
├── alignments/<tipo>/untrimmed/     <gen>.aln.<faa|fna>
├── alignments/<tipo>/trimmed/       <gen>.trimmed.<faa|fna>   ← la etapa 2 lee aquí
├── metrics/        per_gene_metrics.tsv, per_gene_<nivel>_occupancy.tsv
├── rankings/       all_gene_scores.tsv, <perfil>.tsv
├── pca/            pca_scores.tsv, pca_loadings.tsv, pca_explained_variance.tsv
├── panels/<perfil>/n<tamaño>/
│   ├── genes.txt              ← la etapa 2 lee aquí
│   ├── panel_genes.tsv
│   ├── panel_summary.tsv
│   ├── panel.yaml             (no en random_matched)
│   ├── selection_trace.tsv    (no en random_matched)
│   └── alignments/
├── report/index.html
└── provenance/     resolved_config.yaml, software_versions.tsv
```

`genes.txt` es una lista plana sin cabecera. `panel_genes()` verifica que tenga
**exactamente** `n<tamaño>` líneas sin duplicados: un panel más corto de lo nominal
aborta la etapa 2.

---

## 10. Invariantes que los tests bloquean

`tests/select/test_core.py`, 13 tests. Al modificar, comprobar que se mantienen:

1. Los paneles no son top-N — el greedy con penalización debe poder saltarse un gen
   mejor puntuado.
2. La elegibilidad se aplica **antes** del ranking.
3. PC1 nunca es una puntuación de calidad.
4. `diverse_rate` cubre los bins de forma pareja.
5. `occupancy_only`, `information_only` y `random_matched` son controles negativos:
   no "mejorarlos".
6. Determinismo **dado un conjunto de alineamientos**: mismas métricas ⇒ mismo
   panel. Ojo: eso no se extiende al flujo completo, ver abajo.

---

## 10 bis. El flujo completo no es reproducible: MAFFT

`alignment.threads_per_gene: 2` hace que **MAFFT no sea determinista**. Medido
sobre este repositorio en agosto de 2026, con MAFFT v7.526 y `--auto`:

| `--thread` | Tres ejecuciones sobre la misma entrada |
|---|---|
| 1 | idénticas |
| 2 | divergen |
| 4 | divergen |

Consecuencia medida sobre el conjunto Dikarya: **378 de 1.311 genes (29 %)**
recibieron un alineamiento distinto entre dos corridas del *mismo* código con la
*misma* configuración. Eso se propaga a las métricas, al ranking y a la
composición de los paneles: de 16 paneles comparados, 12 cambiaron de contenido.

Implicaciones prácticas:

- Dos corridas del pipeline **no** producen los mismos paneles. Los resultados
  publicados de una corrida concreta no se reproducen volviendo a ejecutar.
- Para comparar dos versiones del código hay que **reutilizar los alineamientos**:
  `align_markers` y `trim_alignments` saltan las salidas existentes no vacías, así
  que basta con copiar `alignments/` al nuevo directorio de salida.
- `threads_per_gene: 1` da reproducibilidad a costa de tiempo de alineamiento.

`provenance/software_versions.tsv` registra las versiones de Python, numpy y
pandas, pero **no las de MAFFT ni trimAl**, que son las que determinan los
alineamientos. Es un hueco de procedencia pendiente.

---

## 11. Trampas

- **Cómputo y escritura están mezclados.** `calculate_metrics` calcula **y** escribe
  dos TSV antes de devolver; igual `create_panels` y `run_exploratory_pca`.
  Separarlos es el siguiente desacoplamiento natural y lo que abriría la puerta a
  ejecución reanudable. Aún sin hacer.
- **Escritura atómica duplicada.** `fasta.write_fasta` y `align.run_to_file`
  implementan cada una por su cuenta el mismo idioma `.tmp` + `.replace()`.
  Candidatas a unificar; hoy pueden divergir en silencio.
- `align.require_executable` falla con mensaje claro si `mafft` o `trimal` no están
  en el PATH.
- Estilo **vertical verboso**: un argumento por línea, comas finales. No arrastrar
  aquí el estilo compacto de la etapa 2.
- Las rutas absolutas de datos están en `configs/select.dikarya.yaml`; los datos
  pesados viven fuera del repo (`data/README.md`).
