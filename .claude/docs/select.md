# Etapa 1 — `phylomarker-select` (referencia de implementación)

Referencia interna para **trabajar sobre el código**. La prosa orientada a usuario
(instalación, interpretación biológica, límites científicos) está en `docs/select.md`;
este documento no la repite.

Todo vive en un único módulo: `src/phylomarker_select/cli.py` (~3.500 líneas).
No hay ejecución parcial ni reanudable: `run_pipeline()` (línea 3047) recorre las
etapas de principio a fin en cada invocación.

---

## 1. Cadena de ejecución

`run_pipeline()` está dirigido íntegramente por el YAML. Orden real de llamadas:

| # | Función | Línea | Escribe en `<out>/` |
|---|---|---|---|
| 0 | `write_provenance` | 2995 | `provenance/resolved_config.yaml`, `provenance/software_versions.tsv` |
| 1 | `load_metadata` | 192 | — |
| 2 | `discover_busco_runs` | 262 | `validation/discovered_runs.tsv` |
| 3 | `validate_runs` | 333 | `validation/validated_runs.tsv`, `validation/warnings.tsv` |
| 4 | `extract_markers` | 451 | `sequences/<seq_type>/per_gene/`, `validation/<seq>_sequence_provenance.tsv` |
| 5 | `align_markers` (MAFFT) | 568 | `alignments/<seq_type>/untrimmed/` |
| 6 | `trim_alignments` (trimAl) | 628 | `alignments/<seq_type>/trimmed/` |
| 7 | `calculate_metrics` | 1023 | `metrics/per_gene_metrics.tsv`, `metrics/per_gene_<level>_occupancy.tsv` |
| 8 | `add_biological_scores` | 1424 | `rankings/all_gene_scores.tsv` |
| 9 | `run_exploratory_pca` | 1674 | `pca/pca_{scores,loadings,explained_variance}.tsv` |
| 10 | `create_panels` | 2569 | `rankings/<profile>.tsv`, `panels/<profile>/n<size>/` |
| 11 | `create_html_report` | 2809 | `report/index.html` |

**Punto de parada duro**: tras el paso 3, si `warnings` contiene alguna fila con
`severity == "error"`, `run_pipeline` lanza `RuntimeError` y no continúa
(cli.py:3172). Los avisos no-error se registran y se sigue.

**Qué alineamientos se analizan.** `trimming.enabled: true` (por defecto) hace que
`analysis_alignment_directory` sea `trimmed/`; con `false` se usa `untrimmed/`.
Los sin recortar **siempre** se conservan, porque `calculate_metrics` los recibe
aparte como `raw_alignment_directory` para calcular la estabilidad al recorte.

---

## 2. Métricas por gen (`calculate_metrics`, 1023)

Escribe una fila por gen en `metrics/per_gene_metrics.tsv`. Familias de columnas:

- **Cobertura** — `n_taxa`, `n_taxa_total`, `taxon_occupancy`,
  `sequence_completeness`, `cell_occupancy`.
- **Alineamiento** — `alignment_length`, `gap_fraction`, `ambiguous_fraction`.
- **Información** — `variable_sites`, `parsimony_informative_sites`,
  `variable_sites_per_length`, `pis_per_length`, `mean_entropy`
  (helpers: `site_entropy` 717, `variable_sites_and_pis` 739).
- **Tasa** — `mean_pairwise_distance` (`mean_pairwise_p_distance`, 928).
- **Sesgo** — `composition_variability` (`composition_variability`, 971).
- **Balance de clados** — agregados por `taxonomy.balance_level` (por defecto
  `order`): `mean_group_occupancy`, `min_replicated_group_occupancy`,
  `min_group_sequence_completeness`, `sd_group_sequence_completeness`,
  `worst_group_gap_fraction`, `n_groups_present`. El desglose por grupo va a
  `metrics/per_gene_<level>_occupancy.tsv`.
- **Estabilidad al recorte** — `raw_*` vs `trimmed_*`, `retained_length_fraction`,
  `pis_retained_fraction`, `variable_sites_retained_fraction`, más
  `trimming_class` y `trimming_stability_score`.

### Clasificación del recorte (`classify_trimming`, 790)

Cuatro clases, evaluadas en este orden:

| Clase | Condición |
|---|---|
| `extreme` | `retained_length_fraction < 0.25` **o** `pis_retained_fraction < 0.50` |
| `stable` | `retained_length_fraction ≥ 0.75` y `pis_retained_fraction ≥ 0.80` |
| `signal_preserved` | `pis_retained_fraction ≥ 0.80` (pero longitud entre 0,25 y 0,75) |
| `signal_sensitive` | resto |

Sin recorte se fuerza `trimming_class = "not_trimmed"` y
`trimming_stability_score = 1.0`.

`calculate_trimming_stability_score` (822) da **más peso a conservar PIS que
longitud**: `0.40 * length_score + 0.60 * information_score`, con penalizaciones
multiplicativas (`×0.25` si la longitud cae por debajo de 0,25; `×0.50` si los PIS
caen por debajo de 0,50).

---

## 3. Puntuación biológica (`add_biological_scores`, 1424)

Dos transformaciones, con papeles distintos:

- **`robust_standardize`** (1370) — `(x − mediana) / (1.4826 · MAD)`, recortado a
  `[−4, 4]`. Si MAD es 0 o NaN devuelve ceros. Se usa para **combinar métricas
  crudas** dentro de cada dimensión.
- **`percentile_score`** (1401) — rango percentil en `[0, 1]`. Se usa para
  **combinar dimensiones entre sí** en el ranking por perfil. Con ≤1 valor
  observado o ≤1 valor único devuelve 0,5 constante.

### Las seis dimensiones y sus pesos internos

| Dimensión | Composición |
|---|---|
| `coverage_score` | 0,45·`taxon_occupancy` + 0,35·`sequence_completeness` + 0,20·`cell_occupancy` |
| `clade_balance_score` | 0,50·`min_group_sequence_completeness` + 0,30·(−`sd_group_sequence_completeness`) + 0,20·(−`worst_group_gap_fraction`) |
| `alignment_score` | 0,45·`alignment_length` + 0,30·(−`gap_fraction`) + 0,25·(−`ambiguous_fraction`) |
| `information_score` | 0,45·`pis_per_length` + 0,25·`variable_sites_per_length` + 0,20·`mean_entropy` + 0,10·`parsimony_informative_sites` |
| `rate_score` | −\|`mean_pairwise_distance` − mediana\| estandarizado: premia la **cercanía a la tasa mediana**, no la tasa alta ni la baja |
| `bias_penalty` | 0,55·`composition_variability` + 0,30·`gap_fraction` + 0,15·`ambiguous_fraction` (mayor = peor) |

`clade_balance_score` tiene **fallback**: si faltan las tres columnas de grupo
finas, usa 0,60·`mean_group_occupancy` + 0,40·`min_replicated_group_occupancy`
(1466).

Después se derivan siete percentiles: `coverage_percentile`,
`clade_balance_percentile`, `alignment_percentile`, `information_percentile`,
`rate_percentile`, `low_bias_percentile` (= percentil **invertido** de
`bias_penalty`) y `trimming_percentile`.

### Elegibilidad — puerta dura previa al ranking

Se calcula aquí (1618), no en el optimizador. Seis condiciones en AND, todas bajo
la clave `eligibility` del YAML:

| Campo | Defecto en código | En `configs/select.dikarya.yaml` |
|---|---|---|
| `min_taxa` | 4 | 7 |
| `min_taxon_occupancy` | 0.5 | 0.50 |
| `min_alignment_length` | 80 | 100 |
| `min_informative_sites` | 2 | 3 |
| `max_gap_fraction` | 0.6 | 0.60 |
| `max_ambiguous_fraction` | 0.1 | 0.10 |

La columna `exclusion_reason` acumula, separados por `;`, los nombres de los
filtros incumplidos — es el sitio al que mirar cuando un gen esperado no aparece.

---

## 4. PCA (`run_exploratory_pca`, 1674)

**Exploratorio y nada más.** `pca.use_for_ranking: false` en el config, y PC1 no
entra en ninguna fórmula de puntuación. Escribe tres TSV en `pca/` y alimenta el
informe HTML. Al tocar esta función, no conectar su salida al ranking.

---

## 5. Ranking por perfil (`calculate_profile_ranking`, 1837)

`profile_score = Σ peso_feature · percentil_feature`, con los pesos de
`PROFILE_WEIGHTS` (cli.py:31). Cada perfil suma 1,0:

| Perfil | Pesos |
|---|---|
| `core_complete` | coverage 0,35 · alignment 0,25 · trimming 0,15 · information 0,15 · low_bias 0,10 |
| `backbone_balanced` | coverage 0,30 · clade_balance 0,25 · alignment 0,20 · trimming 0,15 · information 0,10 |
| `deep_robust` | alignment 0,25 · coverage 0,20 · trimming 0,20 · rate 0,20 · low_bias 0,15 |
| `low_bias` | low_bias 0,35 · alignment 0,25 · trimming 0,20 · coverage 0,10 · information 0,10 |
| `diverse_rate` | coverage 0,25 · alignment 0,25 · trimming 0,20 · information 0,20 · low_bias 0,10 |
| `occupancy_only` | coverage 1,00 — **control negativo** |
| `information_only` | information 1,00 — **control negativo** |

`random_matched` no está en `PROFILE_WEIGHTS`: se atiende por una rama aparte en
`create_panels` (2649) vía `random_panel` (2482), sembrada con
`project.random_seed + size`.

**Segundo filtro, por perfil**: `PROFILE_EXCLUDED_TRIMMING_CLASSES` (cli.py:76)
recorta además por `trimming_class`. Todos los perfiles reales excluyen `extreme`;
`deep_robust` excluye también `signal_sensitive`. Los controles
`information_only` y `random_matched` **no excluyen nada** — es deliberado, son
controles.

Los genes no elegibles reciben `profile_score = −inf` (1902) en lugar de ser
eliminados, para que el TSV de ranking siga conteniendo todos los genes.
Ordenación determinista por `(profile_score desc, gene_id asc)`.

---

## 6. Optimización de paneles

Los paneles **no son top-N**. Esto es lo que hay que preservar.

### `optimize_panel_greedily` (1966)

1. **Filtrar** por `profile_eligible`.
2. **Pool de candidatos**:
   `pool_size = max(candidate_minimum_pool_size, ceil(n · candidate_top_fraction), panel_size · 5)`,
   acotado a `len(candidates)`.
3. **Corte por calidad**: descartar los que caen más de `maximum_score_drop` por
   debajo del mejor `profile_score`. Si eso deja menos genes que `panel_size`, se
   revierte al pool sin recortar (2005) — la red de seguridad que evita paneles
   incompletos.
4. **Selección greedy**, `panel_size` pasos, maximizando

   ```
   Δ_g(P) = profile_score(g) − redundancy_penalty · max_{h∈P} similarity(g, h)
   ```

   `similarity` es coseno sobre `feature_vectors` (1926): `taxon_occupancy`,
   `min_replicated_group_occupancy`, `alignment_length`, `pis_per_length`,
   `mean_pairwise_distance`, `composition_variability`,
   `retained_length_fraction` — imputadas con mediana (`SimpleImputer`) y
   escaladas con `RobustScaler`. Se descartan las features con un solo valor único.

Devuelve `(panel, trace)`; la traza es `selection_trace.tsv` y registra por paso
`profile_score`, `maximum_redundancy`, `marginal_gain` y `candidate_pool_size`.

### `optimize_diverse_rate_panel` (2111)

Mismo pool y mismo penalizador, pero reparte las plazas entre **bins de tasa**:

- `pd.qcut` sobre `mean_pairwise_distance` en `panels.diverse_rate_bins` cuantiles
  (5 por defecto), `duplicates="drop"`.
- Cuotas: `base_quota = size // n_bins`, y el resto se asigna a los primeros bins.
- Recorrido round-robin; dentro de cada bin se elige por la misma ganancia
  marginal. Si una pasada completa no progresa, se rellena desde el conjunto
  restante sin restricción de bin (2327).
- **Degradación**: con ≤1 bin efectivo o `mean_pairwise_distance` constante,
  delega en `optimize_panel_greedily` y marca `rate_bin = "single_bin"` (2174).

Diferencia menor entre ambos: `diverse_rate` desempata explícitamente por
`gene_id` (2268); el greedy estándar depende del orden de iteración, que ya es
determinista por la ordenación previa.

---

## 7. Salidas — contrato con la etapa 2

```
<inputs.output>/
├── validation/     discovered_runs.tsv, validated_runs.tsv, warnings.tsv
├── sequences/<seq_type>/per_gene/
├── alignments/<seq_type>/untrimmed/     <gene>.<faa|fna>
├── alignments/<seq_type>/trimmed/       <gene>.trimmed.<faa|fna>   ← la etapa 2 lee aquí
├── metrics/        per_gene_metrics.tsv, per_gene_<level>_occupancy.tsv
├── rankings/       all_gene_scores.tsv, <profile>.tsv
├── pca/            pca_scores.tsv, pca_loadings.tsv, pca_explained_variance.tsv
├── panels/<profile>/n<size>/
│   ├── genes.txt              ← la etapa 2 lee aquí
│   ├── panel_genes.tsv
│   ├── panel_summary.tsv
│   ├── panel.yaml
│   ├── selection_trace.tsv
│   └── alignments/
├── report/index.html
└── provenance/     resolved_config.yaml, software_versions.tsv
```

**Dos rutas son contrato con la etapa 2** y no se pueden cambiar sin romperla en
silencio (ver `.claude/docs/phylogeny.md`, sección 2):

- `panels/<profile>/n<size>/genes.txt` — leída por `core.panel_genes()`
- `alignments/<seq_type>/trimmed/<gene>.trimmed.<ext>` — leída por
  `core.alignment_path()`

`genes.txt` es una lista plana sin cabecera. `panel_genes()` verifica que tenga
**exactamente** `n<size>` líneas y sin duplicados: un panel más corto de lo
nominal aborta la etapa 2.

---

## 8. Invariantes que los tests bloquean

`tests/select/test_core.py` (300 líneas, 13 tests). Al modificar el módulo,
comprobar que se mantienen:

1. Los paneles no son top-N — el greedy con penalización de redundancia debe
   poder saltarse un gen mejor puntuado.
2. La elegibilidad se aplica **antes** del ranking.
3. PC1 nunca es una puntuación de calidad.
4. `diverse_rate` cubre los bins de forma pareja.
5. `occupancy_only`, `information_only` y `random_matched` son controles
   negativos: no "mejorarlos".
6. Determinismo: mismo config y mismos datos ⇒ mismo panel.

---

## 9. Estilo y trampas

- Estilo **vertical verboso**: un argumento por línea, comas finales, mucho aire.
  Coincide con el fichero; no reformatear hacia el estilo compacto de la etapa 2.
- Comprobación rápida habitual: `python -m py_compile src/phylomarker_select/cli.py`.
- Herramientas externas por `require_executable` (536): `mafft`, `trimal`. Fallan
  con mensaje claro si no están en el PATH.
- `run_pipeline` no acepta reanudar. Una corrida interrumpida se repite entera.
- Las rutas absolutas de datos están en `configs/select.dikarya.yaml`; los datos
  pesados viven fuera del repo (`data/README.md`).
