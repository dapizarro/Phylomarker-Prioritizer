# data/

Este directorio sólo contiene metadatos ligeros y versionados.

- `samples.tsv` — metadatos taxonómicos de las 14 muestras Dikarya.
  `sample_ID` debe coincidir con el sufijo de cada carpeta `run_<sample>` de BUSCO.

## Los datos pesados NO están en el repositorio

| Qué | Dónde | Tamaño |
|---|---|---|
| Salidas BUSCO (`run_<sample>/`) | `~/Documentos/PROJECT_2026/Phylomarker_prioritizer/dicarya_BUSCO/` | 1,2 GB |
| Corrida canónica de la etapa 1 | `~/Documentos/PROJECT_2026/Phylomarker_prioritizer/phylomarker-select/results_dicarya_trimmed_v4/` | 63 MB |
| Corrida local de la etapa 2 | `~/Documentos/PROJECT_2026/Phylomarker_prioritizer/phylomarker_phylogeny_v0.3/phylogenetic_results_local_n10/` | 21 MB |

Las rutas absolutas están en `configs/select.dikarya.yaml` y
`configs/phylogeny.local.n10.yaml`. Si mueves los datos, edita esos dos ficheros.

Material histórico (versiones antiguas, parches, zips, comparativas ad-hoc)
en `~/Documentos/PROJECT_2026/_archive_phylomarker/`.
