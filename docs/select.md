# PhyloMarker Select

PhyloMarker Select caracteriza marcadores BUSCO y construye paneles
filogenómicos adaptados a objetivos evolutivos explícitos.

## Principio científico

Un marcador no es universalmente bueno. Su utilidad depende de:

- escala evolutiva;
- muestreo taxonómico;
- representación de clados;
- occupancy;
- completitud de las secuencias;
- calidad del alineamiento;
- tasa evolutiva;
- información filogenética;
- sesgos composicionales;
- redundancia con otros genes del panel.

El PCA se utiliza como análisis exploratorio. PC1 no se interpreta como una
puntuación universal de calidad y no se utiliza directamente para ordenar los
genes.

## Estructuras BUSCO admitidas

La herramienta busca recursivamente estructuras como:

~~~
run_SAMPLE/
└── single_copy_busco_sequences/
    ├── EOG....faa
    └── EOG....fna
~~~

También admite:

~~~
run_SAMPLE/
└── busco_sequences/
    └── single_copy_busco_sequences/
        ├── BUSCO_ID.faa
        └── BUSCO_ID.fna
~~~

Se utilizan directamente los archivos `.faa` o `.fna`. Las predicciones crudas
de Augustus no se interpretan cuando existen los FASTA de BUSCO single-copy.

## Instalación

Con Mamba:

~~~bash
mamba env create -f environment.yml
mamba activate phylomarker-select
~~~

Con pip:

~~~bash
python -m pip install .
~~~

MAFFT y trimAl deben estar disponibles en `PATH`.

## Metadatos

El archivo debe ser TSV e incluir, como mínimo:

~~~text
sample_ID	species	order	class
~~~

`sample_ID` debe coincidir con el sufijo de la carpeta `run_SAMPLE`.

Las columnas taxonómicas adicionales, como `genus`, `family` y `subphylum`,
permiten definir análisis a distintas escalas.

## Validar los datos

~~~bash
phylomarker-select validate \
    --busco-directory /ruta/dicarya_BUSCO \
    --metadata /ruta/samples.tsv \
    --output validation_dicarya
~~~

La validación informa de:

- outputs BUSCO sin metadatos;
- filas de metadatos sin output BUSCO;
- carpetas `single_copy_busco_sequences` ausentes;
- archivos `.faa` o `.fna` ausentes;
- identificadores duplicados.

## Ejecutar el flujo completo

Copiar la configuración:

~~~bash
cp examples/config.dikarya.yaml config.yaml
~~~

Editar las rutas:

~~~yaml
inputs:
  busco_directory: /ruta/dicarya_BUSCO
  metadata: /ruta/samples.tsv
  output: results_dicarya
~~~

Ejecutar:

~~~bash
phylomarker-select run --config config.yaml
~~~

## Resultados

~~~
results_dicarya/
├── provenance/
├── validation/
├── sequences/
│   └── protein/per_gene/
├── alignments/
│   └── protein/
│       ├── untrimmed/
│       └── trimmed/
├── metrics/
├── pca/
├── rankings/
├── panels/
└── report/
    └── index.html
~~~

## Métricas de occupancy

Occupancy taxonómico:

\[
O_{\mathrm{taxon},g}
=
\frac{n_{\mathrm{taxones\ presentes},g}}
{N_{\mathrm{taxones}}}
\]

Completitud de las secuencias presentes:

\[
O_{\mathrm{sequence},g}
=
1-
\frac{\mathrm{celdas\ ausentes}}
{n_{\mathrm{taxones\ presentes},g}L_g}
\]

Occupancy celular:

\[
O_{\mathrm{cell},g}
=
O_{\mathrm{taxon},g}
O_{\mathrm{sequence},g}
\]

También se calcula occupancy por grupo taxonómico configurable, por ejemplo:

~~~yaml
taxonomy:
  balance_level: order
~~~

Los grupos representados por una sola muestra se identifican como singleton.
Su presencia es informativa, pero no permite estimar de forma robusta el
occupancy dentro del grupo.

## Información filogenética observada

Se calculan:

- sitios variables;
- sitios parsimoniosamente informativos;
- PIS por longitud;
- sitios variables por longitud;
- entropía media;
- distancia p media entre pares.

Estas métricas describen el alineamiento observado. No garantizan que el árbol
génico sea correcto.

## Diagnósticos de sesgo

La versión inicial incluye:

- proporción de gaps;
- caracteres ambiguos;
- variabilidad composicional entre taxones;
- desviación de la tasa respecto al conjunto de genes.

Son diagnósticos preliminares, no pruebas completas de adecuación del modelo.

## Elegibilidad

Antes del ranking se aplican restricciones duras:

~~~yaml
eligibility:
  min_taxa: 7
  min_taxon_occupancy: 0.50
  min_alignment_length: 100
  min_informative_sites: 3
  max_gap_fraction: 0.60
  max_ambiguous_fraction: 0.10
~~~

Una puntuación alta de información no compensa el incumplimiento de una
restricción de elegibilidad.

## PCA

El PCA emplea:

- imputación por mediana;
- escalado robusto;
- orientación determinista de componentes.

Genera:

~~~
pca_scores.tsv
pca_loadings.tsv
pca_explained_variance.tsv
~~~

El PCA es descriptivo y no determina el ranking biológico.

## Dimensiones del ranking

Cada gen recibe puntuaciones separadas para:

- cobertura;
- equilibrio entre clados;
- calidad del alineamiento;
- información;
- adecuación de la tasa;
- penalización por sesgo.

La estandarización utiliza mediana y desviación absoluta mediana.

## Paneles

### core_complete

Prioriza occupancy, completitud y calidad del alineamiento.

### backbone_balanced

Prioriza representación equilibrada entre grupos taxonómicos.

### deep_robust

Favorece cobertura amplia, tasas moderadas y menor sesgo.

### low_bias

Penaliza fuertemente heterogeneidad composicional, gaps y caracteres ambiguos.

### diverse_rate

Busca genes informativos sin concentrar el panel en un único régimen de tasas.

### Controles

- occupancy_only
- information_only
- random_matched

Estos controles permiten comprobar posteriormente si el enfoque multicriterio
supera a estrategias simples o aleatorias.

## Optimización conjunta

Los paneles no se generan tomando simplemente los primeros N genes.

El optimizador selecciona genes mediante ganancia marginal:

\[
\Delta_g(P)
=
S_g
-
\lambda
\max_{h\in P}
\operatorname{similitud}(g,h)
\]

La similitud considera cobertura, completitud, gaps, información, tasa y
composición. Así se penalizan paneles formados por genes casi intercambiables.

## Limitaciones

Este módulo no demuestra que un panel:

- produzca el árbol de especies correcto;
- esté libre de paralogía oculta;
- no esté afectado por introgresión;
- distinga ILS de error de estimación;
- sea óptimo fuera de los taxones analizados;
- funcione igual en todas las escalas evolutivas.

Estas cuestiones se evaluarán en el módulo filogenético independiente mediante
IQ-TREE, árboles génicos, concatenación, ASTRAL, concordancia y análisis de
conflicto.

## Pruebas

~~~bash
pytest
~~~

## Estado

La versión 0.1.0 es un prototipo funcional para validación iterativa. Antes de
uso definitivo debe probarse con:

- outputs BUSCO antiguos y modernos;
- cientos de genomas;
- missingness dependiente del clado;
- distintas estrategias de trimming;
- validación filogenética independiente.
