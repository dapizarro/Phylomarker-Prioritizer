# PhyloMarker Select v0.3 — instalación

Este paquete sustituye los archivos de la versión v0.2.

## Cambios principales

- Scores de los perfiles transformados a percentiles comparables entre 0 y 1.
- `occupancy_only` usa únicamente cobertura.
- Penalización de redundancia predeterminada reducida de 0.25 a 0.10.
- El optimizador trabaja sobre candidatos de alta calidad:
  - top 10 %;
  - mínimo 50 candidatos;
  - al menos 5 veces el tamaño del panel;
  - caída máxima de score de 0.20 respecto al mejor candidato.
- La matriz de redundancia usa métricas primarias y elimina columnas constantes.
- `selection_trace.tsv` añade estabilidad al trimming y tamaño del pool candidato.
- Cada panel genera `panel_summary.tsv`.
- Se mantienen las políticas de exclusión por trimming.
- 11 tests automáticos.

## Instalación

Desde la raíz del repositorio:

```bash
source .venv/bin/activate

cp src/phylomarker_select/cli.py src/phylomarker_select/cli.py.pre_v0.3
cp tests/test_core.py tests/test_core.py.pre_v0.3

cp /RUTA/AL/PAQUETE/src/phylomarker_select/cli.py src/phylomarker_select/cli.py
cp /RUTA/AL/PAQUETE/tests/test_core.py tests/test_core.py
cp /RUTA/AL/PAQUETE/config.trimmed_v3.yaml .
```

## Pruebas

```bash
python -m py_compile src/phylomarker_select/cli.py
python -m pytest -v
```

Resultado esperado: `11 passed`.

## Ejecución

```bash
phylomarker-select --verbose run --config config.trimmed_v3.yaml
```

La salida se escribirá en `results_dicarya_trimmed_v3`.

## Nota científica

`clade_balance_score` seguirá siendo poco informativo cuando todos los órdenes tengan occupancy completo. Esta versión lo trata como componente neutro mediante percentiles constantes. Una futura versión podrá incorporar completitud por orden a nivel de residuo; eso requiere ampliar el cálculo de métricas por taxón y no se ha simulado artificialmente aquí.
