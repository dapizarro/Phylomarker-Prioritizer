# PhyloMarker Select — parche trimming-aware v0.2

## Contenido

- `src/phylomarker_select/cli.py`: archivo principal modificado.
- `tests/test_core.py`: pruebas ampliadas.
- `config.trimmed_v2.yaml`: configuración que escribe en `results_dicarya_trimmed_v2`.
- `phylomarker_select_v0.2.patch`: parche unificado alternativo.

## Instalación recomendada

Desde la raíz del repositorio original:

```bash
source .venv/bin/activate

cp src/phylomarker_select/cli.py src/phylomarker_select/cli.py.pre_v0.2
cp tests/test_core.py tests/test_core.py.pre_v0.2

cp /RUTA/AL/PARCHE/src/phylomarker_select/cli.py src/phylomarker_select/cli.py
cp /RUTA/AL/PARCHE/tests/test_core.py tests/test_core.py
cp /RUTA/AL/PARCHE/config.trimmed_v2.yaml .
```

Después:

```bash
python -m py_compile src/phylomarker_select/cli.py
python -m pytest -v
phylomarker-select --verbose run --config config.trimmed_v2.yaml
```

## Cambios principales

1. Calcula métricas raw y trimmed dentro de la misma ejecución.
2. Añade `retained_length_fraction`, `pis_retained_fraction`,
   `trimming_class` y `trimming_stability_score`.
3. Distingue `stable`, `signal_preserved`, `signal_sensitive` y `extreme`.
4. Aplica políticas de exclusión dependientes del perfil.
5. Evita que similitudes negativas bonifiquen candidatos.
6. Añade métricas de trimming al `selection_trace.tsv`.
7. Mantiene `random_matched` sin filtros de trimming.
