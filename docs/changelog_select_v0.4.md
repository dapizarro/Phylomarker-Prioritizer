# PhyloMarker Select v0.4

## Cambios principales

- `diverse_rate` usa cinco cuantiles de `mean_pairwise_distance`.
- Para paneles de 10 genes intenta seleccionar 2 genes por cuantil.
- Para paneles de 5 genes intenta seleccionar 1 gen por cuantil.
- Dentro de cada cuantil prioriza cobertura, alineamiento, estabilidad al trimming,
  información y bajo sesgo.
- `backbone_balanced` incorpora completitud por orden, dispersión entre órdenes y
  peor fracción de gaps por orden.
- Se añaden métricas:
  - `mean_group_sequence_completeness`
  - `min_group_sequence_completeness`
  - `min_replicated_group_sequence_completeness`
  - `sd_group_sequence_completeness`
  - `worst_group_gap_fraction`
- `panel_summary.tsv` incluye rango y desviación de tasas y diagnósticos por grupo.
- Se mantienen `redundancy_penalty: 0.10` y el pool de candidatos de alta calidad.
- 13 tests pasan.

## Instalación

Desde la raíz del repositorio:

```bash
source .venv/bin/activate

cp src/phylomarker_select/cli.py src/phylomarker_select/cli.py.pre_v0.4
cp tests/test_core.py tests/test_core.py.pre_v0.4

unzip phylomarker_select_v0.4_patch.zip
```

Cuando `unzip` pregunte, selecciona `All`.

## Validación

```bash
python -m py_compile src/phylomarker_select/cli.py
python -m pytest -v
```

Resultado esperado:

```text
13 passed
```

## Ejecución

```bash
phylomarker-select --verbose run --config config.trimmed_v4.yaml
```

La salida se escribirá en `results_dicarya_trimmed_v4`.

## Comprobaciones

```bash
python - <<'PY'
from pathlib import Path
import pandas as pd

root = Path("results_dicarya_trimmed_v4")

for profile in ["diverse_rate", "backbone_balanced"]:
    panel = pd.read_csv(
        root / "panels" / profile / "n10" / "panel_genes.tsv",
        sep="\t",
    )
    print(f"\n{profile}")
    columns = [
        "gene_id",
        "mean_pairwise_distance",
        "rate_bin",
        "min_group_sequence_completeness",
        "sd_group_sequence_completeness",
        "worst_group_gap_fraction",
        "profile_score",
    ]
    print(panel[[c for c in columns if c in panel.columns]].to_string(index=False))

summary = pd.read_csv(
    root / "panels" / "diverse_rate" / "n10" / "panel_summary.tsv",
    sep="\t",
)
print("\nResumen diverse_rate:")
print(summary.to_string(index=False))
PY
```

En `diverse_rate/n10`, `rate_bin` debería contener cinco clases con dos genes cada una,
salvo que un cuantil no disponga de suficientes candidatos elegibles.
