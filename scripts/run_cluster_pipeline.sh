#!/usr/bin/env bash
set -euo pipefail
CONFIG="${1:-examples/phylogeny.cluster.yaml}"

phylomarker-phylogeny validate --config "$CONFIG"
phylomarker-phylogeny prepare --config "$CONFIG"
phylomarker-phylogeny gene-trees --config "$CONFIG"
phylomarker-phylogeny concatenated --config "$CONFIG"
phylomarker-phylogeny jackknife --config "$CONFIG"

OUT=$(python - "$CONFIG" <<'PY'
from pathlib import Path
import sys,yaml
p=Path(sys.argv[1]).resolve()
c=yaml.safe_load(p.read_text())
o=Path(c["inputs"]["output"])
print(o if o.is_absolute() else (p.parent/o).resolve())
PY
)

GENE_JOB=$(sbatch --parsable "$OUT/slurm/gene_trees.sbatch")
CONCAT_JOB=$(sbatch --parsable "$OUT/slurm/concatenated.sbatch")
JK_JOB=$(sbatch --parsable "$OUT/slurm/jackknife.sbatch")
echo "gene trees: $GENE_JOB"
echo "concatenated: $CONCAT_JOB"
echo "jackknife: $JK_JOB"

cat > "$OUT/submit_after_gene_trees.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
phylomarker-phylogeny astral --config "$CONFIG"
ASTRAL_JOB=\$(sbatch --parsable --dependency=afterok:$GENE_JOB "$OUT/slurm/astral.sbatch")
echo "ASTRAL: \$ASTRAL_JOB"
phylomarker-phylogeny concordance --config "$CONFIG"
CF_JOB=\$(sbatch --parsable --dependency=afterok:$GENE_JOB:$CONCAT_JOB "$OUT/slurm/concordance.sbatch")
echo "Concordance: \$CF_JOB"
EOF
chmod +x "$OUT/submit_after_gene_trees.sh"

echo "When gene-tree files exist, run: $OUT/submit_after_gene_trees.sh"
echo "Finally run: phylomarker-phylogeny summarize --config $CONFIG"
