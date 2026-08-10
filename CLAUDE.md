# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A two-stage phylogenomic pipeline shipped as **one installable package with two console scripts**:

1. `phylomarker-select` (stage 1, internally v0.4) — reads BUSCO single-copy outputs, computes per-gene metrics, builds marker panels optimized for explicit evolutionary objectives.
2. `phylomarker-phylogeny` (stage 2, internally v0.3) — validates those panels with IQ-TREE gene trees, partitioned concatenation, ASTRAL, gCF/sCF, jackknife and RF comparison.

Prose docs are in Spanish; code and identifiers in English. Keep that split.

## Commands

```bash
mamba env create -f environment.yml && mamba activate phylomarker
python -m pip install -e ".[test]"

pytest                                              # 27 tests (13 select + 14 phylogeny)
pytest tests/select -q
pytest tests/phylogeny/test_v03.py::test_expected_astral_tree -v
python -m py_compile src/phylomarker_select/cli.py  # habitual smoke check on the monolith

phylomarker-select --verbose run --config configs/select.dikarya.yaml
phylomarker-phylogeny local-run --config configs/phylogeny.local.n10.yaml
```

`tests/select/test_core.py` and `tests/phylogeny/test_core.py` share a basename; `pyproject.toml` sets `addopts = "--import-mode=importlib"` so they collect without a module-name clash. Don't remove it.

External tools: `mafft`, `trimal` (stage 1); `iqtree2`/`iqtree3`, `java`, an ASTRAL jar (stage 2).

## Architecture

Deep implementation references — read the relevant one before non-trivial work on a stage:

- **`.claude/docs/select.md`** — stage 1: the module map, the `run_pipeline` chain, every
  metric family and its formula, the six biological dimensions and their internal weights,
  the `PROFILES` registry, both greedy optimizers, the eligibility gate, the output tree,
  and the MAFFT reproducibility limit.
- **`.claude/docs/phylogeny.md`** — stage 2: the manifest-first contract, resumability via
  `command_expected_output`, the three partition strategies and model reuse, concatenation
  rules, DendroPy RF, deterministic jackknife, and the output tree.

The sections below are the summary; those two files are the detail.

### The coupling between stages

Stage 2's config field `inputs.selection_results` points at stage 1's `inputs.output` directory. `core.panel_genes()` and `core.alignment_path()` read that tree directly — specifically `panels/<profile>/n<size>/genes.txt` and the per-gene alignments. Changing stage 1's output layout breaks stage 2 silently.

### `src/phylomarker_select/`

21 modules (was one ~3.500-line `cli.py` until August 2026). `run_pipeline()` in `pipeline.py` is driven end to end by YAML; there is still **no partial or resumable execution**. Stage order:

`discover_busco_runs` → `validate_runs` → `extract_markers` → `align_markers` (MAFFT) → `trim_alignments` (trimAl; both untrimmed and trimmed retained) → `calculate_metrics` (occupancy, PIS/variable sites, entropy, p-distance, gap/ambiguity, composition variability, trimming stability) → eligibility filter → `add_biological_scores` (median/MAD robust standardization over coverage, clade balance, alignment quality, information, rate fit, bias penalty) → `run_exploratory_pca` → `calculate_profile_ranking` → `optimize_panel_greedily` / `optimize_diverse_rate_panel` → `create_panels` → `create_html_report` → `write_provenance`.

**Three modules own what used to be scattered. Extend there, not elsewhere:**
- `config.py` — every YAML default, once, in typed frozen dataclasses. `SelectConfig.load()`. Keeps `.raw` for the provenance/HTML dump.
- `layout.py` — every output path. `panel_genes_file()` and `trimmed_alignment()` are the stage-2 contract.
- `profiles.py` — the `PROFILES` registry. Adding a profile is one entry. **Do not reorder a profile's `weights`**: ranking sums them in insertion order, so reordering perturbs float results and can flip ranking ties.

If you find yourself writing a path literal or a `get(key, default)` outside those three, the change belongs somewhere else.

Invariants encoded in the code and locked by tests — preserve them:
- Panels are **not** top-N. Greedy marginal gain `Δ_g(P) = S_g − λ·max_{h∈P} similarity(g,h)`, restricted to a high-quality candidate pool (`candidate_top_fraction`, `candidate_minimum_pool_size`, `maximum_score_drop`).
- PCA is exploratory only; `pca.use_for_ranking: false`, PC1 is never a quality score.
- Eligibility is a hard gate applied *before* ranking.
- `diverse_rate` bins by `mean_pairwise_distance` quantiles (`diverse_rate_bins: 5`) and fills bins evenly.
- `occupancy_only`, `information_only`, `random_matched` are deliberate negative controls, not profiles to "improve".

**Reproducibility is load-bearing and was hard-won.** Two independent nondeterminisms were fixed in August 2026; don't reintroduce either:
- `configs/select.dikarya.yaml` pins `threads_per_gene: 1` because MAFFT `--auto` is only deterministic single-threaded (25/25 reproducible at 1 thread, 12/25 at 2). Before the fix, 378/1311 genes aligned differently across runs of identical code, changing 12 of 16 panels.
- `calculate_metrics` iterates `sorted(members)`, not the raw `set` — Python randomizes string hashes per process, so `np.mean` summed in a different order each run and the five group columns drifted by 1 ULP (462 cells across two runs of the same code).

To compare two versions of the code cheaply, reuse the alignments: `align_markers`/`trim_alignments` skip existing non-empty outputs, so copying `alignments/` into the new output directory skips MAFFT entirely. See `.claude/docs/select.md` §10 bis.

### `src/phylomarker_phylogeny/`

`cli.py` (one thin function per subcommand) + `core.py` (FASTA/concatenation, NEXUS `charset` partitions, IQ-TREE/ASTRAL command construction, manifest + SLURM array emission, DendroPy RF, deterministic jackknife subsets).

- **Manifest-first.** `gene-trees`, `concatenated`, `astral`, `concordance`, `jackknife` only write command manifests and SLURM arrays; execution needs `--execute`. `local-run` forces it on. This is what makes one config work on laptop and cluster.
- **Resumability** is by expected-output existence: `command_expected_output()` maps a command to its `.treefile` / `.cf.tree` / ASTRAL output; non-empty results are skipped.
- **Model reuse.** With `concat_partition_strategy: gene_models` (recommended), ModelFinder runs once in `gene-trees`; `gene-models` scrapes each `.iqtree` for the best BIC model. Alternatives: `model_finder_merge` (`MFP+MERGE`), `fixed`. `missing_model_policy: error` fails loudly instead of silently using `fallback_model`.
- RF uses DendroPy, not IQ-TREE, to avoid a bug in IQ-TREE 2.0.7.
- `gene_models` results go to parallel directories (`concatenated_gene_models/`, `concordance_gene_models/`) so earlier runs are never overwritten.

## Data lives outside the repository

Only `data/samples.tsv` is versioned. BUSCO inputs (1,2 GB) and run outputs stay under `~/Documentos/PROJECT_2026/Phylomarker_prioritizer/`; absolute paths are in the two files under `configs/`. See `data/README.md`.

Current Dikarya run state: stage 1 complete (`results_dicarya_trimmed_v4`); stage 2 steps 1–8 complete, jackknife interrupted at 31/120 replicates — `local-run` resumes it.

## Code style

`phylomarker_select` uses a verbose vertical style (one argument per line, trailing commas, wide spacing). `phylomarker_phylogeny` uses a terse compact style (semicolons, single-line functions). Match the file you are editing; do not reformat across the boundary.

## History

This repo was consolidated from a workspace that versioned by zip archives and `.patch` files rather than commits. The old `INSTRUCCIONES_V0.N.md` files are preserved as `docs/changelog_select_v0.N.md`; ad-hoc comparison TSVs are in `docs/analysis/`. Everything else (zips, patches, `cli.py.pre_v0.*` backups, superseded phylogeny v0.1/v0.2, and a `build_phylomarker_select.py` generator that would `rmtree` the source tree) was left behind in the original directory. Use git from now on, not patch archives.

## Scientific caution

The docs are explicit that this software does not demonstrate a correct species tree, absence of hidden paralogy or introgression, and does not distinguish ILS from estimation error; low gCF has many possible causes. Keep that framing in user-facing text, reports and summaries — do not upgrade descriptive metrics into quality guarantees.
