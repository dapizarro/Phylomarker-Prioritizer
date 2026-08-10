from pathlib import Path

import numpy as np
import pandas as pd

from phylomarker_select.fasta import (
    FastaRecord,
    first_fasta_record,
    read_fasta,
    write_fasta,
)
from phylomarker_select.metrics.sites import variable_sites_and_pis
from phylomarker_select.metrics.trimming import (
    calculate_trimming_metrics,
    calculate_trimming_stability_score,
    classify_trimming,
)
from phylomarker_select.optimize import (
    optimize_diverse_rate_panel,
    optimize_panel_greedily,
)
from phylomarker_select.profiles import calculate_profile_ranking
from phylomarker_select.scoring import add_biological_scores, percentile_score


def test_fasta_roundtrip(tmp_path: Path) -> None:
    path = tmp_path / "test.faa"

    records = [
        FastaRecord("sample_a", "ACDEFG"),
        FastaRecord("sample_b", "ACD-FG"),
    ]

    write_fasta(records, path)

    observed = list(read_fasta(path))

    assert observed == records
    assert first_fasta_record(path).identifier == "sample_a"


def test_variable_sites_and_pis() -> None:
    matrix = np.asarray(
        [
            list("AACG"),
            list("AACG"),
            list("AGTG"),
            list("AGTG"),
        ],
        dtype="U1",
    )

    variable, pis, entropy = variable_sites_and_pis(matrix)

    assert variable == 2
    assert pis == 2
    assert entropy > 0


def test_eligibility() -> None:
    metrics = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "n_taxa": [10, 8, 2],
            "taxon_occupancy": [1.0, 0.8, 0.2],
            "sequence_completeness": [0.98, 0.90, 0.40],
            "cell_occupancy": [0.98, 0.72, 0.08],
            "mean_group_occupancy": [1.0, 0.8, 0.2],
            "min_replicated_group_occupancy": [1.0, 0.6, 0.0],
            "alignment_length": [500, 400, 40],
            "gap_fraction": [0.01, 0.10, 0.70],
            "ambiguous_fraction": [0.0, 0.01, 0.20],
            "pis_per_length": [0.10, 0.15, 0.0],
            "variable_sites_per_length": [0.20, 0.25, 0.0],
            "mean_entropy": [0.40, 0.50, 0.0],
            "parsimony_informative_sites": [50, 60, 0],
            "mean_pairwise_distance": [0.10, 0.15, 0.01],
            "composition_variability": [0.01, 0.03, 0.20],
        }
    )

    config = {
        "eligibility": {
            "min_taxa": 4,
            "min_taxon_occupancy": 0.5,
            "min_alignment_length": 80,
            "min_informative_sites": 2,
            "max_gap_fraction": 0.6,
            "max_ambiguous_fraction": 0.1,
        }
    }

    scored = add_biological_scores(
        metrics,
        config,
    )

    assert bool(scored.loc[0, "eligible"])
    assert bool(scored.loc[1, "eligible"])
    assert not bool(scored.loc[2, "eligible"])



def test_classify_trimming_classes() -> None:
    assert classify_trimming(0.90, 0.95) == "stable"
    assert classify_trimming(0.60, 0.95) == "signal_preserved"
    assert classify_trimming(0.60, 0.65) == "signal_sensitive"
    assert classify_trimming(0.20, 0.95) == "extreme"
    assert classify_trimming(0.80, 0.40) == "extreme"


def test_classify_trimming_without_raw_pis() -> None:
    assert classify_trimming(0.90, None) == "stable"
    assert classify_trimming(0.60, None) == "signal_preserved"


def test_trimming_stability_score_is_bounded() -> None:
    cases = [
        (1.00, 1.00),
        (0.60, 0.90),
        (0.30, 0.60),
        (0.10, 0.95),
        (0.80, None),
    ]

    for length_fraction, pis_fraction in cases:
        score = calculate_trimming_stability_score(
            length_fraction,
            pis_fraction,
        )
        assert 0.0 <= score <= 1.0


def test_extreme_trimming_is_penalized() -> None:
    stable_score = calculate_trimming_stability_score(0.90, 0.95)
    extreme_score = calculate_trimming_stability_score(0.15, 0.95)
    assert extreme_score < stable_score


def test_calculate_trimming_metrics_for_column_subset(
    tmp_path: Path,
) -> None:
    raw_path = tmp_path / "gene.aln.faa"
    trimmed_path = tmp_path / "gene.trimmed.faa"

    raw_records = [
        FastaRecord("a", "AACGTT"),
        FastaRecord("b", "AACGTT"),
        FastaRecord("c", "AGTGTT"),
        FastaRecord("d", "AGTGTT"),
    ]
    trimmed_records = [
        FastaRecord("a", "AACG"),
        FastaRecord("b", "AACG"),
        FastaRecord("c", "AGTG"),
        FastaRecord("d", "AGTG"),
    ]

    write_fasta(raw_records, raw_path)
    write_fasta(trimmed_records, trimmed_path)

    observed = calculate_trimming_metrics(
        raw_path,
        trimmed_path,
    )

    assert observed["raw_alignment_length"] == 6
    assert observed["trimmed_alignment_length"] == 4
    assert observed["trimmed_pis"] <= observed["raw_pis"]
    assert observed["trimming_class"] == "signal_preserved"



def test_deep_robust_excludes_signal_sensitive_gene() -> None:
    frame = pd.DataFrame(
        {
            "gene_id": ["stable_gene", "sensitive_gene"],
            "eligible": [True, True],
            "trimming_class": ["stable", "signal_sensitive"],
            "coverage_score": [0.0, 10.0],
            "clade_balance_score": [0.0, 10.0],
            "alignment_score": [0.0, 10.0],
            "trimming_stability_score": [1.0, 0.7],
            "information_score": [0.0, 10.0],
            "rate_score": [0.0, 10.0],
            "bias_penalty": [0.0, 0.0],
        }
    )

    ranked = calculate_profile_ranking(frame, "deep_robust")
    sensitive = ranked.loc[
        ranked["gene_id"] == "sensitive_gene"
    ].iloc[0]

    assert not bool(sensitive["profile_eligible"])
    assert np.isneginf(sensitive["profile_score"])


def test_percentile_score_is_bounded_and_ordered() -> None:
    values = pd.Series([10.0, 20.0, 30.0])
    observed = percentile_score(values)
    assert observed.between(0.0, 1.0).all()
    assert observed.iloc[0] < observed.iloc[1] < observed.iloc[2]


def test_candidate_pool_prevents_low_rank_replacement() -> None:
    frame = pd.DataFrame(
        {
            "gene_id": [f"g{i:03d}" for i in range(100)],
            "profile_score": np.linspace(1.0, 0.0, 100),
            "profile_eligible": [True] * 100,
            "eligible": [True] * 100,
            "taxon_occupancy": [1.0] * 100,
            "alignment_length": np.arange(100, 200),
            "pis_per_length": np.linspace(0.1, 0.5, 100),
            "mean_pairwise_distance": np.linspace(0.01, 0.5, 100),
            "composition_variability": np.linspace(0.001, 0.02, 100),
            "retained_length_fraction": np.linspace(0.8, 1.0, 100),
            "trimming_class": ["stable"] * 100,
            "trimming_stability_score": [1.0] * 100,
        }
    )

    panel, trace = optimize_panel_greedily(
        frame,
        panel_size=5,
        redundancy_penalty=0.10,
        candidate_top_fraction=0.10,
        candidate_minimum_pool_size=10,
        maximum_score_drop=0.20,
    )

    assert len(panel) == 5
    assert panel["profile_score"].min() >= 0.80
    assert int(trace["candidate_pool_size"].max()) <= 20



def test_diverse_rate_panel_covers_rate_bins() -> None:
    rates = np.linspace(0.01, 0.50, 100)
    frame = pd.DataFrame(
        {
            "gene_id": [f"g{i:03d}" for i in range(100)],
            "profile_score": np.linspace(1.0, 0.80, 100),
            "profile_eligible": [True] * 100,
            "eligible": [True] * 100,
            "taxon_occupancy": [1.0] * 100,
            "min_replicated_group_occupancy": [1.0] * 100,
            "alignment_length": np.arange(500, 600),
            "pis_per_length": np.linspace(0.1, 0.5, 100),
            "mean_pairwise_distance": rates,
            "composition_variability": np.linspace(0.001, 0.02, 100),
            "retained_length_fraction": np.linspace(0.8, 1.0, 100),
            "trimming_class": ["stable"] * 100,
            "trimming_stability_score": [1.0] * 100,
        }
    )

    panel, trace = optimize_diverse_rate_panel(
        frame,
        panel_size=10,
        redundancy_penalty=0.10,
        candidate_top_fraction=1.0,
        candidate_minimum_pool_size=100,
        maximum_score_drop=1.0,
        number_rate_bins=5,
    )

    assert len(panel) == 10
    assert panel["rate_bin"].nunique() == 5
    assert sorted(panel["rate_bin"].value_counts().tolist()) == [2, 2, 2, 2, 2]
    assert trace["rate_bin"].nunique() == 5


def test_clade_balance_uses_group_completeness() -> None:
    metrics = pd.DataFrame(
        {
            "gene_id": ["balanced", "unbalanced"],
            "n_taxa": [10, 10],
            "taxon_occupancy": [1.0, 1.0],
            "sequence_completeness": [0.9, 0.9],
            "cell_occupancy": [0.9, 0.9],
            "mean_group_occupancy": [1.0, 1.0],
            "min_replicated_group_occupancy": [1.0, 1.0],
            "min_group_sequence_completeness": [0.9, 0.4],
            "sd_group_sequence_completeness": [0.01, 0.25],
            "worst_group_gap_fraction": [0.1, 0.6],
            "alignment_length": [500, 500],
            "gap_fraction": [0.1, 0.1],
            "ambiguous_fraction": [0.0, 0.0],
            "pis_per_length": [0.2, 0.2],
            "variable_sites_per_length": [0.3, 0.3],
            "mean_entropy": [0.5, 0.5],
            "parsimony_informative_sites": [100, 100],
            "mean_pairwise_distance": [0.2, 0.2],
            "composition_variability": [0.01, 0.01],
            "trimming_stability_score": [1.0, 1.0],
        }
    )
    scored = add_biological_scores(metrics, {"eligibility": {}})
    balanced = scored.loc[scored["gene_id"] == "balanced", "clade_balance_score"].iloc[0]
    unbalanced = scored.loc[scored["gene_id"] == "unbalanced", "clade_balance_score"].iloc[0]
    assert balanced > unbalanced
