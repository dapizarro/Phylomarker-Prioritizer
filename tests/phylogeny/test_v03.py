from pathlib import Path
from phylomarker_phylogeny.core import command_expected_output, deterministic_subsets

def test_expected_iqtree_tree():
    p=command_expected_output("iqtree2 -s a.faa --prefix /tmp/x -m LG+G")
    assert p == Path("/tmp/x.treefile")

def test_expected_concordance_tree():
    p=command_expected_output("iqtree2 -t x --gcf genes -s a --prefix /tmp/cf")
    assert p == Path("/tmp/cf.cf.tree")

def test_expected_astral_tree():
    p=command_expected_output("java -jar astral.jar -i genes.tre -o /tmp/species.tre")
    assert p == Path("/tmp/species.tre")

def test_ten_gene_jackknife_has_unique_subsets():
    genes=[f"g{i}" for i in range(10)]
    subsets=deterministic_subsets(genes,20,0.8,123)
    assert len(subsets)==20
    assert all(len(x)==8 for x in subsets)
    assert len({tuple(x) for x in subsets})==20
