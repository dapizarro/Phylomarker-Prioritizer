from pathlib import Path
import pytest
from phylomarker_phylogeny.core import *

def make_fasta(path, records):
    write_fasta(records,path)

def test_read_fasta_rejects_unaligned(tmp_path):
    p=tmp_path/"x.faa"; p.write_text(">a\nAAA\n>b\nAA\n")
    with pytest.raises(ValueError): read_fasta(p)

def test_concatenate_fills_missing_taxa_and_writes_nexus(tmp_path):
    a=tmp_path/"a.faa"; b=tmp_path/"b.faa"
    make_fasta(a,{"t1":"AAA","t2":"A-A"})
    make_fasta(b,{"t2":"CC","t3":"C-"})
    out=tmp_path/"c.faa"; part=tmp_path/"p.nex"
    stats=concatenate([("g1",a),("g2",b)],["t1","t2","t3"],out,part)
    rec=read_fasta(out)
    assert rec["t1"]=="AAA--"
    assert rec["t3"]=="---C-"
    assert stats["alignment_length"]==5
    text=part.read_text()
    assert text.startswith("#nexus")
    assert "charset g1 = 1-3;" in text
    assert "AA, g1" not in text

def test_read_best_model(tmp_path):
    p=tmp_path/"g.iqtree"
    p.write_text("Best-fit model according to BIC: LG+F+R5\n")
    assert read_iqtree_best_model(p)=="LG+F+R5"

def test_write_gene_model_nexus(tmp_path):
    p=tmp_path/"models.nex"
    write_gene_model_nexus([("g1",1,10),("g2",11,20)],
                           {"g1":"LG+F+G4","g2":"WAG+R3"},p)
    text=p.read_text()
    assert "LG+F+G4:g1," in text
    assert "WAG+R3:g2;" in text

def test_compare_unrooted_trees(tmp_path):
    pytest.importorskip("dendropy")
    a=tmp_path/"a.tre"; b=tmp_path/"b.tre"; c=tmp_path/"c.tre"
    a.write_text("((A,B),(C,(D,E)));\n")
    b.write_text("((A,B),(C,(D,E)));\n")
    c.write_text("((A,C),(B,(D,E)));\n")
    assert compare_unrooted_trees(a,b)["rf"]==0
    assert compare_unrooted_trees(a,c)["rf"]>0

def test_deterministic_subsets_reproducible():
    genes=[f"g{i}" for i in range(10)]
    assert deterministic_subsets(genes,20,.8,7)==deterministic_subsets(genes,20,.8,7)
    assert all(len(x)==8 for x in deterministic_subsets(genes,20,.8,7))

def test_iqtree_command_partition_without_global_model():
    c=iqtree_command("iqtree2",Path("a.faa"),Path("x"),None,4,3,0,0,Path("p.nex"))
    assert "-p p.nex" in c
    assert " -m " not in f" {c} "

def test_astral_command():
    c=astral_command("java",Path("a.jar"),Path("g.tre"),Path("s.tre"))
    assert c=="java -jar a.jar -i g.tre -o s.tre"

def test_panel_genes(tmp_path):
    p=tmp_path/"panels/x/n2"; p.mkdir(parents=True); (p/"genes.txt").write_text("a\nb\n")
    assert panel_genes(tmp_path,"x",2)==["a","b"]

def test_write_manifest(tmp_path):
    assert write_manifest(["a","b"],tmp_path/"m")==2
    assert (tmp_path/"m").read_text()=="a\nb\n"
