from __future__ import annotations
from pathlib import Path
from typing import Iterable
from dataclasses import dataclass
import re
import hashlib, json, random, shutil, subprocess, shlex
import pandas as pd
import yaml

def load_yaml(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {path}")
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict):
        raise ValueError("Configuration root must be a mapping")
    return data

def resolve(base: Path, value: str | Path) -> Path:
    p=Path(value)
    return p if p.is_absolute() else (base/p).resolve()

def read_fasta(path: Path) -> dict[str,str]:
    recs={}; current=None; chunks=[]
    for raw in path.read_text().splitlines():
        line=raw.strip()
        if not line: continue
        if line.startswith(">"):
            if current is not None:
                recs[current]="".join(chunks).upper()
            current=line[1:].split()[0]; chunks=[]
            if current in recs: raise ValueError(f"Duplicate FASTA ID {current}: {path}")
        else:
            if current is None: raise ValueError(f"Sequence before header: {path}")
            chunks.append(line)
    if current is not None: recs[current]="".join(chunks).upper()
    if not recs: raise ValueError(f"Empty FASTA: {path}")
    lengths={len(v) for v in recs.values()}
    if len(lengths)!=1: raise ValueError(f"Unaligned FASTA (different lengths): {path}")
    return recs

def write_fasta(records: dict[str,str], path: Path, width: int=80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as h:
        for name,seq in records.items():
            h.write(f">{name}\n")
            for i in range(0,len(seq),width): h.write(seq[i:i+width]+"\n")

def panel_genes(selection_root: Path, profile: str, size: int) -> list[str]:
    p=selection_root/"panels"/profile/f"n{size}"/"genes.txt"
    if not p.is_file(): raise FileNotFoundError(f"Panel gene list not found: {p}")
    genes=[x.strip() for x in p.read_text().splitlines() if x.strip()]
    if len(genes)!=size: raise ValueError(f"{p}: expected {size} genes, found {len(genes)}")
    if len(set(genes))!=len(genes): raise ValueError(f"Duplicate genes in {p}")
    return genes

def alignment_path(selection_root: Path, gene: str, sequence_type: str) -> Path:
    ext="faa" if sequence_type=="protein" else "fna"
    p=selection_root/"alignments"/sequence_type/"trimmed"/f"{gene}.trimmed.{ext}"
    if not p.is_file(): raise FileNotFoundError(f"Alignment not found: {p}")
    return p

def concatenate(alignments: list[tuple[str,Path]], taxa_order: list[str]|None,
                out_fasta: Path, out_partitions: Path) -> dict:
    """Concatenate alignments and write a valid charset-only NEXUS file."""
    parsed=[]; all_taxa=set()
    for gene,path in alignments:
        recs=read_fasta(path); parsed.append((gene,recs)); all_taxa.update(recs)
    taxa=taxa_order or sorted(all_taxa)
    unknown=all_taxa-set(taxa)
    if unknown: raise ValueError(f"Taxa absent from metadata/order: {sorted(unknown)}")
    concat={t:"" for t in taxa}; start=1; parts=[]
    for gene,recs in parsed:
        L=len(next(iter(recs.values())))
        for t in taxa: concat[t]+=recs.get(t,"-"*L)
        end=start+L-1
        parts.append((gene,start,end)); start=end+1
    write_fasta(concat,out_fasta)
    write_charset_nexus(parts, out_partitions)
    return {"n_taxa":len(taxa),"n_genes":len(parsed),
            "alignment_length":start-1,"partitions":parts}

def write_charset_nexus(parts: list[tuple[str,int,int]], path: Path) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    lines=["#nexus","begin sets;"]
    lines += [f"    charset {g} = {s}-{e};" for g,s,e in parts]
    lines += ["end;"]
    path.write_text("\n".join(lines)+"\n")

BEST_MODEL_PATTERNS = (
    re.compile(r"Best-fit model according to BIC:\s*(\S+)", re.I),
    re.compile(r"Best-fit substitution model(?: according to BIC)?\s*:\s*(\S+)", re.I),
    re.compile(r"Model of substitution:\s*(\S+)", re.I),
)

def read_iqtree_best_model(report: Path) -> str:
    if not report.is_file(): raise FileNotFoundError(f"Gene report not found: {report}")
    text=report.read_text(encoding="utf-8",errors="replace")
    for pattern in BEST_MODEL_PATTERNS:
        match=pattern.search(text)
        if match: return match.group(1).strip().rstrip(";,")
    raise ValueError(f"No BIC-selected model found in {report}; run gene trees with -m MFP")

def write_gene_model_nexus(parts: list[tuple[str,int,int]], models: dict[str,str], path: Path) -> None:
    missing=[g for g,_,_ in parts if g not in models]
    if missing: raise ValueError(f"Missing gene models for: {', '.join(missing)}")
    lines=["#nexus","begin sets;"]
    lines += [f"    charset {g} = {s}-{e};" for g,s,e in parts]
    lines += ["","    charpartition gene_models ="]
    assignments=[f"{models[g]}:{g}" for g,_,_ in parts]
    lines += [f"        {x}{';' if i==len(assignments)-1 else ','}"
              for i,x in enumerate(assignments)]
    lines += ["end;"]
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text("\n".join(lines)+"\n")

def compare_unrooted_trees(tree_a: Path, tree_b: Path) -> dict:
    try:
        import dendropy
        from dendropy.calculate import treecompare
    except ImportError as exc:
        raise RuntimeError("DendroPy is required for tree comparison; install the project environment") from exc
    taxa=dendropy.TaxonNamespace()
    kwargs=dict(schema="newick",taxon_namespace=taxa,rooting="force-unrooted",
                preserve_underscores=True,suppress_internal_node_taxa=True)
    a=dendropy.Tree.get(path=str(tree_a),**kwargs)
    b=dendropy.Tree.get(path=str(tree_b),**kwargs)
    leaves_a={x.taxon.label for x in a.leaf_node_iter()}
    leaves_b={x.taxon.label for x in b.leaf_node_iter()}
    if leaves_a != leaves_b:
        raise ValueError(f"Trees have different taxa; only_a={sorted(leaves_a-leaves_b)}, "
                         f"only_b={sorted(leaves_b-leaves_a)}")
    a.encode_bipartitions(); b.encode_bipartitions()
    rf=treecompare.symmetric_difference(a,b)
    fp,fn=treecompare.false_positives_and_negatives(a,b)
    n=len(leaves_a); rf_max=max(0,2*(n-3))
    return {"n_taxa":n,"rf":rf,"splits_only_first":fp,"splits_only_second":fn,
            "rf_max_binary":rf_max,"rf_normalized":rf/rf_max if rf_max else 0.0}

def command_expected_output(command: str) -> Path | None:
    """Infer the principal completion file for an IQ-TREE or ASTRAL command."""
    tokens=shlex.split(command)
    if "--prefix" in tokens:
        prefix=Path(tokens[tokens.index("--prefix")+1])
        suffix=".cf.tree" if "--gcf" in tokens or "--scf" in tokens else ".treefile"
        return Path(str(prefix)+suffix)
    if "-o" in tokens:
        return Path(tokens[tokens.index("-o")+1])
    return None

def run_command(command: str, log_path: Path, dry_run: bool=False) -> int:
    log_path.parent.mkdir(parents=True,exist_ok=True)
    expected=command_expected_output(command)
    if expected is not None and expected.is_file() and expected.stat().st_size>0:
        log_path.write_text(f"SKIPPED completed output: {expected}\n")
        print(f"SKIPPED completed output: {expected}")
        return 0
    if dry_run:
        log_path.write_text(command+"\n"); return 0
    with log_path.open("w") as log:
        proc=subprocess.run(["bash","-lc",command],stdout=log,stderr=subprocess.STDOUT)
    return proc.returncode

def quote(x) -> str: return shlex.quote(str(x))

def iqtree_command(exe: str, alignment: Path, prefix: Path, model: str|None, threads: str|int, seed: int,
                   bootstrap: int=0, alrt: int=0, partitions: Path|None=None, extra: list[str]|None=None) -> str:
    cmd=[exe,"-s",str(alignment),"--prefix",str(prefix)]
    if model: cmd += ["-m",model]
    cmd += ["-T",str(threads),"--seed",str(seed)]
    if partitions: cmd += ["-p",str(partitions)]
    if bootstrap>0: cmd += ["-B",str(bootstrap)]
    if alrt>0: cmd += ["--alrt",str(alrt)]
    if extra: cmd += list(extra)
    return " ".join(quote(x) for x in cmd)

def astral_command(java: str, jar: Path, trees: Path, output: Path, extra: list[str]|None=None) -> str:
    cmd=[java,"-jar",str(jar),"-i",str(trees),"-o",str(output)]
    if extra: cmd+=list(extra)
    return " ".join(quote(x) for x in cmd)

def write_manifest(commands: Iterable[str], path: Path) -> int:
    commands=[c for c in commands if c.strip()]
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text("\n".join(commands)+("\n" if commands else ""))
    return len(commands)

def write_slurm_array(manifest: Path, output: Path, job_name: str, n: int, cpus: int, mem: str, time: str,
                      partition: str|None=None, modules: list[str]|None=None, conda_activate: str|None=None) -> None:
    lines=["#!/usr/bin/env bash",f"#SBATCH --job-name={job_name}",f"#SBATCH --array=1-{max(1,n)}",
           f"#SBATCH --cpus-per-task={cpus}",f"#SBATCH --mem={mem}",f"#SBATCH --time={time}",
           f"#SBATCH --output={output.parent}/slurm-%A_%a.out",f"#SBATCH --error={output.parent}/slurm-%A_%a.err",
           "set -euo pipefail"]
    if partition: lines.insert(3,f"#SBATCH --partition={partition}")
    lines += modules or []
    if conda_activate: lines += [f"source {conda_activate}"]
    lines += [f'COMMAND=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {quote(manifest)})',
              'if [[ -z "${COMMAND}" ]]; then echo "Empty task"; exit 1; fi',
              'echo "${COMMAND}"','bash -lc "${COMMAND}"']
    output.parent.mkdir(parents=True,exist_ok=True); output.write_text("\n".join(lines)+"\n"); output.chmod(0o755)

def deterministic_subsets(genes: list[str], n_replicates: int, keep_fraction: float, seed: int) -> list[list[str]]:
    k=max(2,min(len(genes),round(len(genes)*keep_fraction)))
    rng=random.Random(seed); seen=set(); out=[]; attempts=0
    while len(out)<n_replicates and attempts<n_replicates*100:
        subset=tuple(sorted(rng.sample(genes,k))); attempts+=1
        if subset not in seen: seen.add(subset); out.append(list(subset))
    return out

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

def executable_version(command: str) -> str:
    try:
        p=subprocess.run(["bash","-lc",command],capture_output=True,text=True,timeout=20)
        text=(p.stdout or p.stderr).strip().splitlines()
        return text[0] if text else f"returncode={p.returncode}"
    except Exception as e: return f"unavailable: {e}"
