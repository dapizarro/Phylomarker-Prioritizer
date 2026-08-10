from __future__ import annotations
from pathlib import Path
import argparse, json, shutil, sys
import pandas as pd
from .core import *

def context(config_path: Path):
    cfg=load_yaml(config_path); base=config_path.parent.resolve()
    inp=cfg["inputs"]; out=resolve(base,inp["output"]); sel=resolve(base,inp["selection_results"])
    sequence_type=inp.get("sequence_type","protein")
    profiles=inp["profiles"]; size=int(inp.get("panel_size",10))
    metadata=resolve(base,inp["metadata"]) if inp.get("metadata") else None
    taxa=None
    if metadata:
        df=pd.read_csv(metadata,sep="\t"); col=inp.get("sample_id_column","sample_ID")
        if col not in df.columns: raise ValueError(f"Metadata missing column: {col}")
        taxa=df[col].astype(str).tolist()
        if len(taxa)!=len(set(taxa)): raise ValueError("Duplicate sample IDs in metadata")
    return cfg,base,out,sel,sequence_type,profiles,size,taxa

def validate(args):
    cfg,base,out,sel,seq,profiles,size,taxa=context(Path(args.config))
    rows=[]; warnings=[]
    for profile in profiles:
        genes=panel_genes(sel,profile,size)
        for gene in genes:
            p=alignment_path(sel,gene,seq); recs=read_fasta(p)
            if taxa:
                extra=set(recs)-set(taxa)
                if extra: raise ValueError(f"{gene}: unknown taxa {sorted(extra)}")
            rows.append({"profile":profile,"gene_id":gene,"alignment":str(p),"n_taxa":len(recs),
                         "alignment_length":len(next(iter(recs.values())))})
    out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out/"validation.tsv",sep="\t",index=False)
    (out/"validation_warnings.tsv").write_text("severity\tcode\tmessage\n")
    print(f"Validated {len(rows)} panel-gene entries across {len(profiles)} profiles.")

def prepare(args):
    cfg,base,out,sel,seq,profiles,size,taxa=context(Path(args.config))
    ext="faa" if seq=="protein" else "fna"
    manifest=[]
    for profile in profiles:
        pdir=out/"panels"/profile; pdir.mkdir(parents=True,exist_ok=True)
        genes=panel_genes(sel,profile,size)
        alignments=[(g,alignment_path(sel,g,seq)) for g in genes]
        table=[]
        for g,p in alignments:
            dst=pdir/"genes"/p.name; dst.parent.mkdir(parents=True,exist_ok=True)
            if dst.exists() or dst.is_symlink(): dst.unlink()
            if cfg.get("prepare",{}).get("copy_alignments",False): shutil.copy2(p,dst)
            else: dst.symlink_to(p.resolve())
            table.append({"gene_id":g,"source_alignment":str(p.resolve()),"panel_alignment":str(dst)})
        pd.DataFrame(table).to_csv(pdir/"gene_manifest.tsv",sep="\t",index=False)
        stats=concatenate(alignments,taxa,pdir/f"{profile}.concat.{ext}",pdir/f"{profile}.partitions.nex")
        (pdir/"concatenation_stats.json").write_text(json.dumps(stats,indent=2))
    print(f"Prepared {len(profiles)} panels in {out/'panels'}")

def gene_trees(args):
    cfg,base,out,sel,seq,profiles,size,taxa=context(Path(args.config))
    iq=cfg["iqtree"]; exe=iq.get("executable","iqtree3"); commands=[]
    unique_genes=sorted({gene for profile in profiles for gene in panel_genes(sel,profile,size)})
    for gene in unique_genes:
        aln=alignment_path(sel,gene,seq)
        pref=out/"gene_trees"/"by_gene"/gene/gene
        pref.parent.mkdir(parents=True,exist_ok=True)
        commands.append(iqtree_command(exe,aln,pref,iq.get("gene_model","MFP"),
            iq.get("threads_per_gene",2),int(iq.get("seed",20260803)),
            int(iq.get("gene_bootstrap",1000)),int(iq.get("gene_alrt",1000)),None,
            iq.get("gene_extra",[])))
    manifest=out/"manifests/gene_trees.commands.txt"; n=write_manifest(commands,manifest)
    emit_slurm(cfg,out,manifest,"gene_trees",n,cfg.get("slurm",{}).get("gene_tree_cpus",2),
               cfg.get("slurm",{}).get("gene_tree_mem","4G"),cfg.get("slurm",{}).get("gene_tree_time","08:00:00"))
    if args.execute:
        failures=sum(run_command(c,out/f"logs/gene_tree_{i:04d}.log")!=0 for i,c in enumerate(commands,1))
        if failures: raise SystemExit(f"{failures} gene-tree jobs failed")
    print(f"Wrote {n} gene-tree commands: {manifest}")

def gene_models(args):
    """Extract BIC models from completed MFP gene analyses."""
    cfg,base,out,sel,seq,profiles,size,taxa=context(Path(args.config))
    unique_genes=sorted({g for p in profiles for g in panel_genes(sel,p,size)})
    policy=cfg.get("iqtree",{}).get("missing_model_policy","error")
    fallback=cfg.get("iqtree",{}).get("fallback_model","LG+G")
    models={}; rows=[]; missing=[]
    for gene in unique_genes:
        report=out/"gene_trees"/"by_gene"/gene/f"{gene}.iqtree"
        try:
            model=read_iqtree_best_model(report); status="ok"
        except Exception as exc:
            missing.append(gene)
            model=fallback if policy=="fallback" else ""
            status=f"{'fallback' if model else 'error'}:{type(exc).__name__}"
        if model: models[gene]=model
        rows.append({"gene_id":gene,"model":model,"report_path":str(report),"status":status})
    outsum=out/"summaries"; outsum.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(outsum/"gene_models.tsv",sep="\t",index=False)
    if missing and policy=="error":
        raise ValueError(f"{len(missing)} genes lack selected models; see {outsum/'gene_models.tsv'}")
    for profile in profiles:
        stats=json.loads((out/"panels"/profile/"concatenation_stats.json").read_text())
        parts=[tuple(x) for x in stats["partitions"]]
        write_gene_model_nexus(parts,models,out/"panels"/profile/f"{profile}.gene_models.nex")
    print(f"Extracted {len(models)} gene models and wrote model-aware partitions.")

def concatenated(args):
    cfg,base,out,sel,seq,profiles,size,taxa=context(Path(args.config))
    iq=cfg["iqtree"]; exe=iq.get("executable","iqtree3")
    ext="faa" if seq=="protein" else "fna"; commands=[]
    strategy=iq.get("concat_partition_strategy","gene_models")
    for profile in profiles:
        pdir=out/"panels"/profile; aln=pdir/f"{profile}.concat.{ext}"
        if not aln.is_file(): raise FileNotFoundError(f"Run prepare first: {aln}")
        if strategy=="gene_models":
            part=pdir/f"{profile}.gene_models.nex"; model=None
            if not part.is_file(): raise FileNotFoundError(f"Run gene-models first: {part}")
            root="concatenated_gene_models"
        elif strategy=="model_finder_merge":
            part=pdir/f"{profile}.partitions.nex"; model="MFP+MERGE"; root="concatenated"
        elif strategy=="fixed":
            part=pdir/f"{profile}.partitions.nex"; model=iq.get("concat_fixed_model","LG+G"); root="concatenated"
        else:
            raise ValueError(f"Unknown concat_partition_strategy: {strategy}")
        pref=out/root/profile/profile; pref.parent.mkdir(parents=True,exist_ok=True)
        commands.append(iqtree_command(exe,aln,pref,model,iq.get("threads_concat","AUTO"),
            int(iq.get("seed",20260803)),int(iq.get("concat_bootstrap",1000)),
            int(iq.get("concat_alrt",1000)),part,iq.get("concat_extra",[])))
    manifest=out/"manifests"/f"concatenated.{strategy}.commands.txt"
    n=write_manifest(commands,manifest)
    emit_slurm(cfg,out,manifest,f"concat_{strategy}",n,cfg.get("slurm",{}).get("concat_cpus",8),
               cfg.get("slurm",{}).get("concat_mem","16G"),cfg.get("slurm",{}).get("concat_time","24:00:00"))
    if args.execute:
        failures=sum(run_command(c,out/f"logs/concat_{strategy}_{i:04d}.log")!=0 for i,c in enumerate(commands,1))
        if failures: raise SystemExit(f"{failures} concatenated jobs failed")
    print(f"Wrote {n} concatenated-tree commands: {manifest}")


def collect_gene_trees(out: Path, profile: str, genes: list[str], destination: Path) -> None:
    lines=[]
    for gene in genes:
        p=out/"gene_trees"/"by_gene"/gene/f"{gene}.treefile"
        if not p.is_file(): raise FileNotFoundError(f"Missing gene tree: {p}")
        text=p.read_text().strip()
        if not text.endswith(";"): raise ValueError(f"Invalid Newick: {p}")
        lines.append(text)
    destination.parent.mkdir(parents=True,exist_ok=True); destination.write_text("\n".join(lines)+"\n")

def astral(args):
    cfg,base,out,sel,seq,profiles,size,taxa=context(Path(args.config))
    acfg=cfg["astral"]; jar=resolve(base,acfg["jar"]); java=acfg.get("java","java"); commands=[]
    if not jar.is_file(): raise FileNotFoundError(f"ASTRAL jar not found: {jar}")
    for profile in profiles:
        genes=panel_genes(sel,profile,size); adir=out/"astral"/profile; adir.mkdir(parents=True,exist_ok=True)
        trees=adir/"gene_trees.tre"; collect_gene_trees(out,profile,genes,trees)
        output=adir/"species_tree.tre"
        commands.append(astral_command(java,jar,trees,output,acfg.get("extra",[])))
    manifest=out/"manifests/astral.commands.txt"; n=write_manifest(commands,manifest)
    emit_slurm(cfg,out,manifest,"astral",n,cfg.get("slurm",{}).get("astral_cpus",2),
               cfg.get("slurm",{}).get("astral_mem","8G"),cfg.get("slurm",{}).get("astral_time","08:00:00"))
    if args.execute:
        failures=sum(run_command(c,out/f"logs/astral_{i:04d}.log")!=0 for i,c in enumerate(commands,1))
        if failures: raise SystemExit(f"{failures} ASTRAL jobs failed")
    print(f"Wrote {n} ASTRAL commands: {manifest}")

def concordance(args):
    cfg,base,out,sel,seq,profiles,size,taxa=context(Path(args.config))
    iq=cfg["iqtree"]; exe=iq.get("executable","iqtree3")
    ext="faa" if seq=="protein" else "fna"; commands=[]
    strategy=iq.get("concat_partition_strategy","gene_models")
    concat_root="concatenated_gene_models" if strategy=="gene_models" else "concatenated"
    concord_root="concordance_gene_models" if strategy=="gene_models" else "concordance"
    for profile in profiles:
        ref=out/concat_root/profile/f"{profile}.treefile"
        trees=out/"astral"/profile/"gene_trees.tre"
        aln=out/"panels"/profile/f"{profile}.concat.{ext}"
        part=out/"panels"/profile/(f"{profile}.gene_models.nex" if strategy=="gene_models" else f"{profile}.partitions.nex")
        for p in (ref,trees,aln,part):
            if not p.is_file(): raise FileNotFoundError(f"Missing prerequisite: {p}")
        pref=out/concord_root/profile/profile; pref.parent.mkdir(parents=True,exist_ok=True)
        cmd=[exe,"-t",str(ref),"--gcf",str(trees),"-s",str(aln),"--scf",str(iq.get("scf_quartets",100)),
             "-p",str(part),"--prefix",str(pref),"-T",str(iq.get("threads_concordance","AUTO")),
             "--seed",str(iq.get("seed",20260803))]
        commands.append(" ".join(quote(x) for x in cmd))
    manifest=out/"manifests"/f"concordance.{strategy}.commands.txt"; n=write_manifest(commands,manifest)
    emit_slurm(cfg,out,manifest,f"concordance_{strategy}",n,cfg.get("slurm",{}).get("concordance_cpus",4),
               cfg.get("slurm",{}).get("concordance_mem","8G"),cfg.get("slurm",{}).get("concordance_time","08:00:00"))
    if args.execute:
        failures=sum(run_command(c,out/f"logs/concordance_{strategy}_{i:04d}.log")!=0 for i,c in enumerate(commands,1))
        if failures: raise SystemExit(f"{failures} concordance jobs failed")
    print(f"Wrote {n} concordance commands: {manifest}")

def compare_trees(args):
    cfg,base,out,sel,seq,profiles,size,taxa=context(Path(args.config))
    strategy=cfg["iqtree"].get("concat_partition_strategy","gene_models")
    concat_root="concatenated_gene_models" if strategy=="gene_models" else "concatenated"
    rows=[]
    for profile in profiles:
        concat=out/concat_root/profile/f"{profile}.treefile"
        astral_tree=out/"astral"/profile/"species_tree.tre"
        for p in (concat,astral_tree):
            if not p.is_file(): raise FileNotFoundError(p)
        result=compare_unrooted_trees(concat,astral_tree)
        rows.append({"profile":profile,"first":"concatenated","second":"astral",**result})
    outsum=out/"summaries"; outsum.mkdir(parents=True,exist_ok=True)
    table=outsum/"tree_comparisons.tsv"
    pd.DataFrame(rows).to_csv(table,sep="\t",index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"Wrote {table}")


def jackknife(args):
    cfg,base,out,sel,seq,profiles,size,taxa=context(Path(args.config))
    jcfg=cfg.get("jackknife",{}); iq=cfg["iqtree"]; ext="faa" if seq=="protein" else "fna"; commands=[]; rows=[]
    nrep=int(jcfg.get("replicates",50)); keep=float(jcfg.get("keep_fraction",0.8)); seed=int(jcfg.get("seed",20260803))
    for pi,profile in enumerate(profiles):
        genes=panel_genes(sel,profile,size)
        subsets=deterministic_subsets(genes,nrep,keep,seed+pi)
        for r,subset in enumerate(subsets,1):
            rdir=out/"jackknife"/profile/f"replicate_{r:03d}"; rdir.mkdir(parents=True,exist_ok=True)
            aligns=[(g,alignment_path(sel,g,seq)) for g in subset]
            aln=rdir/f"{profile}.jk{r:03d}.{ext}"; part=rdir/"partitions.nex"
            concatenate(aligns,taxa,aln,part)
            pref=rdir/"tree"
            commands.append(iqtree_command(iq.get("executable","iqtree3"),aln,pref,
                jcfg.get("model","MFP+MERGE"),jcfg.get("threads",4),seed+r,0,0,part,jcfg.get("extra",[])))
            rows.append({"profile":profile,"replicate":r,"n_genes":len(subset),"genes":",".join(subset)})
    pd.DataFrame(rows).to_csv(out/"jackknife/jackknife_manifest.tsv",sep="\t",index=False)
    manifest=out/"manifests/jackknife.commands.txt"; n=write_manifest(commands,manifest)
    emit_slurm(cfg,out,manifest,"jackknife",n,cfg.get("slurm",{}).get("jackknife_cpus",4),
               cfg.get("slurm",{}).get("jackknife_mem","8G"),cfg.get("slurm",{}).get("jackknife_time","12:00:00"))
    if args.execute:
        failures=sum(run_command(c,out/f"logs/jackknife_{i:04d}.log")!=0 for i,c in enumerate(commands,1))
        if failures: raise SystemExit(f"{failures} jackknife jobs failed")
    print(f"Wrote {n} jackknife commands: {manifest}")


def jackknife_summary(args):
    """Compare every completed jackknife tree with the full concatenated tree."""
    cfg,base,out,sel,seq,profiles,size,taxa=context(Path(args.config))
    strategy=cfg["iqtree"].get("concat_partition_strategy","gene_models")
    concat_root="concatenated_gene_models" if strategy=="gene_models" else "concatenated"
    replicate_rows=[]; summary_rows=[]
    for profile in profiles:
        reference=out/concat_root/profile/f"{profile}.treefile"
        if not reference.is_file():
            raise FileNotFoundError(f"Missing full-panel reference tree: {reference}")
        manifest_path=out/"jackknife"/"jackknife_manifest.tsv"
        if not manifest_path.is_file():
            raise FileNotFoundError(f"Run jackknife first: {manifest_path}")
        manifest=pd.read_csv(manifest_path,sep="\t")
        subset=manifest.loc[manifest["profile"]==profile].copy()
        for _,row in subset.iterrows():
            replicate=int(row["replicate"])
            tree=out/"jackknife"/profile/f"replicate_{replicate:03d}"/"tree.treefile"
            if not tree.is_file() or tree.stat().st_size==0:
                replicate_rows.append({
                    "profile":profile,"replicate":replicate,"n_genes":int(row["n_genes"]),
                    "genes":row["genes"],"status":"missing","rf":None,
                    "rf_normalized":None,"identical":None
                })
                continue
            result=compare_unrooted_trees(reference,tree)
            replicate_rows.append({
                "profile":profile,"replicate":replicate,"n_genes":int(row["n_genes"]),
                "genes":row["genes"],"status":"ok","rf":result["rf"],
                "rf_normalized":result["rf_normalized"],
                "identical":result["rf"]==0
            })
        ok=pd.DataFrame([r for r in replicate_rows if r["profile"]==profile and r["status"]=="ok"])
        expected=len(subset)
        if ok.empty:
            summary_rows.append({
                "profile":profile,"expected_replicates":expected,"completed_replicates":0,
                "identical_fraction":None,"rf_normalized_mean":None,
                "rf_normalized_median":None,"rf_normalized_max":None
            })
        else:
            summary_rows.append({
                "profile":profile,"expected_replicates":expected,
                "completed_replicates":len(ok),
                "identical_fraction":float(ok["identical"].mean()),
                "rf_normalized_mean":float(ok["rf_normalized"].mean()),
                "rf_normalized_median":float(ok["rf_normalized"].median()),
                "rf_normalized_max":float(ok["rf_normalized"].max())
            })
    outsum=out/"summaries"; outsum.mkdir(parents=True,exist_ok=True)
    per_rep=outsum/"jackknife_replicates.tsv"
    summary=outsum/"jackknife_summary.tsv"
    pd.DataFrame(replicate_rows).to_csv(per_rep,sep="\t",index=False)
    frame=pd.DataFrame(summary_rows)
    frame.to_csv(summary,sep="\t",index=False)
    print(frame.to_string(index=False))
    print(f"Wrote {per_rep}")
    print(f"Wrote {summary}")

def summarize(args):
    cfg,base,out,sel,seq,profiles,size,taxa=context(Path(args.config))
    rows=[]
    for profile in profiles:
        genes=panel_genes(sel,profile,size)
        present_gene=sum((out/"gene_trees"/"by_gene"/g/f"{g}.treefile").is_file() for g in genes)
        strategy=cfg["iqtree"].get("concat_partition_strategy","gene_models")
        concat_root="concatenated_gene_models" if strategy=="gene_models" else "concatenated"
        concord_root="concordance_gene_models" if strategy=="gene_models" else "concordance"
        concat=out/concat_root/profile/f"{profile}.treefile"
        astr=out/"astral"/profile/"species_tree.tre"
        cf=out/concord_root/profile/f"{profile}.cf.tree"
        jk=list((out/"jackknife"/profile).glob("replicate_*/tree.treefile")) if (out/"jackknife"/profile).exists() else []
        stats=json.loads((out/"panels"/profile/"concatenation_stats.json").read_text()) if (out/"panels"/profile/"concatenation_stats.json").is_file() else {}
        rows.append({"profile":profile,"expected_gene_trees":len(genes),"completed_gene_trees":present_gene,
                     "concatenated_tree_complete":concat.is_file(),"astral_tree_complete":astr.is_file(),
                     "concordance_complete":cf.is_file(),"jackknife_trees_complete":len(jk),
                     "concatenated_length":stats.get("alignment_length"),"n_taxa":stats.get("n_taxa")})
    outsum=out/"summaries"; outsum.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(outsum/"pipeline_status.tsv",sep="\t",index=False)
    versions={"iqtree":executable_version(f"{cfg['iqtree'].get('executable','iqtree3')} --version"),
              "java":executable_version(f"{cfg['astral'].get('java','java')} -version"),
              "config_sha256":sha256(Path(args.config))}
    (outsum/"provenance.json").write_text(json.dumps(versions,indent=2))
    print(pd.DataFrame(rows).to_string(index=False))

def emit_slurm(cfg,out,manifest,name,n,cpus,mem,time):
    s=cfg.get("slurm",{})
    modules=[f"module load {m}" for m in s.get("modules",[])]
    write_slurm_array(manifest,out/"slurm"/f"{name}.sbatch",f"pm_{name}",n,int(cpus),str(mem),str(time),
                      s.get("partition"),modules,s.get("conda_activate"))

def local_run(args):
    """Run the complete local workflow sequentially and safely resume completed work."""
    args.execute=True
    print("[1/10] validate")
    validate(args)
    print("[2/10] prepare")
    prepare(args)
    print("[3/10] gene trees")
    gene_trees(args)
    print("[4/10] gene models")
    gene_models(args)
    print("[5/10] concatenated trees")
    concatenated(args)
    print("[6/10] ASTRAL")
    astral(args)
    print("[7/10] concordance")
    concordance(args)
    print("[8/10] compare concatenated and ASTRAL trees")
    compare_trees(args)
    print("[9/10] jackknife")
    jackknife(args)
    print("[10/10] jackknife evaluation and final summary")
    jackknife_summary(args)
    summarize(args)
    print("Local workflow completed successfully.")

def all_cmd(args):
    local_run(args)

def parser():
    p=argparse.ArgumentParser(prog="phylomarker-phylogeny")
    p.add_argument("--version",action="version",version="0.3.0")
    sub=p.add_subparsers(dest="command",required=True)
    for name,fn in [("validate",validate),("prepare",prepare),("gene-trees",gene_trees),
                    ("gene-models",gene_models),("concatenated",concatenated),("astral",astral),
                    ("concordance",concordance),("compare-trees",compare_trees),
                    ("jackknife",jackknife),("jackknife-summary",jackknife_summary),("summarize",summarize),("local-run",local_run),("all",all_cmd)]:
        q=sub.add_parser(name); q.add_argument("--config",required=True); q.set_defaults(func=fn)
        if name in {"gene-trees","concatenated","astral","concordance","jackknife"}:
            q.add_argument("--execute",action="store_true",help="Execute locally instead of only writing manifests")
    return p

def main():
    args=parser().parse_args()
    try: args.func(args)
    except Exception as e:
        print(f"ERROR: {e}",file=sys.stderr)
        if "--debug" in sys.argv: raise
        raise SystemExit(1)
