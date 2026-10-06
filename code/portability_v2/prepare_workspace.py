#!/usr/bin/env python3
"""Prepare isolated logical code layout only; never run scientific programs."""
import argparse, csv, hashlib, json, os, shutil
from pathlib import Path

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def safe(root, rel):
    q=Path(rel)
    if q.is_absolute() or ".." in q.parts: raise ValueError("Unsafe relative path")
    return root/q

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--package",type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument("--output",type=Path,required=True,help="New isolated directory; must not exist")
    p.add_argument("--external-inputs",type=Path,help="Optional TSV: logical_path, source_path, sha256 (author-supplied frozen inputs)")
    p.add_argument("--runtime-root",type=Path,help="Optional existing root containing the recorded envs layout")
    a=p.parse_args(); root=a.output.resolve(); package=a.package.resolve()
    if root.exists(): p.error("Output already exists; refusing overwrite")
    rows=list(csv.DictReader((package/"provenance/ORIGIN_MANIFEST.tsv").open(),delimiter="\t"))
    plan=[]
    for r in rows:
        src=safe(package,r["relative_path"])
        if digest(src)!=r["packaged_sha256"]: raise ValueError("Packaged integrity mismatch: "+r["relative_path"])
        dst=safe(root,r["source_locator"])
        plan.append((src,dst))
    ext=[]
    if a.external_inputs:
        for r in csv.DictReader(a.external_inputs.open(),delimiter="\t"):
            src=Path(r["source_path"]).resolve(); dst=safe(root,r["logical_path"])
            if not r["sha256"] or digest(src)!=r["sha256"]: raise ValueError("External integrity mismatch")
            # Any replacement for a packaged source must match its original frozen digest.
            match=[x for x in rows if x["source_locator"]==r["logical_path"]]
            if match and any(r["sha256"]!=x["source_sha256"] for x in match): raise ValueError("Replacement does not match frozen original hash")
            if dst in [d for _,d in ext]:raise ValueError("Duplicate external destination")
            ext.append((src,dst))
    root.mkdir(parents=True)
    for src,dst in plan+ext:
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
    # Original script literals resolve relative to their dataset cwd, with no code edits.
    for dataset in [root]+[root/n for n in ("GSE300977","GSE106329","GSE133111","TMD-RiboSeq") if (root/n).exists()]:
        (dataset/"REPRO_WORKSPACE").symlink_to(root,target_is_directory=True)
        if a.runtime_root: (dataset/"REPRO_ENV").symlink_to(a.runtime_root.resolve(),target_is_directory=True)
    if (root/"GSE133111").exists() and not (root/"GSE133111_stage0C").exists():
        (root/"GSE133111_stage0C").symlink_to(root/"GSE133111",target_is_directory=True)
    result={"copied_source_files":len(plan),"verified_external_files":len(ext),"scientific_execution":"NOT_RUN","runtime_bound":bool(a.runtime_root),"note":"Use the dataset-specific working directory. Set PAPER3_WORK_ROOT explicitly. Missing inputs and original hash gates still stop execution."}
    (root/"WORKSPACE_PREPARATION.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result))
if __name__=="__main__":main()
