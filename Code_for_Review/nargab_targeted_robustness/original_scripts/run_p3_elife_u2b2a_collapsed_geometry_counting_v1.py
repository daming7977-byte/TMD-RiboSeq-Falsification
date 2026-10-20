#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import argparse, gc, hashlib, csv

import numpy as np
import pandas as pd
import pysam
from plastid import (
    BAMGenomeArray,
    CenterMapFactory,
    SizeFilterFactory,
    GenomicSegment,
    SegmentChain,
)

LIBRARIES = [f"B{i:02d}" for i in range(1,13)]

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

def parse_geometry(s):
    blocks = []
    for token in str(s).split("^"):
        chrom, coords = token.split(":",1)
        a,b = map(int, coords.split("-"))
        blocks.append((chrom,a,b))
    total = sum(b-a for _,a,b in blocks)
    if total != 31:
        raise RuntimeError(f"GEOMETRY_NOT_31:{s}:{total}")
    return blocks

def geometry_positions(geometry, strand, bits):
    positions = []
    for chrom,a,b in parse_geometry(geometry):
        it = range(a,b) if strand == "+" else range(b-1,a-1,-1)
        for pos in it:
            positions.append((chrom,pos))
    if len(positions) != 31:
        raise RuntimeError("POSITION_LENGTH_FAIL")
    bits = str(bits)
    if len(bits) != 31 or set(bits) - {"0","1"}:
        raise RuntimeError("MASK_BITS_FAIL")
    keep = [p for p,bit in zip(positions,bits) if bit == "0"]
    if not keep:
        raise RuntimeError("FULLY_MASKED_GEOMETRY")
    return keep

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--contract", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    root = Path(args.root).resolve()
    contract = Path(args.contract).resolve()
    outdir = Path(args.outdir).resolve()

    if outdir.exists():
        print("STOP_OUTDIR_EXISTS =", outdir)
        return 2
    outdir.mkdir(parents=True)
    countdir = outdir / "geometry_counts"
    countdir.mkdir()

    catalog_path = root / "tmd_spatial_secondary_scan_v1/output/SECONDARY_11CENTER_GEOMETRY_CATALOG_v1.tsv"
    if not catalog_path.exists():
        print("STOP_MISSING_CATALOG =", catalog_path)
        return 3

    bams = {
        b: root / f"P3_ELIFE_UPGRADE/U2B1_v1/bam/{b}.NH1.unique.exact5umi_collapse.bam"
        for b in LIBRARIES
    }
    missing = [p for p in bams.values() if not p.exists()]
    if missing:
        print("STOP_MISSING_COLLAPSED_BAMS")
        for p in missing: print(p)
        return 4

    log = []
    def emit(*x):
        s = " ".join(map(str,x))
        print(s)
        log.append(s)

    emit("P3_ELIFE_U2B2A_COLLAPSED_GEOMETRY_COUNTING = START")
    emit("ROOT =", root)
    emit("OUTDIR =", outdir)
    emit("CONTRACT =", contract)
    emit("CONTRACT_SHA256 =", sha256(contract))
    emit("CATALOG =", catalog_path)
    emit("CATALOG_SHA256 =", sha256(catalog_path))
    emit("COUNTING_METHOD = plastid.BAMGenomeArray + CenterMapFactory(nibble=13)")
    emit("READ_LENGTH_FILTER = 27-36")
    emit("BIOLOGICAL_RESULT_CALCULATION = NO")
    emit("TOPOLOGY_INTERACTION_CALCULATION = NO")

    cat = pd.read_csv(catalog_path, sep="\t", dtype=str)
    expected_cols = ["geometry_id","strand","geometry","mask_bits","unmasked_nt"]
    if list(cat.columns) != expected_cols:
        emit("STOP_CATALOG_SCHEMA =", ",".join(cat.columns))
        return 5
    if cat["geometry_id"].duplicated().any():
        emit("STOP_DUPLICATE_GEOMETRY_ID = YES")
        return 6

    geom_positions = {}
    wanted = defaultdict(set)
    for _, r in cat.iterrows():
        gid = str(r["geometry_id"])
        strand = str(r["strand"])
        pos = geometry_positions(str(r["geometry"]), strand, str(r["mask_bits"]))
        if len(pos) != int(r["unmasked_nt"]):
            raise RuntimeError(f"UNMASKED_NT_MISMATCH:{gid}")
        geom_positions[gid] = pos
        for chrom,p in pos:
            wanted[(chrom,strand)].add(p)

    merged_intervals = []
    for (chrom,strand), positions in sorted(wanted.items()):
        xs = sorted(positions)
        if not xs:
            continue
        start = prev = xs[0]
        for pos in xs[1:]:
            if pos == prev + 1:
                prev = pos
                continue
            merged_intervals.append((chrom,start,prev+1,strand))
            start = prev = pos
        merged_intervals.append((chrom,start,prev+1,strand))

    emit("GEOMETRY_CATALOG_ROWS =", len(cat))
    emit("MERGED_GENOMIC_INTERVALS =", len(merged_intervals))
    emit("UNIQUE_UNMASKED_GENOMIC_POSITIONS =", sum(len(x) for x in wanted.values()))

    summary_rows = []

    for lib in LIBRARIES:
        bam = bams[lib]
        bai = Path(str(bam)+".bai")
        if not bai.exists():
            emit(f"STOP_MISSING_BAI:{lib} =", bai)
            return 7

        # Input audit only
        with pysam.AlignmentFile(bam, "rb") as bf:
            input_reads = bf.count(until_eof=True)
        emit(f"{lib} COLLAPSED_BAM_READS =", input_reads)
        emit(f"{lib} COLLAPSED_BAM_SHA256 =", sha256(bam))

        ga = BAMGenomeArray(str(bam))
        ga.set_mapping(CenterMapFactory(nibble=13))
        ga.add_filter("length_27_36", SizeFilterFactory(min=27, max=36))

        position_counts = {}
        for i,(chrom,a,b,strand) in enumerate(merged_intervals, start=1):
            seg = GenomicSegment(chrom,a,b,strand)
            chain = SegmentChain(seg)
            counts = np.asarray(chain.get_counts(ga, stranded=True), dtype=float)
            if len(counts) != b-a:
                raise RuntimeError(f"MERGED_COUNT_LENGTH_FAIL:{lib}:{chrom}:{a}-{b}:{strand}")
            positions = range(a,b) if strand == "+" else range(b-1,a-1,-1)
            for pos,val in zip(positions, counts):
                position_counts[(chrom,strand,pos)] = float(val)
            if i % 2000 == 0 or i == len(merged_intervals):
                emit(lib, "MERGED_INTERVALS_COUNTED =", i, "/", len(merged_intervals))

        rows = []
        total_density = 0.0
        for _,r in cat.iterrows():
            gid = str(r["geometry_id"])
            vals = [position_counts[(chrom,str(r["strand"]),pos)] for chrom,pos in geom_positions[gid]]
            density = float(np.sum(vals)/len(vals))
            total_density += density
            rows.append({
                "geometry_id": gid,
                "unmasked_nt": len(vals),
                "density": density,
            })

        out = pd.DataFrame(rows)
        outfile = countdir / f"{lib}.U2B2_EXACT5UMI_COLLAPSED_GEOMETRY_DENSITIES_v1.tsv"
        out.to_csv(outfile, sep="\t", index=False)
        if len(out) != len(cat) or out["geometry_id"].tolist() != cat["geometry_id"].tolist():
            raise RuntimeError(f"OUTPUT_CATALOG_ALIGNMENT_FAIL:{lib}")
        if not np.isfinite(out["density"].to_numpy(dtype=float)).all():
            raise RuntimeError(f"NONFINITE_DENSITY:{lib}")
        if (out["density"].to_numpy(dtype=float) < 0).any():
            raise RuntimeError(f"NEGATIVE_DENSITY:{lib}")

        emit(f"{lib} COUNTING_GATE = PASS")
        emit(f"{lib} OUTPUT_SHA256 =", sha256(outfile))
        summary_rows.append({
            "library": lib,
            "collapsed_bam": str(bam),
            "collapsed_bam_reads": input_reads,
            "geometry_rows": len(out),
            "sum_geometry_density": total_density,
            "output_tsv": str(outfile),
            "output_sha256": sha256(outfile),
        })

        del ga
        del position_counts
        gc.collect()

    summary = outdir / "U2B2A_COUNTING_SUMMARY_v1.tsv"
    pd.DataFrame(summary_rows).to_csv(summary, sep="\t", index=False)

    manifest = outdir / "U2B2A_SHA256_MANIFEST_v1.tsv"
    with manifest.open("w") as f:
        f.write("sha256\tpath\n")
        for p in sorted(countdir.glob("*.tsv")):
            f.write(f"{sha256(p)}\t{p}\n")
        f.write(f"{sha256(summary)}\t{summary}\n")

    emit("U2B2A_12_LIBRARY_COUNTING = PASS")
    emit("BIOLOGICAL_RESULT_CALCULATION = NO")
    emit("TOPOLOGY_INTERACTION_CALCULATION = NO")
    emit("U2B2A_STATUS = TECHNICAL_COUNTING_COMPLETE")
    emit("P3_ELIFE_U2B2A_COLLAPSED_GEOMETRY_COUNTING = END")
    (outdir / "U2B2A_RUN_LOG_v1.txt").write_text("\n".join(log)+"\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
