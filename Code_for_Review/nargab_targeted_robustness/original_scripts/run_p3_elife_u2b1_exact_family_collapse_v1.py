#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, hashlib, re
import pysam

EXPECTED_TOTAL = 215_999_661
UMI_RE = re.compile(r'(?:^|\|)UMI5=([ACGTN]{5})(?:\||$)')

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

def get_umi(qname):
    m = UMI_RE.search(qname or "")
    return m.group(1) if m else None

def bio5p(read):
    return (read.reference_end - 1) if read.is_reverse else read.reference_start

def finalize_group(fam_counts, stats):
    for n in fam_counts.values():
        stats["exact_valid_umi_families"] += 1
        if n >= 2:
            stats["duplicate_families"] += 1
        if n > stats["max_family_size"]:
            stats["max_family_size"] = n

def process_bam(inp_path, out_path):
    stats = {
        "input_reads": 0,
        "valid_acgt_umi_reads": 0,
        "ambiguous_N_umi_reads": 0,
        "missing_or_malformed_umi_reads": 0,
        "exact_valid_umi_families": 0,
        "duplicate_families": 0,
        "max_family_size": 0,
        "collapsed_reads": 0,
        "output_reads": 0,
        "unmapped_reads": 0,
        "secondary_reads": 0,
        "supplementary_reads": 0,
    }

    with pysam.AlignmentFile(inp_path, "rb") as inp:
        so = (inp.header.to_dict().get("HD") or {}).get("SO", "")
        if so != "coordinate":
            raise RuntimeError(f"INPUT_NOT_COORDINATE_SORTED: {inp_path} SO={so!r}")

        out_path.parent.mkdir(parents=True, exist_ok=True)
        with pysam.AlignmentFile(out_path, "wb", template=inp) as out:
            current_locus = None
            seen = set()
            fam_counts = {}
            last_ref = -1
            last_start = -1

            for read in inp.fetch(until_eof=True):
                stats["input_reads"] += 1

                if read.is_unmapped:
                    stats["unmapped_reads"] += 1
                if read.is_secondary:
                    stats["secondary_reads"] += 1
                if read.is_supplementary:
                    stats["supplementary_reads"] += 1

                if not read.is_unmapped:
                    rid = read.reference_id
                    rst = read.reference_start
                    if rid < last_ref or (rid == last_ref and rst < last_start):
                        raise RuntimeError(
                            f"COORDINATE_ORDER_VIOLATION: {inp_path} "
                            f"{rid}:{rst} after {last_ref}:{last_start}"
                        )
                    last_ref, last_start = rid, rst

                umi = get_umi(read.query_name)

                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    if umi is None:
                        stats["missing_or_malformed_umi_reads"] += 1
                    elif "N" in umi:
                        stats["ambiguous_N_umi_reads"] += 1
                    else:
                        stats["valid_acgt_umi_reads"] += 1
                    out.write(read)
                    stats["output_reads"] += 1
                    continue

                locus = (read.reference_id, read.reference_start)
                if current_locus is None:
                    current_locus = locus
                elif locus != current_locus:
                    finalize_group(fam_counts, stats)
                    seen.clear()
                    fam_counts.clear()
                    current_locus = locus

                if umi is None:
                    stats["missing_or_malformed_umi_reads"] += 1
                    out.write(read)
                    stats["output_reads"] += 1
                    continue

                if "N" in umi:
                    stats["ambiguous_N_umi_reads"] += 1
                    out.write(read)
                    stats["output_reads"] += 1
                    continue

                stats["valid_acgt_umi_reads"] += 1
                qlen = read.query_length if read.query_length is not None else -1
                qaln = read.query_alignment_length if read.query_alignment_length is not None else -1

                key = (
                    read.reference_id,
                    bool(read.is_reverse),
                    bio5p(read),
                    read.reference_start,
                    read.reference_end,
                    read.cigarstring or "",
                    qlen,
                    qaln,
                    umi,
                )

                fam_counts[key] = fam_counts.get(key, 0) + 1
                if key not in seen:
                    seen.add(key)
                    out.write(read)
                    stats["output_reads"] += 1
                else:
                    stats["collapsed_reads"] += 1

            finalize_group(fam_counts, stats)

    if stats["output_reads"] + stats["collapsed_reads"] != stats["input_reads"]:
        raise RuntimeError(
            f"ACCOUNTING_FAILURE: {inp_path} "
            f"out={stats['output_reads']} collapsed={stats['collapsed_reads']} "
            f"in={stats['input_reads']}"
        )

    pysam.index(str(out_path))
    return stats

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

    bamdir = root / "unique_mapping_blinded" / "bam"
    inputs = [bamdir / f"B{i:02d}.NH1.unique.bam" for i in range(1, 13)]
    missing = [p for p in inputs if not p.exists()]
    if missing:
        print("STOP_MISSING_INPUT_BAMS")
        for p in missing:
            print(p)
        return 3

    log = []
    def emit(*x):
        s = " ".join(map(str, x))
        print(s)
        log.append(s)

    emit("P3_ELIFE_U2B1_EXACT_FAMILY_COLLAPSE = START")
    emit("ROOT =", root)
    emit("OUTDIR =", outdir)
    emit("CONTRACT =", contract)
    emit("CONTRACT_SHA256 =", sha256(contract))
    emit("EXPECTED_TOTAL_INPUT_READS =", EXPECTED_TOTAL)
    emit("INPUT_BAM_CLASS = NH:i:1 unique BAM")
    emit("UMI_FIELD = UMI5")
    emit("UMI_ERROR_CORRECTION = NO")
    emit("BC4_FILTERING = NO")
    emit("N_CONTAINING_UMI_COLLAPSE = NO")
    emit("BIOLOGICAL_RESULT_CALCULATION = NO")

    rows = []
    total_in = total_out = total_collapsed = 0
    qname_gate = True
    output_gate = True
    out_bamdir = outdir / "bam"
    out_bamdir.mkdir()

    for i, inp in enumerate(inputs, start=1):
        out = out_bamdir / f"B{i:02d}.NH1.unique.exact5umi_collapse.bam"
        stats = process_bam(inp, out)
        total_in += stats["input_reads"]
        total_out += stats["output_reads"]
        total_collapsed += stats["collapsed_reads"]

        if stats["missing_or_malformed_umi_reads"] != 0:
            qname_gate = False
        if stats["output_reads"] + stats["collapsed_reads"] != stats["input_reads"]:
            output_gate = False
        if not Path(str(out) + ".bai").exists():
            output_gate = False

        row = {
            "library": f"B{i:02d}",
            "input_bam": str(inp),
            "output_bam": str(out),
            **stats,
            "retention_fraction": stats["output_reads"] / stats["input_reads"] if stats["input_reads"] else float("nan"),
            "collapsed_fraction": stats["collapsed_reads"] / stats["input_reads"] if stats["input_reads"] else float("nan"),
            "output_bai_exists": "YES" if Path(str(out) + ".bai").exists() else "NO",
        }
        rows.append(row)

        emit(
            f"LIBRARY=B{i:02d}",
            f"INPUT={stats['input_reads']}",
            f"VALID_UMI={stats['valid_acgt_umi_reads']}",
            f"N_UMI={stats['ambiguous_N_umi_reads']}",
            f"MALFORMED_UMI={stats['missing_or_malformed_umi_reads']}",
            f"FAMILIES={stats['exact_valid_umi_families']}",
            f"DUP_FAMILIES={stats['duplicate_families']}",
            f"MAX_FAMILY={stats['max_family_size']}",
            f"COLLAPSED={stats['collapsed_reads']}",
            f"OUTPUT={stats['output_reads']}",
            f"RETENTION={row['retention_fraction']:.6f}",
        )

    count_gate = (total_in == EXPECTED_TOTAL)

    summary = outdir / "U2B1_EXACT_FAMILY_COLLAPSE_SUMMARY_v1.tsv"
    with summary.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    verdict = outdir / "U2B1_TECHNICAL_VERDICT_v1.txt"
    verdict.write_text("\n".join([
        f"U2B1_TOTAL_INPUT_READS = {total_in}",
        f"U2B1_EXPECTED_TOTAL_INPUT_READS = {EXPECTED_TOTAL}",
        f"U2B1_TOTAL_OUTPUT_READS = {total_out}",
        f"U2B1_TOTAL_COLLAPSED_READS = {total_collapsed}",
        f"U2B1_OVERALL_RETENTION_FRACTION = {total_out/total_in:.12f}",
        f"U2B1_OVERALL_COLLAPSED_FRACTION = {total_collapsed/total_in:.12f}",
        f"U2B1_INPUT_COUNT_GATE = {'PASS' if count_gate else 'FAIL'}",
        f"U2B1_QNAME_UMI_GATE = {'PASS' if qname_gate else 'FAIL'}",
        f"U2B1_OUTPUT_ACCOUNTING_GATE = {'PASS' if output_gate else 'FAIL'}",
        "BIOLOGICAL_RESULT_CALCULATION = NO",
        "RESULT_DEPENDENT_RESCUE = NO",
        f"U2B1_TECHNICAL_STATUS = {'PASS' if (count_gate and qname_gate and output_gate) else 'FAIL'}",
    ]) + "\n")

    emit("U2B1_TOTAL_INPUT_READS =", total_in)
    emit("U2B1_TOTAL_OUTPUT_READS =", total_out)
    emit("U2B1_TOTAL_COLLAPSED_READS =", total_collapsed)
    emit("U2B1_INPUT_COUNT_GATE =", "PASS" if count_gate else "FAIL")
    emit("U2B1_QNAME_UMI_GATE =", "PASS" if qname_gate else "FAIL")
    emit("U2B1_OUTPUT_ACCOUNTING_GATE =", "PASS" if output_gate else "FAIL")
    emit("U2B1_TECHNICAL_STATUS =", "PASS" if (count_gate and qname_gate and output_gate) else "FAIL")
    emit("BIOLOGICAL_RESULT_CALCULATION = NO")
    emit("P3_ELIFE_U2B1_EXACT_FAMILY_COLLAPSE = END")

    (outdir / "U2B1_RUN_LOG_v1.txt").write_text("\n".join(log) + "\n")
    return 0 if (count_gate and qname_gate and output_gate) else 5

if __name__ == "__main__":
    raise SystemExit(main())
