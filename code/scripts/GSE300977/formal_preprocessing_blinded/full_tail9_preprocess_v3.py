#!/usr/bin/env python3

import csv
import gzip
import hashlib
import io
from pathlib import Path
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT = Path.cwd()

MANIFEST = ROOT / "raw_architecture_blinded/BLINDED_RUN_MANIFEST.tsv"
PLAN = ROOT / "GSE300977.DOWNLOAD_PLAN.tsv"
RAW_DIR = ROOT / "ribo_fastq"

OUTDIR = ROOT / "formal_preprocessing_blinded/full_fastq"
REPORTDIR = ROOT / "formal_preprocessing_blinded/run_reports"

TAIL = 9
MIN_INSERT = 18
N_WORKERS = 3


def sha256_file(path):
    h = hashlib.sha256()

    with path.open("rb") as fh:
        for chunk in iter(
            lambda: fh.read(16 * 1024 * 1024),
            b""
        ):
            h.update(chunk)

    return h.hexdigest()


def deterministic_gzip_text_writer(path):
    """
    gzip mtime=0 => reproducible compressed output SHA256.
    """
    raw = open(path, "wb")

    gz = gzip.GzipFile(
        filename="",
        mode="wb",
        fileobj=raw,
        compresslevel=6,
        mtime=0
    )

    txt = io.TextIOWrapper(
        gz,
        encoding="ascii",
        newline="\n"
    )

    return raw, gz, txt


def process_one(r, expected):

    bid = r["blind_id"]
    run = r["run_accession"]

    raw = RAW_DIR / r["filename"]

    if not raw.exists():
        hits = list(
            RAW_DIR.glob(f"{run}*.fastq.gz")
        )

        if len(hits) != 1:
            raise RuntimeError(
                f"{bid}/{run}: raw FASTQ not uniquely found; "
                f"hits={len(hits)}"
            )

        raw = hits[0]

    out = (
        OUTDIR /
        f"{bid}.{run}.tail9_m18.fastq.gz"
    )

    raw_n = 0
    passed_n = 0
    removed_lt18 = 0

    lengths = Counter()
    barcodes = Counter()
    umis = Counter()

    umi_with_n = 0
    barcode_with_n = 0

    rawout, gzout, dst = (
        deterministic_gzip_text_writer(out)
    )

    try:

        with gzip.open(
            raw,
            "rt",
            encoding="ascii",
            errors="strict"
        ) as inp:

            while True:

                header = inp.readline()

                if not header:
                    break

                seq = inp.readline().rstrip("\n")
                plus = inp.readline()
                qual = inp.readline().rstrip("\n")

                if not plus or qual == "":
                    raise RuntimeError(
                        f"{bid}/{run}: truncated FASTQ "
                        f"at read {raw_n + 1}"
                    )

                raw_n += 1

                if len(seq) != len(qual):
                    raise RuntimeError(
                        f"{bid}/{run}: sequence/quality "
                        f"length mismatch at read {raw_n}"
                    )

                if len(seq) < TAIL:
                    removed_lt18 += 1
                    continue

                umi = seq[-9:-4].upper()
                bc = seq[-4:].upper()

                insert = seq[:-9]
                insert_qual = qual[:-9]

                umis[umi] += 1
                barcodes[bc] += 1

                if "N" in umi:
                    umi_with_n += 1

                if "N" in bc:
                    barcode_with_n += 1

                if len(insert) < MIN_INSERT:
                    removed_lt18 += 1
                    continue

                passed_n += 1
                lengths[len(insert)] += 1

                qname = header.strip().split()[0]

                if qname.startswith("@"):
                    qname = qname[1:]

                dst.write(
                    f"@{qname}|UMI5={umi}|BC4={bc}\n"
                )
                dst.write(insert + "\n")
                dst.write("+\n")
                dst.write(insert_qual + "\n")

    finally:
        dst.flush()
        dst.close()

    if raw_n != expected:
        raise RuntimeError(
            f"{bid}/{run}: RAW_READS={raw_n} "
            f"!= ENA_EXPECTED={expected}"
        )

    if raw_n != passed_n + removed_lt18:
        raise RuntimeError(
            f"{bid}/{run}: read accounting failure"
        )

    if not lengths:
        raise RuntimeError(
            f"{bid}/{run}: no passing reads"
        )

    mode_len, mode_n = (
        lengths.most_common(1)[0]
    )

    top_bc, top_bc_n = (
        barcodes.most_common(1)[0]
    )

    barcode_total = sum(barcodes.values())
    umi_total = sum(umis.values())

    result = {
        "blind_id": bid,
        "run_accession": run,

        "raw_reads": raw_n,
        "expected_raw_reads": expected,

        "removed_lt18": removed_lt18,
        "passing_reads": passed_n,
        "passing_fraction":
            passed_n / raw_n,

        "top_barcode": top_bc,
        "top_barcode_fraction":
            top_bc_n / barcode_total,

        "unique_umi5_strings":
            len(umis),

        "umi_with_N_fraction":
            umi_with_n / umi_total,

        "barcode_with_N_fraction":
            barcode_with_n / barcode_total,

        "min_insert_length":
            min(lengths),

        "max_insert_length":
            max(lengths),

        "mode_insert_length":
            mode_len,

        "mode_fraction":
            mode_n / passed_n,

        "fraction_27_34":
            sum(
                v for L, v in lengths.items()
                if 27 <= L <= 34
            ) / passed_n,

        "fraction_28_34":
            sum(
                v for L, v in lengths.items()
                if 28 <= L <= 34
            ) / passed_n,

        "fraction_25_35":
            sum(
                v for L, v in lengths.items()
                if 25 <= L <= 35
            ) / passed_n,

        "output_fastq":
            str(out),

        "output_bytes":
            out.stat().st_size,

        "output_sha256":
            sha256_file(out)
    }

    rep = (
        REPORTDIR /
        f"{bid}.{run}.report.tsv"
    )

    with rep.open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=result.keys(),
            delimiter="\t"
        )

        w.writeheader()
        w.writerow(result)

    return result


def main():

    OUTDIR.mkdir(
        parents=True,
        exist_ok=True
    )

    REPORTDIR.mkdir(
        parents=True,
        exist_ok=True
    )

    expected_reads = {}

    with PLAN.open() as fh:

        for r in csv.DictReader(
            fh,
            delimiter="\t"
        ):
            expected_reads[
                r["run_accession"]
            ] = int(r["read_count"])


    with MANIFEST.open() as fh:

        manifest = list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


    print("=" * 76, flush=True)
    print(
        "GSE300977 FULL FORMAL PREPROCESSING v3",
        flush=True
    )
    print("=" * 76, flush=True)

    print(
        "RUNS_EXPECTED =",
        len(manifest),
        flush=True
    )

    print(
        "TAIL_REMOVAL_NT = 9",
        flush=True
    )

    print(
        "MIN_INSERT_NT = 18",
        flush=True
    )

    print(
        "DEDUPLICATION = NONE",
        flush=True
    )

    print(
        "BIOLOGICAL_IDENTITY_REVEAL = LOCKED",
        flush=True
    )

    results = []

    with ProcessPoolExecutor(
        max_workers=N_WORKERS
    ) as ex:

        futures = {}

        for r in manifest:

            run = r["run_accession"]

            if run not in expected_reads:
                raise SystemExit(
                    f"Missing expected read count: {run}"
                )

            fut = ex.submit(
                process_one,
                r,
                expected_reads[run]
            )

            futures[fut] = r


        for fut in as_completed(futures):

            r = futures[fut]

            try:
                x = fut.result()

            except Exception as e:
                raise SystemExit(
                    f"FAIL: "
                    f"{r['blind_id']} "
                    f"{r['run_accession']}: "
                    f"{e}"
                )

            results.append(x)

            print(
                "COMPLETE",
                x["blind_id"],
                x["run_accession"],
                "RAW="+str(x["raw_reads"]),
                "PASS="+str(x["passing_reads"]),
                "PASS_FRAC=%.6f"
                    % x["passing_fraction"],
                "BC="+x["top_barcode"],
                "BC_FRAC=%.6f"
                    % x["top_barcode_fraction"],
                "MODE="+str(
                    x["mode_insert_length"]
                ),
                "F27_34=%.6f"
                    % x["fraction_27_34"],
                flush=True
            )


    results.sort(
        key=lambda x:x["blind_id"]
    )


    if len(results) != 12:
        raise SystemExit(
            f"FAIL: processed {len(results)} runs"
        )


    summary = (
        ROOT /
        "formal_preprocessing_blinded/"
        "FORMAL_PREPROCESSING_SUMMARY_v3.tsv"
    )


    with summary.open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=results[0].keys(),
            delimiter="\t"
        )

        w.writeheader()
        w.writerows(results)


    total_raw = sum(
        x["raw_reads"]
        for x in results
    )

    total_pass = sum(
        x["passing_reads"]
        for x in results
    )

    total_removed = sum(
        x["removed_lt18"]
        for x in results
    )


    print()
    print("=" * 76)
    print(
        "GSE300977 FORMAL PREPROCESSING SUMMARY v3"
    )
    print("=" * 76)

    print(
        "RUNS_PROCESSED =",
        len(results)
    )

    print(
        "TOTAL_RAW_READS =",
        total_raw
    )

    print(
        "TOTAL_PASSING_READS =",
        total_pass
    )

    print(
        "TOTAL_REMOVED_LT18 =",
        total_removed
    )

    print(
        "OVERALL_PASS_FRACTION =",
        total_pass / total_raw
    )

    print()

    for x in results:

        print(
            x["blind_id"],
            "RAW="+str(x["raw_reads"]),
            "PASS="+str(x["passing_reads"]),
            "PASS_FRAC=%.6f"
                % x["passing_fraction"],
            "BC="+x["top_barcode"],
            "BC_FRAC=%.6f"
                % x["top_barcode_fraction"],
            "UMI_N_FRAC=%.8f"
                % x["umi_with_N_fraction"],
            "MODE="+str(
                x["mode_insert_length"]
            ),
            "F27_34=%.6f"
                % x["fraction_27_34"],
            sep="\t"
        )


    accounting = all(
        x["raw_reads"]
        ==
        x["passing_reads"]
        +
        x["removed_lt18"]
        for x in results
    )

    minpass = all(
        x["min_insert_length"]
        >= MIN_INSERT
        for x in results
    )


    print()
    print(
        "READ_ACCOUNTING_PASS =",
        accounting
    )

    print(
        "MIN_LENGTH_PASS =",
        minpass
    )

    print(
        "FORMAL_PREPROCESSING_COMPLETE =",
        accounting and minpass
    )

    print(
        "DEDUPLICATION_PERFORMED = False"
    )

    print(
        "MAPPING_PERFORMED = False"
    )

    print(
        "PERIODICITY_VIEWED = False"
    )

    print(
        "BIOLOGICAL_IDENTITY_REVEAL = LOCKED"
    )


if __name__ == "__main__":
    main()
