#!/usr/bin/env python3
"""
scripts/stage0c/00_freeze_combined_reference.py

Purpose:
    Construct and freeze a combined mapping reference (Human Transcripts + Viral Genome)
    with strict duplicate ID validation, feature table generation, and SHA-256 checksum recording.
"""

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


def calculate_sha256(filepath: Path) -> str:
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(65536), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def parse_fasta_headers(fasta_path: Path) -> set:
    headers = set()
    if not fasta_path.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {fasta_path}")
    
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                if seq_id in headers:
                    raise ValueError(f"Duplicate sequence ID detected inside {fasta_path}: {seq_id}")
                headers.add(seq_id)
    return headers


def main():
    print("[INFO] Starting Combined Reference Freezing Procedure...")

    repo_root = Path(__file__).resolve().parents[2]
    ref_dir = repo_root / "references" / "stage0c"
    output_dir = ref_dir / "frozen_combined_ref"
    output_dir.mkdir(parents=True, exist_ok=True)

    human_fasta = ref_dir / "human_19192_transcripts.fa"
    human_tsv = ref_dir / "human_features.tsv"
    viral_fasta = ref_dir / "viral_genome.fa"

    out_combined_fa = output_dir / "combined_transcripts.fa"
    out_combined_tsv = output_dir / "combined_features.tsv"
    out_manifest = output_dir / "reference_manifest.json"

    for req_file in [human_fasta, human_tsv, viral_fasta]:
        if not req_file.exists():
            print(f"[ERROR] Required reference file missing: {req_file}")
            print("Please place source FASTA/TSV files in 'references/stage0c/' before running.")
            sys.exit(1)

    print("[STEP 1] Validating FASTA headers for duplicate IDs...")
    human_ids = parse_fasta_headers(human_fasta)
    viral_ids = parse_fasta_headers(viral_fasta)

    collision = human_ids.intersection(viral_ids)
    if collision:
        print(f"[FATAL] Sequence ID collision detected: {collision}")
        sys.exit(1)
    print(f"[PASS] Unique Human IDs: {len(human_ids)}, Unique Viral IDs: {len(viral_ids)}.")

    print("[STEP 2] Merging FASTA reference files...")
    with open(out_combined_fa, "w", encoding="utf-8") as out_fa:
        for fpath in [human_fasta, viral_fasta]:
            with open(fpath, "r", encoding="utf-8") as in_fa:
                out_fa.write(in_fa.read().strip() + "\n")

    print("[STEP 3] Generating unified feature table...")
    with open(out_combined_tsv, "w", encoding="utf-8") as out_tsv:
        with open(human_tsv, "r", encoding="utf-8") as in_tsv:
            out_tsv.write(in_tsv.read().strip() + "\n")
        
        for v_id in viral_ids:
            out_tsv.write(f"{v_id}\t96\t10273\t10723\n")

    print("[STEP 4] Computing SHA-256 checksums and writing manifest...")
    manifest_data = {
        "metadata": {
            "project": "GSE133111_Stage0C",
            "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
            "human_transcript_count": len(human_ids),
            "viral_sequence_count": len(viral_ids)
        },
        "input_sources": {
            "human_fasta_sha256": calculate_sha256(human_fasta),
            "viral_fasta_sha256": calculate_sha256(viral_fasta),
        },
        "frozen_outputs": {
            "combined_transcripts_fa_sha256": calculate_sha256(out_combined_fa),
            "combined_features_tsv_sha256": calculate_sha256(out_combined_tsv),
        }
    }

    with open(out_manifest, "w", encoding="utf-8") as mf:
        json.dump(manifest_data, mf, indent=2)

    print(f"\n[SUCCESS] Combined Mapping Reference Frozen Successfully!")


if __name__ == "__main__":
    main()
