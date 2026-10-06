#!/bin/bash
set -e

cd "$PAPER3_WORK_ROOT"

BT112="REPRO_ENV/envs/gse300977-tophat211-x86/bin/bowtie"

IDX="mapping_reference_blinded/bowtie1_index/genome/GRCh38_release104_primary"

TWOBIT="periodicity_blinded/crossmap/GRCh38_release104_primary.2bit"

OUT="periodicity_blinded/crossmap/GRCh38_release104"

LOG="periodicity_blinded/logs/CROSSMAP_K29_MM1.console.txt"

echo "CROSSMAP_START"
date

crossmap \
  "$IDX" \
  "$OUT" \
  --sequence_file "$TWOBIT" \
  --sequence_format twobit \
  -k 29 \
  --mismatches 1 \
  --bowtie "$BT112" \
  -p 1 \
  > "$LOG" 2>&1

RC=$?

echo "CROSSMAP_EXIT = $RC" \
> periodicity_blinded/crossmap/CROSSMAP_EXIT_STATUS.txt

date >> periodicity_blinded/crossmap/CROSSMAP_EXIT_STATUS.txt

exit "$RC"
