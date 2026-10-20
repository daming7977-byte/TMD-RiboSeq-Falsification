#!/bin/bash
set -euo pipefail

BID="$1"

ROOT="$PAPER3_WORK_ROOT"
INDIR="$ROOT/formal_preprocessing_blinded/full_fastq"
REPORT="$ROOT/mapping_full_blinded/reports/${BID}.formal_mapping.txt"

if [ -s "$REPORT" ] && \
   grep -q '^FORMAL_MAPPING_COMPLETE = True$' "$REPORT"
then
    echo "$BID ALREADY_COMPLETE"
    exit 0
fi

N=$(
    find "$INDIR" \
    -maxdepth 1 \
    -type f \
    -name "${BID}.*.fastq.gz" \
    | wc -l \
    | tr -d ' '
)

if [ "$N" -ne 1 ]; then
    echo "$BID ERROR: expected one FASTQ, found $N" >&2
    exit 1
fi

IN=$(
    find "$INDIR" \
    -maxdepth 1 \
    -type f \
    -name "${BID}.*.fastq.gz" \
    | head -n1
)

IN_REL="${IN#$ROOT/}"

IMAGE_ID=$(
    sed 's/^IMAGE_ID=//' \
    "$ROOT/GSE300977_LINUX_AMD64_IMAGE_ID.txt"
)

docker run --rm \
  --platform linux/amd64 \
  -v "$ROOT:/work" \
  -w /work \
  "$IMAGE_ID" \
  bash \
  /work/mapping_full_blinded/container_full_one_v2.sh \
  "$BID" \
  "$IN_REL"

echo "$BID FORMAL_RUN_FINISHED"
