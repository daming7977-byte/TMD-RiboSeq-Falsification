#!/bin/bash
set -euo pipefail

BID="$1"
IN_REL="$2"

ROOT="/work"
REF="$ROOT/mapping_reference_blinded"

IN="$ROOT/$IN_REL"
NON="$ROOT/mapping_full_blinded/noncontam/${BID}.noncontam.fastq"

CLOG="$ROOT/mapping_full_blinded/logs/${BID}.contaminant.log"
TLOG="$ROOT/mapping_full_blinded/logs/${BID}.tophat.log"

OUT="$ROOT/mapping_full_blinded/tophat/${BID}"
TMP="$ROOT/mapping_full_blinded/tmp/${BID}"
REPORT="$ROOT/mapping_full_blinded/reports/${BID}.formal_mapping.txt"

CONTAM_IDX="$REF/bowtie1_index/contaminant/GRCh38p13_rRNA_tRNA"
GENOME_IDX="$REF/bowtie1_index/genome/GRCh38_release104_primary"
GTF="$REF/Homo_sapiens.GRCh38.104.gtf"
TX="$ROOT/mapping_pilot_blinded/transcriptome_index/GRCh38_104_tx"

mkdir -p \
"$ROOT/mapping_full_blinded/noncontam" \
"$ROOT/mapping_full_blinded/logs" \
"$ROOT/mapping_full_blinded/tophat" \
"$ROOT/mapping_full_blinded/reports"

rm -rf "$OUT" "$TMP"
mkdir -p "$TMP"
rm -f "$NON" "${NON}.gz" "$REPORT"

echo "===== ${BID}: FULL CONTAMINANT FILTER ====="

gzip -dc "$IN" \
| bowtie \
    -q \
    -p 4 \
    --un "$NON" \
    "$CONTAM_IDX" \
    - \
    /dev/null \
    2> "$CLOG"

INPUT_N=$(
    sed -n \
    's/^# reads processed:[[:space:]]*\([0-9,][0-9,]*\).*/\1/p' \
    "$CLOG" \
    | tail -n1 \
    | tr -d ','
)

CONTAM_N=$(
    sed -n \
    's/^# reads with at least one reported alignment:[[:space:]]*\([0-9,][0-9,]*\).*/\1/p' \
    "$CLOG" \
    | tail -n1 \
    | tr -d ','
)

LINES=$(wc -l < "$NON" | tr -d ' ')
NON_N=$((LINES / 4))

case "$INPUT_N" in
    ''|*[!0-9]*)
        echo "ERROR: invalid INPUT_N=[$INPUT_N]" >&2
        exit 1
        ;;
esac

case "$CONTAM_N" in
    ''|*[!0-9]*)
        echo "ERROR: invalid CONTAM_N=[$CONTAM_N]" >&2
        exit 1
        ;;
esac

case "$NON_N" in
    ''|*[!0-9]*)
        echo "ERROR: invalid NON_N=[$NON_N]" >&2
        exit 1
        ;;
esac

EXPECTED_NON=$((INPUT_N - CONTAM_N))

echo "INPUT_N = $INPUT_N"
echo "CONTAM_N = $CONTAM_N"
echo "NON_N = $NON_N"
echo "EXPECTED_NON_N = $EXPECTED_NON"

if [ "$NON_N" -ne "$EXPECTED_NON" ]; then
    echo "ERROR: contaminant read accounting failure" >&2
    exit 1
fi

echo "READ_ACCOUNTING = PASS"

echo "===== ${BID}: FULL TOPHAT ====="

tophat \
  -o "$OUT" \
  --tmp-dir "$TMP" \
  --keep-tmp \
  -p 4 \
  --bowtie1 \
  --read-mismatches 1 \
  -g 64 \
  --no-novel-juncs \
  -T \
  -G "$GTF" \
  --transcriptome-index "$TX" \
  "$GENOME_IDX" \
  "$NON" \
  > "$TLOG" 2>&1

BAM="$OUT/accepted_hits.bam"

test -s "$BAM"

/opt/conda/bin/samtools_0.1.18 view -H "$BAM" >/dev/null

ACCEPTED=$(
    /opt/conda/bin/samtools_0.1.18 view -c "$BAM"
)

NH_COUNTS=$(
    /opt/conda/bin/samtools_0.1.18 view "$BAM" \
    | awk '
    {
        nh=-1

        for(i=12;i<=NF;i++){
            if($i ~ /^NH:i:/){
                split($i,a,":")
                nh=a[3]
                break
            }
        }

        if(nh==1) u++
        else if(nh>1) m++
        else x++
    }

    END{
        print u+0, m+0, x+0
    }'
)

set -- $NH_COUNTS

NH1="$1"
NHMULTI="$2"
NHMISS="$3"

if [ "$NHMISS" -ne 0 ]; then
    echo "ERROR: missing NH tags" >&2
    exit 1
fi

gzip -n -f "$NON"

NON_SHA=$(
    sha256sum "${NON}.gz" | awk '{print $1}'
)

BAM_SHA=$(
    sha256sum "$BAM" | awk '{print $1}'
)

CONTAM_FRAC=$(
    awk -v a="$CONTAM_N" -v b="$INPUT_N" \
    'BEGIN{printf "%.12f",a/b}'
)

{
    echo "BLIND_ID = $BID"
    echo "FORMAL_INPUT_READS = $INPUT_N"
    echo "CONTAMINANT_ALIGNED_READS = $CONTAM_N"
    echo "NONCONTAM_READS = $NON_N"
    echo "CONTAMINANT_FRACTION = $CONTAM_FRAC"
    echo "TOPHAT_ACCEPTED_ALIGNMENTS = $ACCEPTED"
    echo "NH1_UNIQUE_ALIGNMENTS = $NH1"
    echo "NH_MULTI_ALIGNMENTS = $NHMULTI"
    echo "NH_MISSING = $NHMISS"
    echo "NONCONTAM_FASTQ_GZ_SHA256 = $NON_SHA"
    echo "ACCEPTED_HITS_BAM_SHA256 = $BAM_SHA"
    echo "BIOLOGICAL_IDENTITY_REVEAL = LOCKED"
    echo "FORMAL_MAPPING_COMPLETE = True"
} > "$REPORT"

rm -rf "$TMP"

echo "===== ${BID}: COMPLETE ====="
