#!/bin/bash

set -uo pipefail
export MPLBACKEND=Agg
export LC_ALL=C

ROOT="$PWD"

METAGENE="REPRO_ENV/envs/gse300977-plastid061-x86/bin/metagene"
PYTHON="REPRO_ENV/envs/gse300977-plastid061-x86/bin/python"

ROI="periodicity_blinded/roi_masked_v2/GRCh38_104_cds_start_crossmap_rois.txt"
OUTDIR="formal_counting_blinded/counts_v1"
LOGDIR="formal_counting_blinded/logs"
AUDIT="formal_counting_blinded/audit"

EXPECTED_ROI_SHA="0f89aee14812df3e4a79c1173c90eeaf82d81b1f5268962f08de9cb98918b0cc"
EXPECTED_FREEZE_SHA="98d790111f7f3984fb2eac4826c54ca30589c4ef97396d296b788b5f8d830e55"
EXPECTED_MASKED_ROI_STAGE_SHA="e6357cc2bb2ccdca0f461c3e527f92b46a76f7bf5429dbdee5cf831ba3c8b271"

mkdir -p "$OUTDIR" "$LOGDIR" "$AUDIT"

STATUS="$AUDIT/FORMAL_BLINDED_COUNTING_STATUS_v1.tsv"
: > "$STATUS"

echo "FORMAL_BLINDED_COUNTING_START"
date

############################################################
# 0. AUDIT THE FAILED SOURCE-INSPECTION ATTEMPT
############################################################

cat > "$AUDIT/CENTERMAPFACTORY_IMPLEMENTATION_v1_SUPERSEDED.txt" <<'EOF'
CENTERMAPFACTORY SOURCE-INSPECTION ATTEMPT

STATUS = SUPERSEDED BEFORE FORMAL COUNTING

The prior inspect.getsource(CenterMapFactory) attempt failed because
the installed implementation is compiled and source code was not
available through Python inspect.

The resulting file was empty and its SHA256 was therefore the SHA256
of an empty file:

e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855

This empty-file SHA is NOT implementation evidence and is not used
for any scientific decision.

No biological result had been viewed.
No formal counting had been executed.
EOF

############################################################
# 1. RECORD LOCAL CENTERMAPFACTORY DOCSTRING
############################################################

"$PYTHON" - <<'PY' \
> "$AUDIT/CENTERMAPFACTORY_LOCAL_DOCSTRING_v1.txt"
from plastid import CenterMapFactory
print(CenterMapFactory.__doc__)
print()
print("OBJECT =", CenterMapFactory)
print("NIBBLE_TEST_OBJECT =", CenterMapFactory(nibble=13))
PY

if [ $? -ne 0 ]; then
    echo "CENTERMAP_DOCSTRING_GATE = FAIL"
    exit 10
fi

echo "CENTERMAP_DOCSTRING_GATE = PASS"

############################################################
# 2. HARD INPUT SHA GATES
############################################################

ACTUAL_ROI_SHA=$(shasum -a 256 "$ROI" | awk '{print $1}')

if [ "$ACTUAL_ROI_SHA" != "$EXPECTED_ROI_SHA" ]; then
    echo "ROI_SHA_GATE = FAIL"
    echo "EXPECTED = $EXPECTED_ROI_SHA"
    echo "ACTUAL   = $ACTUAL_ROI_SHA"
    exit 11
fi

echo "ROI_SHA_GATE = PASS"

ACTUAL_FREEZE_SHA=$(
    shasum -a 256 \
    "$AUDIT/FORMAL_COUNTING_INPUT_FREEZE_v2.txt" \
    | awk '{print $1}'
)

if [ "$ACTUAL_FREEZE_SHA" != "$EXPECTED_FREEZE_SHA" ]; then
    echo "COUNTING_FREEZE_SHA_GATE = FAIL"
    exit 12
fi

echo "COUNTING_FREEZE_SHA_GATE = PASS"

ACTUAL_MASKED_STAGE_SHA=$(
    shasum -a 256 \
    "$AUDIT/MASKED_CDS_START_ROI_AUDIT_v2.txt" \
    | awk '{print $1}'
)

if [ "$ACTUAL_MASKED_STAGE_SHA" != "$EXPECTED_MASKED_ROI_STAGE_SHA" ]; then
    echo "MASKED_ROI_STAGE_SHA_GATE = FAIL"
    exit 13
fi

echo "MASKED_ROI_STAGE_SHA_GATE = PASS"

############################################################
# 3. BAM + INDEX GATES
############################################################

for i in 1 2 3 4 5 6 7 8 9 10 11 12; do

    B=$(printf 'B%02d' "$i")
    BAM="unique_mapping_blinded/bam/${B}.NH1.unique.bam"

    if [ ! -s "$BAM" ]; then
        echo "$B BAM_GATE = FAIL"
        exit 20
    fi

    if [ ! -s "${BAM}.bai" ]; then
        echo "$B BAM_INDEX_GATE = FAIL"
        exit 21
    fi
done

echo "BAM_AND_INDEX_GATE = 12/12 PASS"

############################################################
# 4. FREEZE COUNTING COMMAND CONTRACT
############################################################

cat > "$AUDIT/FORMAL_BLINDED_COUNTING_CONTRACT_v1.txt" <<EOF
GSE300977 FORMAL BLINDED CDS-START COUNTING CONTRACT V1

DATE = $(date)

BIOLOGICAL_IDENTITY = LOCKED

LIBRARIES = B01-B12
INPUT = NH:i:1 unique BAM
ROI = $ROI
ROI_SHA256 = $EXPECTED_ROI_SHA

READ_LENGTH_MIN = 27
READ_LENGTH_MAX = 36

COUNTFILE_FORMAT = BAM
MAPPING_RULE = CenterMapFactory
CENTER = TRUE
NIBBLE = 13

NORMALIZE_OVER = 20 50
MIN_COUNTS = 10

NOTE:
normalize_over and min_counts affect auxiliary normalized/metagene
outputs only. Formal engineering validation is based on raw count and
mask matrices generated with --keep.

EXPECTED_ROIS = 5130
EXPECTED_WINDOW_WIDTH = 100

BIOLOGICAL_REVEAL = PROHIBITED DURING THIS STAGE
EOF

shasum -a 256 \
"$AUDIT/FORMAL_BLINDED_COUNTING_CONTRACT_v1.txt" \
> "$AUDIT/FORMAL_BLINDED_COUNTING_CONTRACT_v1.sha256"

############################################################
# 5. FUNCTION: RUN ONE BLINDED LIBRARY
############################################################

run_one () {

    B="$1"
    BAM="unique_mapping_blinded/bam/${B}.NH1.unique.bam"
    BASE="$OUTDIR/${B}.cds_start_center_n13"
    LOG="$LOGDIR/${B}.formal_counting.log"

    echo
    echo "============================================================"
    echo "START $B"
    date
    echo "============================================================"

    "$METAGENE" count \
        "$ROI" \
        "$BASE" \
        --count_files "$BAM" \
        --countfile_format BAM \
        --min_length 27 \
        --max_length 36 \
        --center \
        --nibble 13 \
        --normalize_over 20 50 \
        --min_counts 10 \
        --keep \
        > "$LOG" 2>&1

    RC=$?

    if [ "$RC" -ne 0 ]; then
        printf "%s\tFAIL\tMETAGENE_RC=%s\n" "$B" "$RC" \
        | tee -a "$STATUS"

        echo "$B COUNTING_GATE = FAIL"
        tail -50 "$LOG"
        return 1
    fi

    ########################################################
    # Find exact output names (.txt or .txt.gz)
    ########################################################

    RAW=$(find "$OUTDIR" -maxdepth 1 -type f \
        \( -name "${B}.cds_start_center_n13_rawcounts.txt" \
        -o -name "${B}.cds_start_center_n13_rawcounts.txt.gz" \) \
        -print -quit)

    NORM=$(find "$OUTDIR" -maxdepth 1 -type f \
        \( -name "${B}.cds_start_center_n13_normcounts.txt" \
        -o -name "${B}.cds_start_center_n13_normcounts.txt.gz" \) \
        -print -quit)

    MASK=$(find "$OUTDIR" -maxdepth 1 -type f \
        \( -name "${B}.cds_start_center_n13_mask.txt" \
        -o -name "${B}.cds_start_center_n13_mask.txt.gz" \) \
        -print -quit)

    PROFILE=$(find "$OUTDIR" -maxdepth 1 -type f \
        -name "${B}.cds_start_center_n13_metagene_profile.txt*" \
        -print -quit)

    for F in "$RAW" "$NORM" "$MASK" "$PROFILE"; do
        if [ -z "$F" ] || [ ! -s "$F" ]; then
            printf "%s\tFAIL\tMISSING_OUTPUT\n" "$B" \
            | tee -a "$STATUS"

            echo "$B OUTPUT_GATE = FAIL"
            return 1
        fi
    done

    ########################################################
    # Structural matrix validation
    ########################################################

    "$PYTHON" - "$RAW" "$NORM" "$MASK" "$B" <<'PY'
import sys
import numpy as np

raw_path, norm_path, mask_path, lib = sys.argv[1:]

def load(path):
    x = np.loadtxt(path)
    if x.ndim == 1:
        x = x.reshape(1, -1)
    return x

raw  = load(raw_path)
norm = load(norm_path)
mask = load(mask_path)

print("%s_RAW_SHAPE =" % lib, raw.shape)
print("%s_NORM_SHAPE =" % lib, norm.shape)
print("%s_MASK_SHAPE =" % lib, mask.shape)

expected = (5130, 100)

if raw.shape != expected:
    raise SystemExit("%s RAW_SHAPE_FAIL" % lib)

if norm.shape != expected:
    raise SystemExit("%s NORM_SHAPE_FAIL" % lib)

if mask.shape != expected:
    raise SystemExit("%s MASK_SHAPE_FAIL" % lib)

if raw.shape != norm.shape or raw.shape != mask.shape:
    raise SystemExit("%s MATRIX_CONCORDANCE_FAIL" % lib)

print("%s_MATRIX_GATE = PASS" % lib)
PY

    VRC=$?

    if [ "$VRC" -ne 0 ]; then
        printf "%s\tFAIL\tMATRIX_GATE\n" "$B" \
        | tee -a "$STATUS"

        return 1
    fi

    ########################################################
    # Seal outputs
    ########################################################

    {
        shasum -a 256 "$RAW"
        shasum -a 256 "$NORM"
        shasum -a 256 "$MASK"
        shasum -a 256 "$PROFILE"
    } > "$AUDIT/${B}.FORMAL_COUNTING_OUTPUT_SHA256.txt"

    printf "%s\tPASS\tRC=0\n" "$B" \
    | tee -a "$STATUS"

    echo "$B FORMAL_COUNTING_GATE = PASS"
    date

    return 0
}

############################################################
# 6. B01 ENGINEERING GATE
############################################################

echo
echo "=== B01 ENGINEERING PILOT ==="

run_one B01

if [ $? -ne 0 ]; then
    echo "B01_ENGINEERING_GATE = FAIL"
    echo "B02-B12 WILL NOT RUN"
    exit 30
fi

echo "B01_ENGINEERING_GATE = PASS"

############################################################
# 7. ONLY AFTER B01 PASS: RUN B02-B12 SEQUENTIALLY
############################################################

for i in 2 3 4 5 6 7 8 9 10 11 12; do

    B=$(printf 'B%02d' "$i")

    run_one "$B"

    if [ $? -ne 0 ]; then
        echo "$B FAILED"
        echo "QUEUE STOPPED BEFORE NEXT LIBRARY"
        exit 31
    fi
done

############################################################
# 8. FINAL GATE
############################################################

PASS_COUNT=$(awk '$2=="PASS"{n++} END{print n+0}' "$STATUS")

echo
echo "PASS_COUNT = $PASS_COUNT"

if [ "$PASS_COUNT" -ne 12 ]; then
    echo "GSE300977_FORMAL_BLINDED_COUNTING_GATE = FAIL"
    exit 40
fi

shasum -a 256 "$STATUS" \
> "$AUDIT/FORMAL_BLINDED_COUNTING_STATUS_v1.sha256"

{
    echo "GSE300977 FORMAL BLINDED COUNTING FINAL AUDIT"
    echo "DATE = $(date)"
    echo
    echo "B01_ENGINEERING_GATE = PASS"
    echo "LIBRARIES_PASS = 12"
    echo "LIBRARIES_EXPECTED = 12"
    echo
    echo "BIOLOGICAL_IDENTITY = LOCKED"
    echo
    echo "COUNTING_CONTRACT_SHA256:"
    cat "$AUDIT/FORMAL_BLINDED_COUNTING_CONTRACT_v1.sha256"
    echo
    echo "STATUS_SHA256:"
    cat "$AUDIT/FORMAL_BLINDED_COUNTING_STATUS_v1.sha256"
    echo
    echo "GSE300977_FORMAL_BLINDED_COUNTING_GATE = PASS"
} > "$AUDIT/FORMAL_BLINDED_COUNTING_FINAL_AUDIT_v1.txt"

shasum -a 256 \
"$AUDIT/FORMAL_BLINDED_COUNTING_FINAL_AUDIT_v1.txt" \
> "$AUDIT/FORMAL_BLINDED_COUNTING_STAGE_SHA256.txt"

echo
echo "============================================================"
echo "GSE300977_FORMAL_BLINDED_COUNTING_GATE = PASS"
echo "============================================================"

cat "$STATUS"

echo
cat "$AUDIT/FORMAL_BLINDED_COUNTING_STAGE_SHA256.txt"

echo
date
