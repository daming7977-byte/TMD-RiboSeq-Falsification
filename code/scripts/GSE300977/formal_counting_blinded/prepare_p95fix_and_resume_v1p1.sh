#!/bin/bash

set -uo pipefail
export MPLBACKEND=Agg
export LC_ALL=C

ROOT="$PWD"

PYTHON="REPRO_ENV/envs/gse300977-plastid061-x86/bin/python"

ROI="periodicity_blinded/roi_masked_v2/GRCh38_104_cds_start_crossmap_rois.txt"

OUTDIR="formal_counting_blinded/counts_v1"
LOGDIR="formal_counting_blinded/logs"
AUDIT="formal_counting_blinded/audit"
PATCHDIR="formal_counting_blinded/compat_patch"

OLDSTATUS="$AUDIT/FORMAL_BLINDED_COUNTING_STATUS_v1.tsv"
NEWSTATUS="$AUDIT/FORMAL_BLINDED_COUNTING_STATUS_v1p1.tsv"

EXPECTED_ROI_SHA="0f89aee14812df3e4a79c1173c90eeaf82d81b1f5268962f08de9cb98918b0cc"

mkdir -p "$OUTDIR" "$LOGDIR" "$AUDIT" "$PATCHDIR"

echo "============================================================"
echo "GSE300977 P95 COMPATIBILITY FIX + COUNTING RESUME"
date
echo "============================================================"

############################################################
# 1. VERIFY OLD RUN STOPPED EXACTLY AS EXPECTED
############################################################

if [ ! -s "$OLDSTATUS" ]; then
    echo "OLD_STATUS_GATE = FAIL"
    exit 10
fi

cp "$OLDSTATUS" \
"$AUDIT/FORMAL_BLINDED_COUNTING_STATUS_v1_PRE_P95FIX.tsv"

echo "=== PRE-FIX STATUS ==="
cat "$OLDSTATUS"

PASS123=$(
    awk '
    ($1=="B01" || $1=="B02" || $1=="B03") && $2=="PASS" {n++}
    END {print n+0}
    ' "$OLDSTATUS"
)

if [ "$PASS123" -ne 3 ]; then
    echo "B01_B03_PREVIOUS_PASS_GATE = FAIL"
    exit 11
fi

echo "B01_B03_PREVIOUS_PASS_GATE = PASS"

if ! awk '$1=="B04" && $2=="FAIL"{found=1} END{exit(found?0:1)}' "$OLDSTATUS"; then
    echo "B04_EXPECTED_FAILURE_GATE = FAIL"
    exit 12
fi

echo "B04_EXPECTED_FAILURE_GATE = PASS"

############################################################
# 2. VERIFY B01-B03 SEALED OUTPUTS EXIST
############################################################

for i in 1 2 3; do
    B=$(printf 'B%02d' "$i")

    if [ ! -s "$AUDIT/${B}.FORMAL_COUNTING_OUTPUT_SHA256.txt" ]; then
        echo "$B PRIOR_OUTPUT_SEAL_GATE = FAIL"
        exit 13
    fi

    echo "$B PRIOR_OUTPUT_SEAL_GATE = PASS"
done

############################################################
# 3. PRESERVE B04 FAILED-RUN ARTIFACTS
############################################################

FAILDIR="$AUDIT/B04_FAILED_P95_PLOT_RUN"

mkdir -p "$FAILDIR"

for F in "$OUTDIR"/B04.cds_start_center_n13*; do
    if [ -e "$F" ]; then
        mv "$F" "$FAILDIR/"
    fi
done

if [ -s "$LOGDIR/B04.formal_counting.log" ]; then
    cp "$LOGDIR/B04.formal_counting.log" \
       "$FAILDIR/B04.formal_counting.PRE_P95FIX.log"
fi

{
    echo "GSE300977 B04 P95 PLOTTING FAILURE AUDIT"
    echo "DATE = $(date)"
    echo
    echo "B01-B03 = PREVIOUSLY PASS"
    echo "B04 = METAGENE_EXIT_NONZERO"
    echo
    echo "FAILURE_LOCATION = plastid.bin.metagene.do_count overview heatmap p95 calculation"
    echo "ERROR = ValueError: output array is read-only"
    echo
    echo "SCIENTIFIC_PARAMETERS_CHANGED = NO"
    echo "BIOLOGICAL_IDENTITY_VIEWED = NO"
    echo "B05-B12_EXECUTED = NO"
    echo
    echo "ACTION = preserve partial B04 artifacts and apply plotting-only compatibility fix"
} > "$FAILDIR/B04_P95_FAILURE_AUDIT.txt"

############################################################
# 4. LOCATE ORIGINAL PLASTID METAGENE SOURCE
############################################################

ORIG=$(
"$PYTHON" - <<'PY'
import plastid.bin.metagene as m
print(m.__file__)
PY
)

echo "ORIGINAL_METAGENE = $ORIG"

if [ ! -s "$ORIG" ]; then
    echo "ORIGINAL_METAGENE_SOURCE_GATE = FAIL"
    exit 20
fi

cp "$ORIG" "$PATCHDIR/metagene.original.py"

shasum -a 256 "$PATCHDIR/metagene.original.py" \
> "$PATCHDIR/metagene.original.py.sha256"

############################################################
# 5. MAKE PROJECT-LOCAL PLOTTING-ONLY PATCH
############################################################

"$PYTHON" - "$PATCHDIR/metagene.original.py" \
              "$PATCHDIR/metagene.p95fix.py" <<'PY'
from pathlib import Path
import sys

src = Path(sys.argv[1])
dst = Path(sys.argv[2])

text = src.read_text()

old = "p95 = numpy.nanpercentile(rs[rs > 0], 95)"

new = (
    "p95 = numpy.nanpercentile("
    "numpy.ma.filled(rs[rs > 0], numpy.nan).copy(), 95)"
)

n = text.count(old)

print("TARGET_OCCURRENCES =", n)

if n != 1:
    raise SystemExit("PATCH_TARGET_GATE_FAIL")

dst.write_text(text.replace(old, new))

print("PATCH_TARGET_GATE = PASS")
print("PATCH_TYPE = PLOTTING_ONLY_P95_READONLY_COMPATIBILITY")
PY

if [ $? -ne 0 ]; then
    echo "P95_PATCH_GATE = FAIL"
    exit 21
fi

shasum -a 256 "$PATCHDIR/metagene.p95fix.py" \
> "$PATCHDIR/metagene.p95fix.py.sha256"

echo "P95_PATCH_GATE = PASS"

############################################################
# 6. RECORD EXACT PATCH DIFF
############################################################

diff -u \
"$PATCHDIR/metagene.original.py" \
"$PATCHDIR/metagene.p95fix.py" \
> "$PATCHDIR/metagene.p95fix.diff" || true

echo "=== PATCH DIFF ==="
cat "$PATCHDIR/metagene.p95fix.diff"

############################################################
# 7. CREATE PROJECT-LOCAL LAUNCHER
############################################################

cat > "$PATCHDIR/run_metagene_p95fix.py" <<'PY'
import importlib.util
import sys
from pathlib import Path

patched = Path(__file__).with_name("metagene.p95fix.py")

spec = importlib.util.spec_from_file_location(
    "plastid.bin.metagene_p95fix",
    str(patched)
)

mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

ret = mod.main(sys.argv[1:])

if ret is None:
    ret = 0

raise SystemExit(ret)
PY

############################################################
# 8. FREEZE COMPATIBILITY PATCH AUDIT
############################################################

{
    echo "GSE300977 PLASTID P95 COMPATIBILITY PATCH V1P1"
    echo "DATE = $(date)"
    echo
    echo "PLASTID_VERSION = 0.6.1"
    echo
    echo "FAILURE_TRIGGER:"
    echo "numpy.nanpercentile(rs[rs > 0], 95)"
    echo
    echo "PATCH:"
    echo "numpy.nanpercentile(numpy.ma.filled(rs[rs > 0], numpy.nan).copy(), 95)"
    echo
    echo "SCOPE = OVERVIEW HEATMAP P95 COLOR SCALING ONLY"
    echo
    echo "COUNTING_PARAMETERS_CHANGED = NO"
    echo "READ_LENGTH_MIN = 27"
    echo "READ_LENGTH_MAX = 36"
    echo "CENTER = TRUE"
    echo "NIBBLE = 13"
    echo "ROI_CHANGED = NO"
    echo "BAM_CHANGED = NO"
    echo "MASK_CHANGED = NO"
    echo "BIOLOGICAL_IDENTITY = LOCKED"
    echo
    echo "ORIGINAL_SOURCE_SHA256:"
    cat "$PATCHDIR/metagene.original.py.sha256"
    echo
    echo "PATCHED_SOURCE_SHA256:"
    cat "$PATCHDIR/metagene.p95fix.py.sha256"
} > "$AUDIT/PLASTID_P95_COMPATIBILITY_PATCH_v1p1.txt"

shasum -a 256 \
"$AUDIT/PLASTID_P95_COMPATIBILITY_PATCH_v1p1.txt" \
> "$AUDIT/PLASTID_P95_COMPATIBILITY_PATCH_v1p1.sha256"

############################################################
# 9. VERIFY FROZEN ROI
############################################################

ACTUAL_ROI_SHA=$(
    shasum -a 256 "$ROI" | awk '{print $1}'
)

if [ "$ACTUAL_ROI_SHA" != "$EXPECTED_ROI_SHA" ]; then
    echo "ROI_SHA_GATE = FAIL"
    exit 30
fi

echo "ROI_SHA_GATE = PASS"

############################################################
# 10. INITIALIZE CLEAN RESUME STATUS WITH B01-B03 ONLY
############################################################

awk '
($1=="B01" || $1=="B02" || $1=="B03") && $2=="PASS"
' "$OLDSTATUS" > "$NEWSTATUS"

if [ "$(wc -l < "$NEWSTATUS" | tr -d ' ')" -ne 3 ]; then
    echo "RESUME_STATUS_INITIALIZATION_GATE = FAIL"
    exit 31
fi

echo "RESUME_STATUS_INITIALIZATION_GATE = PASS"

############################################################
# 11. FUNCTION: RUN B04-B12 USING PROJECT-LOCAL PATCH
############################################################

run_one () {

    B="$1"

    BAM="unique_mapping_blinded/bam/${B}.NH1.unique.bam"
    BASE="$OUTDIR/${B}.cds_start_center_n13"
    LOG="$LOGDIR/${B}.formal_counting.v1p1.log"

    echo
    echo "============================================================"
    echo "START $B"
    date
    echo "============================================================"

    "$PYTHON" "$PATCHDIR/run_metagene_p95fix.py" \
        count \
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
        | tee -a "$NEWSTATUS"

        echo "$B COUNTING_GATE = FAIL"

        tail -60 "$LOG"

        return 1
    fi

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
            | tee -a "$NEWSTATUS"

            echo "$B OUTPUT_GATE = FAIL"

            return 1
        fi
    done

    ########################################################
    # STRUCTURE ONLY — NO BIOLOGICAL INTERPRETATION
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

expected = (5130, 100)

print("%s_RAW_SHAPE =" % lib, raw.shape)
print("%s_NORM_SHAPE =" % lib, norm.shape)
print("%s_MASK_SHAPE =" % lib, mask.shape)

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
        | tee -a "$NEWSTATUS"

        return 1
    fi

    {
        shasum -a 256 "$RAW"
        shasum -a 256 "$NORM"
        shasum -a 256 "$MASK"
        shasum -a 256 "$PROFILE"
    } > "$AUDIT/${B}.FORMAL_COUNTING_OUTPUT_SHA256.txt"

    printf "%s\tPASS\tRC=0\n" "$B" \
    | tee -a "$NEWSTATUS"

    echo "$B FORMAL_COUNTING_GATE = PASS"
    date

    return 0
}

############################################################
# 12. RESUME EXACTLY AT B04
############################################################

for i in 4 5 6 7 8 9 10 11 12; do

    B=$(printf 'B%02d' "$i")

    run_one "$B"

    if [ $? -ne 0 ]; then

        echo "$B FAILED"
        echo "QUEUE STOPPED BEFORE NEXT LIBRARY"

        exit 40
    fi
done

############################################################
# 13. FINAL 12/12 GATE
############################################################

PASS_COUNT=$(
    awk '$2=="PASS"{n++} END{print n+0}' "$NEWSTATUS"
)

echo
echo "PASS_COUNT = $PASS_COUNT"

if [ "$PASS_COUNT" -ne 12 ]; then
    echo "GSE300977_FORMAL_BLINDED_COUNTING_GATE = FAIL"
    exit 50
fi

shasum -a 256 "$NEWSTATUS" \
> "$AUDIT/FORMAL_BLINDED_COUNTING_STATUS_v1p1.sha256"

{
    echo "GSE300977 FORMAL BLINDED COUNTING FINAL AUDIT V1P1"
    echo "DATE = $(date)"
    echo
    echo "B01-B03 = PASS UNDER ORIGINAL PLASTID 0.6.1"
    echo "B04-B12 = PASS UNDER PROJECT-LOCAL P95 PLOTTING COMPATIBILITY FIX"
    echo
    echo "P95_PATCH_SCOPE = OVERVIEW HEATMAP COLOR SCALING ONLY"
    echo "SCIENTIFIC_COUNTING_PARAMETERS_CHANGED = NO"
    echo
    echo "LIBRARIES_PASS = 12"
    echo "LIBRARIES_EXPECTED = 12"
    echo
    echo "BIOLOGICAL_IDENTITY = LOCKED"
    echo
    echo "PATCH_AUDIT_SHA256:"
    cat "$AUDIT/PLASTID_P95_COMPATIBILITY_PATCH_v1p1.sha256"
    echo
    echo "STATUS_SHA256:"
    cat "$AUDIT/FORMAL_BLINDED_COUNTING_STATUS_v1p1.sha256"
    echo
    echo "GSE300977_FORMAL_BLINDED_COUNTING_GATE = PASS"
} > "$AUDIT/FORMAL_BLINDED_COUNTING_FINAL_AUDIT_v1p1.txt"

shasum -a 256 \
"$AUDIT/FORMAL_BLINDED_COUNTING_FINAL_AUDIT_v1p1.txt" \
> "$AUDIT/FORMAL_BLINDED_COUNTING_STAGE_SHA256_v1p1.txt"

echo
echo "============================================================"
echo "GSE300977_FORMAL_BLINDED_COUNTING_GATE = PASS"
echo "============================================================"

cat "$NEWSTATUS"

echo
echo "FINAL_STAGE_SEAL:"
cat "$AUDIT/FORMAL_BLINDED_COUNTING_STAGE_SHA256_v1p1.txt"

echo
date

