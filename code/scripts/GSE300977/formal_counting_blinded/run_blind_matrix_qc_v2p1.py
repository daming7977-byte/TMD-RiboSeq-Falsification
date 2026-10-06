from pathlib import Path
import numpy as np
import pandas as pd
from plastid import SegmentChain

ROOT = Path("formal_counting_blinded/counts_v1")
AUDIT = Path("formal_counting_blinded/audit")

ROI = Path(
    "periodicity_blinded/roi_masked_v2/"
    "GRCh38_104_cds_start_crossmap_rois.txt"
)

OUT = AUDIT / "BLIND_MATRIX_QC_v2.tsv"

libs = [f"B{i:02d}" for i in range(1,13)]

EXPECTED_SHAPE = (5130,100)

# ---------------------------------------------------------
# Locate outputs
# ---------------------------------------------------------

def find_matrix(lib, kind):
    hits = list(ROOT.glob(
        f"{lib}.cds_start_center_n13_{kind}.txt"
    ))
    hits += list(ROOT.glob(
        f"{lib}.cds_start_center_n13_{kind}.txt.gz"
    ))

    if len(hits) != 1:
        raise SystemExit(
            f"{lib} {kind}: expected 1 file, found {len(hits)}"
        )

    return hits[0]

def find_profile(lib):
    hits = list(ROOT.glob(
        f"{lib}.cds_start_center_n13_metagene_profile.txt"
    ))
    hits += list(ROOT.glob(
        f"{lib}.cds_start_center_n13_metagene_profile.txt.gz"
    ))

    if len(hits) != 1:
        raise SystemExit(
            f"{lib} profile: expected 1 file, found {len(hits)}"
        )

    return hits[0]

def load_matrix(path):
    x = np.loadtxt(path)

    if x.ndim == 1:
        x = x.reshape(1,-1)

    return x

# ---------------------------------------------------------
# Independently reconstruct STATIC ROI/Crossmap mask
#
# plastid's count routine starts fully masked, unmasks the
# actual ROI window, and leaves ROI mask regions masked.
# ---------------------------------------------------------

roi_table = pd.read_table(
    ROI,
    sep="\t",
    comment="#",
    header=0
)

if len(roi_table) != 5130:
    raise SystemExit(
        f"ROI_COUNT_FAIL: {len(roi_table)}"
    )

width = int(roi_table["window_size"].iloc[0])

if width != 100:
    raise SystemExit(
        f"ROI_WIDTH_FAIL: {width}"
    )

static_mask = np.ones(EXPECTED_SHAPE, dtype=bool)

real_mask_rois = 0
mask_positions_not_in_roi = 0

for i, row in roi_table.iterrows():

    roi = SegmentChain.from_str(str(row["region"]))

    offset = int(round(float(row["alignment_offset"])))

    if offset < 0 or offset + roi.length > 100:
        raise SystemExit(
            f"ROI_OFFSET_FAIL row={i}"
        )

    # actual ROI positions become available
    static_mask[
        i,
        offset:offset + roi.length
    ] = False

    v = row["masked"]

    if pd.isna(v):
        continue

    v = str(v).strip()

    if v.lower() in ("", "na", "nan"):
        continue

    real_mask_rois += 1

    m = SegmentChain.from_str(v)

    roi_positions = list(roi.get_position_list())

    # get_position_list is genomic ascending order.
    # metagene vector follows transcript orientation.
    if roi.spanning_segment.strand == "-":
        roi_positions = roi_positions[::-1]

    position_to_index = {
        pos:j for j,pos in enumerate(roi_positions)
    }

    mask_positions = set(m.get_position_list())

    for pos in mask_positions:

        if pos not in position_to_index:
            mask_positions_not_in_roi += 1
            continue

        j = position_to_index[pos]

        static_mask[i, offset+j] = True

if real_mask_rois != 1676:
    raise SystemExit(
        f"STATIC_REAL_MASK_ROI_FAIL: {real_mask_rois}"
    )

if mask_positions_not_in_roi != 0:
    raise SystemExit(
        "STATIC_MASK_COORDINATE_FAIL: "
        f"{mask_positions_not_in_roi}"
    )

static_true_cells = int(static_mask.sum())

print("STATIC_MASK_RECONSTRUCTION = PASS")
print("STATIC_REAL_MASK_ROIS =", real_mask_rois)
print("STATIC_MASK_TRUE_CELLS =", static_true_cells)

# ---------------------------------------------------------
# Library-level blind checks
# ---------------------------------------------------------

rows = []

for lib in libs:

    raw_path  = find_matrix(lib, "rawcounts")
    norm_path = find_matrix(lib, "normcounts")
    mask_path = find_matrix(lib, "mask")
    prof_path = find_profile(lib)

    raw  = load_matrix(raw_path)
    norm = load_matrix(norm_path)
    mask = load_matrix(mask_path)

    if raw.shape != EXPECTED_SHAPE:
        raise SystemExit(
            f"{lib} RAW_SHAPE_FAIL {raw.shape}"
        )

    if norm.shape != EXPECTED_SHAPE:
        raise SystemExit(
            f"{lib} NORM_SHAPE_FAIL {norm.shape}"
        )

    if mask.shape != EXPECTED_SHAPE:
        raise SystemExit(
            f"{lib} MASK_SHAPE_FAIL {mask.shape}"
        )

    # mask file must contain only 0/1
    finite_mask = mask[np.isfinite(mask)]
    unique_mask_values = np.unique(finite_mask)

    if not np.all(
        np.isin(unique_mask_values, [0.0,1.0])
    ):
        raise SystemExit(
            f"{lib} MASK_BINARY_GATE_FAIL "
            f"{unique_mask_values}"
        )

    saved_mask = mask.astype(bool)

    # THE IMPORTANT GATE:
    # Every frozen structural-mask cell MUST be masked
    # in every sample.
    structural_violations = int(
        np.sum(static_mask & ~saved_mask)
    )

    if structural_violations != 0:
        raise SystemExit(
            f"{lib} STATIC_MASK_SUBSET_FAIL "
            f"{structural_violations}"
        )

    # additional masks are permitted and sample-dependent
    extra_dynamic_mask = int(
        np.sum(saved_mask & ~static_mask)
    )

    # raw counts must be nonnegative wherever finite
    finite_raw = raw[np.isfinite(raw)]

    negative_raw = int(
        np.sum(finite_raw < -1e-12)
    )

    if negative_raw != 0:
        raise SystemExit(
            f"{lib} NEGATIVE_RAW_FAIL "
            f"{negative_raw}"
        )

    finite_norm = norm[np.isfinite(norm)]

    negative_norm = int(
        np.sum(finite_norm < -1e-12)
    )

    if negative_norm != 0:
        raise SystemExit(
            f"{lib} NEGATIVE_NORM_FAIL "
            f"{negative_norm}"
        )

    row_total = np.nansum(raw, axis=1)

    active_rois = int(
        np.sum(row_total > 0)
    )

    raw_total = float(
        np.nansum(raw)
    )

    # profile structural gate
    profile = pd.read_table(
        prof_path,
        sep="\t",
        comment="#"
    )

    if len(profile) != 100:
        raise SystemExit(
            f"{lib} PROFILE_LENGTH_FAIL "
            f"{len(profile)}"
        )

    expected_x = np.arange(-50,50)

    if not np.array_equal(
        profile["x"].to_numpy(),
        expected_x
    ):
        raise SystemExit(
            f"{lib} PROFILE_X_FAIL"
        )

    region_min = int(
        profile["regions_counted"].min()
    )

    region_max = int(
        profile["regions_counted"].max()
    )

    rows.append([
        lib,
        5130,
        100,
        static_true_cells,
        int(saved_mask.sum()),
        extra_dynamic_mask,
        structural_violations,
        negative_raw,
        negative_norm,
        active_rois,
        raw_total,
        region_min,
        region_max,
    ])

# ---------------------------------------------------------
# Save descriptive blind QC
# No post-hoc library exclusion thresholds.
# ---------------------------------------------------------

header = [
    "library",
    "n_roi",
    "width",
    "static_mask_cells",
    "saved_mask_cells",
    "extra_dynamic_mask_cells",
    "static_mask_violations",
    "negative_raw_cells",
    "negative_norm_cells",
    "active_rois",
    "raw_total",
    "regions_counted_min",
    "regions_counted_max",
]

with OUT.open("w") as fh:

    fh.write(
        "\t".join(header) + "\n"
    )

    for row in rows:
        fh.write(
            "\t".join(map(str,row)) + "\n"
        )

print()
print("LIBRARIES = 12")
print("EXPECTED_SHAPE = 5130 x 100")
print("ALL_MATRIX_SHAPES = PASS")
print("STATIC_MASK_SUBSET_GATE = PASS")
print("MASK_BINARY_GATE = PASS")
print("NEGATIVE_RAW_COUNT_GATE = PASS")
print("NEGATIVE_NORM_COUNT_GATE = PASS")
print("PROFILE_STRUCTURE_GATE = PASS")

print()
print(
    "NOTE = extra_dynamic_mask_cells may differ "
    "between libraries by design"
)

print()
print(OUT.read_text())
