from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("formal_counting_blinded/counts_v1")
OUTDIR = Path("formal_counting_blinded/audit")
libs = [f"B{i:02d}" for i in range(1,13)]

profiles = {}

for lib in libs:
    f = ROOT / f"{lib}.cds_start_center_n13_metagene_profile.txt"

    if not f.exists():
        raise SystemExit(f"MISSING_PROFILE: {lib}")

    df = pd.read_table(
        f,
        sep="\t",
        comment="#"
    )

    if list(df.columns) != [
        "x",
        "metagene_average",
        "regions_counted"
    ]:
        raise SystemExit(
            f"{lib} PROFILE_COLUMNS_FAIL: {list(df.columns)}"
        )

    if len(df) != 100:
        raise SystemExit(
            f"{lib} PROFILE_LENGTH_FAIL: {len(df)}"
        )

    if not np.array_equal(
        df["x"].to_numpy(),
        np.arange(-50,50)
    ):
        raise SystemExit(f"{lib} PROFILE_X_FAIL")

    y = df["metagene_average"].to_numpy(dtype=float)

    if not np.all(np.isfinite(y)):
        raise SystemExit(f"{lib} PROFILE_NONFINITE_FAIL")

    if np.any(y < -1e-12):
        raise SystemExit(f"{lib} PROFILE_NEGATIVE_FAIL")

    profiles[lib] = y

M = np.vstack([profiles[x] for x in libs])

pearson = np.corrcoef(M)

spearman = (
    pd.DataFrame(M.T, columns=libs)
      .corr(method="spearman")
      .to_numpy()
)

pd.DataFrame(
    pearson,
    index=libs,
    columns=libs
).to_csv(
    OUTDIR / "BLIND_PROFILE_PEARSON_v1.tsv",
    sep="\t"
)

pd.DataFrame(
    spearman,
    index=libs,
    columns=libs
).to_csv(
    OUTDIR / "BLIND_PROFILE_SPEARMAN_v1.tsv",
    sep="\t"
)

summary = []

for i,lib in enumerate(libs):

    other_p = np.delete(pearson[i,:], i)
    other_s = np.delete(spearman[i,:], i)

    summary.append([
        lib,
        float(np.sum(M[i])),
        float(np.linalg.norm(M[i])),
        float(np.min(other_p)),
        float(np.median(other_p)),
        float(np.max(other_p)),
        float(np.min(other_s)),
        float(np.median(other_s)),
        float(np.max(other_s)),
    ])

cols = [
    "library",
    "profile_sum",
    "profile_l2_norm",
    "pearson_min_other",
    "pearson_median_other",
    "pearson_max_other",
    "spearman_min_other",
    "spearman_median_other",
    "spearman_max_other"
]

summary_df = pd.DataFrame(summary, columns=cols)

summary_df.to_csv(
    OUTDIR / "BLIND_PROFILE_CONCORDANCE_SUMMARY_v1.tsv",
    sep="\t",
    index=False
)

print("LIBRARIES = 12")
print("PROFILE_LENGTH = 100")
print("PROFILE_STRUCTURE_GATE = PASS")
print("PROFILE_FINITE_GATE = PASS")
print("PROFILE_NONNEGATIVE_GATE = PASS")
print()
print(summary_df.to_string(index=False))
