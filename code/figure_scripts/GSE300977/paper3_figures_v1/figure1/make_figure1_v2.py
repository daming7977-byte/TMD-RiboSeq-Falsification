from pathlib import Path
import hashlib
import json
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle


# ============================================================
# AUTHORITATIVE SOURCES
# ============================================================

REPORT106 = Path(
    "REPRO_WORKSPACE/"
    "GSE106329_SRP_TMD_Spatial/"
    "PHASE_B_FIRST_BIOLOGICAL_REVEAL_REPORT.md"
)

SENS106 = Path(
    "REPRO_WORKSPACE/"
    "GSE106329_SRP_TMD_Spatial/"
    "phase_b_first_reveal/"
    "sensitivity_summary.tsv"
)

SUMMARY106 = Path(
    "REPRO_WORKSPACE/"
    "GSE106329_SRP_TMD_Spatial/"
    "phase_b_first_reveal/"
    "first_reveal_summary.json"
)

VERDICT133 = Path(
    "REPRO_WORKSPACE/GSE133111_stage0C/"
    "provenance/"
    "GSE133111_FINAL_FORMAL_VERDICT_2026-08-29.md"
)

OUTDIR = Path(
    "paper3_figures_v1/figure1"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

for p in [
    REPORT106,
    SENS106,
    SUMMARY106,
    VERDICT133,
]:
    if not p.exists():
        raise SystemExit(
            f"MISSING_AUTHORITATIVE_SOURCE:{p}"
        )

print(
    "AUTHORITATIVE_SOURCE_GATE = PASS"
)


# ============================================================
# READ SOURCES
# ============================================================

report106 = REPORT106.read_text(
    encoding="utf-8",
    errors="ignore"
)

sens106 = pd.read_csv(
    SENS106,
    sep="\t"
)

summary106 = json.loads(
    SUMMARY106.read_text(
        encoding="utf-8",
        errors="ignore"
    )
)

verdict133 = VERDICT133.read_text(
    encoding="utf-8",
    errors="ignore"
)


# ============================================================
# EXTRACT GSE106329 VALUES
# ============================================================

needed106 = [
    "primary",
    "umi_primary",
    "Sec63_noCHX_30_minus_0",
]

missing = [
    x
    for x in needed106
    if x not in set(
        sens106["analysis"]
        .astype(str)
    )
]

if missing:
    raise SystemExit(
        f"GSE106329_SENSITIVITY_ROWS_MISSING:{missing}"
    )


def sens_value(
    analysis,
    column
):
    return float(
        sens106.loc[
            sens106["analysis"].eq(
                analysis
            ),
            column
        ].iloc[0]
    )


primary_mean = sens_value(
    "primary",
    "domain_mean_effect"
)

primary_peak = sens_value(
    "primary",
    "absolute_peak_effect"
)

umi_mean = sens_value(
    "umi_primary",
    "domain_mean_effect"
)

umi_peak = sens_value(
    "umi_primary",
    "absolute_peak_effect"
)

umi_attenuation = sens_value(
    "umi_primary",
    "primary_mean_attenuation_fraction"
)

nochx_mean = sens_value(
    "Sec63_noCHX_30_minus_0",
    "domain_mean_effect"
)


# pseudo domain mean from authoritative report
sec8 = report106.split(
    "## 8. Matched pseudo-anchor diagnostic",
    1
)

if len(sec8) != 2:
    raise SystemExit(
        "GSE106329_PSEUDO_SECTION_GATE_FAIL"
    )

m = re.search(
    r"domain mean effect[^0-9+\-]*"
    r"([+\-]?[0-9]+\.[0-9]+)",
    sec8[1],
)

if not m:
    raise SystemExit(
        "GSE106329_PSEUDO_MEAN_PARSE_FAIL"
    )

pseudo_mean = float(
    m.group(1)
)

pseudo_peak = float(
    summary106[
        "pseudo_anchor"
    ][
        "maximum_absolute_effect"
    ]
)

null_p_max = float(
    summary106[
        "coordinate_null"
    ][
        "null_max_absolute_effect_p"
    ]
)

null_p_mean = float(
    summary106[
        "coordinate_null"
    ][
        "null_absolute_locked_band_mean_p"
    ]
)


# Frozen consistency gates
if not np.isclose(
    primary_mean,
    0.5562068559743386,
    atol=1e-6
):
    raise SystemExit(
        f"GSE106329_PRIMARY_MEAN_GATE_FAIL:{primary_mean}"
    )

if not np.isclose(
    umi_attenuation,
    0.8903519435825291,
    atol=1e-6
):
    raise SystemExit(
        f"GSE106329_UMI_ATTENUATION_GATE_FAIL:{umi_attenuation}"
    )

if not np.isclose(
    pseudo_peak,
    0.7615380774031452,
    atol=1e-6
):
    raise SystemExit(
        "GSE106329_PSEUDO_PEAK_GATE_FAIL"
    )

print(
    "GSE106329_FROZEN_VALUE_GATE = PASS"
)


# ============================================================
# EXTRACT GSE133111 VALUES
# ============================================================

section133 = verdict133.split(
    "## Primary all-frozen-length analysis",
    1
)

if len(section133) != 2:
    raise SystemExit(
        "GSE133111_PRIMARY_SECTION_GATE_FAIL"
    )

primary133 = section133[1]

mD = re.search(
    r"median D\s*=\s*"
    r"([+\-]?[0-9]+\.[0-9]+)",
    primary133
)

mP = re.search(
    r"one-sided Wilcoxon P\s*=\s*"
    r"([0-9]+\.[0-9]+)",
    primary133
)

mCI = re.search(
    r"90% Hodges-Lehmann interval:\s*"
    r"\[\s*([+\-]?[0-9]+\.[0-9]+)"
    r"\s*,\s*"
    r"([+\-]?[0-9]+\.[0-9]+)\s*\]",
    primary133
)

mNull = re.search(
    r"empirical P\s*=\s*"
    r"([0-9]+\.[0-9]+)",
    primary133
)

if not all([
    mD,
    mP,
    mCI,
    mNull,
]):
    raise SystemExit(
        "GSE133111_VALUE_PARSE_GATE_FAIL"
    )

d133 = float(
    mD.group(1)
)

p133 = float(
    mP.group(1)
)

ci133_lo = float(
    mCI.group(1)
)

ci133_hi = float(
    mCI.group(2)
)

null133_p = float(
    mNull.group(1)
)

if not np.isclose(
    d133,
    0.0197015,
    atol=1e-7
):
    raise SystemExit(
        "GSE133111_D_GATE_FAIL"
    )

print(
    "GSE133111_FROZEN_VALUE_GATE = PASS"
)


# ============================================================
# SAVE FIGURE-DERIVED INPUT TABLE
# ============================================================

derived = pd.DataFrame([
    {
        "dataset":
            "GSE106329",

        "metric":
            "primary_domain_mean",

        "value":
            primary_mean,

        "source":
            str(SENS106),
    },
    {
        "dataset":
            "GSE106329",

        "metric":
            "umi_domain_mean",

        "value":
            umi_mean,

        "source":
            str(SENS106),
    },
    {
        "dataset":
            "GSE106329",

        "metric":
            "umi_attenuation_fraction",

        "value":
            umi_attenuation,

        "source":
            str(SENS106),
    },
    {
        "dataset":
            "GSE106329",

        "metric":
            "pseudo_domain_mean",

        "value":
            pseudo_mean,

        "source":
            str(REPORT106),
    },
    {
        "dataset":
            "GSE106329",

        "metric":
            "pseudo_peak_absolute_effect",

        "value":
            pseudo_peak,

        "source":
            str(SUMMARY106),
    },
    {
        "dataset":
            "GSE106329",

        "metric":
            "coordinate_null_max_effect_p",

        "value":
            null_p_max,

        "source":
            str(SUMMARY106),
    },
    {
        "dataset":
            "GSE106329",

        "metric":
            "coordinate_null_domain_mean_p",

        "value":
            null_p_mean,

        "source":
            str(SUMMARY106),
    },
    {
        "dataset":
            "GSE133111",

        "metric":
            "real_vs_pseudo_median_D",

        "value":
            d133,

        "source":
            str(VERDICT133),
    },
    {
        "dataset":
            "GSE133111",

        "metric":
            "real_vs_pseudo_p",

        "value":
            p133,

        "source":
            str(VERDICT133),
    },
    {
        "dataset":
            "GSE133111",

        "metric":
            "coordinate_null_empirical_p",

        "value":
            null133_p,

        "source":
            str(VERDICT133),
    },
])

DERIVED = (
    OUTDIR /
    "FIGURE1_DERIVED_INPUTS_v1.tsv"
)

derived.to_csv(
    DERIVED,
    sep="\t",
    index=False
)

print(
    "FIGURE1_DERIVED_INPUT_GATE = PASS"
)


# ============================================================
# STYLE
# ============================================================

plt.rcParams.update({
    "font.family":
        "Arial",

    "font.size":
        7.5,

    "axes.labelsize":
        8,

    "axes.titlesize":
        8.5,

    "xtick.labelsize":
        7,

    "ytick.labelsize":
        7,

    "legend.fontsize":
        7,

    "axes.linewidth":
        0.8,

    "figure.facecolor":
        "white",

    "axes.facecolor":
        "white",

    "savefig.facecolor":
        "white",

    "savefig.transparent":
        False,

    "pdf.fonttype":
        42,

    "ps.fonttype":
        42,

    "svg.fonttype":
        "none",
})


COL_PRIMARY = "#4D4D4D"
COL_UMI = "#0072B2"
COL_PSEUDO = "#D55E00"
COL_NOCHX = "#009E73"

COL_PASS = "#4DAF4A"
COL_FAIL = "#E41A1C"
COL_NEUTRAL = "#BDBDBD"


# ============================================================
# HELPERS
# ============================================================

def panel_label(
    ax,
    letter
):
    ax.text(
        -0.16,
        1.07,
        letter,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
        ha="left"
    )


def status_box(
    ax,
    x,
    y,
    label,
    status,
    detail="",
):

    color = (
        COL_PASS
        if status == "PASS"
        else COL_FAIL
    )

    box = FancyBboxPatch(
        (x,y),
        0.82,
        0.115,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        transform=ax.transAxes,
        facecolor="white",
        edgecolor="0.78",
        linewidth=0.7,
    )

    ax.add_patch(
        box
    )

    ax.add_patch(
        Rectangle(
            (x+0.02,y+0.025),
            0.055,
            0.065,
            transform=ax.transAxes,
            facecolor=color,
            edgecolor="none",
        )
    )

    ax.text(
        x+0.095,
        y+0.072,
        label,
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=7.2,
        fontweight="bold",
    )

    ax.text(
        x+0.095,
        y+0.035,
        detail,
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=6.1,
        color="0.35",
    )


# ============================================================
# FIGURE LAYOUT
# ============================================================

fig = plt.figure(
    figsize=(7.15,6.5),
    facecolor="white"
)

fig.patch.set_facecolor(
    "white"
)

fig.patch.set_alpha(
    1.0
)

gs = fig.add_gridspec(
    2,
    3,
    left=0.07,
    right=0.985,
    bottom=0.10,
    top=0.955,
    wspace=0.42,
    hspace=0.48,
)

axA = fig.add_subplot(
    gs[0,0]
)

axB = fig.add_subplot(
    gs[0,1]
)

axC = fig.add_subplot(
    gs[0,2]
)

axD = fig.add_subplot(
    gs[1,0]
)

axE = fig.add_subplot(
    gs[1,1]
)

axF = fig.add_subplot(
    gs[1,2]
)

for ax in [
    axA,axB,axC,
    axD,axE,axF
]:
    ax.set_facecolor(
        "white"
    )


# ============================================================
# A — FALSIFICATION FRAMEWORK
# ============================================================

axA.axis(
    "off"
)

panel_label(
    axA,
    "A"
)

axA.text(
    0.5,
    0.95,
    "Falsification before inference",
    transform=axA.transAxes,
    ha="center",
    va="top",
    fontsize=9,
    fontweight="bold"
)

steps = [
    (
        "Apparent\nTMD-relative\npattern",
        COL_NEUTRAL
    ),
    (
        "Matched\npseudo\ncontrol",
        COL_PSEUDO
    ),
    (
        "Coordinate\nnull",
        "#9E9AC8"
    ),
    (
        "Technical\nrobustness",
        COL_UMI
    ),
    (
        "Biological\ninference",
        COL_PASS
    ),
]

ys = [
    0.77,
    0.61,
    0.45,
    0.29,
    0.10,
]

for i,((text,color),y) in enumerate(
    zip(
        steps,
        ys
    )
):

    box = FancyBboxPatch(
        (0.18,y),
        0.64,
        0.105,
        boxstyle=(
            "round,pad=0.012,"
            "rounding_size=0.02"
        ),
        transform=axA.transAxes,
        facecolor="white",
        edgecolor=color,
        linewidth=1.5,
    )

    axA.add_patch(
        box
    )

    axA.text(
        0.50,
        y+0.052,
        text,
        transform=axA.transAxes,
        ha="center",
        va="center",
        fontsize=7.3,
        fontweight=(
            "bold"
            if i == 4
            else "normal"
        )
    )

    if i < len(steps)-1:

        axA.annotate(
            "",
            xy=(0.50,ys[i+1]+0.115),
            xytext=(0.50,y-0.008),
            xycoords=axA.transAxes,
            textcoords=axA.transAxes,
            arrowprops=dict(
                arrowstyle="-|>",
                color="0.45",
                lw=0.9
            )
        )

axA.text(
    0.5,
    0.01,
    "A strong signal is not sufficient.",
    transform=axA.transAxes,
    ha="center",
    va="bottom",
    fontsize=6.7,
    color="0.3",
)


# ============================================================
# B — GSE106329 EFFECT SUMMARY
# ============================================================

vals = [
    primary_mean,
    umi_mean,
    pseudo_mean,
    nochx_mean,
]

labels = [
    "Primary\nCHX",
    "UMI\ndedup.",
    "Matched\npseudo",
    "No-CHX",
]

colors = [
    COL_PRIMARY,
    COL_UMI,
    COL_PSEUDO,
    COL_NOCHX,
]

x = np.arange(
    len(vals)
)

bars = axB.bar(
    x,
    vals,
    color=colors,
    width=0.68,
    edgecolor="none",
)

axB.axhline(
    0,
    color="0.65",
    lw=0.8
)

axB.set_xticks(
    x
)

axB.set_xticklabels(
    labels
)

axB.set_ylabel(
    r"Domain mean effect $\Delta$log$_2$"
)

axB.set_title(
    "GSE106329: apparent signal fails controls",
    pad=8
)

axB.spines["top"].set_visible(
    False
)

axB.spines["right"].set_visible(
    False
)

for bar,val in zip(
    bars,
    vals
):
    axB.text(
        bar.get_x()
        +
        bar.get_width()/2,
        val + (
            0.025
            if val >= 0
            else -0.045
        ),
        f"{val:+.3f}",
        ha="center",
        va=(
            "bottom"
            if val >= 0
            else "top"
        ),
        fontsize=6.5,
    )

axB.annotate(
    "89% attenuation",
    xy=(
        1,
        umi_mean
    ),
    xytext=(
        0.65,
        0.34
    ),
    textcoords="data",
    ha="center",
    fontsize=6.6,
    color=COL_UMI,
    arrowprops=dict(
        arrowstyle="->",
        lw=0.8,
        color=COL_UMI
    )
)

axB.text(
    0.97,
    0.84,
    "pseudo > real",
    transform=axB.transAxes,
    ha="right",
    va="top",
    fontsize=6.7,
    color=COL_PSEUDO,
    fontweight="bold"
)

panel_label(
    axB,
    "B"
)


# ============================================================
# C — GSE106329 GATE AUDIT
# ============================================================

axC.axis(
    "off"
)

panel_label(
    axC,
    "C"
)

axC.text(
    0.5,
    0.96,
    "GSE106329 frozen-criterion audit",
    transform=axC.transAxes,
    ha="center",
    va="top",
    fontsize=8.5,
    fontweight="bold",
)

status_box(
    axC,
    0.08,
    0.76,
    "Nominal broad band",
    "PASS",
    "131/131 centers"
)

status_box(
    axC,
    0.08,
    0.62,
    "UMI robustness",
    "FAIL",
    "89% attenuation"
)

status_box(
    axC,
    0.08,
    0.48,
    "Matched pseudo",
    "FAIL",
    "pseudo peak > real"
)

status_box(
    axC,
    0.08,
    0.34,
    "Coordinate null",
    "FAIL",
    f"P={null_p_max:.3f} / {null_p_mean:.3f}"
)

status_box(
    axC,
    0.08,
    0.20,
    "No-CHX support",
    "FAIL",
    "no spatial band"
)

status_box(
    axC,
    0.08,
    0.06,
    "Coverage / LOO stability",
    "PASS",
    "broad shift is real, not specific"
)


# ============================================================
# D — GSE133111 REAL-vs-PSEUDO EFFECT
# ============================================================

panel_label(
    axD,
    "D"
)

axD.axvline(
    0,
    color="0.65",
    lw=0.9,
    ls="--"
)

axD.errorbar(
    d133,
    0,
    xerr=np.array([
        [d133-ci133_lo],
        [ci133_hi-d133]
    ]),
    fmt="o",
    markersize=6,
    color=COL_PRIMARY,
    ecolor=COL_PRIMARY,
    capsize=4,
    lw=1.5,
)

axD.set_xlim(
    -0.12,
    0.15
)

axD.set_ylim(
    -0.48,
    0.48
)

axD.set_yticks(
    []
)

axD.set_xlabel(
    r"Median real $-$ pseudo $D$"
)

axD.set_title(
    "GSE133111: matched-pseudo test",
    pad=8
)

axD.spines["top"].set_visible(
    False
)

axD.spines["right"].set_visible(
    False
)

axD.spines["left"].set_visible(
    False
)

axD.text(
    d133,
    0.25,
    f"D = {d133:+.3f}",
    ha="center",
    va="bottom",
    fontsize=7,
    fontweight="bold"
)

axD.text(
    0.5,
    0.15,
    "90% Hodges–Lehmann interval",
    transform=axD.transAxes,
    ha="center",
    fontsize=6.6,
    color="0.35"
)

axD.text(
    0.5,
    0.91,
    f"one-sided Wilcoxon P = {p133:.3f}",
    transform=axD.transAxes,
    ha="center",
    fontsize=7
)

axD.text(
    0.5,
    0.80,
    "REAL-vs-PSEUDO: FAIL",
    transform=axD.transAxes,
    ha="center",
    fontsize=7.5,
    color=COL_FAIL,
    fontweight="bold"
)


# ============================================================
# E — GSE133111 GATE AUDIT
# ============================================================

axE.axis(
    "off"
)

panel_label(
    axE,
    "E"
)

axE.text(
    0.5,
    0.96,
    "GSE133111 formal Stage 1 verdict",
    transform=axE.transAxes,
    ha="center",
    va="top",
    fontsize=8.5,
    fontweight="bold"
)

status_box(
    axE,
    0.08,
    0.76,
    "Real vs pseudo",
    "FAIL",
    f"D={d133:+.020f}"[:9]
)

status_box(
    axE,
    0.08,
    0.62,
    "Coordinate-null V2",
    "PASS",
    f"empirical P={null133_p:.5f}"
)

status_box(
    axE,
    0.08,
    0.48,
    "Localized downstream band",
    "FAIL",
    "Wilcoxon P=0.738"
)

status_box(
    axE,
    0.08,
    0.34,
    "No-global-shift",
    "PASS",
    "Q interval includes 0"
)

status_box(
    axE,
    0.08,
    0.20,
    "Transcript LOO",
    "FAIL",
    "sign instability"
)

status_box(
    axE,
    0.08,
    0.06,
    "Gene LOO",
    "FAIL",
    "formal gate failed"
)


# ============================================================
# F — CROSS-BENCHMARK LESSON
# ============================================================

axF.axis(
    "off"
)

panel_label(
    axF,
    "F"
)

axF.text(
    0.5,
    0.95,
    "Cross-benchmark inference rule",
    transform=axF.transAxes,
    ha="center",
    va="top",
    fontsize=9,
    fontweight="bold"
)

# left source
box1 = FancyBboxPatch(
    (0.04,0.64),
    0.40,
    0.19,
    boxstyle="round,pad=0.015",
    transform=axF.transAxes,
    facecolor="#F7F7F7",
    edgecolor="0.7",
    linewidth=0.8,
)

axF.add_patch(
    box1
)

axF.text(
    0.24,
    0.735,
    "GSE106329",
    transform=axF.transAxes,
    ha="center",
    va="center",
    fontsize=8,
    fontweight="bold"
)

axF.text(
    0.24,
    0.675,
    "strong broad signal\n≠ TMD-specific",
    transform=axF.transAxes,
    ha="center",
    va="center",
    fontsize=6.8,
)

# right source
box2 = FancyBboxPatch(
    (0.56,0.64),
    0.40,
    0.19,
    boxstyle="round,pad=0.015",
    transform=axF.transAxes,
    facecolor="#F7F7F7",
    edgecolor="0.7",
    linewidth=0.8,
)

axF.add_patch(
    box2
)

axF.text(
    0.76,
    0.735,
    "GSE133111",
    transform=axF.transAxes,
    ha="center",
    va="center",
    fontsize=8,
    fontweight="bold"
)

axF.text(
    0.76,
    0.675,
    "coordinate-null PASS\n≠ TMD validation",
    transform=axF.transAxes,
    ha="center",
    va="center",
    fontsize=6.8,
)

# arrows down
for xpos in [
    0.24,
    0.76
]:
    axF.annotate(
        "",
        xy=(0.50,0.50),
        xytext=(xpos,0.61),
        xycoords=axF.transAxes,
        textcoords=axF.transAxes,
        arrowprops=dict(
            arrowstyle="-|>",
            lw=1.0,
            color="0.45"
        )
    )

# conclusion
final_box = FancyBboxPatch(
    (0.08,0.20),
    0.84,
    0.29,
    boxstyle="round,pad=0.018",
    transform=axF.transAxes,
    facecolor="#EEF5F9",
    edgecolor=COL_UMI,
    linewidth=1.1,
)

axF.add_patch(
    final_box
)

axF.text(
    0.50,
    0.405,
    "Credible TMD-specific inference requires",
    transform=axF.transAxes,
    ha="center",
    va="center",
    fontsize=7.5,
    fontweight="bold"
)

axF.text(
    0.50,
    0.315,
    "matched controls + technical robustness\n"
    "+ positional specificity\n"
    "+ preserved inference hierarchy",
    transform=axF.transAxes,
    ha="center",
    va="center",
    fontsize=6.4,
)

axF.text(
    0.5,
    0.07,
    "These rules motivate the independent\nGSE300977 spatial-state analysis.",
    transform=axF.transAxes,
    ha="center",
    va="center",
    fontsize=6.8,
    color="0.3",
)


# ============================================================
# GLOBAL FOOTNOTE
# ============================================================

fig.text(
    0.5,
    0.022,
    "Values are reproduced from frozen authoritative benchmark outputs; "
    "no benchmark re-analysis was performed.",
    ha="center",
    va="bottom",
    fontsize=6.1,
)


# ============================================================
# SAVE
# ============================================================

pdf = (
    OUTDIR /
    "FIGURE1_FALSIFICATION_FRAMEWORK_v2.pdf"
)

svg = (
    OUTDIR /
    "FIGURE1_FALSIFICATION_FRAMEWORK_v2.svg"
)

png = (
    OUTDIR /
    "FIGURE1_FALSIFICATION_FRAMEWORK_v2.png"
)

fig.savefig(
    pdf,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="white",
    transparent=False,
)

fig.savefig(
    svg,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="white",
    transparent=False,
)

fig.savefig(
    png,
    dpi=450,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="white",
    transparent=False,
)

plt.close(
    fig
)


# ============================================================
# OUTPUT SHA
# ============================================================

for p in [
    DERIVED,
    pdf,
    svg,
    png,
]:

    sha = hashlib.sha256(
        p.read_bytes()
    ).hexdigest()

    print(
        "SHA256",
        sha,
        p
    )

print()
print(
    "FIGURE1_GENERATION_GATE = PASS"
)
