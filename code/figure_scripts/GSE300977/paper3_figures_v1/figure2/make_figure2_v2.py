from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

# ============================================================
# SOURCES
# ============================================================

REP = Path(
    "tmd_spatial_primary_reveal_v1/output/"
    "PRIMARY_PLUS75_REPLICATE_CONTRASTS_v1.tsv"
)

LOC = Path(
    "tmd_spatial_primary_reveal_v1/output/"
    "PRIMARY_PLUS75_LOCATION_SUMMARY_v1.tsv"
)

OUTDIR = Path(
    "paper3_figures_v1/figure2"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

for p in [REP, LOC]:
    if not p.exists():
        raise SystemExit(
            f"MISSING_PRIMARY_SOURCE:{p}"
        )

print("PRIMARY_SOURCE_GATE = PASS")

# ============================================================
# LOAD
# ============================================================

rep = pd.read_csv(
    REP,
    sep="\t"
)

loc = pd.read_csv(
    LOC,
    sep="\t"
)

print(
    "REPLICATE_COLUMNS =",
    list(rep.columns)
)

print(
    "LOCATION_COLUMNS =",
    list(loc.columns)
)

# ============================================================
# FLEXIBLE COLUMN RESOLUTION
# ============================================================

def pick(df, candidates, label):
    for c in candidates:
        if c in df.columns:
            return c
    raise SystemExit(
        f"COLUMN_GATE_FAIL:{label}:"
        f"{candidates}:"
        f"OBSERVED={list(df.columns)}"
    )

REP_LOC = pick(
    rep,
    [
        "location",
        "compartment",
        "spatial_state",
    ],
    "REP_LOCATION"
)

REP_REP = pick(
    rep,
    [
        "replicate",
        "rep",
        "biological_replicate",
    ],
    "REP_REPLICATE"
)

REP_D = pick(
    rep,
    [
        "D",
        "D_value",
        "tmd_specific_D",
        "delta_D",
    ],
    "REP_D"
)

REP_EREAL = pick(
    rep,
    [
        "E_real",
        "real_enrichment",
        "Ereal",
    ],
    "REP_E_REAL"
)

REP_EPSEUDO = pick(
    rep,
    [
        "E_pseudo",
        "pseudo_enrichment",
        "Epseudo",
    ],
    "REP_E_PSEUDO"
)

LOC_LOC = pick(
    loc,
    [
        "location",
        "compartment",
        "spatial_state",
    ],
    "LOC_LOCATION"
)

LOC_MEAN = pick(
    loc,
    [
        "D_location_mean",
        "D_mean",
        "mean_D",
        "location_mean_D",
    ],
    "LOC_MEAN"
)

LOC_LO = pick(
    loc,
    [
        "bootstrap_ci_2.5",
        "CI_2.5",
        "CI_low",
        "CI_lower",
        "bootstrap_CI_low",
    ],
    "LOC_CI_LOW"
)

LOC_HI = pick(
    loc,
    [
        "bootstrap_ci_97.5",
        "CI_97.5",
        "CI_high",
        "CI_upper",
        "bootstrap_CI_high",
    ],
    "LOC_CI_HIGH"
)

print("PRIMARY_SCHEMA_GATE = PASS")

# ============================================================
# NORMALIZE LOCATION NAMES
# ============================================================

def normalize_location(x):
    s = str(x).strip().lower()

    if s in {
        "er",
        "endoplasmic reticulum",
    }:
        return "ER"

    if s in {
        "mito",
        "mitochondria",
        "mitochondrial",
    }:
        return "Mitochondria"

    if s in {
        "cytosol",
        "cytosolic",
    }:
        return "Cytosol"

    return str(x)

rep["_location"] = (
    rep[REP_LOC]
    .map(normalize_location)
)

loc["_location"] = (
    loc[LOC_LOC]
    .map(normalize_location)
)

ORDER = [
    "ER",
    "Mitochondria",
    "Cytosol",
]

# ============================================================
# NUMERIC GATES
# ============================================================

for c in [
    REP_D,
    REP_EREAL,
    REP_EPSEUDO,
]:
    rep[c] = pd.to_numeric(
        rep[c]
    )

rep[REP_REP] = pd.to_numeric(
    rep[REP_REP]
).astype(int)

for c in [
    LOC_MEAN,
    LOC_LO,
    LOC_HI,
]:
    loc[c] = pd.to_numeric(
        loc[c]
    )

# exactly 2 replicate rows per location
for location in ORDER:
    z = rep[
        rep["_location"].eq(
            location
        )
    ]

    if len(z) != 2:
        raise SystemExit(
            f"REPLICATE_ACCOUNTING_GATE_FAIL:"
            f"{location}:{len(z)}"
        )

# exactly one location summary row
for location in ORDER:
    z = loc[
        loc["_location"].eq(
            location
        )
    ]

    if len(z) != 1:
        raise SystemExit(
            f"LOCATION_ACCOUNTING_GATE_FAIL:"
            f"{location}:{len(z)}"
        )

print("PRIMARY_ACCOUNTING_GATE = PASS")

# ============================================================
# FROZEN RESULT CONSISTENCY
# ============================================================

expected = {
    ("ER",1): 0.082889583159,
    ("ER",2): 0.095748374758,

    ("Mitochondria",1): 0.124859849072,
    ("Mitochondria",2): -0.039749358711,

    ("Cytosol",1): 0.100608092444,
    ("Cytosol",2): 0.067930823332,
}

for key,val in expected.items():

    location,r = key

    z = rep[
        rep["_location"].eq(location)
        &
        rep[REP_REP].eq(r)
    ]

    obs = float(
        z[REP_D].iloc[0]
    )

    if format(obs, ".12f") != format(val, ".12f"):
        raise SystemExit(
            f"FROZEN_D_GATE_FAIL:"
            f"{location}:rep{r}:"
            f"{obs}:{val}"
        )

expected_loc = {
    "ER":
        (
            0.089318978958,
            -0.018090483015,
            0.221254999092,
        ),

    "Mitochondria":
        (
            0.042555245181,
            -0.168335631922,
            0.261755265636,
        ),

    "Cytosol":
        (
            0.084269457888,
            -0.010151258823,
            0.171142696192,
        ),
}

for location,(m,lo,hi) in expected_loc.items():

    z = loc[
        loc["_location"].eq(
            location
        )
    ].iloc[0]

    vals = [
        float(z[LOC_MEAN]),
        float(z[LOC_LO]),
        float(z[LOC_HI]),
    ]

    exp = [
        m,lo,hi
    ]

    if [format(v, ".12f") for v in vals] != [format(v, ".12f") for v in exp]:
        raise SystemExit(
            f"FROZEN_LOCATION_GATE_FAIL:"
            f"{location}:{vals}:{exp}"
        )

assert loc["strong_support"].eq("NO").all(), "FROZEN_SUPPORT_GATE_FAIL"
print("PRIMARY_FROZEN_VALUE_GATE = PASS (exact to 12 decimal places)")

# ============================================================
# DERIVED DISPLAY TABLE
# ============================================================

display_rows = []

for location in ORDER:

    for r in [1,2]:

        z = rep[
            rep["_location"].eq(
                location
            )
            &
            rep[REP_REP].eq(r)
        ].iloc[0]

        display_rows.append({
            "location":
                location,

            "replicate":
                r,

            "E_real":
                float(
                    z[REP_EREAL]
                ),

            "E_pseudo":
                float(
                    z[REP_EPSEUDO]
                ),

            "D":
                float(
                    z[REP_D]
                ),
        })

DISPLAY = (
    OUTDIR /
    "FIGURE2_PRIMARY_DISPLAY_VALUES_v1.tsv"
)

pd.DataFrame(
    display_rows
).to_csv(
    DISPLAY,
    sep="\t",
    index=False
)

print(
    "FIGURE2_DISPLAY_TABLE_GATE = PASS"
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

COL_ER = "#009E73"
COL_MITO = "#D55E00"
COL_CYTO = "#0072B2"

COL_REAL = "#4D4D4D"
COL_PSEUDO = "#CC79A7"

COLORS = {
    "ER": COL_ER,
    "Mitochondria": COL_MITO,
    "Cytosol": COL_CYTO,
}

# ============================================================
# HELPERS
# ============================================================

def panel_label(
    ax,
    letter
):
    ax.text(
        -0.13,
        1.07,
        letter,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
        ha="left",
    )

def clean_axis(ax):
    ax.spines["top"].set_visible(
        False
    )
    ax.spines["right"].set_visible(
        False
    )
    ax.axhline(
        0,
        color="0.72",
        lw=0.8,
        zorder=0,
    )

# ============================================================
# FIGURE
# ============================================================

fig = plt.figure(
    figsize=(7.15,5.9),
    facecolor="white",
)

fig.patch.set_facecolor(
    "white"
)

fig.patch.set_alpha(
    1.0
)

gs = fig.add_gridspec(
    2,
    2,
    left=0.075,
    right=0.985,
    bottom=0.105,
    top=0.95,
    wspace=0.34,
    hspace=0.44,
)

axA = fig.add_subplot(
    gs[0,0]
)

axB = fig.add_subplot(
    gs[0,1]
)

axC = fig.add_subplot(
    gs[1,0]
)

axD = fig.add_subplot(
    gs[1,1]
)

for ax in [
    axA,
    axB,
    axC,
    axD,
]:
    ax.set_facecolor(
        "white"
    )

# ============================================================
# A — 12-LIBRARY DESIGN
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
    0.96,
    "GSE300977 spatial-state design",
    transform=axA.transAxes,
    ha="center",
    va="top",
    fontsize=9,
    fontweight="bold",
)

design = [
    (
        "ER",
        COL_ER,
        [
            ("Rep 1","B01 PD","B02 input"),
            ("Rep 2","B11 PD","B12 input"),
        ]
    ),
    (
        "Mitochondria",
        COL_MITO,
        [
            ("Rep 1","B09 PD","B10 input"),
            ("Rep 2","B07 PD","B08 input"),
        ]
    ),
    (
        "Cytosol",
        COL_CYTO,
        [
            ("Rep 1","B05 PD","B06 input"),
            ("Rep 2","B03 PD","B04 input"),
        ]
    ),
]

ys = [
    0.72,
    0.45,
    0.18,
]

for (
    location,
    color,
    reps
),y in zip(
    design,
    ys
):

    axA.text(
        0.19,
        y+0.085,
        location,
        transform=axA.transAxes,
        color=color,
        fontsize=6.5,
        fontweight="bold",
        ha="right",
    )

    for i,(
        repname,
        pdlib,
        inlib
    ) in enumerate(
        reps
    ):

        x = (
            0.22
            +
            i*0.39
        )

        box = FancyBboxPatch(
            (x,y),
            0.34,
            0.15,
            boxstyle=(
                "round,pad=0.012,"
                "rounding_size=0.015"
            ),
            transform=axA.transAxes,
            facecolor="white",
            edgecolor=color,
            linewidth=1.1,
        )

        axA.add_patch(
            box
        )

        axA.text(
            x+0.17,
            y+0.115,
            repname,
            transform=axA.transAxes,
            ha="center",
            va="center",
            fontsize=6.8,
            fontweight="bold",
        )

        axA.text(
            x+0.17,
            y+0.067,
            pdlib,
            transform=axA.transAxes,
            ha="center",
            va="center",
            fontsize=6.5,
        )

        axA.text(
            x+0.17,
            y+0.029,
            inlib,
            transform=axA.transAxes,
            ha="center",
            va="center",
            fontsize=6.5,
            color="0.35",
        )

axA.text(
    0.5,
    0.03,
    "3 spatial states × pulldown/input × 2 biological replicates",
    transform=axA.transAxes,
    ha="center",
    va="bottom",
    fontsize=6.8,
    color="0.3",
)

# ============================================================
# B — FROZEN PRIMARY ESTIMAND
# ============================================================

axB.axis(
    "off"
)

panel_label(
    axB,
    "B"
)

axB.text(
    0.5,
    0.96,
    "Preregistered +75-nt primary estimand",
    transform=axB.transAxes,
    ha="center",
    va="top",
    fontsize=9,
    fontweight="bold",
)

# Equivalent anchor-relative geometry on separate coordinate tracks.
# Identical mapping for both anchors: x = 0.22 + relative_nt / 150.
# The 31-nt half-open window [anchor+60, anchor+91) contains center +75.
for label, color, y in [
    ("Real TMD anchor (0)", COL_REAL, 0.73),
    ("Matched pseudo anchor (0)", COL_PSEUDO, 0.52),
]:
    axB.text(0.02, y+0.085, label, transform=axB.transAxes,
             fontsize=6.7, fontweight="bold", color=color, ha="left")
    axB.plot([0.16, 0.91], [y, y], transform=axB.transAxes,
             color="0.65", lw=1.0)
    axB.plot([0.22, 0.22], [y-0.025, y+0.035],
             transform=axB.transAxes, color=color, lw=1.5)
    axB.add_patch(Rectangle((0.62, y-0.025), 31/150, 0.05,
                            transform=axB.transAxes,
                            facecolor="#DCEEF8", edgecolor=COL_CYTO,
                            linewidth=0.8))
    axB.plot([0.72, 0.72], [y-0.035, y+0.048],
             transform=axB.transAxes, color=COL_CYTO, lw=1.0, ls=":")
    axB.text(0.72, y+0.067, "Center: +75 nt", transform=axB.transAxes,
             ha="center", fontsize=6.5, fontweight="bold", color=COL_CYTO)
    for x, text in [(0.62, "+60"), (0.22+91/150, "+91")]:
        axB.text(x, y-0.035, text, transform=axB.transAxes,
                 ha="center", va="top", fontsize=6)

axB.text(0.50, 0.38, "31-nt window: [anchor +60, anchor +91)",
         transform=axB.transAxes, ha="center", fontsize=6.5)
axB.text(0.50, 0.30, r"$E_{real}=\log_2(R_{PD}/R_{IN})$",
         transform=axB.transAxes, ha="center", fontsize=8)
axB.text(0.50, 0.22, r"$E_{pseudo}=\log_2(P_{PD}/P_{IN})$",
         transform=axB.transAxes, ha="center", fontsize=8)
axB.text(0.50, 0.12, r"$D=E_{real}-E_{pseudo}$",
         transform=axB.transAxes, ha="center", fontsize=9, fontweight="bold")
axB.text(0.50, 0.015,
         "Primary window, controls, read lengths and mapping rules\n"
         "were frozen before biological reveal.",
         transform=axB.transAxes, ha="center", va="bottom",
         fontsize=6.5, color="0.3")

# ============================================================
# C — REPLICATE-SPECIFIC PRIMARY D
# ============================================================

panel_label(
    axC,
    "C"
)

xbase = np.arange(
    len(ORDER)
)

offsets = [
    -0.10,
    +0.10,
]

markers = [
    "o",
    "s",
]

for i,location in enumerate(
    ORDER
):

    z = (
        rep[
            rep["_location"].eq(
                location
            )
        ]
        .sort_values(
            REP_REP
        )
    )

    for j,(_,r) in enumerate(
        z.iterrows()
    ):

        axC.scatter(
            i+offsets[j],
            float(
                r[REP_D]
            ),
            s=42,
            marker=markers[j],
            facecolors=(
                COLORS[location]
                if j == 0
                else "white"
            ),
            edgecolors=COLORS[
                location
            ],
            linewidths=1.2,
            zorder=4,
        )

    # connect biological replicates visually
    ys = z[
        REP_D
    ].to_numpy(float)

    axC.plot(
        [
            i+offsets[0],
            i+offsets[1],
        ],
        ys,
        color=COLORS[
            location
        ],
        lw=1.0,
        alpha=0.55,
        zorder=2,
    )

clean_axis(
    axC
)

axC.set_xticks(
    xbase
)

axC.set_xticklabels(
    ORDER
)

axC.set_ylabel(
    r"Primary $D$ at +75 nt (log$_2$)"
)

axC.set_title(
    "Replicate-specific primary effects"
)

axC.set_ylim(
    -0.10,
    0.18
)

# legend
axC.scatter(
    [],
    [],
    s=35,
    marker="o",
    facecolors="0.35",
    edgecolors="0.35",
    label="Replicate 1",
)

axC.scatter(
    [],
    [],
    s=35,
    marker="s",
    facecolors="white",
    edgecolors="0.35",
    label="Replicate 2",
)

axC.legend(
    frameon=False,
    loc="lower left",
)

# ============================================================
# D — LOCATION MEAN + BOOTSTRAP CI
# ============================================================

panel_label(
    axD,
    "D"
)

for i,location in enumerate(
    ORDER
):

    z = loc[
        loc["_location"].eq(
            location
        )
    ].iloc[0]

    mean = float(
        z[LOC_MEAN]
    )

    lo = float(
        z[LOC_LO]
    )

    hi = float(
        z[LOC_HI]
    )

    axD.errorbar(
        i,
        mean,
        yerr=np.array([
            [mean-lo],
            [hi-mean],
        ]),
        fmt="o",
        markersize=6,
        color=COLORS[
            location
        ],
        ecolor=COLORS[
            location
        ],
        elinewidth=1.5,
        capsize=4,
        capthick=1.2,
        zorder=4,
    )

    axD.text(
        i,
        hi+0.025,
        f"{mean:+.3f}",
        ha="center",
        va="bottom",
        fontsize=6.8,
        fontweight="bold",
        color=COLORS[
            location
        ],
    )

clean_axis(
    axD
)

axD.set_xlim(-0.4, 2.4)

axD.set_xticks(
    xbase
)

axD.set_xticklabels(
    ORDER
)

axD.set_ylabel(
    r"Mean primary $D$ (log$_2$)"
)

axD.set_title(
    "Preregistered primary conclusion"
)

axD.set_ylim(
    -0.36,
    0.32
)

axD.text(
    0.5,
    0.12,
    "No compartment met the frozen\n"
    "strong-support criterion",
    transform=axD.transAxes,
    ha="center",
    va="center",
    fontsize=8,
    fontweight="bold",
)

axD.text(
    0.5,
    0.025,
    "Primary conclusion retained unchanged",
    transform=axD.transAxes,
    ha="center",
    va="bottom",
    fontsize=6.5,
    color="0.35",
)

# ============================================================
# GLOBAL FOOTNOTE
# ============================================================

fig.text(
    0.5,
    0.022,
    "Primary +75-nt estimand was frozen before reveal; "
    "bootstrap intervals are unit-level and do not replace biological-replicate uncertainty.",
    ha="center",
    va="bottom",
    fontsize=6.1,
)

# ============================================================
# SAVE
# ============================================================

pdf = (
    OUTDIR /
    "FIGURE2_GSE300977_PRIMARY_v2.pdf"
)

svg = (
    OUTDIR /
    "FIGURE2_GSE300977_PRIMARY_v2.svg"
)

png = (
    OUTDIR /
    "FIGURE2_GSE300977_PRIMARY_v2.png"
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
# SHA
# ============================================================

for p in [
    DISPLAY,
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
    "FIGURE2_GENERATION_GATE = PASS"
)
