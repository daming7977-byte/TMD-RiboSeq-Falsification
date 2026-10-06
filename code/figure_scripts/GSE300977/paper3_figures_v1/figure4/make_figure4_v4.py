from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ============================================================
# PATHS
# ============================================================

PROFILE = Path(
    "tmd_spatial_first_internal_v1/output/"
    "FIRST_INTERNAL_11CENTER_PROFILE_v1.tsv"
)

MP_CONTRAST = Path(
    "tmd_spatial_multipass_first_internal_v1/output/"
    "MULTIPASS_FIRST_MINUS_INTERNAL_11CENTER_v1.tsv"
)

INTERACTION = Path(
    "tmd_spatial_topology_compartment_interaction_v1/output/"
    "MULTIPASS_TOPOLOGY_COMPARTMENT_INTERACTIONS_v1.tsv"
)

OUTDIR = Path("paper3_figures_v1/figure4")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# LOAD + SCHEMA GATES
# ============================================================

prof = pd.read_csv(PROFILE, sep="\t")
mp = pd.read_csv(MP_CONTRAST, sep="\t")
inter = pd.read_csv(INTERACTION, sep="\t")

req_prof = {
    "center","TMD_class","N_units","location",
    "D_rep1","D_rep2","D_mean",
    "CI_2.5","CI_97.5"
}

req_mp = {
    "center","location",
    "FIRST_minus_INTERNAL_mean",
    "CI_2.5","CI_97.5"
}

req_inter = {
    "center","interaction",
    "rep1","rep2","mean",
    "CI_2.5","CI_97.5"
}

for name,df,req in [
    ("PROFILE",prof,req_prof),
    ("MULTIPASS",mp,req_mp),
    ("INTERACTION",inter,req_inter),
]:
    missing = req - set(df.columns)
    if missing:
        raise SystemExit(
            f"{name}_SCHEMA_GATE_FAIL:{sorted(missing)}"
        )
    print(name, "SCHEMA_GATE = PASS")

CENTERS = [-30,-15,0,15,30,45,60,75,90,105,120]

# ============================================================
# STYLE
# ============================================================

plt.rcParams.update({
    "font.family": "Arial",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.transparent": False,
    "font.size": 7.5,
    "axes.labelsize": 8,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

FIRST_COLOR = "#7B3294"
INTERNAL_COLOR = "#4D4D4D"

ER_COLOR = "#009E73"
MITO_COLOR = "#D55E00"
CYTO_COLOR = "#0072B2"

# ============================================================
# HELPERS
# ============================================================

def clean_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.axhline(0, color="0.72", lw=0.8, zorder=0)
    ax.axvline(0, color="0.88", lw=0.7, ls="--", zorder=0)
    ax.set_xticks(CENTERS)
    ax.set_xticklabels(
        [f"{x:+d}" if x != 0 else "0" for x in CENTERS],
        rotation=45,
        ha="right"
    )
    ax.set_xlim(-36,126)

def panel_label(ax, letter):
    ax.text(
        -0.13, 1.08, letter,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
        ha="left"
    )

def plot_class_profile(ax, location):

    z = prof[prof["location"].eq(location)].copy()

    for cls,color,label in [
        ("FIRST",FIRST_COLOR,"First TMD"),
        ("INTERNAL",INTERNAL_COLOR,"Internal TMD"),
    ]:

        q = (
            z[z["TMD_class"].eq(cls)]
            .sort_values("center")
        )

        x = q["center"].to_numpy(float)
        y = q["D_mean"].to_numpy(float)
        lo = q["CI_2.5"].to_numpy(float)
        hi = q["CI_97.5"].to_numpy(float)

        ax.fill_between(
            x,lo,hi,
            color=color,
            alpha=0.14,
            linewidth=0
        )

        ax.plot(
            x,y,
            color=color,
            marker="o",
            markersize=3.2,
            label=label,
            zorder=3
        )

        # biological replicate values
        ax.scatter(
            x,
            q["D_rep1"],
            s=8,
            color=color,
            alpha=0.38,
            linewidths=0,
            zorder=2
        )

        ax.scatter(
            x,
            q["D_rep2"],
            s=8,
            facecolors="white",
            edgecolors=color,
            linewidths=0.55,
            alpha=0.65,
            zorder=2
        )

    clean_axis(ax)

    # Shared scale across ER / mitochondria / cytosol
    # for direct visual comparability.
    ax.set_ylim(-0.80, 0.90)

    ax.set_xlabel("Position relative to TMD start (nt)")
    ax.set_ylabel(
        r"TMD-specific effect $D$ (log$_2$)"
    )

# ============================================================
# FIGURE
# ============================================================

fig = plt.figure(
    figsize=(7.15,7.25),
    facecolor="white"
)

fig.patch.set_facecolor("white")
fig.patch.set_alpha(1.0)

gs = fig.add_gridspec(
    2,3,
    left=0.08,
    right=0.985,
    bottom=0.125,
    top=0.965,
    wspace=0.38,
    hspace=0.42
)

axA = fig.add_subplot(gs[0,0])
axB = fig.add_subplot(gs[0,1])
axC = fig.add_subplot(gs[0,2])
axD = fig.add_subplot(gs[1,0])
axE = fig.add_subplot(gs[1,1])
axF = fig.add_subplot(gs[1,2])

for _ax in [axA, axB, axC, axD, axE, axF]:
    _ax.set_facecolor("white")
    _ax.patch.set_alpha(1.0)

# ============================================================
# A — SCHEMATIC
# ============================================================

axA.set_xlim(0,10)
axA.set_ylim(0,10)
axA.axis("off")

# membrane bands
# single membrane schematic
axA.add_patch(
    Rectangle(
        (0.8,6.87),
        8.3,0.26,
        facecolor="0.86",
        edgecolor="none"
    )
)

# top multipass protein
xs = [2.0,4.4,6.7,8.0]

for i,x in enumerate(xs):
    color = FIRST_COLOR if i == 0 else INTERNAL_COLOR
    axA.add_patch(
        Rectangle(
            (x-0.24,5.9),
            0.48,2.2,
            facecolor=color,
            edgecolor="none",
            alpha=0.95
        )
    )

axA.plot(
    [2.0,4.4,6.7,8.0],
    [8.35,8.35,8.35,8.35],
    color="0.35",
    lw=1.1
)

axA.text(
    2.0,5.48,
    "FIRST",
    color=FIRST_COLOR,
    ha="center",
    va="top",
    fontweight="bold"
)

axA.text(
    6.7,5.48,
    "INTERNAL",
    color=INTERNAL_COLOR,
    ha="center",
    va="top",
    fontweight="bold"
)

axA.annotate(
    "",
    xy=(2.0,5.95),
    xytext=(2.0,4.65),
    arrowprops=dict(
        arrowstyle="-|>",
        lw=1,
        color=FIRST_COLOR
    )
)

axA.annotate(
    "",
    xy=(6.7,5.95),
    xytext=(6.7,4.65),
    arrowprops=dict(
        arrowstyle="-|>",
        lw=1,
        color=INTERNAL_COLOR
    )
)

axA.text(
    4.9,9.25,
    "TMD order in a multi-pass protein",
    ha="center",
    va="center",
    fontsize=8.5,
    fontweight="bold"
)

axA.text(
    4.9,1.8,
    "Question:\nDoes TMD order interact with\nribosome spatial state?",
    ha="center",
    va="center",
    fontsize=8
)

panel_label(axA,"A")

# ============================================================
# B–D — PREREGISTERED FIRST vs INTERNAL
# ============================================================

plot_class_profile(axB,"ER")
axB.set_title("ER")
panel_label(axB,"B")

plot_class_profile(axC,"mitochondria")
axC.set_title("Mitochondria")
panel_label(axC,"C")

plot_class_profile(axD,"cytosol")
axD.set_title("Cytosol")
panel_label(axD,"D")

axB.legend(
    frameon=False,
    loc="upper left"
)

# ============================================================
# E — MULTI-PASS-ONLY FIRST − INTERNAL
# ============================================================

for loc,color,label in [
    ("ER",ER_COLOR,"ER"),
    ("mitochondria",MITO_COLOR,"Mitochondria"),
    ("cytosol",CYTO_COLOR,"Cytosol"),
]:

    z = (
        mp[mp["location"].eq(loc)]
        .sort_values("center")
    )

    x = z["center"].to_numpy(float)
    y = z["FIRST_minus_INTERNAL_mean"].to_numpy(float)
    lo = z["CI_2.5"].to_numpy(float)
    hi = z["CI_97.5"].to_numpy(float)

    axE.fill_between(
        x,lo,hi,
        color=color,
        alpha=0.10,
        linewidth=0
    )

    axE.plot(
        x,y,
        color=color,
        marker="o",
        markersize=3.0,
        label=label
    )

clean_axis(axE)

axE.set_title(
    "Multi-pass-only topology contrast"
)

axE.set_xlabel(
    "Position relative to TMD start (nt)"
)

axE.set_ylabel(
    r"First $-$ internal $D$ (log$_2$)"
)

axE.legend(
    frameon=False,
    loc="lower right"
)

panel_label(axE,"E")

# ============================================================
# F — FOCAL DIRECT TOPOLOGY × COMPARTMENT INTERACTION
# ============================================================

z = (
    inter[
        inter["interaction"].eq(
            "MITO_MINUS_ER"
        )
    ]
    .sort_values("center")
)

x = z["center"].to_numpy(float)
y = z["mean"].to_numpy(float)
lo = z["CI_2.5"].to_numpy(float)
hi = z["CI_97.5"].to_numpy(float)

axF.fill_between(
    x,
    lo,
    hi,
    color=MITO_COLOR,
    alpha=0.15,
    linewidth=0,
)

axF.plot(
    x,
    y,
    color=MITO_COLOR,
    marker="o",
    markersize=3.8,
    linewidth=1.8,
    label="Mitochondria − ER",
    zorder=4,
)

# Biological-replicate interaction values
axF.scatter(
    x,
    z["rep1"],
    s=10,
    color=MITO_COLOR,
    alpha=0.32,
    linewidths=0,
    zorder=3,
)

axF.scatter(
    x,
    z["rep2"],
    s=10,
    facecolors="white",
    edgecolors=MITO_COLOR,
    linewidths=0.65,
    alpha=0.75,
    zorder=3,
)

# Emphasize the downstream estimates that constitute
# the clearest interaction pattern, without significance stars.
for center in [90,105,120]:

    q = z[
        z["center"].eq(center)
    ].iloc[0]

    axF.scatter(
        [center],
        [q["mean"]],
        s=30,
        facecolors="white",
        edgecolors=MITO_COLOR,
        linewidths=1.25,
        zorder=6,
    )

clean_axis(axF)

axF.set_ylim(
    -1.05,
    1.48
)

axF.set_title(
    "Direct topology × compartment interaction"
)

axF.set_xlabel(
    "Position relative to TMD start (nt)"
)

axF.set_ylabel(
    "Mito − ER topology interaction\n"
    r"(log$_2$)"
)

axF.legend(
    frameon=False,
    loc="upper left"
)

# Compact numerical annotation
axF.text(
    0.97,
    0.05,
    "+90   +0.77\n"
    "+105  +0.97\n"
    "+120  +1.06",
    transform=axF.transAxes,
    ha="right",
    va="bottom",
    fontsize=6.7,
)

axF.text(
    -0.16,
    1.09,
    "F",
    transform=axF.transAxes,
    fontsize=11,
    fontweight="bold",
    va="top",
    ha="left"
)

# ============================================================
# GLOBAL NOTE
# ============================================================

fig.text(
    0.5,
    0.025,
    "Shaded regions: pointwise 95% bootstrap intervals; "
    "positional intervals are not multiplicity-corrected.",
    ha="center",
    va="bottom",
    fontsize=6.2
)

# ============================================================
# SAVE
# ============================================================

pdf = OUTDIR / "FIGURE4_TOPOLOGY_SPATIAL_STATE_v4.pdf"
svg = OUTDIR / "FIGURE4_TOPOLOGY_SPATIAL_STATE_v4.svg"
png = OUTDIR / "FIGURE4_TOPOLOGY_SPATIAL_STATE_v4.png"

fig.savefig(
    pdf,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="white",
    transparent=False
)

fig.savefig(
    svg,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="white",
    transparent=False
)

fig.savefig(
    png,
    dpi=450,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="white",
    transparent=False
)

plt.close(fig)

print()
print("FIGURE4_PDF =", pdf)
print("FIGURE4_SVG =", svg)
print("FIGURE4_PNG =", png)

# ============================================================
# SHA
# ============================================================

for p in [pdf,svg,png]:

    sha = hashlib.sha256(
        p.read_bytes()
    ).hexdigest()

    print(
        "SHA256",
        sha,
        p
    )

print()
print("FIGURE4_GENERATION_GATE = PASS")
