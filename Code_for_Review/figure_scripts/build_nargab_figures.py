#!/usr/bin/env python3
"""Rebuild NARGAB Figures 1, 2 and 4 from packaged frozen source tables.

Figure 3 is the unchanged frozen GSE300977 primary figure and is copied into the
submission package. This script performs presentation-only transformations; it
does not fit models, choose windows or recompute scientific source tables.
"""
from pathlib import Path
import json

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "Source_Data"
if not SRC.exists():
    SRC = ROOT / "source_data"
OUT = ROOT / "Figures"
EDIT = OUT / "Editable_and_PNG"
OUT.mkdir(exist_ok=True)
EDIT.mkdir(exist_ok=True)

COL = {"ER": "#009E73", "mitochondria": "#D55E00", "cytosol": "#0072B2",
       "FIRST": "#7B3294", "INTERNAL": "#555555"}


def save(fig, number):
    fig.savefig(OUT / f"Figure_{number}.pdf", bbox_inches="tight")
    fig.savefig(EDIT / f"Figure_{number}.svg", bbox_inches="tight")
    fig.savefig(EDIT / f"Figure_{number}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def panel(ax, letter):
    ax.text(-0.08, 1.04, letter, transform=ax.transAxes, fontsize=16,
            fontweight="bold", va="bottom")


def box(ax, xy, wh, text, edge="#777777", fill="white", fontsize=10,
        weight="normal"):
    x, y = xy; w, h = wh
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02",
                       linewidth=1.6, edgecolor=edge, facecolor=fill)
    ax.add_patch(p)
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fontsize,
            fontweight=weight)
    return p


def figure1():
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4))
    for ax in axes.ravel():
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax = axes[0, 0]; panel(ax, "A")
    ax.set_title("Equivalent real and pseudo estimands", fontweight="bold", pad=10)
    for y, c, lab in [(0.68, "#333333", "Real TMD anchor (0)"),
                      (0.34, "#CC79A7", "Matched pseudo anchor (0)")]:
        ax.plot([0.10, 0.92], [y, y], color="#BBBBBB", lw=2)
        ax.plot([0.24, 0.24], [y-0.07, y+0.07], color=c, lw=4)
        ax.add_patch(Rectangle((0.60, y-0.045), 0.22, 0.09,
                               facecolor="#E6F2FA", edgecolor="#0072B2", lw=1.8))
        ax.plot([0.71, 0.71], [y-0.07, y+0.07], color="#0072B2", ls=":", lw=2)
        ax.text(0.10, y+0.11, lab, color=c, fontweight="bold", fontsize=10)
        ax.text(0.60, y-0.12, "+60", ha="center", fontsize=9)
        ax.text(0.82, y-0.12, "+91", ha="center", fontsize=9)
    ax.text(0.5, 0.08, "D = log₂(real pulldown/input) − log₂(pseudo pulldown/input)",
            ha="center", fontsize=10)

    ax = axes[0, 1]; panel(ax, "B")
    ax.set_title("Falsification before interpretation", fontweight="bold", pad=10)
    items = [("Apparent TMD-relative profile", "#999999"),
             ("Matched pseudo anchors", "#D55E00"),
             ("Coordinate null", "#7B6DB0"),
             ("Deduplication / read-length sensitivity", "#0072B2"),
             ("Calibrated biological interpretation", "#009E73")]
    ys = [0.80, 0.62, 0.44, 0.26, 0.08]
    for i, ((txt, ec), y) in enumerate(zip(items, ys)):
        box(ax, (0.17, y), (0.66, 0.11), txt, edge=ec,
            fill="#F7F7F7" if i < 4 else "#EAF7F2", fontsize=9.5,
            weight="bold" if i == 4 else "normal")
        if i < 4:
            ax.annotate("", xy=(0.50, y-0.06), xytext=(0.50, y-0.005),
                        arrowprops=dict(arrowstyle="-|>", color="#777777", lw=1.5))

    ax = axes[1, 0]; panel(ax, "C")
    ax.set_title("Frozen evidence hierarchy", fontweight="bold", pad=10)
    rows = [("PRESPECIFIED PRIMARY", "+75 nt universal TMD-specific enrichment", "#B2182B"),
            ("PRESPECIFIED SECONDARY", "11-center profiles; FIRST vs INTERNAL", "#2166AC"),
            ("POST-REVEAL", "Common anchors; multi-pass; direct interaction", "#7B3294"),
            ("POST-REVEAL REGION", "+90/+105/+120 composite; +60…+120 sensitivity", "#7B3294")]
    for (head, txt, ec), y in zip(rows, [0.76, 0.55, 0.34, 0.13]):
        box(ax, (0.07, y), (0.86, 0.14), f"{head}\n{txt}", edge=ec,
            fill="#FAFAFA", fontsize=9, weight="normal")
    ax.text(0.5, 0.03, "Later analyses do not replace the negative primary result.",
            ha="center", fontsize=9, color="#555555")

    ax = axes[1, 1]; panel(ax, "D")
    ax.set_title("Interpretation boundary", fontweight="bold", pad=10)
    box(ax, (0.08, 0.58), (0.84, 0.24),
        "Permitted\npost-reveal, replicate-concordant association\ninternally stable under specified sensitivities",
        edge="#009E73", fill="#EAF7F2", fontsize=10)
    box(ax, (0.08, 0.18), (0.84, 0.24),
        "Not established\ntargeting • insertion • causality • universal mechanism\nvalidation • independent replication",
        edge="#B2182B", fill="#FBEDEE", fontsize=10)
    ax.annotate("evidence boundary", xy=(0.50, 0.46), xytext=(0.50, 0.52),
                ha="center", fontsize=9, color="#555555",
                arrowprops=dict(arrowstyle="-|>", color="#777777"))
    fig.suptitle("TMD-relative ribosome profiles: estimand, controls and hierarchy",
                 fontsize=16, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save(fig, 1)


def audit_boxes(ax, title, rows):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_title(title, fontweight="bold", pad=8)
    h = 0.12; y = 0.82
    for good, head, sub in rows:
        ec = "#4DAF4A" if good else "#E41A1C"
        ax.add_patch(FancyBboxPatch((0.08, y), 0.84, h,
                                    boxstyle="round,pad=0.012", facecolor="white",
                                    edgecolor="#BBBBBB", lw=1.2))
        ax.add_patch(Rectangle((0.11, y+0.026), 0.045, 0.068,
                               facecolor=ec, edgecolor=ec))
        ax.text(0.18, y+0.073, head, fontsize=9.5, fontweight="bold", va="center")
        ax.text(0.18, y+0.032, sub, fontsize=8.5, color="#555555", va="center")
        y -= 0.15


def figure2():
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.3))
    bsrc = SRC / "Figure2_Benchmarks"
    ss = pd.read_csv(bsrc / "sensitivity_summary.tsv", sep="\t")
    vals = [float(ss.loc[ss.analysis == "primary", "domain_mean_effect"].iloc[0]),
            float(ss.loc[ss.analysis == "umi_primary", "domain_mean_effect"].iloc[0]),
            0.579, -0.025]
    ax = axes[0, 0]; panel(ax, "A")
    labs = ["Primary\nCHX", "UMI\ndedup.", "Matched\npseudo", "No-CHX"]
    colors = ["#555555", "#0072B2", "#D55E00", "#009E73"]
    ax.bar(range(4), vals, color=colors, width=0.68)
    ax.axhline(0, color="#999999", lw=1)
    ax.set_xticks(range(4), labs); ax.set_ylabel("Domain-mean effect (Δlog₂)")
    ax.set_title("GSE106329: apparent signal fails controls", fontweight="bold")
    for i, v in enumerate(vals):
        ax.text(i, v + (0.018 if v >= 0 else 0.09), f"{v:+.3f}", ha="center",
                va="bottom", fontsize=9)
    ax.annotate("89% attenuation", xy=(1, vals[1]), xytext=(0.72, 0.32),
                color="#0072B2", ha="center",
                arrowprops=dict(arrowstyle="->", color="#0072B2"))

    ax = axes[0, 1]; panel(ax, "B")
    audit_boxes(ax, "GSE106329 frozen-criterion audit", [
        (True, "Nominal broad band", "131/131 centers"),
        (False, "UMI sensitivity", "89% attenuation"),
        (False, "Matched pseudo", "pseudo effect > real effect"),
        (False, "Coordinate null", "empirical P = 0.222 / 0.207"),
        (False, "No-CHX comparison", "no spatial band")])

    ax = axes[1, 0]; panel(ax, "C")
    pdv = pd.read_csv(bsrc / "Panel_D_frozen_values.tsv", sep="\t").iloc[0]
    med = float(str(pdv["median_D"])); lo = float(str(pdv["HL_90_CI_lower"])); hi = float(str(pdv["HL_90_CI_upper"]))
    hl = json.loads((bsrc / "Complete_GSE133111_Summaries/primary_all_frozen_lengths/real_pseudo_summary.json").read_text())["hodges_lehmann"]
    ax.errorbar([hl], [0], xerr=[[hl-lo], [hi-hl]], fmt="o", color="#555555",
                markersize=8, capsize=5, lw=2)
    ax.axvline(0, color="#999999", ls="--")
    ax.set_xlim(-0.12, 0.15); ax.set_ylim(-0.6, 0.6); ax.set_yticks([])
    ax.set_xlabel("Hodges–Lehmann real − pseudo location")
    ax.set_title("GSE133111: matched-pseudo test", fontweight="bold")
    ax.text(0.02, 0.38, "REAL vs PSEUDO: FAIL", transform=ax.transAxes,
            color="#E41A1C", fontweight="bold")
    ax.text(0.5, 0.12, "median D = +0.020\none-sided Wilcoxon P = 0.153",
            transform=ax.transAxes, ha="center", fontsize=9)

    ax = axes[1, 1]; panel(ax, "D")
    audit_boxes(ax, "GSE133111 formal Stage 1 audit", [
        (False, "Real vs pseudo", "median D = +0.0197"),
        (True, "Coordinate null", "empirical P = 0.00010"),
        (False, "Localized downstream band", "Wilcoxon P = 0.738"),
        (True, "No global shift", "90% Q interval within ±log₂(1.20)"),
        (False, "Transcript / gene LOO", "mandatory gate failed")])
    fig.text(0.5, 0.02,
             "Strong profiles or coordinate-null separation do not by themselves establish TMD specificity.",
             ha="center", fontsize=10, color="#444444")
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    save(fig, 2)


def line_with_band(ax, d, group_col, groups, title, colors):
    for g in groups:
        z = d[d[group_col] == g].sort_values("center")
        ax.plot(z.center, z.D_mean, marker="o", lw=2, color=colors[g], label=g.title())
        ax.fill_between(z.center, z["CI_2.5"], z["CI_97.5"], color=colors[g], alpha=.14)
    ax.axhline(0, color="#999999", lw=1); ax.axvline(0, color="#BBBBBB", ls="--", lw=1)
    ax.set_title(title); ax.set_xlabel("Position relative to TMD start (nt)")
    ax.set_ylabel("D (log₂)"); ax.legend(frameon=False, fontsize=8)


def figure4():
    tsrc = SRC / "Figure4_Topology"
    fi = pd.read_csv(tsrc / "FIRST_INTERNAL_11CENTER_PROFILE_v1.tsv", sep="\t")
    mp = pd.read_csv(tsrc / "MULTIPASS_FIRST_MINUS_INTERNAL_11CENTER_v1.tsv", sep="\t")
    inter = pd.read_csv(tsrc / "MULTIPASS_TOPOLOGY_COMPARTMENT_INTERACTIONS_v1.tsv", sep="\t")
    av = pd.read_csv(tsrc / "Region_Composite_Available_Case/REGION_COMPOSITE_SUMMARY_v1.tsv", sep="\t")
    cc = pd.read_csv(tsrc / "Region_Composite_Complete_Case/REGION_COMPOSITE_SUMMARY_v1.tsv", sep="\t")
    infl = pd.read_csv(tsrc / "Region_Composite_Available_Case/REGION_COMPOSITE_INFLUENCE_SUMMARY_v1.tsv", sep="\t")

    fig = plt.figure(figsize=(14, 10.5))
    gs = fig.add_gridspec(2, 3, wspace=.60, hspace=.48)

    ax = fig.add_subplot(gs[0, 0]); panel(ax, "A"); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("Topology class and evidence tier", fontweight="bold")
    ax.plot([.14,.88],[.73,.73],color="#777777",lw=2)
    for x,c in [(.26,COL["FIRST"]),(.50,"#555555"),(.68,"#555555"),(.82,"#555555")]:
        ax.add_patch(Rectangle((x-.025,.60),.05,.25,color=c))
    ax.text(.26,.54,"FIRST",ha="center",color=COL["FIRST"],fontweight="bold")
    ax.text(.68,.54,"INTERNAL",ha="center",color="#555555",fontweight="bold")
    box(ax,(.08,.25),(.84,.15),"Prespecified stratification",edge="#2166AC",fill="#F7FAFD",fontsize=9)
    box(ax,(.08,.04),(.84,.15),"Post-reveal\nmulti-pass and interaction\nregion summaries",edge="#7B3294",fill="#FAF5FB",fontsize=8.2)

    ax = fig.add_subplot(gs[0, 1]); panel(ax, "B")
    for loc, ls in [("ER","-"),("mitochondria","--")]:
        for cls in ["FIRST","INTERNAL"]:
            z=fi[(fi.location==loc)&(fi.TMD_class==cls)].sort_values("center")
            ax.plot(z.center,z.D_mean,color=COL[cls],ls=ls,lw=2,
                    label=f"{loc if loc == 'ER' else loc.title()} {cls.title()}")
    ax.axhline(0,color="#999999",lw=1); ax.axvline(75,color="#999999",ls=":")
    ax.set_title("Prespecified class profiles")
    ax.set_xlabel("Position (nt)"); ax.set_ylabel("D (log₂)")
    ax.legend(frameon=True,facecolor="white",framealpha=.9,fontsize=8,ncol=1)

    ax = fig.add_subplot(gs[0, 2]); panel(ax, "C")
    for loc in ["ER","mitochondria"]:
        z=mp[mp.location==loc].sort_values("center")
        label = "ER" if loc == "ER" else "Mitochondria"
        ax.plot(z.center,z.FIRST_minus_INTERNAL_mean,marker="o",lw=2,color=COL[loc],label=label)
        ax.fill_between(z.center,z["CI_2.5"],z["CI_97.5"],color=COL[loc],alpha=.15)
    ax.axhline(0,color="#999999",lw=1); ax.axvline(75,color="#999999",ls=":")
    ax.set_title("Post-reveal multi-pass contrast")
    ax.set_xlabel("Position (nt)"); ax.set_ylabel("FIRST − INTERNAL D")
    ax.legend(frameon=False,fontsize=8)

    ax = fig.add_subplot(gs[1, 0]); panel(ax, "D")
    z=inter[inter.interaction=="MITO_MINUS_ER"].sort_values("center")
    ax.plot(z.center,z["mean"],marker="o",lw=2.3,color=COL["mitochondria"],label="Mean")
    ax.fill_between(z.center,z["CI_2.5"],z["CI_97.5"],color=COL["mitochondria"],alpha=.16)
    ax.scatter(z.center,z.rep1,s=18,color=COL["mitochondria"],alpha=.45,label="Replicate 1")
    ax.scatter(z.center,z.rep2,s=24,facecolor="white",edgecolor=COL["mitochondria"],label="Replicate 2")
    ax.axhline(0,color="#999999",lw=1); ax.axvline(75,color="#999999",ls=":")
    ax.axvspan(90,120,color="#F4A261",alpha=.10)
    ax.set_title("Post-reveal MITO − ER interaction")
    ax.set_xlabel("Position (nt)"); ax.set_ylabel("Interaction (log₂)")
    ax.legend(frameon=False,fontsize=7)

    ax = fig.add_subplot(gs[1, 1]); panel(ax, "E")
    names=["+90/+105/+120\navailable", "+60…+120\navailable", "+90/+105/+120\ncomplete", "+60…+120\ncomplete"]
    rows=[av.iloc[0],av.iloc[1],cc.iloc[0],cc.iloc[1]]
    for i,r in enumerate(rows):
        c=COL["mitochondria"] if i<2 else "#888888"
        ax.errorbar(i,r.MITO_minus_ER_mean,
                    yerr=[[r.MITO_minus_ER_mean-r["interaction_CI_2.5"]],[r["interaction_CI_97.5"]-r.MITO_minus_ER_mean]],
                    fmt="o",ms=7,color=c,capsize=4,lw=2)
        ax.scatter([i-.08,i+.08],[r.MITO_minus_ER_rep1,r.MITO_minus_ER_rep2],s=22,
                   facecolor=[c,"white"],edgecolor=c,zorder=3)
    ax.axhline(0,color="#999999",lw=1)
    ax.set_xticks(range(4),names,rotation=25,ha="right",fontsize=8)
    ax.set_ylabel("Region interaction (log₂)")
    ax.set_title("Gene-clustered region summaries")
    ax.text(.02,.96,"Available-case is the main estimand",transform=ax.transAxes,va="top",fontsize=7.5)

    ax = fig.add_subplot(gs[1, 2]); panel(ax, "F")
    y=np.arange(len(infl))
    ax.barh(y,infl.max_relative_change*100,color=["#7B3294","#B07CC6"])
    ax.set_yticks(y,["+90/+105/+120","+60…+120"])
    ax.invert_yaxis(); ax.set_xlabel("Maximum relative LOO change (%)")
    ax.set_title("Leave-one-gene-out influence")
    for i,r in infl.iterrows():
        ax.text(r.max_relative_change*100+.6,i,f"{r.max_relative_change*100:.1f}%\nno sign reversal",va="center",fontsize=8)
    ax.set_xlim(0,45)
    ax.text(.02,.05,"Broader-region maximum deletion\nincreased the estimate.",transform=ax.transAxes,fontsize=8,color="#555555",bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=2))

    fig.suptitle("FIRST-versus-INTERNAL topology and compartment-specific ribosome organization",
                 fontsize=15,fontweight="bold",y=.995)
    fig.text(.5,.012,"All interaction and region panels are post-reveal; gene bootstrap does not create biological replication.",
             ha="center",fontsize=9,color="#444444")
    fig.tight_layout(rect=[0,.035,1,.965])
    save(fig,4)


if __name__ == "__main__":
    plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,
                         "axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    figure1(); figure2(); figure4()
    print("Wrote Figures 1, 2 and 4 from packaged frozen source tables.")
