#!/usr/bin/env python3
"""Verify frozen NARGAB statistics and deterministic region TSV outputs."""
from pathlib import Path
import argparse
import hashlib
import math
import pandas as pd

ROOT = Path(__file__).resolve().parent
JOURNAL = ROOT.parent
SOURCE = JOURNAL / "Source_Data"
if not SOURCE.exists():
    SOURCE = JOURNAL / "source_data"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a, b, tol=1e-12):
    if not math.isclose(float(a), float(b), rel_tol=0, abs_tol=tol):
        raise AssertionError(f"numeric mismatch: {a} != {b}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--rerun-available", type=Path)
    p.add_argument("--rerun-complete", type=Path)
    args = p.parse_args()

    primary = pd.read_csv(SOURCE / "Figure3_Primary/PRIMARY_PLUS75_LOCATION_SUMMARY_v1.tsv", sep="\t")
    expected_primary = {
        "ER": (0.089319, -0.018090, 0.221255),
        "mitochondria": (0.042555, -0.168336, 0.261755),
        "cytosol": (0.084269, -0.010151, 0.171143),
    }
    for loc, vals in expected_primary.items():
        r = primary[primary.location == loc].iloc[0]
        for got, exp in zip((r.D_location_mean, r["bootstrap_ci_2.5"], r["bootstrap_ci_97.5"]), vals):
            close(got, exp, 1e-6)

    av = pd.read_csv(SOURCE / "Figure4_Topology/Region_Composite_Available_Case/REGION_COMPOSITE_SUMMARY_v1.tsv", sep="\t")
    exp = {
        "PLATEAU_90_120": (0.9325467977159678, 0.712214155436259, 1.1328308818519723, 1.073626423160913, 0.7914671722710228),
        "BROAD_60_120": (0.6518202467809296, 0.35735786295986155, 1.0554072349773076, 0.7840297741332183, 0.519610719428641),
    }
    for region, vals in exp.items():
        r = av[av.region == region].iloc[0]
        for got, target in zip((r.MITO_minus_ER_mean, r["interaction_CI_2.5"], r["interaction_CI_97.5"], r.MITO_minus_ER_rep1, r.MITO_minus_ER_rep2), vals):
            close(got, target)

    for rerun, frozen in [
        (args.rerun_available, SOURCE / "Figure4_Topology/Region_Composite_Available_Case"),
        (args.rerun_complete, SOURCE / "Figure4_Topology/Region_Composite_Complete_Case"),
    ]:
        if rerun:
            for name in ["REGION_CENTER_ESTIMATES_v1.tsv", "REGION_COMPOSITE_GENE_INFLUENCE_v1.tsv", "REGION_COMPOSITE_INFLUENCE_SUMMARY_v1.tsv", "REGION_COMPOSITE_SUMMARY_v1.tsv", "REGION_COMPOSITE_UNIT_INVENTORY_v1.tsv"]:
                if sha(rerun / name) != sha(frozen / name):
                    raise AssertionError(f"byte mismatch: {rerun/name}")
    print("PASS: frozen primary and region statistics verified; supplied rerun TSVs are byte-identical.")


if __name__ == "__main__":
    main()
