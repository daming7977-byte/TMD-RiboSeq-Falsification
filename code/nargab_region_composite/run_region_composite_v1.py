from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
import sys

import numpy as np
import pandas as pd


CENTERS_BY_REGION = {
    "PLATEAU_90_120": [90, 105, 120],
    "BROAD_60_120": [60, 75, 90, 105, 120],
}

SEEDS = {
    "PLATEAU_90_120": 30097890,
    "BROAD_60_120": 30097860,
}

LIBRARIES = [f"B{i:02d}" for i in range(1, 13)]

CONTRASTS = [
    ("ER", 1, "B01", "B02"),
    ("ER", 2, "B11", "B12"),
    ("mitochondria", 1, "B09", "B10"),
    ("mitochondria", 2, "B07", "B08"),
]

CLASSES = ["FIRST", "INTERNAL"]
BOOT_N = 10_000
BOOT_BATCH = 200


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--frozen-project",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "frozen_intermediates"
        / "region_composite_project",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "output",
    )
    parser.add_argument(
        "--population-mode",
        choices=["COMPLETE_CASE", "AVAILABLE_CASE"],
        default="COMPLETE_CASE",
    )
    return parser.parse_args()


def split_delimited(value: object, delimiter: str) -> list[str]:
    return [item for item in str(value).split(delimiter) if item]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def percentile_interval(values: np.ndarray) -> tuple[float, float]:
    return (
        float(np.percentile(values, 2.5)),
        float(np.percentile(values, 97.5)),
    )


def main() -> None:
    args = parse_args()
    root = args.frozen_project.resolve()
    outdir = args.output.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    formal_path = root / "tmd_anchor_preregister/bridge/FORMAL_TMD_ANCHOR_CANDIDATES_v1p1.tsv"
    pair_path = root / "crossmap_window_preregister/output/GSE300977_PAIRED_REAL_PSEUDO_CROSSMAP_MASK_AUDIT_v1.tsv"
    catalog_path = root / "tmd_spatial_secondary_scan_v1/output/SECONDARY_11CENTER_GEOMETRY_CATALOG_v1.tsv"
    count_dir = root / "tmd_spatial_secondary_scan_v1/geometry_counts"

    required_paths = [formal_path, pair_path, catalog_path]
    required_paths.extend(
        count_dir / f"{lib}.SECONDARY_GEOMETRY_DENSITIES_v1.tsv"
        for lib in LIBRARIES
    )
    missing = [str(path) for path in required_paths if not path.exists()]
    if missing:
        raise SystemExit(f"MISSING_INPUTS:{missing}")

    formal = pd.read_csv(formal_path, sep="\t", dtype=str)
    transcript_map = (
        formal[["target_transcript", "pass_class"]]
        .dropna()
        .drop_duplicates()
    )
    conflict = transcript_map.groupby("target_transcript")["pass_class"].nunique()
    if bool((conflict > 1).any()):
        raise SystemExit("TRANSCRIPT_PASS_CLASS_CONFLICT")
    tx_to_pass = transcript_map.set_index("target_transcript")["pass_class"].to_dict()

    pairs = pd.read_csv(pair_path, sep="\t", dtype=str)
    pairs["center"] = pd.to_numeric(pairs["center"]).astype(int)
    pairs["real_crossmap_unmasked_nt"] = pd.to_numeric(
        pairs["real_crossmap_unmasked_nt"]
    )
    pairs["pseudo_min_unmasked_nt"] = pd.to_numeric(
        pairs["pseudo_min_unmasked_nt"]
    )
    pairs["n_gene_ids"] = pd.to_numeric(pairs["n_gene_ids"]).astype(int)
    pairs = pairs[
        (pairs["real_crossmap_unmasked_nt"] >= 1)
        & (pairs["pseudo_min_unmasked_nt"] >= 1)
    ].copy()

    def topology_class(row: pd.Series) -> str:
        transcripts = split_delimited(row["target_transcripts"], ",")
        classes = sorted({tx_to_pass[tx] for tx in transcripts})
        if classes == ["MULTI_PASS"]:
            return "MULTI_PASS"
        if classes == ["SINGLE_PASS"]:
            return "SINGLE_PASS"
        return "MIXED_TOPOLOGY"

    pairs["topology_class"] = pairs.apply(topology_class, axis=1)
    pairs = pairs[
        pairs["topology_class"].eq("MULTI_PASS")
        & pairs["position_class_consensus"].isin(CLASSES)
        & pairs["unique_gene_assignment"].eq("YES")
        & pairs["n_gene_ids"].eq(1)
    ].copy()
    pairs["gene_id"] = pairs["gene_ids"]

    key_columns = [
        "chrom",
        "strand",
        "real_anchor_genomic_start0",
        "target_transcripts",
        "gene_ids",
        "position_class_consensus",
        "tmd_indices",
    ]
    pairs["unit_key"] = pairs[key_columns].astype(str).agg("|".join, axis=1)
    duplicated = pairs.duplicated(["center", "unit_key"], keep=False)
    if bool(duplicated.any()):
        examples = pairs.loc[duplicated, ["center", "unit_key"]].head().to_dict("records")
        raise SystemExit(f"CENTER_UNIT_KEY_DUPLICATION:{examples}")

    catalog = pd.read_csv(catalog_path, sep="\t", dtype=str)
    catalog_keys = list(
        zip(catalog["strand"], catalog["geometry"], catalog["mask_bits"])
    )
    if len(catalog_keys) != len(set(catalog_keys)):
        raise SystemExit("GEOMETRY_CATALOG_KEY_DUPLICATION")
    key_to_gid = dict(zip(catalog_keys, catalog["geometry_id"].astype(str)))

    densities: dict[str, dict[str, float]] = {}
    for lib in LIBRARIES:
        frame = pd.read_csv(
            count_dir / f"{lib}.SECONDARY_GEOMETRY_DENSITIES_v1.tsv",
            sep="\t",
        )
        densities[lib] = dict(
            zip(frame["geometry_id"].astype(str), frame["density"].astype(float))
        )

    matrix_columns: list[str] = []
    for lib in LIBRARIES:
        matrix_columns.extend([f"{lib}_real", f"{lib}_pseudo"])
    column_index = {name: index for index, name in enumerate(matrix_columns)}

    def build_matrix(frame: pd.DataFrame) -> np.ndarray:
        rows: list[list[float]] = []
        for _, row in frame.iterrows():
            real_key = (
                str(row["strand"]),
                str(row["real_window_geometry"]),
                str(row["real_crossmap_mask_bits"]),
            )
            if real_key not in key_to_gid:
                raise SystemExit("REAL_GEOMETRY_NOT_IN_CATALOG")
            real_gid = key_to_gid[real_key]
            pseudo_geometries = split_delimited(row["pseudo_window_geometries"], ";")
            pseudo_masks = split_delimited(row["pseudo_crossmap_mask_bits"], ";")
            if len(pseudo_geometries) != len(pseudo_masks):
                raise SystemExit("PSEUDO_GEOMETRY_MASK_CARDINALITY_FAIL")
            pseudo_gids = []
            for geometry, mask in zip(pseudo_geometries, pseudo_masks):
                pseudo_key = (str(row["strand"]), geometry, mask)
                if pseudo_key not in key_to_gid:
                    raise SystemExit("PSEUDO_GEOMETRY_NOT_IN_CATALOG")
                pseudo_gids.append(key_to_gid[pseudo_key])
            values: list[float] = []
            for lib in LIBRARIES:
                values.append(float(densities[lib][real_gid]))
                values.append(
                    float(np.mean([densities[lib][gid] for gid in pseudo_gids]))
                )
            rows.append(values)
        matrix = np.asarray(rows, dtype=float)
        if not np.isfinite(matrix).all() or bool((matrix < 0).any()):
            raise SystemExit("MATRIX_VALUE_GATE_FAIL")
        return matrix

    def d_vector(means: np.ndarray, pd_lib: str, input_lib: str) -> np.ndarray:
        real_pd = means[:, column_index[f"{pd_lib}_real"]]
        real_input = means[:, column_index[f"{input_lib}_real"]]
        pseudo_pd = means[:, column_index[f"{pd_lib}_pseudo"]]
        pseudo_input = means[:, column_index[f"{input_lib}_pseudo"]]
        if bool(
            (real_pd <= 0).any()
            or (real_input <= 0).any()
            or (pseudo_pd <= 0).any()
            or (pseudo_input <= 0).any()
        ):
            raise SystemExit("ZERO_AGGREGATE_GATE_FAIL")
        return np.log2(real_pd / real_input) - np.log2(pseudo_pd / pseudo_input)

    summary_rows: list[dict[str, object]] = []
    center_rows: list[dict[str, object]] = []
    influence_frames: list[pd.DataFrame] = []
    inventory_frames: list[pd.DataFrame] = []
    distribution_rows: list[dict[str, object]] = []

    for region_name, centers in CENTERS_BY_REGION.items():
        regional = pairs[pairs["center"].isin(centers)].copy()
        if args.population_mode == "COMPLETE_CASE":
            center_counts = regional.groupby("unit_key")["center"].nunique()
            complete_keys = set(center_counts[center_counts.eq(len(centers))].index)
            regional = regional[regional["unit_key"].isin(complete_keys)].copy()
        if regional.empty:
            raise SystemExit(f"EMPTY_REGION:{args.population_mode}:{region_name}")

        unit_meta = (
            regional[
                ["unit_key", "gene_id", "gene_names", "position_class_consensus"]
            ]
            .drop_duplicates()
            .sort_values("unit_key")
            .reset_index(drop=True)
        )
        if bool(unit_meta["unit_key"].duplicated().any()):
            raise SystemExit(f"UNIT_METADATA_CONFLICT:{region_name}")
        unit_meta.insert(0, "region", region_name)
        inventory_frames.append(unit_meta)

        genes = sorted(unit_meta["gene_id"].unique())
        gene_index = {gene: idx for idx, gene in enumerate(genes)}
        n_genes = len(genes)

        gene_sums: dict[tuple[int, str], np.ndarray] = {}
        gene_counts: dict[tuple[int, str], np.ndarray] = {}
        observed_meta: dict[str, pd.DataFrame] = {}

        for cls in CLASSES:
            cls_meta = (
                unit_meta[unit_meta["position_class_consensus"].eq(cls)]
                .sort_values("unit_key")
                .reset_index(drop=True)
            )
            if cls_meta.empty:
                raise SystemExit(f"EMPTY_CLASS:{region_name}:{cls}")
            observed_meta[cls] = cls_meta
            expected_keys = cls_meta["unit_key"].tolist()
            for center in centers:
                frame = regional[
                    regional["center"].eq(center)
                    & regional["position_class_consensus"].eq(cls)
                ].sort_values("unit_key")
                if (
                    args.population_mode == "COMPLETE_CASE"
                    and frame["unit_key"].tolist() != expected_keys
                ):
                    raise SystemExit(
                        f"COMPLETE_CASE_ORDER_GATE_FAIL:{region_name}:{center}:{cls}"
                    )
                matrix = build_matrix(frame)
                sums = np.zeros((n_genes, matrix.shape[1]), dtype=float)
                counts = np.zeros(n_genes, dtype=float)
                for row_index, gene in enumerate(frame["gene_id"].tolist()):
                    idx = gene_index[gene]
                    sums[idx] += matrix[row_index]
                    counts[idx] += 1.0
                gene_sums[(center, cls)] = sums
                gene_counts[(center, cls)] = counts

        def evaluate(gene_weights: np.ndarray) -> dict[str, object]:
            effects: dict[tuple[int, str, int, str], float] = {}
            for center in centers:
                for cls in CLASSES:
                    sums = gene_sums[(center, cls)]
                    counts = gene_counts[(center, cls)]
                    denominator = gene_weights @ counts
                    if bool((denominator <= 0).any()):
                        raise SystemExit(
                            f"EMPTY_WEIGHTED_CLASS:{region_name}:{center}:{cls}"
                        )
                    means = (gene_weights @ sums) / denominator[:, None]
                    for location, replicate, pd_lib, input_lib in CONTRASTS:
                        effects[(center, location, replicate, cls)] = d_vector(
                            means, pd_lib, input_lib
                        )

            fi_region: dict[tuple[str, int], np.ndarray] = {}
            for location in ["ER", "mitochondria"]:
                for replicate in [1, 2]:
                    per_center = []
                    for center in centers:
                        first = effects[(center, location, replicate, "FIRST")]
                        internal = effects[(center, location, replicate, "INTERNAL")]
                        per_center.append(first - internal)
                    fi_region[(location, replicate)] = np.mean(
                        np.vstack(per_center), axis=0
                    )
            interaction_rep1 = (
                fi_region[("mitochondria", 1)] - fi_region[("ER", 1)]
            )
            interaction_rep2 = (
                fi_region[("mitochondria", 2)] - fi_region[("ER", 2)]
            )
            interaction_mean = (interaction_rep1 + interaction_rep2) / 2.0
            return {
                "effects": effects,
                "fi_region": fi_region,
                "interaction_rep1": interaction_rep1,
                "interaction_rep2": interaction_rep2,
                "interaction_mean": interaction_mean,
            }

        observed_weights = np.ones((1, n_genes), dtype=float)
        observed = evaluate(observed_weights)

        for center in centers:
            for location in ["ER", "mitochondria"]:
                for replicate in [1, 2]:
                    first = float(
                        observed["effects"][(center, location, replicate, "FIRST")][0]
                    )
                    internal = float(
                        observed["effects"][(center, location, replicate, "INTERNAL")][0]
                    )
                    center_rows.append(
                        {
                            "region": region_name,
                            "center": center,
                            "location": location,
                            "replicate": replicate,
                            "FIRST_D": first,
                            "INTERNAL_D": internal,
                            "FIRST_minus_INTERNAL": first - internal,
                        }
                    )

        rng = np.random.default_rng(SEEDS[region_name])
        boot_interaction: list[np.ndarray] = []
        boot_er: list[np.ndarray] = []
        boot_mito: list[np.ndarray] = []
        complete = 0
        probabilities = np.full(n_genes, 1.0 / n_genes)
        while complete < BOOT_N:
            batch = min(BOOT_BATCH, BOOT_N - complete)
            weights = rng.multinomial(n_genes, probabilities, size=batch)
            values = evaluate(weights)
            boot_interaction.append(values["interaction_mean"])
            boot_er.append(
                (
                    values["fi_region"][("ER", 1)]
                    + values["fi_region"][("ER", 2)]
                )
                / 2.0
            )
            boot_mito.append(
                (
                    values["fi_region"][("mitochondria", 1)]
                    + values["fi_region"][("mitochondria", 2)]
                )
                / 2.0
            )
            complete += batch

        boot_interaction_array = np.concatenate(boot_interaction)
        boot_er_array = np.concatenate(boot_er)
        boot_mito_array = np.concatenate(boot_mito)
        interaction_ci = percentile_interval(boot_interaction_array)
        er_ci = percentile_interval(boot_er_array)
        mito_ci = percentile_interval(boot_mito_array)

        observed_interaction_rep1 = float(observed["interaction_rep1"][0])
        observed_interaction_rep2 = float(observed["interaction_rep2"][0])
        observed_interaction = float(observed["interaction_mean"][0])
        observed_er = float(
            (
                observed["fi_region"][("ER", 1)][0]
                + observed["fi_region"][("ER", 2)][0]
            )
            / 2.0
        )
        observed_mito = float(
            (
                observed["fi_region"][("mitochondria", 1)][0]
                + observed["fi_region"][("mitochondria", 2)][0]
            )
            / 2.0
        )

        first_units = int(
            unit_meta["position_class_consensus"].eq("FIRST").sum()
        )
        internal_units = int(
            unit_meta["position_class_consensus"].eq("INTERNAL").sum()
        )
        first_per_center = [
            int(gene_counts[(center, "FIRST")].sum()) for center in centers
        ]
        internal_per_center = [
            int(gene_counts[(center, "INTERNAL")].sum()) for center in centers
        ]
        summary_rows.append(
            {
                "region": region_name,
                "population_mode": args.population_mode,
                "centers": ",".join(map(str, centers)),
                "N_genes": n_genes,
                "N_FIRST_units_union": first_units,
                "N_INTERNAL_units_union": internal_units,
                "N_total_units_union": len(unit_meta),
                "N_FIRST_units_per_center_min": min(first_per_center),
                "N_FIRST_units_per_center_max": max(first_per_center),
                "N_INTERNAL_units_per_center_min": min(internal_per_center),
                "N_INTERNAL_units_per_center_max": max(internal_per_center),
                "ER_FIRST_minus_INTERNAL_mean": observed_er,
                "ER_CI_2.5": er_ci[0],
                "ER_CI_97.5": er_ci[1],
                "MITO_FIRST_minus_INTERNAL_mean": observed_mito,
                "MITO_CI_2.5": mito_ci[0],
                "MITO_CI_97.5": mito_ci[1],
                "MITO_minus_ER_rep1": observed_interaction_rep1,
                "MITO_minus_ER_rep2": observed_interaction_rep2,
                "MITO_minus_ER_mean": observed_interaction,
                "interaction_CI_2.5": interaction_ci[0],
                "interaction_CI_97.5": interaction_ci[1],
                "replicate_direction": (
                    "SAME_SIGN"
                    if np.sign(observed_interaction_rep1)
                    == np.sign(observed_interaction_rep2)
                    else "OPPOSITE_SIGN"
                ),
                "gene_clustered_CI_excludes_zero": (
                    "YES"
                    if interaction_ci[0] > 0 or interaction_ci[1] < 0
                    else "NO"
                ),
                "bootstrap_resamples": BOOT_N,
                "bootstrap_seed": SEEDS[region_name],
            }
        )

        leave_one_out_weights = np.ones((n_genes, n_genes), dtype=float)
        np.fill_diagonal(leave_one_out_weights, 0.0)
        loo = evaluate(leave_one_out_weights)["interaction_mean"]
        deletion_change = loo - observed_interaction
        relative_change = np.abs(deletion_change) / abs(observed_interaction)

        unit_counts = (
            unit_meta.groupby(["gene_id", "position_class_consensus"])
            .size()
            .unstack(fill_value=0)
            .reindex(genes, fill_value=0)
        )
        for cls in CLASSES:
            if cls not in unit_counts.columns:
                unit_counts[cls] = 0
        gene_names = (
            unit_meta.drop_duplicates("gene_id").set_index("gene_id")["gene_names"]
        )
        influence = pd.DataFrame(
            {
                "region": region_name,
                "gene_id": genes,
                "gene_name": [gene_names.get(gene, "") for gene in genes],
                "FIRST_units": unit_counts["FIRST"].astype(int).values,
                "INTERNAL_units": unit_counts["INTERNAL"].astype(int).values,
                "total_units": unit_counts[CLASSES].sum(axis=1).astype(int).values,
                "full_interaction": observed_interaction,
                "leave_one_gene_out_interaction": loo,
                "signed_change": deletion_change,
                "absolute_change": np.abs(deletion_change),
                "relative_change": relative_change,
                "sign_reversal": np.sign(loo) != np.sign(observed_interaction),
            }
        ).sort_values("absolute_change", ascending=False)
        influence_frames.append(influence)

        ranks_x = influence["total_units"].rank(method="average")
        ranks_y = influence["absolute_change"].rank(method="average")
        spearman = float(np.corrcoef(ranks_x, ranks_y)[0, 1])
        top_decile_n = max(1, int(np.ceil(0.10 * n_genes)))
        high_count = influence.sort_values("total_units", ascending=False).head(top_decile_n)
        total_absolute = float(influence["absolute_change"].sum())
        top_decile_share = (
            float(high_count["absolute_change"].sum()) / total_absolute
            if total_absolute > 0
            else 0.0
        )
        distribution_rows.append(
            {
                "region": region_name,
                "N_genes": n_genes,
                "max_absolute_change": float(influence["absolute_change"].max()),
                "max_relative_change": float(influence["relative_change"].max()),
                "any_sign_reversal": bool(influence["sign_reversal"].any()),
                "signed_change_q025": float(influence["signed_change"].quantile(0.025)),
                "signed_change_q50": float(influence["signed_change"].quantile(0.5)),
                "signed_change_q975": float(influence["signed_change"].quantile(0.975)),
                "absolute_change_q50": float(influence["absolute_change"].quantile(0.5)),
                "absolute_change_q90": float(influence["absolute_change"].quantile(0.9)),
                "absolute_change_q99": float(influence["absolute_change"].quantile(0.99)),
                "spearman_total_units_vs_absolute_change": spearman,
                "top_TMD_count_decile_share_of_total_absolute_change": top_decile_share,
            }
        )

    summary = pd.DataFrame(summary_rows)
    centers = pd.DataFrame(center_rows)
    influence_all = pd.concat(influence_frames, ignore_index=True)
    inventory = pd.concat(inventory_frames, ignore_index=True)
    distributions = pd.DataFrame(distribution_rows)

    output_files = {
        "summary": outdir / "REGION_COMPOSITE_SUMMARY_v1.tsv",
        "center_estimates": outdir / "REGION_CENTER_ESTIMATES_v1.tsv",
        "influence": outdir / "REGION_COMPOSITE_GENE_INFLUENCE_v1.tsv",
        "inventory": outdir / "REGION_COMPOSITE_UNIT_INVENTORY_v1.tsv",
        "influence_summary": outdir / "REGION_COMPOSITE_INFLUENCE_SUMMARY_v1.tsv",
    }
    summary.to_csv(output_files["summary"], sep="\t", index=False)
    centers.to_csv(output_files["center_estimates"], sep="\t", index=False)
    influence_all.to_csv(output_files["influence"], sep="\t", index=False)
    inventory.to_csv(output_files["inventory"], sep="\t", index=False)
    distributions.to_csv(output_files["influence_summary"], sep="\t", index=False)

    run_record = {
        "status": "PASS",
        "analysis_role": "outcome-informed post-reveal region-level robustness",
        "population_mode": args.population_mode,
        "frozen_project": str(root),
        "input_sha256": {str(path): sha256(path) for path in required_paths},
        "outputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in output_files.items()},
        "bootstrap_resamples": BOOT_N,
        "seeds": SEEDS,
        "regions": CENTERS_BY_REGION,
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "analysis_files_sha256": {
            "script": sha256(Path(__file__).resolve()),
            "contract": sha256(
                Path(__file__).resolve().parent
                / "REGION_COMPOSITE_ANALYSIS_CONTRACT_v1.md"
            ),
            "amendment": sha256(
                Path(__file__).resolve().parent
                / "REGION_COMPOSITE_ANALYSIS_AMENDMENT_v1p1.md"
            ),
        },
    }
    run_record_path = outdir / "REGION_COMPOSITE_RUN_RECORD_v1.json"
    run_record_path.write_text(json.dumps(run_record, indent=2) + "\n")

    print(summary.to_string(index=False))
    print()
    print(distributions.to_string(index=False))
    print()
    for path in [*output_files.values(), run_record_path]:
        print("SHA256", sha256(path), path)


if __name__ == "__main__":
    main()
