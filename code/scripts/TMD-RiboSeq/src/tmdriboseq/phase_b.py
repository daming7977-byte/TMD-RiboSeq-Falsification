"""Gene-clustered spatial inference for Paper 3 Phase B.

This module is additive to the v1.0 workflow.  It does not change the released
profile, differential, or paper-reproduction paths.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class SpatialBand:
    start_index: int
    end_index: int
    start_nt: int
    end_nt: int
    peak_index: int
    peak_nt: int
    peak_effect: float
    direction: int


@dataclass(frozen=True)
class ClusterScanResult:
    centers: np.ndarray
    estimate: np.ndarray
    standard_error: np.ndarray
    pointwise_low: np.ndarray
    pointwise_high: np.ndarray
    simultaneous_low: np.ndarray
    simultaneous_high: np.ndarray
    maxT_adjusted_p: np.ndarray
    n_tmds: np.ndarray
    n_genes: np.ndarray
    maxT_critical: float
    bands: tuple[SpatialBand, ...]


def _aggregate_by_gene(
    effects: np.ndarray, gene_ids: Sequence[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if effects.ndim != 2:
        raise ValueError("effects must be a two-dimensional TMD-by-position array")
    if effects.shape[0] != len(gene_ids):
        raise ValueError("gene_ids length must equal the number of effect rows")
    genes, inverse = np.unique(np.asarray(gene_ids, dtype=str), return_inverse=True)
    sums = np.zeros((len(genes), effects.shape[1]), dtype=float)
    counts = np.zeros((len(genes), effects.shape[1]), dtype=float)
    finite = np.isfinite(effects)
    np.add.at(sums, inverse, np.where(finite, effects, 0.0))
    np.add.at(counts, inverse, finite.astype(float))
    return genes, sums, counts


def _cluster_standard_error(
    gene_sums: np.ndarray,
    gene_counts: np.ndarray,
    estimate: np.ndarray,
) -> np.ndarray:
    n = gene_counts.sum(axis=0)
    active = (gene_counts > 0).sum(axis=0)
    residual_cluster_sums = gene_sums - gene_counts * estimate
    variance = np.full(estimate.shape, np.nan, dtype=float)
    valid = (n > 0) & (active > 1)
    variance[valid] = (
        active[valid]
        / (active[valid] - 1.0)
        * np.square(residual_cluster_sums[:, valid]).sum(axis=0)
        / np.square(n[valid])
    )
    return np.sqrt(variance)


def find_spatial_bands(
    centers: Sequence[int],
    estimates: np.ndarray,
    simultaneous_low: np.ndarray,
    simultaneous_high: np.ndarray,
    minimum_consecutive: int = 3,
) -> tuple[SpatialBand, ...]:
    centers_array = np.asarray(centers, dtype=int)
    direction = np.zeros(len(centers_array), dtype=int)
    direction[simultaneous_low > 0] = 1
    direction[simultaneous_high < 0] = -1
    bands: list[SpatialBand] = []
    index = 0
    while index < len(direction):
        sign = int(direction[index])
        if sign == 0:
            index += 1
            continue
        end = index + 1
        while end < len(direction) and int(direction[end]) == sign:
            end += 1
        if end - index >= minimum_consecutive:
            segment = np.abs(estimates[index:end])
            peak = index + int(np.nanargmax(segment))
            bands.append(
                SpatialBand(
                    start_index=index,
                    end_index=end - 1,
                    start_nt=int(centers_array[index]),
                    end_nt=int(centers_array[end - 1]),
                    peak_index=peak,
                    peak_nt=int(centers_array[peak]),
                    peak_effect=float(estimates[peak]),
                    direction=sign,
                )
            )
        index = end
    return tuple(bands)


def gene_cluster_maxT_scan(
    effects: np.ndarray,
    gene_ids: Sequence[str],
    centers: Sequence[int],
    *,
    iterations: int = 10_000,
    seed: int = 106_329,
    alpha: float = 0.05,
    minimum_consecutive: int = 3,
    batch_size: int = 200,
) -> ClusterScanResult:
    """Compute a studentized gene-cluster bootstrap maxT spatial scan."""

    effects = np.asarray(effects, dtype=float)
    centers_array = np.asarray(centers, dtype=int)
    if effects.shape[1] != len(centers_array):
        raise ValueError("centers length must equal the number of effect columns")
    if iterations <= 0:
        raise ValueError("iterations must be positive")
    genes, gene_sums, gene_counts = _aggregate_by_gene(effects, gene_ids)
    total_count = gene_counts.sum(axis=0)
    estimate = np.divide(
        gene_sums.sum(axis=0),
        total_count,
        out=np.full(effects.shape[1], np.nan),
        where=total_count > 0,
    )
    standard_error = _cluster_standard_error(gene_sums, gene_counts, estimate)
    observed_t = np.divide(
        estimate,
        standard_error,
        out=np.full_like(estimate, np.nan),
        where=standard_error > 0,
    )

    rng = np.random.default_rng(seed)
    maximum_t = np.full(iterations, np.nan, dtype=float)
    probability = np.full(len(genes), 1.0 / len(genes))
    completed = 0
    while completed < iterations:
        size = min(batch_size, iterations - completed)
        multiplicity = rng.multinomial(len(genes), probability, size=size).astype(float)
        bootstrap_sum = multiplicity @ gene_sums
        bootstrap_count = multiplicity @ gene_counts
        bootstrap_mean = np.divide(
            bootstrap_sum,
            bootstrap_count,
            out=np.full_like(bootstrap_sum, np.nan),
            where=bootstrap_count > 0,
        )
        term1 = multiplicity @ np.square(gene_sums)
        term2 = multiplicity @ (gene_sums * gene_counts)
        term3 = multiplicity @ np.square(gene_counts)
        residual_square_sum = term1 - 2.0 * bootstrap_mean * term2 + np.square(
            bootstrap_mean
        ) * term3
        active = multiplicity @ (gene_counts > 0).astype(float)
        bootstrap_variance = np.divide(
            active * residual_square_sum,
            (active - 1.0) * np.square(bootstrap_count),
            out=np.full_like(bootstrap_mean, np.nan),
            where=(active > 1) & (bootstrap_count > 0),
        )
        bootstrap_se = np.sqrt(np.maximum(bootstrap_variance, 0.0))
        bootstrap_t = np.divide(
            bootstrap_mean - estimate,
            bootstrap_se,
            out=np.full_like(bootstrap_mean, np.nan),
            where=bootstrap_se > 0,
        )
        maximum_t[completed : completed + size] = np.nanmax(
            np.abs(bootstrap_t), axis=1
        )
        completed += size

    finite_maximum = maximum_t[np.isfinite(maximum_t)]
    if not len(finite_maximum):
        raise ValueError("bootstrap produced no finite maximum-t statistics")
    try:
        critical = float(np.quantile(finite_maximum, 1.0 - alpha, method="higher"))
    except TypeError:  # NumPy < 1.22 compatibility.
        critical = float(
            np.quantile(finite_maximum, 1.0 - alpha, interpolation="higher")
        )
    adjusted_p = np.full_like(estimate, np.nan)
    for position, statistic in enumerate(np.abs(observed_t)):
        if np.isfinite(statistic):
            adjusted_p[position] = (
                1.0 + np.count_nonzero(finite_maximum >= statistic)
            ) / (1.0 + len(finite_maximum))

    simultaneous_low = estimate - critical * standard_error
    simultaneous_high = estimate + critical * standard_error
    pointwise_low = estimate - 1.96 * standard_error
    pointwise_high = estimate + 1.96 * standard_error
    bands = find_spatial_bands(
        centers_array,
        estimate,
        simultaneous_low,
        simultaneous_high,
        minimum_consecutive=minimum_consecutive,
    )
    return ClusterScanResult(
        centers=centers_array,
        estimate=estimate,
        standard_error=standard_error,
        pointwise_low=pointwise_low,
        pointwise_high=pointwise_high,
        simultaneous_low=simultaneous_low,
        simultaneous_high=simultaneous_high,
        maxT_adjusted_p=adjusted_p,
        n_tmds=np.isfinite(effects).sum(axis=0),
        n_genes=(gene_counts > 0).sum(axis=0),
        maxT_critical=critical,
        bands=bands,
    )


def leave_one_gene_out_locked_band(
    effects: np.ndarray,
    gene_ids: Sequence[str],
    band_indices: Iterable[int],
    peak_index: int,
) -> list[dict[str, float | int | str | bool]]:
    """Recompute locked-band summaries after removing one complete gene."""

    effects = np.asarray(effects, dtype=float)
    indices = np.asarray(list(band_indices), dtype=int)
    full_curve = np.nanmean(effects, axis=0)
    full_band = float(np.nanmean(full_curve[indices]))
    full_peak = float(full_curve[peak_index])
    output: list[dict[str, float | int | str | bool]] = []
    gene_array = np.asarray(gene_ids, dtype=str)
    for gene in np.unique(gene_array):
        keep = gene_array != gene
        curve = np.nanmean(effects[keep], axis=0)
        band = float(np.nanmean(curve[indices]))
        peak = float(curve[peak_index])
        attenuation = (
            (abs(full_band) - abs(band)) / abs(full_band)
            if full_band != 0
            else float("nan")
        )
        reversal = bool(full_band != 0 and band != 0 and np.sign(full_band) != np.sign(band))
        output.append(
            {
                "gene_id": str(gene),
                "n_tmds": int(np.count_nonzero(gene_array == gene)),
                "full_band_mean": full_band,
                "leave_one_out_band_mean": band,
                "absolute_change": abs(band - full_band),
                "proportional_attenuation": attenuation,
                "direction_reversal": reversal,
                "full_peak_effect": full_peak,
                "leave_one_out_peak_effect": peak,
            }
        )
    return output

