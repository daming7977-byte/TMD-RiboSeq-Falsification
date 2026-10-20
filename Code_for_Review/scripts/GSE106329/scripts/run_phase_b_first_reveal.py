#!/usr/bin/env python3
"""Execute the checksum-frozen Paper 3 first biological reveal.

The script is deliberately limited to spatial discovery, prespecified
sensitivities, and artifact diagnostics.  It does not calculate targeting fate,
sequence/topology associations, or candidate rankings.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import pysam
import yaml

from tmdriboseq.coordinates import TranscriptModel, load_gtf_transcripts
from tmdriboseq.phase_b import (
    ClusterScanResult,
    gene_cluster_maxT_scan,
    leave_one_gene_out_locked_band,
)
from tmdriboseq.psite import aligned_five_prime_reference_positions
from tmdriboseq.tmds import TMDRecord, load_paper_legacy_tmds


PREREG_SHA256 = "c722911d86a17702b843906f827b5e87a0f1cdaa4e2bb26b8e7e5e175d1fbe76"
EXECUTION_LOCK_SHA256 = "c604d66e6eb527b44c6c3cbe7a821268e70a53bf39cd42697bbd849b699c94b9"
CENTERS = np.arange(-90, 301, 3, dtype=int)
PRIMARY_GSMS = ("GSM2836139", "GSM2836143")


@dataclass(frozen=True)
class Anchor:
    anchor_id: str
    source_tmd_id: str
    gene_id: str
    transcript_id: str
    start_cds_nt: int
    end_cds_nt: int
    tmd_index: int
    tmd_count_gene: int
    is_real_tmd: bool


@dataclass
class PairProfile:
    effects: np.ndarray
    eligible: np.ndarray
    local_0: np.ndarray
    local_1: np.ndarray
    background_0: np.ndarray
    background_1: np.ndarray
    background_length: np.ndarray
    neighboring_tmd: np.ndarray


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_freeze(root: Path) -> None:
    prereg = root / "frozen/phase_b/phase_b_preregistration_amended_v1.yaml"
    execution = root / "config/phase_b_first_reveal_execution_lock.yaml"
    if sha256(prereg) != PREREG_SHA256:
        raise RuntimeError("amended Phase-B preregistration checksum mismatch")
    if sha256(execution) != EXECUTION_LOCK_SHA256:
        raise RuntimeError("first-reveal execution-lock checksum mismatch")
    data = yaml.safe_load(prereg.read_text())
    if data["positional_scan"]["centers_nt"] != {"start": -90, "end": 300, "step": 3}:
        raise RuntimeError("frozen scan domain differs from executable constants")
    if data["positional_scan"]["primary_window_width_nt"] != 21:
        raise RuntimeError("frozen primary window differs from executable constants")


def load_offsets(root: Path) -> dict[tuple[str, str], dict[int, int]]:
    consensus = pd.read_csv(
        root / "frozen/psite/psite_series_consensus_frozen.tsv", sep="\t"
    )
    sensitivity = pd.read_csv(
        root / "frozen/psite/psite_sensitivity_schemes_frozen.tsv", sep="\t"
    )
    manifest = pd.read_csv(root / "metadata/gsm_manifest.tsv", sep="\t")
    series_by_gsm = manifest.set_index("gsm")["series_id"].to_dict()
    output: dict[tuple[str, str], dict[int, int]] = {}
    for gsm, series in series_by_gsm.items():
        rows = consensus.loc[consensus.series_id == series]
        output[(gsm, "primary")] = dict(
            zip(rows.read_length.astype(int), rows.selected_offset_5prime_0based.astype(int))
        )
    for (gsm, scheme), rows in sensitivity.groupby(["sample_id", "scheme"]):
        key = {
            "sample_specific": "sample_specific",
            "author_L_minus_15": "author_L_minus_15",
            "adjacent_L_minus_16": "adjacent_L_minus_16",
        }.get(str(scheme))
        if key:
            output[(str(gsm), key)] = dict(
                zip(rows.read_length.astype(int), rows.offset_5prime_0based.astype(int))
            )
    return output


def build_layout(
    transcripts: Mapping[str, TranscriptModel], tmds: Sequence[TMDRecord]
) -> tuple[list[str], dict[str, slice], int]:
    transcript_ids = sorted({tmd.transcript_id for tmd in tmds})
    slices: dict[str, slice] = {}
    cursor = 0
    for transcript_id in transcript_ids:
        length = transcripts[transcript_id].cds_length_nt
        slices[transcript_id] = slice(cursor, cursor + length)
        cursor += length
    return transcript_ids, slices, cursor


def project_counter(
    counter: Mapping[tuple[str, str, int], int],
    transcripts: Mapping[str, TranscriptModel],
    transcript_ids: Sequence[str],
    slices: Mapping[str, slice],
    total_length: int,
) -> np.ndarray:
    output = np.zeros(total_length, dtype=np.int32)
    for transcript_id in transcript_ids:
        model = transcripts[transcript_id]
        output[slices[transcript_id]] = np.fromiter(
            (
                counter.get((model.contig, model.strand, position), 0)
                for position in model.genomic_positions0
            ),
            dtype=np.int32,
            count=model.cds_length_nt,
        )
    return output


def count_sample_tracks(
    bam_path: Path,
    track_offsets: Mapping[str, Mapping[int, int]],
    nm0_tracks: set[str],
) -> tuple[dict[str, Counter[tuple[str, str, int]]], dict[str, int]]:
    counters = {name: Counter() for name in track_offsets}
    metrics: Counter[str] = Counter()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            metrics["alignments_seen"] += 1
            if (
                read.is_unmapped
                or read.mapping_quality < 20
                or read.is_secondary
                or read.is_supplementary
            ):
                continue
            read_length = read.query_length
            if read_length is None or int(read_length) not in range(26, 33):
                continue
            metrics["eligible_26_32"] += 1
            positions = aligned_five_prime_reference_positions(read)
            strand = "-" if read.is_reverse else "+"
            for name, offsets in track_offsets.items():
                if name in nm0_tracks and (not read.has_tag("NM") or read.get_tag("NM") != 0):
                    continue
                offset = offsets.get(int(read_length))
                if offset is None or offset >= len(positions):
                    continue
                counters[name][(read.reference_name, strand, positions[offset])] += 1
                metrics[f"assigned__{name}"] += 1
    return counters, dict(metrics)


def ensure_sample_arrays(
    root: Path,
    gsm: str,
    bam_kind: str,
    tracks: Sequence[str],
    offsets: Mapping[tuple[str, str], Mapping[int, int]],
    transcripts: Mapping[str, TranscriptModel],
    transcript_ids: Sequence[str],
    slices: Mapping[str, slice],
    total_length: int,
    cache_dir: Path,
) -> dict[str, np.ndarray]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    output: dict[str, np.ndarray] = {}
    missing: list[str] = []
    for track in tracks:
        cache = cache_dir / f"{gsm}.{bam_kind}.{track}.npz"
        if cache.is_file():
            output[track] = np.load(cache)["counts"]
        else:
            missing.append(track)
    if not missing:
        return output
    if bam_kind == "full":
        bam = root / f"derived/alignments/{gsm}.sorted.bam"
    elif bam_kind == "umi":
        bam = root / f"derived/umi_deduplicated/{gsm}.umi_dedup.sorted.bam"
    else:
        raise ValueError(bam_kind)
    if not bam.is_file():
        raise FileNotFoundError(bam)
    track_offsets = {
        name: offsets[(gsm, "primary" if name in {"primary", "primary_nm0", "umi_primary"} else name)]
        for name in missing
    }
    start = time.time()
    counters, metrics = count_sample_tracks(
        bam, track_offsets, {name for name in missing if name.endswith("nm0")}
    )
    for name, counter in counters.items():
        array = project_counter(
            counter, transcripts, transcript_ids, slices, total_length
        )
        cache = cache_dir / f"{gsm}.{bam_kind}.{name}.npz"
        np.savez_compressed(cache, counts=array)
        output[name] = array
    metrics.update(
        {
            "gsm": gsm,
            "bam_kind": bam_kind,
            "bam_path": str(bam),
            "bam_size": bam.stat().st_size,
            "elapsed_seconds": time.time() - start,
            "tracks": missing,
        }
    )
    (cache_dir / f"{gsm}.{bam_kind}.metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True)
    )
    return output


def real_anchors(tmds: Sequence[TMDRecord]) -> list[Anchor]:
    counts = Counter(tmd.gene_id for tmd in tmds)
    output: list[Anchor] = []
    for tmd in tmds:
        try:
            index = int(tmd.tmd_id.rsplit("TMD", 1)[1])
        except (IndexError, ValueError):
            index = 0
        output.append(
            Anchor(
                anchor_id=tmd.tmd_id,
                source_tmd_id=tmd.tmd_id,
                gene_id=tmd.gene_id,
                transcript_id=tmd.transcript_id,
                start_cds_nt=tmd.start_cds_nt,
                end_cds_nt=tmd.end_cds_nt,
                tmd_index=index,
                tmd_count_gene=counts[tmd.gene_id],
                is_real_tmd=True,
            )
        )
    return output


def background_masks(
    transcripts: Mapping[str, TranscriptModel],
    tmds: Sequence[TMDRecord],
    width: int,
) -> dict[str, np.ndarray]:
    half = width // 2
    by_transcript: dict[str, list[TMDRecord]] = defaultdict(list)
    for tmd in tmds:
        by_transcript[tmd.transcript_id].append(tmd)
    masks: dict[str, np.ndarray] = {}
    for transcript_id, records in by_transcript.items():
        length = transcripts[transcript_id].cds_length_nt
        mask = np.ones(length, dtype=bool)
        for tmd in records:
            left = max(1, tmd.start_cds_nt - 90 - half)
            right = min(length, max(tmd.end_cds_nt, tmd.start_cds_nt + 300 + half))
            mask[left - 1 : right] = False
        masks[transcript_id] = mask
    return masks


def pair_profile(
    counts0: np.ndarray,
    counts1: np.ndarray,
    anchors: Sequence[Anchor],
    transcripts: Mapping[str, TranscriptModel],
    tmds: Sequence[TMDRecord],
    slices: Mapping[str, slice],
    width: int,
) -> PairProfile:
    n, p = len(anchors), len(CENTERS)
    local0 = np.zeros((n, p), dtype=np.int32)
    local1 = np.zeros((n, p), dtype=np.int32)
    eligible = np.zeros((n, p), dtype=bool)
    neighboring = np.zeros((n, p), dtype=bool)
    effects = np.full((n, p), np.nan, dtype=float)
    bg0 = np.zeros(n, dtype=np.int64)
    bg1 = np.zeros(n, dtype=np.int64)
    bg_length = np.zeros(n, dtype=np.int32)
    masks = background_masks(transcripts, tmds, width)
    by_transcript: dict[str, list[TMDRecord]] = defaultdict(list)
    for tmd in tmds:
        by_transcript[tmd.transcript_id].append(tmd)
    background_cache: dict[str, tuple[int, int, int]] = {}
    for transcript_id, mask in masks.items():
        slc = slices[transcript_id]
        background_cache[transcript_id] = (
            int(mask.sum()),
            int(counts0[slc][mask].sum()),
            int(counts1[slc][mask].sum()),
        )
    half = width // 2
    for row, anchor in enumerate(anchors):
        model = transcripts[anchor.transcript_id]
        slc = slices[anchor.transcript_id]
        array0 = counts0[slc]
        array1 = counts1[slc]
        cumulative0 = np.concatenate(([0], np.cumsum(array0, dtype=np.int64)))
        cumulative1 = np.concatenate(([0], np.cumsum(array1, dtype=np.int64)))
        length_bg, count_bg0, count_bg1 = background_cache[anchor.transcript_id]
        bg_length[row] = length_bg
        bg0[row] = count_bg0
        bg1[row] = count_bg1
        starts = anchor.start_cds_nt + CENTERS - half
        ends = anchor.start_cds_nt + CENTERS + half
        valid = (starts >= 1) & (ends <= model.cds_length_nt)
        valid_indices = np.flatnonzero(valid)
        start0 = starts[valid_indices] - 1
        end_exclusive = ends[valid_indices]
        local0[row, valid_indices] = cumulative0[end_exclusive] - cumulative0[start0]
        local1[row, valid_indices] = cumulative1[end_exclusive] - cumulative1[start0]
        if anchor.is_real_tmd:
            for column in valid_indices:
                left, right = int(starts[column]), int(ends[column])
                neighboring[row, column] = any(
                    other.tmd_id != anchor.anchor_id
                    and max(left, other.start_cds_nt) <= min(right, other.end_cds_nt)
                    for other in by_transcript[anchor.transcript_id]
                )
        coverage = (
            valid
            & (length_bg >= 300)
            & (count_bg0 >= 20)
            & (count_bg1 >= 20)
            & ((local0[row] + local1[row]) >= 5)
        )
        eligible[row] = coverage
        if np.any(coverage):
            density0 = ((local0[row, coverage] + 0.5) / width) / (
                (count_bg0 + 0.5) / length_bg
            )
            density1 = ((local1[row, coverage] + 0.5) / width) / (
                (count_bg1 + 0.5) / length_bg
            )
            effects[row, coverage] = np.log2(density1) - np.log2(density0)
    return PairProfile(
        effects=effects,
        eligible=eligible,
        local_0=local0,
        local_1=local1,
        background_0=bg0,
        background_1=bg1,
        background_length=bg_length,
        neighboring_tmd=neighboring,
    )


def scan_rows(name: str, result: ClusterScanResult) -> list[dict[str, object]]:
    return [
        {
            "analysis": name,
            "center_nt": int(result.centers[i]),
            "effect": result.estimate[i],
            "cluster_se": result.standard_error[i],
            "pointwise_95_low": result.pointwise_low[i],
            "pointwise_95_high": result.pointwise_high[i],
            "simultaneous_95_low": result.simultaneous_low[i],
            "simultaneous_95_high": result.simultaneous_high[i],
            "maxT_FWER_p": result.maxT_adjusted_p[i],
            "n_tmds": int(result.n_tmds[i]),
            "n_genes": int(result.n_genes[i]),
            "maxT_critical": result.maxT_critical,
        }
        for i in range(len(result.centers))
    ]


def band_rows(name: str, result: ClusterScanResult) -> list[dict[str, object]]:
    return [{"analysis": name, **asdict(band)} for band in result.bands]


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def coverage_rows(
    anchors: Sequence[Anchor], profile: PairProfile, primary: ClusterScanResult
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    peak_index = (
        primary.bands[0].peak_index
        if primary.bands
        else int(np.nanargmax(np.abs(primary.estimate)))
    )
    peak_members = profile.eligible[:, peak_index]
    for i, center in enumerate(CENTERS):
        members = profile.eligible[:, i]
        genes = {anchor.gene_id for anchor, keep in zip(anchors, members) if keep}
        union = np.count_nonzero(members | peak_members)
        peak_jaccard = np.count_nonzero(members & peak_members) / union if union else np.nan
        if i:
            previous = profile.eligible[:, i - 1]
            adjacent_union = np.count_nonzero(members | previous)
            adjacent_jaccard = (
                np.count_nonzero(members & previous) / adjacent_union
                if adjacent_union
                else np.nan
            )
            shared = members & previous
            shared_step = (
                float(np.nanmean(profile.effects[shared, i] - profile.effects[shared, i - 1]))
                if np.any(shared)
                else np.nan
            )
        else:
            adjacent_jaccard = np.nan
            shared_step = np.nan
        rows.append(
            {
                "center_nt": int(center),
                "n_tmds": int(members.sum()),
                "n_genes": len(genes),
                "adjacent_membership_jaccard": adjacent_jaccard,
                "peak_membership_jaccard": peak_jaccard,
                "neighboring_tmd_fraction": (
                    float(profile.neighboring_tmd[members, i].mean()) if np.any(members) else np.nan
                ),
                "median_local_pooled": (
                    float(np.median((profile.local_0 + profile.local_1)[members, i]))
                    if np.any(members)
                    else np.nan
                ),
                "median_background_0": (
                    float(np.median(profile.background_0[members])) if np.any(members) else np.nan
                ),
                "median_background_1": (
                    float(np.median(profile.background_1[members])) if np.any(members) else np.nan
                ),
                "observed_adjacent_effect_step": (
                    float(primary.estimate[i] - primary.estimate[i - 1]) if i else np.nan
                ),
                "shared_members_adjacent_effect_step": shared_step,
            }
        )
    return rows


def build_pseudo_anchors(
    anchors: Sequence[Anchor], transcripts: Mapping[str, TranscriptModel]
) -> tuple[list[Anchor], dict[str, np.ndarray]]:
    by_transcript: dict[str, list[Anchor]] = defaultdict(list)
    for anchor in anchors:
        by_transcript[anchor.transcript_id].append(anchor)
    output: list[Anchor] = []
    candidate_pools: dict[str, np.ndarray] = {}
    for transcript_id, records in by_transcript.items():
        length = transcripts[transcript_id].cds_length_nt
        candidates = np.arange(1, length + 1, 3, dtype=int)
        candidates = candidates[(candidates - 105 >= 1) & (candidates + 315 <= length)]
        clean = np.ones(len(candidates), dtype=bool)
        for record in records:
            clean &= ~(
                (candidates - 105 <= record.end_cds_nt)
                & (candidates + 315 >= record.start_cds_nt)
            )
        candidates = candidates[clean]
        candidate_pools[transcript_id] = candidates
        used: set[int] = set()
        for record in sorted(records, key=lambda value: (value.start_cds_nt, value.anchor_id)):
            for side, pool in (
                ("up", candidates[candidates < record.start_cds_nt]),
                ("down", candidates[candidates > record.start_cds_nt]),
            ):
                available = [int(value) for value in pool if int(value) not in used]
                if not available:
                    continue
                coordinate = (
                    max(available) if side == "up" else min(available)
                )
                used.add(coordinate)
                output.append(
                    Anchor(
                        anchor_id=f"{record.anchor_id}__pseudo_{side}_{coordinate}",
                        source_tmd_id=record.anchor_id,
                        gene_id=record.gene_id,
                        transcript_id=transcript_id,
                        start_cds_nt=coordinate,
                        end_cds_nt=coordinate,
                        tmd_index=record.tmd_index,
                        tmd_count_gene=record.tmd_count_gene,
                        is_real_tmd=False,
                    )
                )
    return output, candidate_pools


def longest_sign_run(values: np.ndarray) -> int:
    best = current = 0
    previous = 0
    for value in values:
        sign = int(np.sign(value)) if np.isfinite(value) else 0
        if sign and sign == previous:
            current += 1
        elif sign:
            current = 1
        else:
            current = 0
        previous = sign
        best = max(best, current)
    return best


def coordinate_null(
    counts0: np.ndarray,
    counts1: np.ndarray,
    real: Sequence[Anchor],
    transcripts: Mapping[str, TranscriptModel],
    tmds: Sequence[TMDRecord],
    slices: Mapping[str, slice],
    candidate_pools: Mapping[str, np.ndarray],
    observed: ClusterScanResult,
    primary_eligible_tmd_ids: set[str],
    iterations: int = 10_000,
    seed: int = 106_330,
) -> tuple[pd.DataFrame, dict[str, object]]:
    by_transcript: dict[str, list[Anchor]] = defaultdict(list)
    for anchor in real:
        by_transcript[anchor.transcript_id].append(anchor)
    total_sum = np.zeros((iterations, len(CENTERS)), dtype=float)
    total_count = np.zeros((iterations, len(CENTERS)), dtype=np.int32)
    rng = np.random.default_rng(seed)
    eligible_genes: set[str] = set()
    retained_primary_tmds = 0
    pool_rows: list[dict[str, object]] = []
    for transcript_id, records in sorted(by_transcript.items()):
        candidates = candidate_pools[transcript_id]
        k = len(records)
        if len(candidates) < k:
            pool_rows.append(
                {"gene_id": records[0].gene_id, "tmd_count": k, "candidate_count": len(candidates), "eligible": False}
            )
            continue
        candidate_anchors = [
            Anchor(
                anchor_id=f"null_{transcript_id}_{coordinate}",
                source_tmd_id="",
                gene_id=records[0].gene_id,
                transcript_id=transcript_id,
                start_cds_nt=int(coordinate),
                end_cds_nt=int(coordinate),
                tmd_index=0,
                tmd_count_gene=k,
                is_real_tmd=False,
            )
            for coordinate in candidates
        ]
        profile = pair_profile(
            counts0, counts1, candidate_anchors, transcripts, tmds, slices, 21
        )
        effects = profile.effects
        eligible_genes.add(records[0].gene_id)
        retained_primary_tmds += sum(
            record.anchor_id in primary_eligible_tmd_ids for record in records
        )
        pool_rows.append(
            {"gene_id": records[0].gene_id, "tmd_count": k, "candidate_count": len(candidates), "eligible": True}
        )
        for start in range(0, iterations, 250):
            stop = min(iterations, start + 250)
            size = stop - start
            if k == 1:
                chosen = rng.integers(0, len(candidates), size=(size, 1))
            else:
                chosen = np.empty((size, k), dtype=int)
                for row in range(size):
                    chosen[row] = rng.choice(len(candidates), size=k, replace=False)
            selected = effects[chosen]
            finite = np.isfinite(selected)
            total_sum[start:stop] += np.where(finite, selected, 0.0).sum(axis=1)
            total_count[start:stop] += finite.sum(axis=1)
    curves = np.divide(
        total_sum,
        total_count,
        out=np.full_like(total_sum, np.nan),
        where=total_count > 0,
    )
    max_abs = np.nanmax(np.abs(curves), axis=1)
    run = np.asarray([longest_sign_run(curve) for curve in curves], dtype=int)
    if observed.bands:
        locked = observed.bands[0]
        band_indices = np.arange(locked.start_index, locked.end_index + 1)
    else:
        peak = int(np.nanargmax(np.abs(observed.estimate)))
        band_indices = np.asarray([peak])
    locked_mean = np.nanmean(curves[:, band_indices], axis=1)
    observed_max = float(np.nanmax(np.abs(observed.estimate)))
    observed_run = longest_sign_run(observed.estimate)
    observed_locked = float(np.nanmean(observed.estimate[band_indices]))
    summary = {
        "eligible_gene_count": len(eligible_genes),
        "retained_primary_tmd_count": retained_primary_tmds,
        "primary_eligible_tmd_count": len(primary_eligible_tmd_ids),
        "retained_primary_tmd_fraction": retained_primary_tmds / len(primary_eligible_tmd_ids),
        "observed_max_absolute_effect": observed_max,
        "null_max_absolute_effect_p": (1 + int(np.count_nonzero(max_abs >= observed_max))) / (iterations + 1),
        "observed_longest_sign_run": observed_run,
        "null_longest_sign_run_p": (1 + int(np.count_nonzero(run >= observed_run))) / (iterations + 1),
        "observed_absolute_locked_band_mean": abs(observed_locked),
        "null_absolute_locked_band_mean_p": (1 + int(np.count_nonzero(np.abs(locked_mean) >= abs(observed_locked)))) / (iterations + 1),
        "iterations": iterations,
        "seed": seed,
    }
    return pd.DataFrame(pool_rows), summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    verify_freeze(root)
    output = root / "phase_b_first_reveal"
    output.mkdir(parents=True, exist_ok=True)
    cache = output / "cache"
    repo = root / "software/TMD-RiboSeq"
    gtf = repo / "revision/input/Saccharomyces_cerevisiae.R64-1-1.115.gtf"
    tmd_path = repo / "revision/input/yeast_tmd_genomic_table.tsv"
    transcripts = load_gtf_transcripts(gtf)
    tmds = load_paper_legacy_tmds(tmd_path, gtf)
    anchors = real_anchors(tmds)
    transcript_ids, slices, total_length = build_layout(transcripts, tmds)
    offsets = load_offsets(root)

    full_samples = {
        "GSM2836139": ["primary", "sample_specific", "author_L_minus_15", "adjacent_L_minus_16", "primary_nm0"],
        "GSM2836143": ["primary", "sample_specific", "author_L_minus_15", "adjacent_L_minus_16", "primary_nm0"],
        "GSM2836141": ["primary"],
        "GSM2836145": ["primary"],
        "GSM2836155": ["primary"],
        "GSM2836157": ["primary"],
        "GSM2836159": ["primary"],
        "GSM2836147": ["primary"],
        "GSM2836149": ["primary"],
        "GSM2836151": ["primary"],
        "GSM2836153": ["primary"],
    }
    umi_samples = {
        gsm: ["umi_primary"]
        for gsm in [
            "GSM2836139", "GSM2836141", "GSM2836143", "GSM2836145",
            "GSM2836155", "GSM2836157", "GSM2836159",
        ]
    }
    arrays: dict[tuple[str, str, str], np.ndarray] = {}
    for bam_kind, sample_map in (("full", full_samples), ("umi", umi_samples)):
        for gsm, tracks in sample_map.items():
            print(f"COUNT {bam_kind} {gsm} {','.join(tracks)}", flush=True)
            found = ensure_sample_arrays(
                root, gsm, bam_kind, tracks, offsets, transcripts, transcript_ids,
                slices, total_length, cache
            )
            for track, array in found.items():
                arrays[(gsm, bam_kind, track)] = array

    profile_specs = [
        ("primary", "GSM2836139", "GSM2836143", "full", "primary", 21),
        ("umi_primary", "GSM2836139", "GSM2836143", "umi", "umi_primary", 21),
        ("sample_specific_psite", "GSM2836139", "GSM2836143", "full", "sample_specific", 21),
        ("author_L_minus_15", "GSM2836139", "GSM2836143", "full", "author_L_minus_15", 21),
        ("adjacent_L_minus_16", "GSM2836139", "GSM2836143", "full", "adjacent_L_minus_16", 21),
        ("author_zero_mismatch_alignment", "GSM2836139", "GSM2836143", "full", "primary_nm0", 21),
        ("window_15nt", "GSM2836139", "GSM2836143", "full", "primary", 15),
        ("window_31nt", "GSM2836139", "GSM2836143", "full", "primary", 31),
        ("Sec63_noCHX_30_minus_0", "GSM2836155", "GSM2836157", "full", "primary", 21),
        ("Sec63_noCHX_UMI_30_minus_0", "GSM2836155", "GSM2836157", "umi", "umi_primary", 21),
        ("Ubc6_descriptive_30_minus_0", "GSM2836147", "GSM2836151", "full", "primary", 21),
    ]
    profiles: dict[str, PairProfile] = {}
    scans: dict[str, ClusterScanResult] = {}
    all_scan_rows: list[dict[str, object]] = []
    all_band_rows: list[dict[str, object]] = []
    genes = [anchor.gene_id for anchor in anchors]
    for name, gsm0, gsm1, kind, track, width in profile_specs:
        print(f"PROFILE {name}", flush=True)
        profile = pair_profile(
            arrays[(gsm0, kind, track)], arrays[(gsm1, kind, track)], anchors,
            transcripts, tmds, slices, width
        )
        profiles[name] = profile
        print(f"INFERENCE {name}", flush=True)
        result = gene_cluster_maxT_scan(
            profile.effects, genes, CENTERS, iterations=10_000, seed=106_329
        )
        scans[name] = result
        all_scan_rows.extend(scan_rows(name, result))
        all_band_rows.extend(band_rows(name, result))
    neighbor_effects = profiles["primary"].effects.copy()
    neighbor_effects[profiles["primary"].neighboring_tmd] = np.nan
    neighbor_scan = gene_cluster_maxT_scan(
        neighbor_effects, genes, CENTERS, iterations=10_000, seed=106_329
    )
    scans["neighboring_TMD_excluded"] = neighbor_scan
    all_scan_rows.extend(scan_rows("neighboring_TMD_excluded", neighbor_scan))
    all_band_rows.extend(band_rows("neighboring_TMD_excluded", neighbor_scan))
    write_tsv(output / "spatial_scans.tsv", all_scan_rows)
    write_tsv(output / "spatial_bands.tsv", all_band_rows)

    primary = scans["primary"]
    primary_profile = profiles["primary"]
    membership_rows = []
    for row, anchor in enumerate(anchors):
        for column, center in enumerate(CENTERS):
            membership_rows.append(
                {
                    "tmd_id": anchor.anchor_id,
                    "gene_id": anchor.gene_id,
                    "tmd_index": anchor.tmd_index,
                    "tmd_count_gene": anchor.tmd_count_gene,
                    "center_nt": int(center),
                    "eligible": bool(primary_profile.eligible[row, column]),
                    "neighboring_tmd_overlap": bool(primary_profile.neighboring_tmd[row, column]),
                    "local_0": int(primary_profile.local_0[row, column]),
                    "local_30": int(primary_profile.local_1[row, column]),
                    "background_0": int(primary_profile.background_0[row]),
                    "background_30": int(primary_profile.background_1[row]),
                    "background_length_nt": int(primary_profile.background_length[row]),
                    "effect": primary_profile.effects[row, column],
                }
            )
    pd.DataFrame(membership_rows).to_csv(
        output / "primary_coverage_membership.tsv.gz", sep="\t", index=False
    )
    write_tsv(output / "primary_coverage_by_position.tsv", coverage_rows(anchors, primary_profile, primary))
    fixed_members = primary_profile.eligible.all(axis=1)
    fixed_rows = [
        {
            "center_nt": int(center),
            "effect": float(np.nanmean(primary_profile.effects[fixed_members, column])),
            "n_tmds": int(fixed_members.sum()),
            "n_genes": len(
                {
                    anchor.gene_id
                    for anchor, keep in zip(anchors, fixed_members)
                    if keep
                }
            ),
        }
        for column, center in enumerate(CENTERS)
    ]
    write_tsv(output / "fixed_all_positions_population_curve.tsv", fixed_rows)

    pseudo, candidate_pools = build_pseudo_anchors(anchors, transcripts)
    print(f"PSEUDO anchors={len(pseudo)}", flush=True)
    pseudo_profile = pair_profile(
        arrays[("GSM2836139", "full", "primary")],
        arrays[("GSM2836143", "full", "primary")],
        pseudo, transcripts, tmds, slices, 21
    )
    pseudo_scan = gene_cluster_maxT_scan(
        pseudo_profile.effects,
        [anchor.gene_id for anchor in pseudo],
        CENTERS,
        iterations=10_000,
        seed=106_329,
    )
    write_tsv(output / "pseudo_anchor_scan.tsv", scan_rows("pseudo_anchor", pseudo_scan))
    write_tsv(output / "pseudo_anchor_bands.tsv", band_rows("pseudo_anchor", pseudo_scan))
    primary_eligible_tmd_ids = {
        anchor.anchor_id
        for row, anchor in enumerate(anchors)
        if primary_profile.eligible[row].any()
    }
    source_with_eligible_pseudo = {
        anchor.source_tmd_id
        for row, anchor in enumerate(pseudo)
        if pseudo_profile.eligible[row].any()
        and anchor.source_tmd_id in primary_eligible_tmd_ids
    }
    pseudo_summary = {
        "pseudo_anchor_count": len(pseudo),
        "eligible_gene_count_any_position": len({anchor.gene_id for i, anchor in enumerate(pseudo) if pseudo_profile.eligible[i].any()}),
        "source_primary_tmd_with_eligible_pseudo_count": len(source_with_eligible_pseudo),
        "primary_eligible_tmd_count": len(primary_eligible_tmd_ids),
        "source_primary_tmd_with_eligible_pseudo_fraction": len(source_with_eligible_pseudo) / len(primary_eligible_tmd_ids),
        "maximum_absolute_effect": float(np.nanmax(np.abs(pseudo_scan.estimate))),
        "band_count": len(pseudo_scan.bands),
    }
    (output / "pseudo_anchor_summary.json").write_text(json.dumps(pseudo_summary, indent=2))

    print("COORDINATE NULL", flush=True)
    null_pools, null_summary = coordinate_null(
        arrays[("GSM2836139", "full", "primary")],
        arrays[("GSM2836143", "full", "primary")],
        anchors, transcripts, tmds, slices, candidate_pools, primary,
        primary_eligible_tmd_ids
    )
    null_pools.to_csv(output / "coordinate_null_gene_pools.tsv", sep="\t", index=False)
    (output / "coordinate_null_summary.json").write_text(
        json.dumps(null_summary, indent=2, sort_keys=True)
    )

    if primary.bands:
        locked = primary.bands[0]
        band_indices = range(locked.start_index, locked.end_index + 1)
        loo = leave_one_gene_out_locked_band(
            primary_profile.effects, genes, band_indices, locked.peak_index
        )
        write_tsv(output / "leave_one_gene_out.tsv", loo)
    else:
        (output / "leave_one_gene_out.tsv").write_text("")

    trajectory_specs = [
        ("Sec63_CHX_15_minus_0", "GSM2836139", "GSM2836141", "full", "primary"),
        ("Sec63_CHX_60_minus_0", "GSM2836139", "GSM2836145", "full", "primary"),
        ("Sec63_CHX_UMI_15_minus_0", "GSM2836139", "GSM2836141", "umi", "umi_primary"),
        ("Sec63_CHX_UMI_60_minus_0", "GSM2836139", "GSM2836145", "umi", "umi_primary"),
        ("Sec63_noCHX_60_minus_0", "GSM2836155", "GSM2836159", "full", "primary"),
        ("Sec63_noCHX_UMI_60_minus_0", "GSM2836155", "GSM2836159", "umi", "umi_primary"),
        ("Ubc6_descriptive_15_minus_0", "GSM2836147", "GSM2836149", "full", "primary"),
        ("Ubc6_descriptive_60_minus_0", "GSM2836147", "GSM2836153", "full", "primary"),
    ]
    trajectory_rows: list[dict[str, object]] = []
    if primary.bands:
        locked = primary.bands[0]
        indices = np.arange(locked.start_index, locked.end_index + 1)
    else:
        indices = np.asarray([int(np.nanargmax(np.abs(primary.estimate)))])
    for name, gsm0, gsm1, kind, track in trajectory_specs:
        profile = pair_profile(
            arrays[(gsm0, kind, track)], arrays[(gsm1, kind, track)], anchors,
            transcripts, tmds, slices, 21
        )
        curve = np.nanmean(profile.effects, axis=0)
        trajectory_rows.append(
            {
                "analysis": name,
                "locked_band_mean_effect": float(np.nanmean(curve[indices])),
                "locked_peak_effect": float(curve[indices[np.nanargmax(np.abs(primary.estimate[indices]))]]),
                "eligible_tmds_at_primary_peak": int(np.isfinite(profile.effects[:, indices[0]]).sum()),
            }
        )
    write_tsv(output / "locked_band_time_trajectories.tsv", trajectory_rows)

    summary = {
        "preregistration_sha256": PREREG_SHA256,
        "execution_lock_sha256": EXECUTION_LOCK_SHA256,
        "tmd_count": len(anchors),
        "gene_count": len(set(genes)),
        "centers": CENTERS.tolist(),
        "primary_band_count": len(primary.bands),
        "primary_bands": [asdict(band) for band in primary.bands],
        "primary_max_absolute_effect": float(np.nanmax(np.abs(primary.estimate))),
        "primary_min_maxT_FWER_p": float(np.nanmin(primary.maxT_adjusted_p)),
        "pseudo_anchor": pseudo_summary,
        "coordinate_null": null_summary,
        "analysis_names": list(scans),
    }
    (output / "first_reveal_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True)
    )


if __name__ == "__main__":
    main()
