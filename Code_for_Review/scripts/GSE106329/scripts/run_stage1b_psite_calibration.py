#!/usr/bin/env python3
"""Run frozen v1 coding-frame calibration on all Stage 1B series."""

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "frozen/v1.0.0_owner_review_candidate/wheel_extracted"))
from tmdriboseq.qc import run_calibration_stage  # noqa: E402


config = {
    "reference": {"annotation": str(ROOT / "references/genome/Saccharomyces_cerevisiae.R64-1-1.115.gtf")},
    "samples": {
        "metadata": str(ROOT / "metadata/stage1b_calibration_samples.tsv"),
        "columns": {
            "sample_id": "sample_id", "bam": "bam", "replicate": "replicate",
            "role": "role", "series": "series", "fraction": "fraction",
        },
    },
    "design": {"factors": ["series", "fraction"]},
    "reads": {
        "allowed_lengths": [26, 27, 28, 29, 30, 31, 32],
        "minimum_mapping_quality": 20, "include_secondary": False,
        "include_supplementary": False, "include_duplicates": False,
    },
    "psite": {
        "mode": "calibrate", "candidate_offsets": list(range(10, 18)),
        "canonical_offset_by_phase": {0: 12, 1: 13, 2: 14},
        "minimum_reads_per_length": 10000, "warning_frame0_fraction": 0.50,
        "start_codon_metagene_start_nt": 0, "start_codon_metagene_end_nt": 120,
        "maximum_alignments_per_sample": None,
    },
    "output": {"directory": str(ROOT / "qc/stage1b_calibration")},
}
summary = run_calibration_stage(config)
(ROOT / "qc/stage1b_calibration/calibration_summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(summary, indent=2, sort_keys=True))
