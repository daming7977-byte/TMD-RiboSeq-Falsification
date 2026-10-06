# Stage 0C Module 04 Blinded Technical-QC Amendment

**Date:** 2026-08-26  
**Module:** `scripts/stage0c/04_apply_psite_offsets.py`  
**Amended module version:** `0.2.0`  
**Authoritative P-site freeze:** `stage0c_psite_v1_2026-08-24`  
**Blinding state:** TMD-blind, benchmark-blind, condition-blind, biological-result-blind

## Scope and rationale

The P-site offsets were calibrated independently in Pilot A and Pilot B and
were frozen before full-library application. Both pilots supported the frozen
length-28 and length-29 offsets `28 -> 9` and `29 -> 9`. Full-library phase
dominance is therefore a
post-application technical diagnostic, not a recalibration or sample/read-
length selection rule.

This amendment changes only the handling of the full-library status
`PHASE0_NOT_DOMINANT`: Module 04 records it as a technical warning and still
writes valid production outputs. It does not change an offset, remove a sample,
or remove a calibrated read length. `INSUFFICIENT_VALIDATION_EVIDENCE`, count
conservation failures, per-length accounting failures, phase-count accounting
failures, validation-region accounting failures, non-frozen offset use, and
all other integrity failures remain hard failures.

## Blinded technical observation

The full-library length-28 phase QC supplied for this amendment was:

| Library | Phase 0 | Phase 1 | Phase 2 | Status |
|---|---:|---:|---:|---|
| GSM3900208 | 46.04% | 11.47% | 42.48% | `PHASE0_DOMINANT` |
| GSM3900209 | 44.95% | 11.30% | 43.74% | `PHASE0_DOMINANT` |
| GSM3900210 | 43.41% | 11.57% | 45.02% | `PHASE0_NOT_DOMINANT` |
| GSM3900211 | 48.35% | 10.26% | 41.39% | `PHASE0_DOMINANT` |
| GSM3900212 | 45.54% | 12.55% | 41.91% | `PHASE0_DOMINANT` |
| GSM3900213 | 48.59% | 11.18% | 40.24% | `PHASE0_DOMINANT` |

The full-library length-29 phase QC subsequently supplied for this amendment
was:

| Library | Phase 0 | Phase 1 | Phase 2 | Status |
|---|---:|---:|---:|---|
| GSM3900208 | 43.06% | 14.93% | 42.01% | `PHASE0_DOMINANT` |
| GSM3900209 | 42.64% | 14.85% | 42.51% | `PHASE0_DOMINANT` |
| GSM3900210 | 40.16% | 15.59% | 44.25% | `PHASE0_NOT_DOMINANT` |
| GSM3900211 | 50.51% | 13.21% | 36.29% | `PHASE0_DOMINANT` |
| GSM3900212 | 46.32% | 14.61% | 39.08% | `PHASE0_DOMINANT` |
| GSM3900213 | 49.25% | 13.03% | 37.72% | `PHASE0_DOMINANT` |

The discrepancy is sample-specific to GSM3900210 and affects lengths 28 and
29. At both lengths it is a phase-0/phase-2 ambiguity while phase 1 remains
low. This technical observation does not justify recalibration: frozen offsets
`28 -> 9` and `29 -> 9` remain unchanged. Each discrepancy is retained exactly
in the generated `phase_check.tsv`; Module 04 does not rewrite the phase
counts, percentages, dominant phase, or `PHASE0_NOT_DOMINANT` status. No
biological, TMD, benchmark, condition, or group-level result was inspected in
making this amendment.

## Production policy

1. The frozen mapping remains exactly `27->9, 28->9, 29->9, 30->12,
   32->13, 33->13, 34->13, 35->13`; unresolved lengths remain unresolved.
2. Primary production output retains all samples and every read length allowed
   by the authoritative freeze. A full-library phase warning causes no new
   exclusion.
3. `phase_check.tsv` remains the exact technical record of phase counts,
   percentages, dominant phase, and status.
4. Run statistics explicitly state that full-library phase dominance is
   diagnostic after independent frozen calibration and record warning count
   and warning read length(s).

## Pre-registered sensitivity analysis

Before any biological result is inspected, the mandatory blinded sensitivity
analysis is pre-registered as follows: repeat the core analysis with read
lengths 28 **and** 29 excluded from **all six libraries** (GSM3900208,
GSM3900209, GSM3900210, GSM3900211, GSM3900212, and GSM3900213). Do not single
out GSM3900210 for read removal in either the primary or sensitivity datasets.
The primary production output remains unchanged and includes both lengths 28
and 29. The joint 28+29 sensitivity analysis is a separate, explicitly labeled
analysis and must not be used to revise either frozen offset.

Optional 28-only and 29-only exclusions may be run only as descriptive
technical diagnostics. They are not substitutes for the mandatory joint 28+29
sensitivity analysis and must not be used for inferential selection or
recalibration.

This document pre-registers the sensitivity design only. No biological
analysis was run for this amendment.
