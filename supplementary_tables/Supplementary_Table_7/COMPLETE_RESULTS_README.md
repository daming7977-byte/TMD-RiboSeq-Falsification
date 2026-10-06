# Complete GSE133111 benchmark results

The Complete_GSE133111_Results folder supplies the previously completed, unrounded summary JSON and full edgeR result tables for both frozen variants. These files accompany Supplementary Table 7; no additional test or analysis was performed to prepare them.

| Frozen variant | Matched real–pseudo transcripts | Real–pseudo pairs | Localization transcripts |
|---|---:|---:|---:|
| All frozen lengths | 899 | 1,903 | 460 |
| Excluding lengths 28 and 29 in all six libraries | 810 | 1,664 | 396 |

The counts are transcribed from the summary JSON files. They refer to analysis units, not biological replicates. The localization and matched real–pseudo tests have different eligibility populations. The complete TSV files retain their original column names, precision and NA entries. PValue/FDR entries intentionally suppressed for descriptive reagent-omission results have not been filled. The original exports did not retain W/z or edgeR degrees-of-freedom fields; no values for those fields are supplied or inferred here. The formal benchmark classification remains FAIL, and coordinate-null PASS does not establish TMD-specific biological validation.
