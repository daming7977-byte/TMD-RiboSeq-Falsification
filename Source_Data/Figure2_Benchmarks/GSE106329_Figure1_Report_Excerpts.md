# Figure 1B–C: GSE106329 recorded values and categorical outcomes

The lines and gate table below are verbatim excerpts from the frozen benchmark report, at its recorded display precision. They are not recalculated estimates. The unrounded sensitivity table supplies the primary, UMI and no-CHX bars; first_reveal_summary.json supplies the coordinate-null probabilities.

## Primary

- domain mean effect：`+0.556 Δlog2(local/background)`；
- simultaneous 95% interval excludes zero：131/131 centers；

## Matched pseudo anchors

- domain mean effect：`+0.579`；
- simultaneously significant centers：131/131。

## Recorded gate outcomes

| Criterion | Result |
|---|---|
| Nominal scan-global band | PASS, but full-domain rather than localized |
| Effect approximately ≥0.25 | PASS only in non-deduplicated primary |
| UMI technical robustness | FAIL: 89% attenuation, no spatial band |
| P-site/alignment sensitivity | PASS for broad shift, not for specificity |
| no-CHX not significantly reversed | PASS |
| no-CHX support | FAIL / absent |
| Matched pseudo-anchor falsification | FAIL |
| Coordinate-null falsification | FAIL |
| Leave-one-gene-out stability | PASS |
| Coverage/population stability | PASS |
| Ubc6 directional support | DESCRIPTIVE POSITIVE ONLY |

The decisive failures are UMI robustness and TMD specificity against both
prespecified negative controls. They are not marginal threshold failures.


