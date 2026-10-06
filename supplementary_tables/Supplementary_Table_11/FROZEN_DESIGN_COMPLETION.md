# Task16 completion and gate

{
  "technical_gate": {
    "accepted": {
      "27": 14,
      "28": 15,
      "31": 16
    },
    "pass_gate": false
  },
  "annotation_genes": "NOT_EVALUATED",
  "annotation_anchors": "NOT_EVALUATED",
  "common_genes": "NOT_EVALUATED",
  "common_anchors": "NOT_EVALUATED",
  "executable": false,
  "stop_reason": "CROSS_LIBRARY_TECHNICAL_GATE_FAILED"
}

Exact accepted offsets are in GSE106329_ASITE_CALIBRATION_v0.tsv. Failed/ambiguous lengths are not rescued. Annotation and support outputs are only realized if technical gate passed. No TMD profile, fraction ratio, effect or bootstrap has been computed. The original pre-calibration contracts remain unchanged. Stop for scientific adjudication.

The 0-minute input calibrated 27 nt at A+14, 28 nt at A+15, and 31 nt at A+16. These are control-derived assignments, not an executable common accepted population. In GSM2836157, 28-nt fixed-A-site frame0 = 0.508144 and frame gap = 0.100964. In GSM2836158, corresponding values = 0.536754 and 0.144752. Both fail the pre-frozen frame0 >= 0.60 and gap >= 0.20 criteria. All other accepted-length/library checks pass. The frozen all-accepted-lengths common gate therefore fails. No reduced subset was substituted.

TMD/reference identity mapping, real/pseudo construction and common-support enumeration were not executed. Their placeholder tables report NOT_EXECUTED; these are not evidence of annotation incompatibility or zero support. No executable biological contract is asserted.
