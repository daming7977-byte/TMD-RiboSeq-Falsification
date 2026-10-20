Current mapping: Figure 2C–D. The legacy Figure 1D–E label below identifies the historical report; frozen scientific values are unchanged.

# Figure 1D–E: GSE133111 recorded benchmark outcomes

The following numerical and gate records are transcribed verbatim from the frozen formal report. The complete saved summary JSON and model-result tables accompany Source Data and Supplementary Table 7. No values are recomputed or substituted for the plotted report precision.


GSE133111 Stage 1:
FAIL


## Primary all-frozen-length analysis

REAL-vs-PSEUDO:
FAIL

median D = +0.0197015
one-sided Wilcoxon P = 0.153195
90% Hodges-Lehmann interval:
[-0.022692, +0.099421]

COORDINATE-NULL V2:
PASS

S_real = 1.904227
null median = 1.567389
empirical P = 0.00009999

LOCALIZED DOWNSTREAM BAND:
FAIL

median P = +0.024279
median(P-Q) = -0.033197
one-sided Wilcoxon P = 0.737748

NO-GLOBAL-SHIFT:
PASS

90% Q interval:
[-0.025435,+0.042487]

epsilon:
+/-0.263034

REAGENT OMISSION:
sign reversal observed.

Full model median logFC:
+0.012408

omit siP1_1:
-0.032442

omit siP2_1/siP2_4:
-0.159371

## Mandatory read-length sensitivity

Excluding 28/29 nt footprints did not rescue the benchmark.

Decision:
FAIL

Required gates:

A real-vs-pseudo: FAIL
B coordinate-null V2: PASS
C localized downstream band: FAIL
D no-global-shift: PASS
E transcript LOO: FAIL
E gene LOO: FAIL

Real-vs-pseudo:
median D = -0.001935
P = 0.330939

Localized downstream band:
median(P-Q) = -0.036566
P = 0.744564

Thus the formal failure is not attributable to inclusion of 28/29 nt
footprint classes.

