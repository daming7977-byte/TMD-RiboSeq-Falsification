# GSE133111 Stage 0C final Paper 3 analysis-contract freeze

Date: 2026-08-25  
Freeze version: `stage0c_paper3_analysis_contract_v1_2026-08-25`  
Decision: **PASS**

## 1. Scope and reveal boundary

This is the final pre-reveal contract for the first GSE133111 Paper 3 biological analysis. It freezes all target/background, design, estimand, model, multiplicity, control, robustness, and benchmark-gate decisions before any result is inspected.

```text
DECISION TABLE: 25/25 RESOLVED
BIOLOGICAL RESULTS INSPECTED: NO
BENCHMARK OUTCOMES INSPECTED: NO
STAGE 1 ANALYSIS RUN IN THIS TASK: NO
```

No GSE133111 biological/profile matrix, depletion-versus-control effect, real-versus-pseudo result, coordinate-null statistic, benchmark outcome, or individual-gene response was opened or calculated during this freeze.

## 2. Primary design and estimand

The primary model contains six Ribo-seq biological libraries as separate response columns:

- Control: GSM3900208, GSM3900209, GSM3900210;
- Knockdown: GSM3900211 (`siP1_1`), GSM3900212 (`siP2_1`), GSM3900213 (`siP2_4`).

```text
DESIGN: ~ condition
REFERENCE: Control
CONTRAST: Knockdown - Control
PRIMARY COEFFICIENT: conditionKnockdown
POSITIVE BETA: greater exposure-adjusted local target occupancy in Knockdown
```

RNA-seq is orthogonal/secondary and is excluded from the local TMD Ribo-seq response, offset, design, dispersion, and contrast. The authoritative design matrix is `docs/STAGE0C_GSE133111_DESIGN_MATRIX_2026-08-25.tsv`, SHA256 `86679bd767900a7d705867dcbfedb293f4520676c85ea854ddfc8d3b6102821d`.

## 3. Primary universe, target, and background

The primary real-TMD starting universe is exactly 15,597 primary-window-valid alpha-helical `TMhelix` anchors. Sample-dependent inferential eligibility is then mechanical under the inherited Paper 1 rules.

```text
PRIMARY CENTER: +75
TARGET: [+60,+91)
TARGET LENGTH: 31 nt
BACKGROUND THRESHOLD: background_count >= 20 in every one of six Ribo libraries
ELIGIBILITY: strict six-library intersection
INVALID VALUES: missing, never zero-filled
ALL-ZERO TARGET ROWS: retained in accounting, excluded from inference
```

The ordinary background is translated CDS after removal of every frozen DeepTMHMM `TMhelix`, `signal`, and `Beta sheet` interval under the inherited exclusion:

```text
[A-30,max(E,A+136))
```

clipped to CDS bounds, followed by removal of the complete tested target. The primary TMD universe remains `TMhelix` only; excluding `signal` and `Beta sheet` from background does not relabel them as primary TMDs.

## 4. Primary model and multiplicity

The primary inference inherits Paper 1's edgeR negative-binomial quasi-likelihood framework:

```text
response: integer target-window P-site count
offset: log((background_count+0.5)*31/background_length)
estimateDisp(..., robust=TRUE)
glmQLFit(..., robust=TRUE)
QL F-test
BH FDR < 0.05
```

At `+75`, every inferential real anchor belongs to the single primary confirmatory BH family. Each anchor remains a separate primary response row; there is no gene/transcript cluster term in the primary edgeR model.

The ten other center-specific BH families are mandatory predeclared secondary analyses. The complete anchor-by-center family across all 11 frozen centers is also mandatory secondary. Every positional family uses FDR 0.05; no center may replace `+75` as primary.

## 5. Matched real-versus-pseudo contract

The annotation-matched universe contains 9,829 pairs. A pair enters inference only if real and pseudo loci independently pass identical frozen structural validity, six-library `background_count >= 20`, strict-intersection, and model-estimability rules.

For the pseudo segment, with `A` equal to pseudo-segment start/pseudo anchor and `E` its end-exclusive coordinate, apply the full inherited focal background exclusion:

```text
[A-30,max(E,A+136))
```

then remove the 31-nt pseudo target. Retain the ordinary P3-08 real-topology masks. The length-matched pseudo segment is a geometric control only and must not be labeled or interpreted as a predicted TMD.

For eligible pair `i` and transcript `t`:

```text
delta_i = beta_real_i - beta_pseudo_i
D_t = mean(delta_i for eligible pairs in transcript t)
```

Only transcript-level `D_t` values are inferential units. The primary global test is a one-sided Wilcoxon signed-rank test with `alternative="greater"`. Core PASS requires `median(D_t) > 0` and `P < 0.05`. Pair-level deltas are descriptive only.

## 6. Coordinate-null v2 contract

Coordinate-null v2 is the sole authoritative coordinate-null version:

```text
POLICY: stage0c_coordinate_null_policy_v2_2026-08-25
PRODUCTION: stage0c_coordinate_null_realizations_v2_2026-08-25
UNIVERSE: 13,060 anchors / 3,997 transcripts
REALIZATIONS: 10,000
REPRODUCTION MISMATCHES: 0
CANONICAL CONTENT SHA256: c98f4b899a3f01398000b246b04c22002560b4b6f787f87801673b61e287aaf4
MATRIX FILE SHA256: 7983d6dc1c382e57521e676669891a03271a4f1ffda0dda8d9ecf101e01071b3
V1: SUPERSEDED BEFORE UNBLINDING; PROHIBITED FOR INFERENCE
```

The authoritative downstream addendum is `stage0c_coordinate_null_statistical_addendum_v2_2026-08-25`, recorded in `docs/STAGE0C_COORDINATE_NULL_STATISTICAL_BACKGROUND_ADDENDUM_V2_2026-08-25.md`, SHA256 `43b5a9340a5fb2f304d51b899a2c496c954b1c91158f1ca28860d56b78dc3e8e`.

For null coordinate `q`, remove `[q-30,q+136)` in addition to the frozen ordinary topology masks, then explicitly remove `[q+60,q+91)`. For anchor `i` and library `s`:

```text
l_is = log2((((target_count_is+0.5)/31))/((background_count_is+0.5)/background_length_is))
```

First average over all 13,060 anchors within each library, then average the six library means to obtain `S_real`. Apply the identical statistic to each complete realization `b` to obtain `S_null,b`.

```text
p_empirical = (count(S_null,b >= S_real)+1)/10001
```

Coordinate-null core PASS requires `S_real > median(S_null)` and right-tailed `p_empirical < 0.05`. This is a coordinate-localization control, not a Knockdown-versus-Control test. The primary-model count threshold does not filter this fixed v2 universe, and realization-specific filtering is prohibited.

## 7. Localized-band and no-global-shift gates

Use the common all-center eligibility set `I*`. For transcript `t` and frozen center `c`:

```text
R_t(c) = mean beta_real over eligible anchors in transcript t
B = {+60,+75,+90}
O = the other eight frozen centers
P_t = mean R_t(c) over B
Q_t = mean R_t(c) over O
```

LDB-1 PASS requires all of:

- `median(P_t) > 0`;
- `median(P_t-Q_t) > 0`;
- one-sided Wilcoxon signed-rank `P < 0.05` for `P_t-Q_t > 0`.

NGS-1 freezes `epsilon = log2(1.20)`. Calculate a 90% signed-rank/Hodges-Lehmann-compatible confidence interval for the transcript-level location of `Q_t`. NGS PASS requires the entire interval to lie strictly inside `[-log2(1.20),+log2(1.20)]`. A nonsignificant test against zero is not an equivalence substitute.

## 8. Robustness and benchmark gate

Reagent-omission refits are descriptive:

1. omit Treatment_A (`siP1_1`);
2. omit Treatment_B (`siP2_1`, `siP2_4`).

Recompute eligibility in each refit and report direction, median real-anchor beta, Spearman correlation with the full model, and sign concordance. No reagent-specific P/FDR discoveries are permitted.

Transcript LOO is mandatory. Gene LOO is mandatory using the authoritative unambiguous one-to-one transcript-to-gene mapping covering 19,192/19,192 frozen transcripts. Each mandatory cluster omission must retain positive direction and the corresponding `<0.05` core conclusion. Anchor LOO is descriptive only.

BENCHMARK PASS requires the intersection of:

1. real-versus-pseudo core PASS;
2. coordinate-null core PASS;
3. LDB-1 PASS;
4. NGS-1 PASS;
5. every mandatory transcript/gene LOO retaining required direction and core conclusion.

Failure of any required component is `FAIL`. If all required components pass but either reagent-omission summary reverses the full-model median real-anchor beta sign, the decision is `PASS WITH CAUTION`; otherwise it is `PASS`. No minimum number of FDR-significant primary TMD anchors is required. A null primary biological-result set is scientifically permissible.

## 9. Runtime and executable provenance

```text
R: 4.6.0 (2026-04-24)
R executable: /usr/local/bin/R
Rscript executable: /usr/local/bin/Rscript
R_HOME: /Library/Frameworks/R.framework/Resources
edgeR: 4.10.3
edgeR path: REPRO_WORKSPACE/Library/R/arm64/4.6/library/edgeR
limma: 3.68.5
Platform: aarch64-apple-darwin23
OS: macOS Tahoe 26.5.2
BLAS: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib; LAPACK 3.12.1
Locale: C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8
Time zone: Asia/Tokyo
Library paths: REPRO_WORKSPACE/Library/R/arm64/4.6/library; /Library/Frameworks/R.framework/Versions/4.6/Resources/library
```

Full loaded session:

```text
attached base packages: stats, graphics, grDevices, utils, datasets, methods, base
other attached packages: edgeR_4.10.3, limma_3.68.5
loaded namespaces: compiler_4.6.0, grid_4.6.0, locfit_1.5-9.12, lattice_0.22-9, statmod_1.5.2
```

The complete captured output and executable/library paths are frozen in `docs/STAGE0C_PAPER3_R_SESSIONINFO_2026-08-25.txt`, SHA256 `fbbdbfc873eba33f89a2f176f55927cc5f7f61bb3486f7401e7574bb81bd8ab4`.

The installed edgeR and limma packages report that they were built under R 4.6.1. They loaded successfully under the human-approved R 4.6.0 runtime; the warning is preserved as provenance and no dependency was installed or changed.

Paper 3 adapter provenance at contract freeze:

```text
Repository HEAD: e76c51768864c3b8b111cfdf3eed4f2e5f9411a8
Stage 1 Paper 3 adapter executed: NO
Stage 1 Paper 3 adapter frozen: NO — implementation must cite this contract and record its own version/hash before reading biological outputs
Paper 1 backend may not be run unchanged because its 2x2 design is incompatible
```

## 10. Paper 1 clean-release authority

| Artifact | SHA256 |
|---|---|
| `docs/STAGE0C_PAPER1_CLEAN_RELEASE_ANALYSIS_CONTRACT_AUDIT_2026-08-25.md` | `865152609651f3a1a8991e98879f50b4e9dc10a64469582699c536132b1a46c1` |
| `configs/paper_v1.0.yaml` | `e101330f566a5545489f18c9f3e31ed449a80675cd89c6cd6d824946b7fa28b6` |
| `docs/methodology.md` | `75f335f7adc494f6f0898ecf7e1252f68d413168843a721d7721ca0173a21389` |
| `src/tmdriboseq/profiles.py` | `7dd268704dba98fd5739ff3064a905b529cc56007c06c25feb3e10bf0eb36eb2` |
| `src/tmdriboseq/coverage.py` | `5be9cd83a9707bc936a1acf6d7d3a18649d73e64fb2d81f987f673a41a938b53` |
| `R/fit_edger_offset.R` | `ed1e258e7d5cfbce7a0b939fb054bea06879922faabdc298e9c6a16aec050631` |
| `R/fit_edger_positional_scan.R` | `02002b79ec2d5e385fe492d7ef52575c020edd73d5a7e7284ea48aaddea5f7f7` |
| `R/fit_target_background_sensitivity.R` | `363a99e390dbf5d2771087d4d56b5c8e2cea85c33520a8dfd9bd431fea9ab55d` |

These files provide method provenance. The frozen Paper 3 decisions in this record override nontransferable Paper 1 design fields.

## 11. Upstream freeze manifest

| Component | Freeze version | Primary record SHA256 |
|---|---|---|
| P-site | `stage0c_psite_v1_2026-08-24` | `77aa0d5e624833abc57ffdeae43314bc82c2e3e37e1b332ab0b96cd768ad3ce9` |
| DeepTMHMM source | `stage0c_deeptmhmm_v1_2026-08-24` | `579d8d9c62142527f6dab1bb121c8da1ca0e6afc990081ffd9f4bbdeae57ab79` |
| Full raw topology | `stage0c_deeptmhmm_raw_v1_2026-08-24` | `922bcc9857c816eaf1184e92faf62b7ad71f0b7d118c64b0d34929a800162c63` |
| TMD label | `stage0c_tmd_label_policy_v1_2026-08-24` | `f0830549719e9da5bad559303b88d6b0b53cdaf7e8cf9843cdd7b44324082a78` |
| TMD anchor | `stage0c_tmd_anchor_policy_v1_2026-08-24` | `e97592a9ffd832ac199f1776fcbfe835488d9619294214f0a9b820aaeec62b87` |
| Primary anchor table | `stage0c_primary_tmd_anchors_v1_2026-08-24` | `957720c5d7fbf07e19a7c5b7d42ae77837432367313296fe6695f3ed8389c29e` |
| Metagene/window | `stage0c_metagene_window_policy_v1_2026-08-24` | `2e6082eb50f2b3e7a92fbc982ebcc550dc5fddd372486638f5a7e57e1fa73b58` |
| Pseudo policy | `stage0c_pseudo_anchor_policy_v1_2026-08-24` | `ab38f4a7ffc10141539cd9cab922dcab5c6eb17e279273026338b080f7b242fe` |
| Pseudo production | `stage0c_pseudo_anchor_table_v1_2026-08-25` | `e86d6e08c1577ff36d56e9e6ca6f4884933e2f3217e22ddf1682144fb63c0a88` |
| Coordinate-null policy | `stage0c_coordinate_null_policy_v2_2026-08-25` | `6aaec94cb79864543398b41eb65e97eaba9d0ff42f21b8823488850e3187a66b` |
| Coordinate-null production | `stage0c_coordinate_null_realizations_v2_2026-08-25` | `fc9e7d65a1364f1a0a8023dfc9f1d3d1460016d4939a9b4e787fca112e22db08` |
| Coordinate-null statistical addendum | `stage0c_coordinate_null_statistical_addendum_v2_2026-08-25` | `43b5a9340a5fb2f304d51b899a2c496c954b1c91158f1ca28860d56b78dc3e8e` |

Additional authoritative data hashes:

| Artifact | SHA256 |
|---|---|
| Frozen human FASTA | `ecda9babdf126cffa6b91e5731c355db9e40db96e2c5c5f737125342961dbf56` |
| Frozen human features | `c54965e9bca3cca852f422402180b8f0811a9ed58e8c9fd7a18b92b0e3099f36` |
| Full DeepTMHMM input manifest | `e08a6caa640c57474199319f8edba10bd797cb5c006f5ee03a2b67dc73d72dc8` |
| Normalized topology | `ca60babf576a8547077d68be1f52d29835ce8829a5d2b671c4ef3fd175c94040` |
| Primary TMD anchors | `7257a07c51158ea4235c1c45a01ada899fcd78bcaddba68518bbb9a2146a92df` |
| Pseudo pairs | `1f6150d9985c82460d9f12ccc565edb416a1539ef70c761149f35b8fefba60a5` |
| Coordinate-null v2 anchor universe | `963daf18ce5a789de2c354e4742dc355a2f40a9ab09ae78bba28ce11b1c39fdd` |
| Coordinate-null v2 validation JSON | `68d7b1a089b470079a678a2567fd60e5a386ae0536241796502cea7b7aa36e47` |
| Metadata manifest | `dbf9d950cc10e46c4bf20ce1aa6edc2cb00e4981ed73794d6a30bd419a3c1f9d` |
| Final statistical decision table | `5009f1e2eccdbefaf7f86f3f909154014d590fb07970afa4682f52a35673c99a` |
| R/edgeR sessionInfo | `fbbdbfc873eba33f89a2f176f55927cc5f7f61bb3486f7401e7574bb81bd8ab4` |

## 12. Freeze decision

```text
FINAL PAPER3 CONTRACT AUDIT: PASS
DECISION TABLE: 25/25 RESOLVED
ANALYSIS CONTRACT FREEZE: PASS
FREEZE VERSION: stage0c_paper3_analysis_contract_v1_2026-08-25
PRIMARY DESIGN MATRIX: FROZEN
COORDINATE-NULL AUTHORITATIVE VERSION: V2
P3-08: FROZEN — TMhelix + signal + Beta sheet excluded from background
P3-14: FROZEN — full [A-30,max(E,A+136)) pseudo focal exclusion
BIOLOGICAL RESULTS INSPECTED: NO
STAGE 1 ANALYSIS RUN: NO
```
