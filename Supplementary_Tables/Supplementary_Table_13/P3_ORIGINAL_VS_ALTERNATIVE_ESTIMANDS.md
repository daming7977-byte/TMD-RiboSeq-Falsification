# Original versus alternative estimands

**Frozen estimates independently verified.** Both available-case point estimates, both replicate-specific estimates and both frozen percentile endpoints match their approved source table at absolute tolerance 1e-10 (10/10 comparisons). A second independently written groupby implementation confirms all 30 branch point metrics. No original or approved scientific value is replaced.

The original statistic is a log contrast of class density means, equally averaged over fixed centers and two replicate labels. It is not the mean of a gene's FIRST-minus-INTERNAL log contrast. Resampling genes changes uncertainty, not the point weighting. Canonical TMD units are geometry-resolved and need not correspond one-to-one to distinct protein molecules.

Equal-gene density means first average all units within each gene/class/center, then average genes before taking logs. This prevents a gene with many INTERNAL units from automatically receiving many unit weights, but higher-density genes can still dominate a class density total. A mean within-gene log contrast removes that intensity aggregation only on a strictly positive-density subset where the per-gene logs exist. No pseudocount or selective positive window was used.

## All prespecified results

All intervals below are conditional gene-resampling diagnostics, not biological confidence-coverage validation. The17-gene rows fail the prespecified n >= 30 robustness standard and are reported as sparse descriptive estimates. Common-transcript and within-gene mean points are identical because their genes are identical; small interval differences reflect independent prespecified random seeds, not new evidence.

| Region | Estimand | Genes | Estimate | 95% percentile | 95% basic | Scope |
| --- | --- | --- | --- | --- | --- | --- |
| Plateau | Frozen equal-unit | 1750 | 0.932547 | [0.712214, 1.132831] | [0.732263, 1.152879] | Conditional |
| Plateau | Equal-gene density, all | 1750 | 0.910209 | [0.681452, 1.110474] | [0.709945, 1.138967] | Conditional |
| Plateau | Unit weighting, shared genes | 410 | 1.262429 | [1.022943, 1.503658] | [1.021200, 1.501915] | Conditional |
| Plateau | Equal-gene density, shared genes | 410 | 1.185824 | [0.963258, 1.412960] | [0.958688, 1.408389] | Conditional |
| Plateau | Unit weighting, positive shared | 17 | 1.374591 | [1.034473, 1.747444] | [1.001737, 1.714708] | Sparse: 17 genes |
| Plateau | Equal-gene density, positive shared | 17 | 1.258929 | [0.964158, 1.640055] | [0.877803, 1.553700] | Sparse: 17 genes |
| Plateau | Within-gene log mean | 17 | 1.215265 | [0.734448, 1.699407] | [0.731123, 1.696082] | Sparse: 17 genes |
| Plateau | Within-gene log median | 17 | 1.315505 | [0.840655, 1.929129] | [0.701880, 1.790355] | Sparse: 17 genes |
| Plateau | Within-gene, common transcript | 17 | 1.215265 | [0.728330, 1.696865] | [0.733665, 1.702200] | Sparse: 17 genes |
| Plateau | Core support, unit weighting | 1749 | 0.929540 | [0.707400, 1.132856] | [0.726224, 1.151681] | Conditional |
| Plateau | Core support, equal genes | 1749 | 0.908495 | [0.681428, 1.112948] | [0.704043, 1.135563] | Conditional |
| Plateau | Core coarsened standardization | 1749 | 0.920510 | [0.681289, 1.124660] | [0.716360, 1.159732] | Conditional |
| Plateau | Extended support, unit weighting | 575 | 0.881041 | [0.440027, 1.291033] | [0.471049, 1.322055] | Conditional |
| Plateau | Extended support, equal genes | 575 | 0.688740 | [0.273592, 1.109450] | [0.268031, 1.103889] | Conditional |
| Plateau | Extended coarsened standardization | 575 | 0.590465 | [0.119587, 1.042098] | [0.138833, 1.061344] | Conditional |
| Broad | Frozen equal-unit | 1848 | 0.651820 | [0.357358, 1.055407] | [0.248233, 0.946283] | Conditional |
| Broad | Equal-gene density, all | 1848 | 0.675684 | [0.389574, 1.021114] | [0.330253, 0.961793] | Conditional |
| Broad | Unit weighting, shared genes | 407 | 1.264204 | [1.044947, 1.449914] | [1.078495, 1.483461] | Conditional |
| Broad | Equal-gene density, shared genes | 407 | 1.163146 | [0.950503, 1.363106] | [0.963185, 1.375788] | Conditional |
| Broad | Unit weighting, positive shared | 17 | 1.365322 | [0.993718, 1.646117] | [1.084527, 1.736925] | Sparse: 17 genes |
| Broad | Equal-gene density, positive shared | 17 | 1.227367 | [0.927181, 1.541251] | [0.913482, 1.527552] | Sparse: 17 genes |
| Broad | Within-gene log mean | 17 | 1.256206 | [0.850619, 1.636645] | [0.875767, 1.661793] | Sparse: 17 genes |
| Broad | Within-gene log median | 17 | 1.287536 | [0.815645, 1.675791] | [0.899281, 1.759427] | Sparse: 17 genes |
| Broad | Within-gene, common transcript | 17 | 1.256206 | [0.849438, 1.640033] | [0.872379, 1.662974] | Sparse: 17 genes |
| Broad | Core support, unit weighting | 1847 | 0.648369 | [0.358843, 1.047198] | [0.249541, 0.937896] | Conditional |
| Broad | Core support, equal genes | 1847 | 0.673460 | [0.382261, 1.021666] | [0.325255, 0.964659] | Conditional |
| Broad | Core coarsened standardization | 1847 | 0.683759 | [0.398387, 1.015769] | [0.351748, 0.969131] | Conditional |
| Broad | Extended support, unit weighting | 966 | 0.624611 | [0.234288, 1.004013] | [0.245209, 1.014934] | Conditional |
| Broad | Extended support, equal genes | 966 | 0.463166 | [0.100687, 0.820703] | [0.105628, 0.825645] | Conditional |
| Broad | Extended coarsened standardization | 966 | 0.408269 | [0.034515, 0.750890] | [0.065647, 0.782022] | Conditional |

## Population and weighting decomposition

The full machine-readable decomposition is `evidence/estimand_decomposition.tsv`. On identical extended-support populations, plateau estimates move 0.881041 ->0.688740 ->0.590465 when changing unit weights to gene weights and then common-support standardization; broad estimates move0.624611 ->0.463166 ->0.408269. Starting from the full frozen populations, restriction alone changes plateau by -0.051506 and broad by -0.027209. These are algebraic decompositions of chosen estimands, not causal estimates of how much biology is explained by composition.

Restricting to shared genes increases the unit-weighted estimates to 1.262429 and 1.264204. Equal-gene weighting on those same populations yields1.185824 and 1.163146. Both differ from the frozen population questions. On the identical 17 positive-support genes, equal-gene density means yield1.258929/1.227367, whereas means of gene log contrasts yield1.215265/1.256206. Thus the large shift from the frozen statistic cannot be attributed to aggregation order alone.

## Non-equivalence and limits

FIRST and INTERNAL differ structurally in ordinal position; balancing ordinal index across them is impossible. The extended design standardizes coarsened relative CDS position and other measured quantities only within observed overlap. Gene IDs are shared across compartments, but sample-level culture/batch pairing remains unverified. None of these contrasts randomizes topology, creates biological replicates, establishes causal insertion/pausing, or supplies independent validation of the selected downstream regions. Both regions were selected after the original pointwise results were known.
