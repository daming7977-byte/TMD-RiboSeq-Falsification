# Route A prespecified ordinal-dose secondary v0

DESIGN ONLY; never replace the primary endpoint with this secondary.

Use SEVEN_DOSE_COMMON_SUPPORT, the fixed intersection of positive CDS-wide accepted-A-site totals in all eight libraries and both footprint populations. Same frozen real/pseudo anchors within retained genes, normalization, bands and gene-level S as in the primary contract. Collapse100x only by equal gene-level averaging of rep1/rep2, requiring both. Do not pool reads or claim replication at other dose levels.

Seven states in exact order: noCHX,1/64x,1/16x,1/4x,1x,8x,100x. Their ordinal ranks r are0,1,2,3,4,5,6. No numeric CHX concentration is assigned to noCHX. Ordinary least squares with an intercept yields b(g)=sum_{r=0..6}[(r-3)*S(g,r)]/28, since sum(r-3)^2=28. Do not transform ranks, fit nonlinear terms or select doses. Secondary statistic beta=equal-gene mean b(g).

Use10000 whole-gene bootstrap draws with NumPy Generator(PCG64(5957315)), ASCII-sorted genes, rng.integers(0,n,size=n) for each draw, mean of sampled slopes, float64 and two-sided95% percentile interval using numpy.quantile([.025,.975],method='linear'). Preserve all seven states together in each gene record. Reset the same seed/order for28nt sensitivity so both lengths use identical resampling indices. Fewer than2 common genes means TECHNICALLY_NON_EVALUABLE. No secondary p-value, replacement primary verdict or additional hypothesis test is introduced.

Report beta and interval for each footprint. Report sign/interval-zero-exclusion agreement or disagreement with the endpoint as a prespecified descriptive comparison; retain the endpoint hierarchy even when they differ. No post-result rescue or alternative dose encoding.
