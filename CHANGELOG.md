# Changelog

## 0.5.0

- Made the backward structured-coalescent simulation piecewise exact across
  generation-by-generation Wright–Fisher frequency boundaries.
- Added topology counts (`n1`, `n2`, `n3`) to replicate output.
- Added nearest-MSC-arm fitting, off-arm distance, minor-topology contrast,
  standard error, z score, p value, and 95% confidence interval.
- Added maximum-likelihood fits for all three single-tree MSC quartet models.
- Added constrained fits for all three pairs of two-tree quartet mixtures.
- Added nondegenerate-network and model-classification diagnostics.
- Added `msrc-sim-compare` for analyzing existing replicate CSV files.
- Separated mechanistically discordant 2:2 histories, statistically supported
  off-arm histories, and network-interior fits in prevalence summaries.
- Removed obsolete pre-0.3 tests that targeted APIs no longer shipped by the
  repository.
