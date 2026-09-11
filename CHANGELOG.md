# Changelog

## 0.8.6
- Add exact finite-state Wright-Fisher DP utilities for rearrangement count
  transitions, persistent-polymorphism probabilities, and four-tip arrangement
  pattern probabilities.
- Add `msrc-sim-pattern-probabilities` to compare theoretical pattern classes
  and quartet-partition weights against Monte Carlo simulations.

## 0.8.5
- Add regression-based threshold validation with bootstrap intervals and an
  independent-block statistical-consistency benchmark.

## 0.8.4
- Add a backward-compatible parameter-grid validation layer for the v0.8.3
  linked-spatial robustness benchmark, including threshold, kappa calibration,
  correction-strategy CSV summaries and PDF diagnostics.

## 0.8.3
- Replace deterministic benchmark genealogy-block lengths with stochastic
  Poisson/exponential breakpoint sampling while preserving
  `rho_inside = kappa * rho_background`.
- Keep rearrangement interval starts and ends as structural boundaries without
  forcing internal one-window blocks.
- Add per-replicate and aggregate linkage diagnostics for block counts,
  mean/median block length, breakpoint density, inside/outside density ratio,
  expected `kappa`, and rearrangement intervals with 1, 2, or 3+ genealogy
  blocks.
- Preserve MSRC/MSC quartet marginal draws, robustness strategies, and existing
  CLI defaults.

## 0.8.2
- Fix the linked spatial benchmark block process so `kappa` is interpreted as a
  multiplier on the rearranged genealogy-breakpoint rate. Expected rearranged
  block length now scales as approximately background block length divided by
  `kappa`.
- Preserve structural boundaries at rearrangement interval entry and exit
  without forcing every rearranged window into a separate genealogy block.
- Keep quartet marginal generators and species-tree robustness strategies
  unchanged.
- Add linkage diagnostic output with mean block counts, mean block lengths, and
  breakpoint densities inside and outside rearrangements.

## 0.8.1
- Refine the species-tree robustness benchmark with proper replicated sweeps
  over rearrangement fractions from 0.00 to 0.60 in 0.05 increments.
- Separate genealogy-block collapse from rearrangement-interval collapse:
  `genealogy_block_collapse` gives each `block_id` total weight 1, while
  `rearrangement_interval_collapse` gives each rearrangement interval bounded
  total weight 1 even when it contains multiple genealogy blocks.
- Add genuinely soft MSRC-aware weights using `w_l = 1 - P_l(MSRC)`, with
  oracle 0/1 probabilities for debugging and noisy probability simulation via
  sensitivity, specificity, and Gaussian noise controls.
- Add paired benchmark mode that keeps one baseline MSC chromosome realization
  fixed while increasing central fractions are replaced by MSRC signal.
- Add recovery summaries and PDF figures for weighted quartet support and
  `P(inferred quartet = T1)` with constrained binomial confidence intervals and
  the theoretical flip threshold when available.

## 0.8.0
- Add opt-in linked spatial genealogies with real bp coordinates, ordered windows,
  contiguous block IDs, rearrangement intervals, and rearrangement-specific
  genealogy-breakpoint suppression via `kappa`.
- Preserve existing MSRC and MSC marginal genealogy simulators: linked mode
  resamples only when entering a new piecewise-correlated block and is not a
  full ARG.
- Add MSRC-aware four-taxon species-tree robustness utilities and
  `msrc-sim-species-tree-robustness` for all-window, oracle-filtered,
  block-collapsed, and soft-weighted quartet support comparisons.
- Export coordinate-level `spatial_genealogies.csv` and Newick gene-tree files
  for later empirical/SBI features or external ASTRAL runs.

## 0.7.0
- Add `msrc-sim-plot-history` for schematic Wright-Fisher frequency-history
  trees from existing `frequency_history.csv` output.
- Add opt-in `output.make_history_plot` support for generating that figure
  directly during mechanistic runs.
- Keep the history-tree visualizer separate from replicate/model-comparison
  plotting.
- Add spatially ordered MSRC profiles and `msrc-sim-hybridization`, a
  pulse-hybridization ancestry-tract comparator with moving-window quartet
  summaries and an automated chromosome-wide profile figure.

## 0.6.0
- Freeze and replay complete Wright–Fisher rearrangement histories.
- Export branch-level trajectory summaries and integrated arrangement exposure.
- Add Benjamini–Hochberg FDR-adjusted off-arm tests.
- Add automated quartet-simplex, off-arm, pattern-prevalence, and model-comparison figures.
- Preserve all v0.5.1 single-run, replicate, grid, and comparison commands.

# Changelog

## 0.5.1

- Preserve and normalize four-character terminal arrangement patterns (for example, `0101`).
- Separate network geometric representability from statistical model preference.
- Add `network_representable`, `off_arm_supported`, `network_aic_preferred`, and `network_strongly_preferred`.
- Add practical network-boundary diagnostics for gamma and branch lengths.
- Split model output into geometry and evidence classifications.
- Expand replicate and prevalence outputs with these diagnostics.


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
