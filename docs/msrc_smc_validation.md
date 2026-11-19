# Corrected multispecies SMC′ validation

## Decision

The root-age defect in the first transition implementation was corrected and
the neutral kernel now passes the recorded transition tests and the finite
parameterized marginal checks. This is a **conditional validation**, not a
proof that every implementation detail matches the complete McKenzie–Eaton
kernel. Waiting-distance comparisons against the published `ipcoal`
implementation were not run because that package was unavailable. Do not infer
that this validation applies to backgrounds or the planned MSRC-SMC extension.

## Kernel correction and derivation

At a recombination, the cut lineage can re-coalesce with guide lineages only
while the lineages occupy a compatible species-tree population. The original
edge target is eligible only until its parent coalescence time. Previously,
the event-time partition included guide-edge and species-population
boundaries, but omitted that parent time. If an exponential wait crossed it,
the target set stayed stale and allowed an impossible original-edge event,
which appeared as an invisible recombination and biased root ages downward.

The transition now adds the original edge's upper time to the next event
boundary. The process advances to that boundary, expires the original target,
and continues with the updated compatible target set and pair hazards
`1/(2Ne)`. This is a renewal of the event simulation at a real change in the
hazard set; no rate adjustment or post-hoc tree correction is used.

For a two-tip guide tree with root time `T`, a cut at `tr`, and equal pair
hazard `c=1/(2Ne)`, two targets are available up to `T` (the original edge and
the guide stem). The probability of re-coalescing to the original edge before
`T` is `0.5 * (1 - exp(-2c(T-tr)))`. Its expiration at `T` leaves only the
guide-stem target. Thus for `tr <= t <= T`,
`P(t_recoal <= t)=1-exp(-2c(t-tr))`; for `t >= T`, the CDF is
`1-exp(-2c(T-tr))*exp(-c(t-T))`. The unit tests check target expiry and the
resulting conditional time distribution against these analytic probabilities.

This is the multispecies SMC′ approximation: recombination cuts a current
genealogy edge and its detached lineage re-coalesces to an eligible lineage in
the guide genealogy, retaining the same-edge invisible outcome. The exact ARG
tracks ancestral material and can re-coalesce with lineages outside the
current marginal genealogy; the SMC′ state omits that full ancestral graph.
This distinction is substantive. The implementation does not treat arbitrary
subtree-prune-and-regraft operations as a valid multispecies SMC transition.

## Validation performed

The saved run uses 40 independently seeded linked chromosomes and 1,000
independent MSC reference trees in each of seven settings (2, 3, 4, and 8
taxa; short and long internal branches; constant and branch-specific sizes).
Marginal local-tree metrics are weighted by genomic span. For each method, 439
root-age, pair-coalescence-time, total-tree-length, and quartet-frequency
comparisons were assessed against their independent-MSC references.

| Method | Holm significant (0.05) | Nominal 95% intervals excluding zero |
|---|---:|---:|
| This SMC′ | 0 / 439 | 16 / 439 |
| msprime SMC′ | 0 / 439 | 18 / 439 |
| msprime Hudson ARG | 0 / 439 | 2 / 439 |
| msrc-arg | 0 / 439 | 26 / 439 |

The nominal counts are consistent with multiple-comparison noise; none remains
significant under Holm correction. Root-age spot checks include a four-taxon
long-branch case (reference 1008.0 generations, SMC′ 1026.1; chromosome-level
95% CI for the difference [-2.4, 38.6]) and an eight-taxon long-branch case
(142516.4 versus 142619.9; CI [-2228.1, 2435.1]). These are finite simulation
checks, not exact equality claims.

Across 23,821 integrated genomic-hazard samples, the mean was 0.9974 (expected
1) and the exponential goodness-of-fit test had p=0.731. At 1,000 bp
separation, topology agreement and pair-time correlation declined as rate
increased, in the same direction and similar range as msprime SMC′ and Hudson
ARG references. In the suppressed-region experiment, inside-region
recombination density rose monotonically with multiplier κ: 0 at κ=0,
0.000003 at κ=0.01, 0.000073 at κ=0.1, 0.000362 at κ=0.5, and 0.000827 at
κ=1 events/bp in the tested run. The finite inside/outside marginal checks
showed no Holm-significant quartet or tree-metric differences; lowering the
map rate was evaluated as a spatial-dependence change, not assumed to alter
marginal quartet frequencies.

The deterministic demo figure with tree tracks, map, change positions, and
topology autocorrelation is `validation/msrc_smc_corrected/demo_smc_map_comparison.png`.
Other compact result tables and plots are in the same directory. Raw
chromosome-level genealogies and event-level waiting samples were retained in
the temporary run directory for this audit and are intentionally not part of
the repository artifact set.

## Performance and scope

The diagnostic benchmark covers 4, 8, 16, and 32 taxa, two chromosome lengths,
and two rate levels. Median SMC′ chromosome runtime rose from about 0.004 s at
4 taxa to 0.86 s at 32 taxa; `msrc-arg` was about 0.006 s and 0.09 s at those
sizes. The current Python SMC′ implementation is slower at larger taxon counts
in this small benchmark, so no performance advantage is claimed. Maximum
traced Python allocation for SMC′ was about 0.66 MiB at 32 taxa. Timings are
seeded microbenchmarks, not production estimates.

The corrected full repository test suite and exact command are recorded in
`PROJECT_STATUS.md`. The transition-specific analytic tests include original
target expiration and a 12,000-replicate fixed-guide conditional-time check.
The marginal suite adds four- and eight-taxon high-recombination comparisons.

## Reproduction

From the repository root (with the project virtual environment and msprime
installed):

```bash
MPLCONFIGDIR=/private/tmp/mpl .venv/bin/python scripts/validate_msrc_smc.py \
  --output /private/tmp/msrc-smc-validation \
  --replicates 40 --reference-trees 1000 --length 10000 \
  --spatial-replicates 40 --spatial-length 20000
MPLCONFIGDIR=/private/tmp/mpl .venv/bin/python scripts/demo_msrc_smc.py \
  --output /private/tmp/msrc-smc-demo
```

`--benchmark-only`, `--suppression-only`, and `--resummarize-only` run the
corresponding audit sections independently. Validation results are compacted
under `validation/msrc_smc_corrected/`; regenerate raw results in a writable
output directory when needed.

## Remaining limits and next work

- Run the paper-associated conditional waiting-distance checks directly using
  `ipcoal` or the authors' released calculations and expand the tree/Ne grid.
- Profile and optimize transition tree reconstruction at higher taxon counts.
- Keep background histories out of state until a joint backward model has been
  derived for genomic recombination, background switching, selection, and
  shared demographic trajectories.
- Retain `msrc-arg` unchanged as an independent benchmark.
