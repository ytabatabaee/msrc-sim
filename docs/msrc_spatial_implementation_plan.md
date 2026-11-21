# Spatial MSRC implementation and validation plan

## Architecture and completed prototype

Keep the production ARG, legacy single-locus MSRC, and neutral SMC′ unchanged.
The isolated `msrcsim.spatial_reservoir` package exposes
`simulate_reservoir_arg(tree, rec_map, history, marker, suppressed_region,
kappa, seed, sampled_backgrounds)`. It returns the existing ARGResult, which
can be converted to a tskit tree sequence with
`msrcsim.arg.output.tree_sequence.build_tree_sequence`. The API is executable;
this research phase does not add a production CLI.

`ReservoirBackground` supplies a versioned recombination domain (active hull
extended to the marker unconditionally), validates explicit origin and
frequencies, samples marks using the inherited integrated measure, and
replaces the inherited floor with positive-frequency rates plus an explicit
error on absent occupied backgrounds. Existing Moran/WF histories satisfy the
API, but generated trajectories hitting incompatible boundary states are
outside the validated domain. Never discard a failed chromosome and retry
without reporting it: that would condition and potentially bias the law.

State storage reuses `ARGLineage`, active intervals, `Segment`, event logs,
and nodes/edges. A carrier split creates two ancestral carriers; a transfer
between marker and active hull updates only the background. Coalescence
combines overlapping and disjoint material through the existing union
routine. This retains full active segment ancestry; output is a succinct
collection of local trees rather than an expanded tree per genomic base.

## Minimal execution

From the repository root:

```bash
MPLCONFIGDIR=/private/tmp/mpl .venv/bin/python -m pytest -q -p no:cacheprovider \
  tests/test_spatial_reservoir.py tests/test_msrc_arg.py \
  tests/test_msrc_smc.py tests/test_msrc_smc_marginals.py
MPLCONFIGDIR=/private/tmp/mpl .venv/bin/python \
  validation/spatial_reservoir/validate.py \
  --output /private/tmp/spatial-reservoir --replicates 1000
MPLCONFIGDIR=/private/tmp/mpl .venv/bin/python \
  validation/spatial_reservoir/benchmark.py /private/tmp/spatial-benchmark
```

The validation program is a reproducible two-/four-taxon runnable example.
It writes raw chromosome spatial summaries, marginal comparisons with standard
errors/95% intervals and Holm p-values, diagnostic plots, runtime data, and a
machine-readable summary. Its raw data are never window-level independent
replicates. The benchmark varies 4/8/16/32 taxa and 100/1000-bp chromosomes
under a fixed rate, with two seeds; it is a diagnostic scaling check.

## Acceptance criteria, set before interpreting results

For the simple-case prototype, use 1,000 independent chromosomes per setting,
fixed histories `p=0.2,0.5,0.8`, cross retention `k=0,0.1,1`, explicit origin
100 generations, and a 100-bp map at `r0=0.0002` per bp per generation. The
primary focal coordinate is 75, marker 0, and the whole chromosome is the
suppressed region. This isolates the cross-retention mechanism.

The primary family includes root ages, mean pair coalescence times, and
four-taxon quartet frequencies compared with independent legacy draws. Use
Holm family correction at 0.05. Also compare the two-tip root-age mean to an
independent finite-state CTMC calculation. Do not interpret lack of
significance as equivalence. Report practical intervals against absolute
quartet error 0.05 and relative root/pair mean error 10%; if the 95% difference
interval extends outside a tolerance, classify that comparison as unresolved
for practical equivalence and increase replication in a future run. D=q2-q1
and the internal-branch distribution should be added in the next validation
stage; the current figure focuses on q1 and mean root times.

Direct mark tests have 12,000 proposals per p setting and compare the
cross-parent/focal-separating probability to the integrated measure within
4.5 binomial standard errors. The neutral single-present-background test uses
1,200 chromosomes and a four-standard-error comparison of the survival
probability at `2Ne` to `exp(-1)`. These tests assess distributions, not
identical RNG streams between algorithms.

No acceptance criterion assumes that lower cross recombination must leave
structured marginal quartets unchanged: here the map-integrated switching
rate also changes. The neutral map experiment remains a separate control.

## Independent two-tip CTMC

For contemporaneous tips in one species population, transient states before
coalescence are `00,01,11`. For constant p and focal switching scalar rho:

```
Q = [ -(2*rho*p+c0),       2*rho*p,                    0;
          rho*(1-p),          -rho,                rho*p;
                  0, 2*rho*(1-p), -(2*rho*(1-p)+c1) ],
c0=1/[2Ne(1-p)], c1=1/[2Ne*p].
```

The row deficit is the coalescence hazard. Starting from opposite states,
`S(t)=e_01 exp(Qt) 1` for `t<T0`. Beyond the origin,
`S(t)=S(T0) exp[-(t-T0)/(2Ne)]`.
The expected root time is `integral_0^T0 S(t)dt + S(T0)*2Ne`.
The validation uses matrix exponentiation and quadrature independently of the
ARG event sampler. For rho=0, opposite backgrounds cannot coalesce until the
origin, so `E[T]=T0+2Ne`, another exact check.

## Next validation stages

1. Expand marginals to root/internal/pair time CDFs, all quartet frequencies,
   D, background dwell times and per-lineage switching intensities. Compute
   empirical switch counts divided by the matching integrated background
   exposure, rather than by all recombination events.
2. Use the same frozen piecewise history in every method, including branch
   size differences and changes in p. Exercise all frequency/species/origin
   boundaries. Add randomized event-generator checks on restricted sites.
3. Test explicit structure-off dispatch against Hudson/msprime and compare
   neutral spatial correlations with corrected SMC′. A constant single
   background has the neutral ARG law, not the exact SMC′ joint law.
4. Extend the spatial design to distances spanning the suppression boundaries
   and k=0,0.01,0.1,0.5,1. Measure exact tree identity, topology agreement,
   rooted RF, pair-time correlation, waiting and tract distributions, and
   background switches. Bootstrap independent chromosomes; do not treat
   adjacent windows as independent.
5. Verify restriction consistency numerically: simulate a full chromosome and
   compare its focal-site distributions to a restricted-material sampler
   using the same original map/marker domain. The generator proof gives the
   target; an independent implementation catches data-structure mistakes.
6. Compare to an MSMS implementation in compatible constant/interior-frequency
   settings or a finite forward pedigree. Do not require finite pedigree
   agreement at founder/rare-copy limits where the model deliberately differs.
7. Specify founder/bridge behavior for backward loss/fixation jumps before
   enabling automatic general Moran/WF simulations. Benchmark conditioning
   failure rates explicitly. Preserve the no-floor v1 behavior in a versioned
   mode if a new boundary rule is added.

## Complexity and performance expectations

With A active carriers, S stored interval pieces, and K rate/history
boundaries, the reused scheduler enumerates O(A²) candidate pairs and scans
recombination integration pieces each epoch. Interval union/retirement can
cost O(S*A) or worse in the current Python representation. Memory depends on
active pieces and recorded graph events; total recorded ancestry can grow
with recombination map length and total ancestral time. There is no promised
bound linear in taxon count or genome length.

Optimize after profiling: cache population/background pair buckets, maintain
incremental interval coverage counts, integrate map prefixes, thin null
outside-hull same-parent marks, index frequency boundaries, and store graph
edges directly in compact tables. Each optimization needs regression against
the unchanged generator and restriction properties. Do not parallelize windows
independently when dependence is a scientific output; chromosome replicates
can be parallelized.

A future SMC can be constructed against this precise reference law and tested
for marginal and spatial errors. Structured prune/re-coalescence proposals
must include background paths and rates at the current coordinate. With
position-dependent rho, reversibility at one coordinate does not by itself
transport the correct distribution to the next coordinate; that is an
additional derivation problem. Use bounded lookback or retained ancestral
material only with a stated discarded-state approximation and quantified
errors on the chosen observables.
