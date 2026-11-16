# Neutral MSC-SMC validation audit

## Scope and environment

This audit validates the neutral SMC′ prototype against independent MSC
genealogies, msprime's SMC′ and Hudson ARG models, and the retained
`msrc-arg` implementation. It also checks integrated recombination hazard,
spatial dependence, a region of suppressed recombination, tree validity,
reproducibility, and runtime. It does not validate an MSRC background model.

The full run used 20 independent chromosomes and 1,000 independent MSC trees
per parameter case, with seven cases spanning 2, 3, 4, and 8 taxa, short and
long internal branches, small and large population sizes, and branch-specific
population sizes. Linked-tree uncertainty is summarized across chromosomes,
not across adjacent trees. msprime 1.4.5, tskit 1.0.4, NumPy 2.5.3, and SciPy
1.18.1 were used. For the msprime reference, population sizes are doubled when
using one haploid sample per species, so its pairwise rate `1/N` matches this
project's `1/(2Ne)` convention. Full-ARG recombination flags count the two
ancestral products; that count is divided by two. Marginal trees are simplified
before branch summaries so unary recombination nodes are not counted as
genealogy branches.

## Results

- **Spatial recombination hazard:** 11,087 integrated event hazards had mean
  0.982 (expected 1.0); the exponential goodness-of-fit test gave p = 0.112.
  The piecewise map's exact boundary-integration tests also pass.
- **Independent MSC marginal comparison: FAIL for this implementation.** The
  harness performs 395 metric comparisons for the prototype over the seven
  cases. Nine remain significant after Holm correction across the prototype's
  metrics. Deviations are concentrated in higher-recombination cases. In the
  8-taxon long-branch case, the chromosome-weighted root-age mean was 132,245
  generations versus 142,516 in 1,000 independent MSC trees (difference about
  -9.3 standard errors). In the 4-taxon long-branch case, it was 924 versus
  1,008 generations (about -6.4 standard errors).
- **Reference calibration:** msprime SMC′ and Hudson ARG had zero Holm-significant
  metric comparisons in this run (417 and 426 comparisons respectively).
  `msrc-arg` had one of 422 comparisons remain significant. This supports the
  reference scaling and analysis pipeline, while showing the prototype's
  marginal discrepancy is not explained by finite independent-tree noise
  alone.
- **Recombination map:** across κ = 0, 0.01, 0.1, 0.5, and 1, recombination
  density inside the region increased from zero at κ = 0 to approximately
  0.00001, 0.000095, 0.00041, and 0.00063 events per bp, respectively. The
  region's rate therefore controls spatial persistence. The κ = 1 case uses
  the same constant rate inside and outside the region. The region-specific
  quartet comparisons did not flag a significant difference in this finite
  sample, but have only 20 chromosome replicates and do not establish equality
  of marginal distributions.
- **Spatial dependence:** the diagnostic suite compares exact-genealogy,
  topology, rooted RF, pair-time, and total-length similarity at five genomic
  distances and three constant rates, against both msprime SMC′ and Hudson ARG.
  Machine-readable chromosome-level outputs and confidence intervals are
  emitted by the validation script. These descriptive dependence results do
  not override the failed marginal test above.
- **Tree validity and event semantics:** focused automated tests verify valid
  tip labels, nonnegative branches, ultrametricity at large generation scales,
  compatible-population re-coalescence, explicit original-edge identity for
  invisible recombination, ancestral-stem re-coalescence, and preservation of
  consecutive identical local trees.
- **Benchmark:** over 4–32 taxa, chromosome lengths 5,000 and 20,000 bp, and
  two recombination settings, median runtimes per seeded chromosome ranged
  from about 0.004 seconds (4 taxa) to 2.464 seconds (32 taxa). Peak traced
  Python memory was at most about 0.5 MiB in this small benchmark. These are
  prototype timings, not production performance guarantees.

## Interpretation and status

The map inversion and spatial event process are verified against their stated
hazard convention. The transition kernel is **not validated as an
MSC-marginal-preserving SMC′ kernel**. The failures are recorded rather than
hidden by loosening statistical thresholds. Do not use this implementation for
scientific linked-locus inference or as the basis of Phase 2 until the
conditional transition is derived against the multispecies SMC′ event
probabilities and the marginal comparisons pass over the stated parameter
grid. The two-background MSRC extension remains out of scope.

The code in `scripts/validate_msrc_smc.py` writes CSV diagnostics, a JSON
summary, and spatial/suppression figures. Reproduce the full suite with:

```bash
MPLCONFIGDIR=/tmp/msrc-mpl .venv/bin/python scripts/validate_msrc_smc.py \
  --output validation/msrc_smc --replicates 20 --reference-trees 1000 \
  --length 10000 --spatial-replicates 20 --spatial-length 20000
```

Run only the suppressed-region suite with:

```bash
MPLCONFIGDIR=/tmp/msrc-mpl .venv/bin/python scripts/validate_msrc_smc.py \
  --output validation/msrc_smc_suppression --suppression-only \
  --replicates 20 --reference-trees 1000 --spatial-replicates 20 \
  --spatial-length 20000
```

Run the reproducible four-taxon demonstration with:

```bash
MPLCONFIGDIR=/tmp/msrc-mpl .venv/bin/python scripts/demo_msrc_smc.py \
  --output smc_demo
```

The associated papers and code repositories are listed in
[`msrc_smc.md`](msrc_smc.md). The `ipcoal` package was not installed in this
environment; its published waiting-distance code was not called. The
independent implementation checks here use msprime and `msrc-arg` instead.
