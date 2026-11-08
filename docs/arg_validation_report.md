# Neutral ARG validation report

## Scope

This report covers the corrected implementation committed after the initial
neutral ARG prototype. The validation engine is repository code; msprime is
used only as an independent reference in
`validation/arg/validate_neutral.py`.

The corrected process uses pairwise rate `1/(2N)`, map-weighted crossover rate
`integral_L^R r(x) dx` over the active ancestral span, common-ancestor events
for disjoint carriers, and explicit retirement of intervals whose active
carrier count is one. Completed local roots remain in the edge table but are
not exposed to later events.

## Reproducible validation

Install the optional reference and run:

```text
uv pip install msprime --python .venv/bin/python
.venv/bin/python validation/arg/validate_neutral.py --replicates 100
```

The 100-replicate run used seed `20261108`. The resulting small summary was:

| Check | msrcsim | msprime/reference |
|---|---:|---:|
| Two-sample mean MRCA time | 184.7 | 199.2 |
| Two-sample theoretical mean | 200.0 | 200.0 |
| Two-locus MRCA correlation, 25 kb vs 75 kb | 0.930 | 0.941 |
| Heterogeneous-map mean tree count | 1.22 | 1.15 |
| Four-taxon species-tree concordant topology | 0.590 | 0.610 |
| Expected concordant topology | 0.595 | 0.595 |

The expected quartet value uses `t = 100/(2N) = 0.5` coalescent units for
the combined two 50-generation ancestral branches, so
`1 - (2/3) exp(-t) = 0.595`. This avoids applying the formula to only one of
the two branch lengths.

Mean total branch lengths for the single-population 100-generation checks
were, for 4/8/16 taxa respectively, `760/1044/1345` generations for msrcsim
and `750/1104/1296` for msprime. These are small exploratory replicate counts,
not a claim of a formal equivalence test; increase `--replicates` for tighter
confidence intervals.

## Diagnostics and performance

Run the fixed-position spatial diagnostic with:

```text
.venv/bin/python scripts/arg_diagnostics.py arg_output/simulated_arg.trees \
  --output arg_linkage_diagnostics.png --csv arg_linkage_diagnostics.csv
```

It bins actual pairwise genomic separations and reports topology agreement,
Pearson MRCA-time correlation, and Pearson internal-branch correlation. The
chromosome is one linked replicate; uncertainty must come from independently
seeded chromosomes.

The benchmark harness is `validation/arg/benchmark.py`. It records runtime,
peak Python-traced memory, event count, and peak active lineages for 4/8/16
taxa and three length/rate combinations, alongside msprime runtime. The main
cost in this version is Python event enumeration over all same-population pairs
and interval-boundary bookkeeping; tskit construction is performed only after
simulation.

## Remaining limitations and readiness

The implemented validation is deliberately small and does not yet cover a
large 8/16-taxon topology-frequency grid, formal confidence-interval reports,
or every arbitrary species-tree shape. Those are recommended before a release
claim. The neutral process itself now has the required segment-aware Hudson
behavior for constant and piecewise-constant maps, valid rooted marginal
trees, and no independent-block resampling.

The engine is suitable as the base for the two-background MSRC phase, but the
MSRC extension should wait until a larger replicate grid is run and archived.
No background switching, inversion suppression, or frequency-dependent
coalescence is implemented here.
