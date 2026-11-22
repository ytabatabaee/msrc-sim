# Exact multispecies theory versus spatial simulators

This validation package imports the original theory from an external checkout;
it contains no replacement quartet formula. It changes no production simulator.

Set `MSRC_THEORY_SOURCE` to the theory repository (locally `msrc-paper`), or
pass `--theory-source` explicitly. Both a repository root and its
`simulations/quartet/moran_theory` directory are accepted. No user-specific
path is embedded in this workflow.

```sh
export MSRC_THEORY_SOURCE=/path/to/theory-repository
python validation/moran_spatial/numerics.py --theory-source "$MSRC_THEORY_SOURCE" --config validation/moran_spatial/cases.json --output /tmp/numerics.json
python validation/moran_spatial/run.py --theory-source "$MSRC_THEORY_SOURCE" --config validation/moran_spatial/cases.json --output /tmp/moran-spatial --replicates 4000
python validation/moran_spatial/run.py --theory-source "$MSRC_THEORY_SOURCE" --config validation/moran_spatial/additional_cases.json --output /tmp/moran-additional --replicates 4000
python validation/moran_spatial/collect.py --inputs /tmp/moran-spatial /tmp/moran-additional --output /tmp/moran-combined
python validation/moran_spatial/aggregate.py --output /tmp/moran-combined
python validation/moran_spatial/plots.py --config validation/moran_spatial/cases.json --output /tmp/moran-combined
python validation/moran_spatial/conditional_moran.py --theory-source "$MSRC_THEORY_SOURCE" --config validation/moran_spatial/cases.json --output /tmp/moran-conditional --attempts 160000
python validation/moran_spatial/early_topology.py --theory-source "$MSRC_THEORY_SOURCE" --config validation/moran_spatial/cases.json --output /tmp/moran-early --replicates 4000
python validation/moran_spatial/history_audit.py --config validation/moran_spatial/cases.json validation/moran_spatial/additional_cases.json --output /tmp/history_audit.json --replicates 4000
MSRC_THEORY_SOURCE="$MSRC_THEORY_SOURCE" python -m pytest tests/test_moran_spatial_validation.py
```

Use the project's environment (`.venv/bin/python`) when necessary. For an
individual case add `--case CASE_ID` to `run.py`. A fresh output directory is
recommended for each run; each case is checkpointed after it completes.
`aggregate.py` reconstructs the combined summary from completed case files.

## Input and prediction API

`theory.load_theory(source)` loads the existing five modules under an isolated
package namespace, preserving relative imports, and returns their SHA256 hashes.
`predict(modules, case, rho, history=None)` calls the original full balanced
or pectinate recursion, with joint and conditional terminal-pattern output.
A `PrescribedHistory` argument selects the original prescribed-frequency
operator, not the unconditional Moran operator.

JSON cases declare the rooted dated species tree, root and branch-specific
M=2Ne, root duration/origin, initial root count law, native Moran selection,
prescribed or unconditional history mode, sampling, base recombination map,
marker, region of suppressed recombination, cross retention, reference and
positions. Positive even M is required by the simulator's integer Ne interface.
Map theorem taxa 1,2,3,4 to A,B,C,D and name internal branches AB/CD/ROOT or
AB/ABC/ROOT. These mappings are checked. For these comparisons the structured
origin is at the older root boundary: `root.age + root_duration`; an optional
`origin_age` must agree. Origins within a descendant population are outside
this adapter's present scope, rather than silently coerced.

The physical map and cross-retention map are separate: base-rate variation
changes all physical crossovers; kappa changes cross-background marks only.
`rho_effective(case,x)` integrates their product between marker and position,
splitting at every boundary. Theory rho has units inverse generations and is
not a per-base rate. The SMC frozen target is evaluated separately at x*.

## Sampling and interpretation

- `fixed`: externally fixed frequencies on every branch; no binomial daughter
  draw is performed. This is a prescribed-history model, not finite-M Moran
  evolution with its mean frequency substituted.
- `realized_moran`: one declared seed produces one full history. The original
  conditional operator evaluates that history; chromosomes are independent
  conditional replicates. No history search, survival filter or retry occurs.
- `unconditional_moran`: independently draw the root count from pi_old, evolve
  each branch with the production Moran generator/binomial inheritance, then
  sample terminal backgrounds. Compare to the original integrated theorem.
- Explicit four-bit sampling is used only for prescribed histories. Under an
  unconditional history law, imposing fixed tip bits on prior histories would
  give the wrong posterior. `conditional_moran.py` instead uses a fixed budget
  of prior-history/terminal draws, retaining only declared terminal patterns
  before any genealogy simulation. The pattern-selection probability is
  checked against the theorem. This is explicit conditioning on Z, not filtering
  trajectories by fixation/loss or simulator success. Failed selected simulations
  are recorded and never replaced.

One independent chromosome is one replicate; positions on that chromosome
are never pooled as independent observations. Practical tolerances, fixed
before the run, are ±0.05 for each q and D. PASS requires the entire pointwise
95% interval of the difference to fit the tolerance. Intervals are Wilson for
q and multinomial normal for D; they are not simultaneous bands across the
whole grid. A mismatch requires the interval to lie outside the practical
band. Failure to establish equivalence due to uncertainty remains a failed
criterion, separately labeled from a demonstrated mismatch.

Any simulator failure invalidates its corresponding whole-chromosome empirical
comparison. `observed`, errors and intervals are omitted rather than reporting
survivor-only estimates. Failure details and raw successful records are retained.
Joint probabilities use the actual pattern frequencies; conditional distributions
are never averaged uniformly over the 16 patterns.

## Outputs and numerical scope

Per case: `theory.json`, `comparisons.csv`, `draws.csv`, `failures.json`,
`history.json`, `summary.json`. Combined: comparison tables and all 16 joint
terminal-pattern tables, source/dependency/seed provenance, PNG and editable
SVG profiles. Committed results omit large raw draws; these regenerate from
configuration and seeds. Figures show pointwise Monte Carlo uncertainty.

Numerical cross-checks compare sparse and dense exponentials in both full
species-tree recursions, sparse/matrix-free shared-population calculations,
semigroup composition, all four generator dimensions, and the existing
fixed-half balanced regression. Matrix-free checks do not substitute a
shared-population result for the multispecies comparison. Python traced
allocation is reported separately from any claim about process RSS; no
large-M performance/convergence claim is made from these modest matrices.

See `docs/MSRC_MORAN_THEORY_ARG_SMC_VALIDATION.md` for results and boundary
reachability analysis. Existing mathematical conventions and production
ARG/SMC implementations are unchanged.
