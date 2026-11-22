# MSRC-ARG consolidation validation (2026-11-22)

## Software and mathematical audit

One ARG scheduler now serves the primary `msrcsim.arg.simulate_msrc_arg`
interface. `legacy` and `reservoir_v1` select background adapters; the new
interface defaults to `reservoir_v1`. Historical low-level calls and
`background_model: msrc_two_state` configurations retain their defaults.
The old `msrcsim.spatial_reservoir` import and simulation wrapper remain.
No structured SMC′ source or report was changed.

The scheduler, segment operations, event measure, recombination inheritance,
species constraints and ordinary ancestral MSC phase are shared. The
reservoir adapter was moved, not rederived: adapter method ASTs and the ARG
stochastic event loop match the starting commit
`c1823e78edf257f2d6c12c86830723ee2e110008` exactly. See
[the identity audit](../validation/msrc_arg_consolidated/results/source_identity.json).

| Policy | legacy | reservoir_v1 |
|---|---|---|
| Marker exchange after marker material disappears | Disabled | Reservoir exchange remains available |
| Occupied background coalescence | Historical finite-copy floor | Exact positive-frequency rate; absent occupied state errors |
| Missing structured origin | Historical indefinitely structured behavior | Explicit origin required |
| At/older than origin | Ordinary MSC when an origin exists | Ordinary MSC |
| Low-level missing sample states | Historical threshold initialization | Existing adapter initialization |
| Unified API missing sample states | One Bernoulli draw per tip | One Bernoulli draw per tip |

Both adapters use the same physical crossover measure and inherited
background split rules. Reservoir exchange removes the legacy marker-presence
restriction; it also assumes an exchangeable unobserved marker-parent
reservoir. Removing that restriction restores the intended focal generator
under this assumption. It does not establish finite-population phase accuracy.
The finite-copy floor and missing-origin behavior remain versioned historical
conventions, not newly endorsed mechanisms. Genuine fixation/loss founder
bridges and arbitrary Moran/Wright–Fisher paths remain unresolved.

For reservoir_v1, with marker z, the focal switching rates are
`q01=p(t)*rho(u,z)` and `q10=(1-p(t))*rho(u,z)`, where rho integrates the base
map times cross-background suppression. Coalescence of a compatible pair in
occupied background h is `1/(2*Ne*p_h(t))`; different backgrounds cannot
coalesce during structure. Within-background recombination remains possible
at kappa=0. No independent background-switching clock was added. See
[the reference specification](msrc_spatial_recommended_model.md).

## Deterministic regression

The baseline fixture records scientific-output hashes from the starting
commit: 24 neutral/legacy/reservoir scenarios, including Bernoulli and
threshold initialization, missing origins, and historically accepted duplicate
map breakpoints/outside-marker legacy inputs. Five historical CLI artifacts
also match byte-for-byte. New additive metadata is excluded from scientific
hashes. Runtime is excluded. Numpy is 2.5.3; seeded identity is asserted in
this environment, not across arbitrary dependency versions.

Tests cover compatibility wrappers, generator continuation, numpy seed
variants, explicit version selection, exact rates, missing marker ancestry,
origin boundaries, rare positive frequencies, fixation/loss errors, and
complete cross-background suppression. Legacy input acceptance is retained;
reservoir-specific guards do not constrain legacy calls.

Every result records the adapter version and policy metadata. CLI runs write
`metadata.json` and include the version in `arg_summary.json`. Historical
configuration chooses legacy unless explicitly overridden; new
`background_model: msrc` chooses reservoir_v1. Unknown model/version names
now raise rather than silently falling back to neutral. This input error
correction does not affect recognized historical configurations.

## Scientific reruns

All existing criteria were retained. Independent chromosomes, not positions,
are the replication unit.

| Audit | Outcome |
|---|---|
| Neutral, 1,000 replicates | Mean root 200.9324 versus msprime 200.8666; difference CI [-17.1404,17.2720], theoretical mean 200 |
| Legacy existing suite, 1,000 | Entire scientific JSON equals the historical result |
| Reservoir existing suite, 18 settings × 1,000 | Marginal table numerically identical to the historical table; 0/63 Holm-significant discrepancies, 61/63 practical intervals pass |
| Reservoir refined four-taxon p=.2/kappa=1, 10,000 | 0/5 Holm-significant; all 5 practical intervals pass |
| Physical event-rate audit, 216,000 marks | Byte-identical to the historical rate-audit CSV |
| Supplemental unified API, 15 settings × 1,000 per version | q1/q2/q3/D, root age and mean pair time: 0/90 Holm-significant per version; practical CI pass 85/90 legacy, 82/90 reservoir |
| Focal switching exposure audit | Legacy 3/15 Holm-significant deficits; reservoir 0/15 |

The neutral scientific summary matches its historical JSON after excluding
software labels and elapsed time. The reservoir suite additionally compares
against the historical ARG (0/63 Holm-significant, 60/63 practical intervals)
and analytic checks (0/9 Holm-significant, 9/9 practical intervals). Results
are stored in [the consolidated results directory](../validation/msrc_arg_consolidated/results).

The supplemental audit uses p=.2/.5/.8 and kappa=0/.01/.1/.5/1, an explicit
origin, fixed tip states, a uniform physical map, and independent single-locus
references at the focal marker-path integral. Mean pair-time autocorrelation
at four genomic separations is reported with 250 chromosome bootstrap
resamples in `autocorrelation.csv`. These are descriptive spatial diagnostics,
not a proof of the full spatial joint law.

The preserved legacy switching deficits occur at (p,kappa)=(.2,1),(.5,.5),
and (.8,.5), with mean count-minus-integrated-rate differences -0.1199,
-0.0870 and -0.0925. These are expected consequences of the historical
marker-presence restriction, not a consolidation regression. They remain
reported failures of the intended switching generator in that mode.

Nonsignificance is not equivalence: practical intervals that fail to fit the
predeclared tolerance remain unresolved at 1,000 replicates. The refined case
resolves its five comparisons only. General fixation/loss paths, founder
conditioning, marker phase correlations and genome-scale performance are
outside this validation. No speed improvement is claimed.

## Repository test history

The first full run returned **1 failed, 313 passed**, with an intermittent
normalization failure in the unchanged analytic `test_Hm`. Isolated reruns on
both current code and the archived starting commit passed; the cause is
unknown, and the failure was not reproduced on the baseline. A second full
run passed **314 tests** with one existing validation-grid warning. After
adding two further legacy-input regression cases, the focused ARG suite
passed **81 tests**. The final full run passed **316 tests, 1 warning in 186.17 s**; see
`regression_summary.json` and `full_pytest.txt` in the results directory. No scientific assertion or tolerance was weakened.

## Reproduction

```sh
.venv/bin/python -m pytest -o addopts='' -q -p no:cacheprovider
.venv/bin/python -m msrcsim.cli --config examples/arg_msrc_reservoir_v1.yaml
PYTHONPATH=. .venv/bin/python validation/arg/validate_neutral.py --replicates 1000 --output /tmp/neutral.json --figures /tmp/neutral_figures
PYTHONPATH=. .venv/bin/python validation/msrc_arg/validate_msrc_arg.py --replicates 1000 --output /tmp/legacy.json --figure /tmp/legacy.png
.venv/bin/python validation/spatial_reservoir/validate.py --replicates 1000 --output /tmp/reservoir
.venv/bin/python validation/spatial_reservoir/validate.py --replicates 10000 --case 4,0.2,1 --output /tmp/reservoir_refined
.venv/bin/python validation/spatial_reservoir/rate_audit.py /tmp/rate_audit.csv
.venv/bin/python validation/msrc_arg_consolidated/validate.py --replicates 1000 --output /tmp/unified
```

Set `MPLCONFIGDIR` to a writable directory if needed. Historical validators
continue to exercise compatibility entry points; the supplemental audit
exercises the new primary API. Raw supplemental chromosome draws can be
regenerated with the recorded seeds; compact comparison tables are committed.
