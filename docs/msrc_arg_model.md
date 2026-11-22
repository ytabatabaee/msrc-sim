# Unified two-background MSRC ancestral recombination graph

The single shared Hudson engine in `msrcsim/arg/engine.py` supports two
versioned background adapters. Use `reservoir_v1` for new simulations within
its documented admissible domain. Select `legacy` explicitly to reproduce
historical MSRC-ARG experiments. This consolidation changes organization,
selection and metadata, not the two event generators. Structured SMC′ is a
separate approximation and is unchanged.

## Public API and CLI

```python
from msrcsim.arg import simulate_msrc_arg

# history is one already generated or frozen population-frequency history.
result = simulate_msrc_arg(
    tree, rec_map, history, marker=0, suppressed_region=(20, 80),
    kappa=0.1, seed=93, sampled_backgrounds=states,
    background_model="reservoir_v1",  # default for this new API
)
# Same interface, original background rules:
old = simulate_msrc_arg(
    tree, rec_map, history, marker=0, suppressed_region=(20, 80),
    kappa=0.1, seed=93, sampled_backgrounds=states,
    background_model="legacy",
)
```

`rec_map` is `msrcsim.arg.recombination_map.PiecewiseConstantRecombinationMap`.
The history supplies `frequency_at`, `next_frequency_boundary`, and
`terminal_frequency`, plus a positive explicit `structured_origin_age` for
reservoir_v1. Ancestral and terminal species branches can have distinct Ne.
The API draws terminal states once by Bernoulli in `tree.taxa` order if they
are omitted. Supplying every tip's state consumes no sampling RNG draws.
An optional `rng=` continues a caller generator (seed is then unused);
`record_events=False` disables logging, not the stochastic process.

The existing command handles both versions:

```bash
msrc-sim --config examples/arg_msrc_reservoir_v1.yaml
# Equivalent without reinstalling console entry points:
.venv/bin/python -m msrcsim.cli --config examples/arg_msrc_reservoir_v1.yaml
```

In YAML, `mode: arg` selects the shared engine. The recommended selection is
`arg.background_model: msrc`, with `arg.background_model_version: reservoir_v1`.
`msrc` defaults to reservoir_v1 if the version is omitted; specify it explicitly
in archived configurations. To use the old rules, set version `legacy`.
Historical `arg.background_model: msrc_two_state` still defaults to legacy.
An explicit version on either MSRC selection opts into that version. Neutral
remains `background_model: neutral`, including the historical no-key default.
Unknown names/versions now raise rather than silently falling back to neutral.

Use `arg.suppressed_region: {start: a, end: b}` for new files. The historical
`arg.rearrangement_interval` key remains accepted and takes precedence if both
are present, preserving its previous effect in historical configurations.
`background_anchor_bp`, `cross_background_suppression`, baseline-map and
sampling keys retain their previous meanings. A `frequency_history_path`
loads the same history for the whole chromosome; it needs no redundant
forward rearrangement block. Otherwise the existing Moran/WF history generator
uses the configuration and the same CLI RNG stream as before.

## State and shared event machinery

Each active carrier holds a species population, one background h in {0,1},
and ancestral segments/active intervals. A single realized frequency path
p_B(t) is shared by all loci. Recombination splits carrier material; a
same-background coalescence unions material, adding local binary ancestry
only for overlaps. Disjoint material can coalesce at the carrier level.
The scheduler, segment splitting/merging, material retirement, output graph,
and species/frequency event boundaries are shared by both adapters.

For a cut y and potential distal-parent background g, both versions use

```
mu_h,g(dy,t) = r0(y) p_B,g(t) eta(h,g,y) dy,
eta(h,g,y) = 1 if h=g, otherwise kappa_cross(y).
```

The sampler integrates and samples the same marked measure. If y is left of
the marker z, the left parent receives g and the right retains h; if y is
right of z, the right receives g and the left retains h. At y=z the current
left-side convention is retained (a zero-measure continuous-map tie).
Within-background recombination remains possible when cross retention is
zero. Map boundaries do not create events. Lower baseline r0 affects both
same- and cross-background recombination; lower cross retention affects only
the cross component. No independent background-switching clock is added.

## The two versioned policies

| Property | legacy | reservoir_v1 |
|---|---|---|
| Adapter | `MSRCTwoBackground` in `arg.backgrounds.msrc_two_state` | `ReservoirBackground` in `arg.backgrounds.reservoir` |
| Marker-to-carrier domain | Extended to marker only if stored segments contain it | Always extended to marker, independent of sampled marker ancestry |
| Same-background pair rate | `1/[2Ne max(p_h,1/(2Ne))]` | `1/(2Ne p_h)` for occupied positive p_h |
| Different-background pair | Zero while structured | Zero while structured |
| Occupied absent background | Finite-copy floor permits a positive same-state pair rate | Explicit error, including stranded lineages at history jumps |
| Out-of-range recorded p | Existing frequency clamping retained | Explicit error |
| Missing origin | Structured forever, as historically implemented | Explicit error |
| Explicit origin | At/older than origin, ordinary MSC | Same ordinary MSC boundary |

The adapters retain the existing 1e-10-generation origin tolerance and all
frequency-boundary scheduling, including otherwise inert boundaries. Those
implementation details can affect fixed-seed streams, so they were not
optimized or redesigned here.

For active hull [l,r), reservoir_v1 uses
`J=[min(l,z),max(r,z))`. If a selected cut lies between the marker and active
material, it changes the carrier's background without inventing a genealogy
edge. Same-background marks there are null state events and remain in the
proposal/event process as before. Legacy restricts J to [l,r) once the marker
is no longer stored. This is the old marker-presence-dependent rate deficit,
retained only for reproduction.

For a focal coordinate u in the reservoir model, cross marks that separate u
from z integrate to

```
q_0_to_1(u,t) = p_B(t) integral_between(u,z) r0(y) kappa_cross(y) dy,
q_1_to_0(u,t) = [1-p_B(t)] integral_between(u,z) r0(y) kappa_cross(y) dy.
```

That interval is always in J, giving the preserved focal-generator argument
in [the reservoir specification](msrc_spatial_recommended_model.md). Legacy
cannot claim this equality in states where marker storage has been lost.
Both versions have ordinary same-population pair rate 1/(2Ne) at and older
than an explicit origin. Origin erasure follows the existing MSRC convention;
it is not a newly implemented unique-founder merger.

## Reproduction and compatibility

The old `simulate_arg(tree, ARGConfig(...), rng)` remains available and
neutral by default. The old `MSRCTwoBackground` constructor and defaults are
unchanged. Direct low-level background objects historically use threshold
sampling (`terminal_frequency >= 0.5`) when no explicit initial state is
provided; that fallback remains, including on a directly constructed
ReservoirBackground. The preferred unified API always makes a Bernoulli draw
or passes the supplied complete assignment. To reproduce an old low-level
threshold run through the unified API, pass those threshold states explicitly
and use its original ancestry seed, avoiding additional RNG draws.

`msrcsim.spatial_reservoir.ReservoirBackground` and `.model.ReservoirBackground`
are aliases of the canonical class. `simulate_reservoir_arg(...)` remains a
thin compatibility call with its original signature, defaults, Bernoulli draw
order and seed behavior. No deprecation warning or extra random draw is added.
The class's historical `inversion_interval` parameter name remains internal
compatibility; public documentation and the unified API use **region of
suppressed recombination**. No inference of inversion biology is intended.

Every `ARGResult` now includes metadata recording model version, policy names,
origin and actual terminal assignments. The unified API adds seed description,
RNG source, marker, suppression interval and retention. CLI output adds
`metadata.json` and the version in `arg_summary.json`. Its resolved YAML,
frequency history, event log, local Newick/CSV and tree-sequence outputs are
retained. Runtime values are inherently variable; additional metadata fields
are not a claim of byte-identical output files. A caller-supplied mutable RNG
or generator used as seed must be archived by the caller for exact replay;
its type descriptor is not a replacement for its initial state.

## Audit and scientific limits

The refactor has no new biological model. Always-available reservoir exchange
repairs the rate deficit under the chosen effective marker-reservoir target;
it also makes the deliberate exchangeability/phase approximation explicit.
Rejecting absent occupied backgrounds and removing the old floor are the
existing reservoir-v1 boundary choices, not newly derived founder rules.
Speciation, frequency jumps, small copy numbers and fixation/loss remain
scientific concerns. General Moran/WF histories crossing incompatible states
are not declared fully supported, and failed paths must not be silently
redrawn. Genome-scale performance is not established.

Pre-refactor fixed-seed fixtures, rerun scientific results and new D/switching/
autocorrelation diagnostics are described in
[the consolidation validation report](msrc_arg_consolidation_validation.md).
The original reports and their compact results remain available. The software
consolidation verdict and the model's limited scientific validation domain
are distinct assessments.
