# Trajectory-conditioned reservoir ARG: mathematical specification v1

## Target and approximation

The target is the legacy MSRC effective structured coalescent conditional on
one realized history `H={p_B(t)}`, a species tree, sample background vector,
and an explicit origin age `T0`. The spatial extension uses Hudson ancestral
material with a marker-background reservoir. It is a continuous-time
structured ARG, not an SMC′ and not an exact reproductive pedigree.

The new assumption is **unobserved marker-parent exchangeability**: each
crossover can draw the background of the parent distal to the marker from
`(1-p_B(t),p_B(t))`, with cross-background exchange attenuated by `k(y)`.
This reservoir remains available even if the sampled marker ancestry has
separated or resolved. No marker genealogy is required to sample this model.
Missing marker-phase correlations are approximated by the reservoir rule.
Only sampled material is tracked; discarded parental material can never
re-enter except through the same background reservoir/coalescent rules.

## State and inputs

Chromosome coordinates are continuous on `[0,L)`. Marker `z` is in `[0,L]`.
Rates `r0(y)` and cross-retention factors `k(y)` are piecewise constant;
`k(y)=kappa` in one suppressed region `[a,b)` and 1 outside. This is
cross-background suppression. A reduced baseline map is a separate mechanism
that also reduces within-background events.

At evolutionary time `t` the state is a list of carriers
`(id, population B, background h, ancestral segments, active intervals)` and
an output ancestry graph. Segments have coordinates and node identifiers;
active intervals have at least two remaining carriers and need further
ancestry. Each carrier has a single background for all its material.
Coalescence merges only matching backgrounds while structure is active;
recombination can split a carrier and assign different backgrounds to its
parents. Thus a focal lineage can change background repeatedly along a tree
branch; the background is not an immutable tip color.

A single shared history object supplies `frequency_at(B,t)`,
`next_frequency_boundary(B,t)`, `terminal_frequency(taxon)`, and finite
positive `structured_origin_age=T0`. Supplied histories may come from Moran,
WF, or frozen records. All runs conditional on a history use the same object;
no per-locus redraw is allowed. The v1 prototype has been validated only with
constant interior frequencies. Boundary compatibility must be checked before
claiming support for arbitrary generated histories.

Terminal states are supplied explicitly or sampled independently once per
taxon as Bernoulli(`p_taxon(0)`). This assignment applies to its entire sampled
chromosome. It is a sampling assumption, including across separated loci.

## Evolutionary-time event measure

For carrier `i`, let `[l_i,r_i)` be the hull of its active intervals. Define
`J_i=[min(l_i,z),max(r_i,z))`. For `t<T0`, candidate crossover/parental-state
marks `(y,g)` have measure

```
mu_i(dy,g;t) = 1_{y in J_i} r0(y) p_(B_i,g)(t) eta(h_i,g,y) dy,
eta(h,g,y) = 1 if h=g, else k(y).
R_i(t) = integral_Ji r0(y) [p_(B_i,h_i)(t) + k(y)p_(B_i,1-h_i)(t)] dy.
```

For `t>=T0`, labels are ignored and the event domain and rate reduce to the
ordinary Hudson active hull and `integral r0(y)dy`. Each same-population pair
has rate `1/(2Ne_B)`.

For `t<T0`, any two same-population, same-background carriers, including those
with disjoint material, coalesce at
`c_(B,h)=1/(2Ne_B p_(B,h)(t))`. Different-background pairs have rate zero.
At a merge, overlapping ancestry receives a new binary parent node;
disjoint ancestry is unioned without inventing a local coalescence. Resolved
material can be retired after its local root is recorded. Retiring it cannot
alter a surviving locus's projected generator under the reservoir rule.

### Inheritance at a crossover

Use the *same* candidate event to split material and assign backgrounds:

- If `y<z`, left parent gets `g`, right parent retains `h`.
- If `y>z`, left parent retains `h`, right parent gets `g`.
- `y=z` has zero probability under a continuous map; a deterministic tie may
  use the left-side rule (the implementation uses `y<=z`).

If both sides have active material, create two parent carriers, add unary
recombination parent nodes for output, and preserve their segment/active
partitions. If only the marker-distal side has active material (a cut between
the marker and hull), update its background to `g` without inserting a
coalescence or splitting empty material. If `g=h`, this outside-hull proposal
is a null state event. The prototype retains these proposals for simple
sampling; they can be thinned later. Cuts in gaps between disjoint active
intervals are legitimate carrier splits even without a local tree edge at
the cut. Background transfer and tree change are separate observables.

The total evolutionary-time event hazard is `sum_i R_i + sum_pairs c_pair`.
There is no additional independent background-switching clock.

## Complete sampler

1. Read and validate tree, maps, marker, origin, and history. Initialize one
   full-chromosome carrier per taxon with the fixed terminal assignments.
2. At time `t`, first move carriers at species boundaries to their parent
   population. At `T0`, ignore all background labels for future rates. Their
   stored integers may remain in output but have no biological rate effect.
3. Compute each carrier's integrated recombination measure and each eligible
   pair's coalescence rate. Partition map integrals at all baseline and
   suppression boundaries.
4. Determine the next deterministic boundary: minimum occupied species
   boundary, next occupied population frequency change, and `T0`.
5. Draw `W ~ Exp(total hazard)`. If `t+W` reaches the boundary, advance to it
   without a stochastic event and recompute. If total hazard is zero, advance
   directly to the next finite boundary; if none exists, report an invalid
   configuration.
6. Otherwise advance by `W` and draw one event proportional to its rate.
   For coalescence, merge the selected pair. For recombination, choose a
   piece with weight `(b-a)r[p_h+k p_cross]`, draw `y` uniformly on it,
   and select `g` with weights `p_h` and `k p_cross`. This is exactly the
   same measure integrated in step 3. Apply inheritance above once.
7. Record the event, retire complete material, update carriers, and repeat
   until every coordinate has one ancestral carrier. Export a tree sequence;
   simplify unary nodes for marginal gene-tree metrics.

The joint next-event distribution is
`P(W in dw,event=e)=exp(-A*w)*rate_e*dw` on a constant epoch; multiplying
survivals across deterministic boundaries gives the inhomogeneous law.
Every random choice is specified by this construction. No genomic scan or
`Exp(r L_G)` clock is used here: breakpoints and genealogy correlations emerge
from ancestral events across the entire chromosome.

## Projection proof for a focal locus

Fix `u` while at least two sampled focal ancestors remain. Project the state
to carriers of `u`, their background/population labels, and their focal
coalescence tree. Each such ancestor is on exactly one carrier.

A candidate changes its background from `h` to `g!=h` exactly when `y` lies
between `u` and `z`. That interval is always contained in `J_i`, independent
of other material or whether the marker is stored. Consequently

```
q_(h->g)(u,t) = p_(B,g)(t) integral_between(u,z) r0(y) k(y) dy.
```

The same-background proposals and cuts outside the marker-to-focal path
leave the projected background unchanged. Recombining one carrier never
creates two ancestors carrying the same point `u`. A pair of distinct focal
ancestors coalesces at exactly `1/(2Ne p_h)` when compatible. Coalescing a
focal carrier with a nonfocal carrier changes no focal state or future focal
rates. All events invisible under this projection have no effect on its
transition rates. Thus the full generator applied to any function of the
projected state equals the legacy one-site generator with
`rho_eff(u)=integral_between(u,z) r0 k`.

This generator closure (lumpability) establishes equality in law of the
one-site effective MSRC process, not merely matching average switches.
It holds conditional on a deterministic admissible history, and hence for
mixtures over histories drawn once with the same distribution, provided the
boundary rules also match. It proves the full local genealogy law, including
quartets, internal times, and background paths, for the stated target.

The legacy finite-copy floor equals the exact formula whenever
`p_h>=1/(2Ne)`. For smaller positive frequencies, or any occupied absent
state, the target differs. Do not claim legacy equivalence there.

## Multi-position consistency

For any finite set of sites, restrict all carrier intervals to those sites
and erase carriers with no retained site. Candidate events that affect
retained sites have the same marked measure and inheritance rule as in the
full process; events on discarded material are invisible. Same-population,
same-background pair rates do not depend on unretained material. Therefore
the generator also closes on the retained-site carrier partition. This gives
consistent two-, three-, and multiple-site laws under this effective model.
It does not determine that these laws equal finite-haplotype phase dynamics.

This is the central advantage over choosing an unspecified local-tree
transition: a complete coupled process is defined first; local trees are
projections of it. A small sufficient state in *evolutionary time* exists,
although its number of carriers/intervals is variable. A small sufficient
state in *genomic scan position* is neither assumed nor needed.

## Boundaries and established limits

- **Zero map:** no recombination or background transfer; all sites share one
  genealogy, conditional on the same terminal assignments.
- **Cross suppression `k=0`:** all cross-background proposals vanish in that
  interval. Same-background rate `r0 p_h` remains; marginal switching across
  paths elsewhere may remain nonzero.
- **No suppression `k=1`:** total crossover rate on `J_i` is `r0`; structured
  coalescence persists while two backgrounds are active.
- **No structure:** dispatch to the ordinary neutral ARG; neutral local trees
  have MSC marginals. With a single present background (`p=0` or 1 and all
  carriers in it), the projected process also has ordinary ARG laws; the
  prototype may generate additional null proposals outside active hulls, so
  identical seed streams are not guaranteed.
- **Origin:** use an explicit `T0`; at and older than it, switching is disabled
  and ordinary MSC pair rates apply. The v1 convention erases the partition;
  it does not force all derived-background copies to coalesce into a unique
  founder. This preserves the legacy boundary, not exact single-mutation
  founder ancestry.
- **Fixation/loss:** entry into an absent parental background has zero marked
  rate. A supplied carrier already in an absent background, including one
  stranded by a backward history jump or species merge, raises a named error.
  Do not substitute a frequency floor or silently resample its chromosome.
  Admissible fixed-boundary examples (all carriers in the present state) are
  well defined. General histories crossing loss/fixation require an explicit
  founder/bridge law. This remains a limitation of the v1 application domain.

Moran/WF trajectories can be supplied through the existing history interface,
but their finite-population ancestry is not automatically generated by these
conditional rates, especially at rare copy numbers and founder boundaries.
A realized path conditions the simulator's hazards; it does not prove that
the chosen conditional genealogy is the exact genealogy of the forward path.

## Prototype and scientific status

`ReservoirBackground` subclasses the old background adapter to replace the
marker-presence condition and boundary policy. `simulate_reservoir_arg`
accepts the species tree, ARG recombination map, history, marker, suppression
interval, retention factor, seed, and optional sample backgrounds. It reuses
the existing scheduler without changing production code. Baseline maps and
species-dependent population sizes are already accepted by that scheduler.
The inherited internal parameter name `inversion_interval` is compatibility
only; the scientific parameter is a region of suppressed recombination.

This v1 is a reference prototype. It carries an ancestral-material graph and
cannot yet claim genome-scale efficiency. Background reservoir/phase,
continuous-time approximation, and origin erasure are scientific assumptions;
event simulation is exact for the resulting piecewise-constant generator up
to floating-point arithmetic. Existing SMC′ original-edge expiration code is
untouched: this model has no guide-tree original-edge approximation.
