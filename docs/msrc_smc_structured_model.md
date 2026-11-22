# Structured SMC′ with a frozen reference marginal, version 1

## Scope and the chosen approximation

This is an opt-in, sequential local-genealogy simulator implemented in
`msrcsim/smc/structured.py`. Its state contains one dated rooted tree and
piecewise constant background histories, with an augmented root-stem history.
It never stores ancestral genomic material or runs the ARG simulator.

The approximation is deliberately **window-local**. At a declared reference
position `u*`, compute `rho*=integral_between(marker,u*) r0(y) k(y) dy`.
The single-locus MSRC distribution at that coordinate is the target
`pi*(G,H)`. Use that same target throughout the simulated window. Spatial
cut hazards still follow the actual baseline map and cross-retention map.
A proposed path/genealogy change is accepted with an extended-state
Metropolis ratio. There are no independent genealogy resets or background
updates alongside a neutral tree process.

Consequently exact reference-marginal preservation is a mathematical target
for this model, whereas correct marker-distance marginals at *other*
coordinates are not guaranteed. It approximates those marginals by freezing
rho at u*. It is suitable for assessing a local-window approximation and its
spatial errors; it is not an exact whole-chromosome marker-mediated SMC′.
An explicit output metadata field records the reference and the limitation.

The conditional detached-path construction inherits the state at the cut;
it integrates out unobserved marker ancestry through the single-locus
switching generator. In particular it does **not** force a new parent color
at exactly the pruning time as the reservoir ARG does. The cut-weight measure
sums candidate parental-state weights, but the immediate inheritance is an
inherited-state closure. This is the main additional spatial approximation.

## Retained state and probability distribution

At genomic coordinate x, retain `S_x=(G,H,P,x,u*)`:

- G: dated binary tree with fixed labelled contemporary tips;
- H: background jump times and states along every edge;
- root-stem H: states from the root to the origin if the root is younger than
  the origin, so a subsequent detached lineage can reattach above the root;
- P: one fixed, population-specific realized history, with an explicit origin;
- maps, marker, reference coordinate, and scan-window bounds.

The finite tree is recombination-eligible; the root stem is a reattachment
opportunity but not part of its finite spatial cut length. Backgrounds are
ignored at/older than the origin. The stored stem extends no biological tree
branch for output metrics. Population assignments are inferred from dated
clades and species-tree divergence times using the existing embedding rule.

For an occupied background h in population B before origin T0:

```
p_B,1(t)=p_B(t), p_B,0(t)=1-p_B(t),
c_B,h(t)=1/[2Ne_B p_B,h(t)],
q_h,1-h(t)=rho* p_B,1-h(t).
```

At and older than T0, Q=0 and every same-population pair rate is 1/(2Ne_B).
Tips have supplied states or one Bernoulli draw each, shared across positions.
The target is conditional on those terminal states and the realized history.

Given all background paths, the augmented density with respect to dated-tree
node ages, topology counting measure, and background jump times is

```
pi*(G,H) = product_coals c_B,h(t)
           * product_background_jumps q_h,1-h(t)
           * exp[-integral (sum_eligible_pairs c + sum_extant_lineages q) dt].
```

The integral includes the single root stem until T0 when necessary; after the
last coalescence and T0 both hazards vanish. Population, frequency, tree-node,
background-jump and origin times partition the integral into exact constant
intervals. Incompatible coalescences have zero density.

Initialization uses the existing reference MSRC event simulator with event
records, builds full-precision dated nodes and paths directly, and extends
the root state independently to T0. The prototype requires histories for
which exact positive-frequency rates match the legacy finite-copy convention
on occupied states. It raises on an absent occupied state or a sub-one-copy
occupied frequency rather than silently changing the initial law.

## Spatial waiting distance and cut selection

For an edge e, a point at age t in population B and background h has effective
candidate cut intensity

```
w_x(e,t)=r0(x) [p_B,h(t)+k(x)p_B,1-h(t)]    for t<T0,
w_x(e,t)=r0(x)                             for t>=T0.
A_x(G,H)=sum_edges integral_edge w_x(e,t) dt.
```

Sample E~Exp(1) and solve `integral_x^xnext A_y(G,H)dy=E` across map and
suppression boundaries, holding the current tree/path state fixed. For a
constant map epoch this is an exponential distance with rate A_x. Map
boundaries themselves are not events. Zero-rate intervals are crossed
without changing the state. Candidate events are retained even when they
are rejected or invisible.

Partition each edge at species, frequency, background, and origin boundaries.
Select one edge-time interval with weight `length * [p_h+k p_other]`, then
sample its cut time uniformly. r0(x) cancels from this conditional selection.
The state at the cut is inherited from H; the topology and H strictly below
the cut remain fixed. Suppress the old parent after detaching the subtree,
concatenating the sibling and parent histories to form the guide backbone.

## Conditional detached-path/re-coalescence proposal

The detached lineage starts at `(descendants,cut_time,inherited_h)`. At age t
its species population follows the species tree. It can switch background
with q, or re-coalesce with a compatible guide lineage at c per target.

Compatible targets have the same population and, before T0, the same
background. The guide includes finite edges, a continuing root stem, and the
original cut edge up to its **old parent coalescence time**. That last target
implements the SMC′ invisible outcome. It expires at the old parent time,
even when no guide edge ends there after pruning. Its original background
path is retained while it remains eligible. This preserves the neutral bug
fix that prevented spurious post-parent invisible events.

If m(t) targets are compatible, detached event rate is `q(t)+m(t)c(t)`.
Within each constant interval sample an exponential wait. If it reaches a
species, guide-node, guide-background, original-background, frequency, or
origin boundary, advance without an event and renew hazards. Otherwise
choose a switch with weight q or each target with weight c. Repeat until
reattachment. Background switching ceases at T0.

If the original-edge target is chosen, preserve dated topology and replace
only the detached path between cut and reconnection; this can be a
background-history change without a genealogy change. Otherwise insert a
new binary coalescence at the target time, splice the detached path into its
edge, and split the target's background history around the new parent.
The rest of the guide histories are inherited. Rebuild dated-tree output
from absolute node ages; no artificial coalescence is inserted to encode a
background change.

## Proposal density and marginal correction

For detached background path gamma from cut a to reattachment b, the density
of the particular target (not the target-summed event) is

```
f_S(gamma,b,target | e,a)
 = exp[-integral_a^b (q(t)+m_S(t)c(t))dt]
   * product_switches q(t_i) * c(b).
```

For an original-edge outcome the reverse auxiliary move has the same cut
and reconnection ages and restores the old path segment. For an ordinary
regraft the reverse move cuts the same retained subtree at the same age,
then restores its old path and old parent time. Pruning the proposed state
produces the same backbone; its temporary original-edge target is the new
edge, with its own path and expiration. Concatenation and swapping the old
and new coalescence/path-time coordinates have Jacobian one.

Accept with

```
alpha=min(1, pi*(S') f_S'(gamma_reverse) /
             [pi*(S) f_S(gamma_forward)]).
```

These are **continuous-position generator rates**. The unnormalised marked
cut measure w_x(e,a) is the same forward and backward because the subtree,
state and population below the cut are unchanged. It cancels in the ratio.
Do not include a ratio of total tree cut lengths from a discrete embedded
proposal chain; that would produce the wrong spatial stationary law.
For each paired auxiliary move, detailed balance is
`pi*(S) w f_forward alpha = pi*(S') w f_reverse alpha_reverse`.
Rejections retain the entire previous state.

With a shared Q and c in the target and proposal, the target and detached-path
terms cancel analytically: guide events cancel, detached migration terms
cancel, and the two temporary original-edge survival terms have symmetric
same-state overlap. The expected log ratio is zero for this inherited-state
construction. The default fast mode uses this analytic cancellation and
retains the acceptance uniform draw for seed compatibility. Setting
`audit_density: true` evaluates every target/forward/reverse density explicitly.
The primary scientific grid used that audit mode. Explicit calculations are
also available for future proposals; the fast-path cancellation must not be
reused if the inheritance or detached-path law changes. Tests check the cancellation on actual
multi-background moves; unexpected nonzero values must be investigated.

Detailed balance implies that if the initial augmented state has law pi*,
every genomic position has that law, even under spatially varying cut-map
intensities. This is a property of the stated frozen-target Markov process,
not a claim that it reproduces the reservoir ARG's joint law. Code validation
is still required; an algebraic assertion does not prove the implementation.

## What the events mean

`recombination_events.csv` contains **candidate SMC cut/path proposals** with
acceptance, cut/re-coalescence times, target kind, genealogy/topology/background
history flags, and forward/reverse density diagnostics. They are not a census
of physical crossovers in an ancestral graph. `proposed_path_switches`
counts switches generated while proposing a detached path. These are
conditional ancestral-history marks from the map-derived Q, not additional
spatial events. Do not add them to the spatial proposal count or interpret
repeated path draws as independent physical crossovers.

Only committed path histories describe the current local ancestral state;
proposal path marks rejected by MH are not committed background transitions.
`background_transitions.csv` indexes accepted changes to stochastic background
jump marks (`background_changed`), excluding edge-birth/time bookkeeping.
The separate `background_history_changed` flag also includes changes in path
representation caused by changed coalescence times or clades, even when no
background switch occurs. Neither flag is a backward-time switching-rate table. The full
paths in `background_histories.jsonl` supply jump times and dwell durations.
The old/new background columns are endpoint states of the proposed lineage,
not necessarily a transfer at the cut. These distinctions avoid counting one
auxiliary event twice or falsely labelling the approximation a physical ARG.

## Boundaries, histories and limitations

- `neutral: true` delegates directly to the corrected neutral SMC′, retaining
  its RNG semantics and original-edge behavior.
- At zero baseline recombination both the spatial proposal rate and map-
  derived switching vanish. One genealogy is shared across the scan.
- k=0 removes the cross contribution to cuts within the region; within-
  background recombination remains. Q may still receive contributions from
  unsuppressed parts of the marker path.
- k=1 does not erase structured coalescence. Holding rho fixed does not turn
  a two-background model into the MSC.
- `strict_v1` rejects occupied absent backgrounds and sub-one-copy positive
  frequencies, and requires an explicit finite origin. It does not alter the
  legacy floor. No failed genealogy/history is silently redrawn.
- Frequency changes, including time-varying Moran/WF or frozen histories, are
  shared by every local tree and included as deterministic rate boundaries.
  Arbitrary forward histories can strand a backward lineage in an absent
  state; a general founder/bridge treatment remains unimplemented.
- Origin comparisons use the existing 1e-10-generation numerical tolerance.
  This matters for a Moran origin record stored at `nextafter(T0, younger)`:
  the process is already treated as unstructured at that numerical boundary,
  matching the legacy simulator rather than querying an absent pre-origin
  state a floating-point step too early.
- The whole target is frozen at u*. Away from u*, the switching-rate error is
  `[rho*-rho_eff(x)] p_other(t)`. It vanishes where the marker-path integral
  is constant, e.g. through a fully cross-suppressed interval with no added
  marker-path exchange. For a varying integral, reduce the window or report
  the discrepancy. Invariance of pi_x separately at each x would not alone
  transport the correct marginal law along a nonstationary scan.
- Rejected proposals and same-edge events can increase persistence. Tract
  lengths and multi-distance correlations must be compared to the reservoir
  reference; marginal correctness does not establish spatial accuracy.

## Sources and implementation status

Structured path-density evaluation and simulated coalescence/migration SPR
moves have established MCMC precedents; this implementation's frozen-window
spatial law and cut measure are explicit choices here, not attributed as an
existing MSRC derivation:

- [Efficient Bayesian inference under the multispecies coalescent with migration](https://pmc.ncbi.nlm.nih.gov/articles/PMC10622872/).
- [Efficient Bayesian inference under the structured coalescent](https://academic.oup.com/bioinformatics/article/30/16/2272/2748160).
- [The Structured Coalescent and Its Approximations](https://pmc.ncbi.nlm.nih.gov/articles/PMC5850743/).
- [SMC′ waiting distances](https://pmc.ncbi.nlm.nih.gov/articles/PMC8556732/).

The validated neutral engine and both ARG implementations are unchanged.
Results, failure cases, performance, reproduction and readiness are recorded
in `msrc_smc_structured_validation.md` rather than inferred from these formulas.
