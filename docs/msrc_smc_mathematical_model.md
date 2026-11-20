# Two-background MSRC-SMC′: Phase A mathematical audit

**Decision: do not implement the requested transition kernel yet.** The
current local-tree SMC′ state, even when augmented with one background label
per branch, is insufficient to determine a physically coupled transition
along the chromosome. The single-locus switching intensity does not uniquely
specify the joint ancestry of the focal locus and the background-defining
marker, or the conditional recombination/re-coalescence distribution between
neighboring positions. A kernel built from the supplied marginal rates alone
would be an additional approximation, not a derived MSRC-SMC′ process.

This document records the code audit, the equations that can be derived, a
counterexample to the proposed state, and the information needed to resume
design. No simulator code or production path is changed in this phase.

## 1. What the implementations share and what differs

| Component | Legacy single-locus MSRC | Current MSRC-ARG | Neutral SMC′ |
|---|---|---|---|
| Time and pair rate | Generations; species branch `Ne`; pair rate `1/(2Ne)` in the ordinary ancestor | Same | Same |
| Structured pair rate | Same-background only; `1/(2Ne * max(p_h,1/(2Ne)))` | Same function and floor in `MSRCTwoBackground` | No background state; ordinary MSC rate |
| Background changes | Independent effective backward clock: `q01=rho*p`, `q10=rho*(1-p)` | Parental background sampled at physical recombination points between active material and anchor | Not represented |
| Physical map | Effective scalar `rho = recombination_rate * effective_fraction` | Piecewise map integrated over the lineage's active span and, when the marker is retained, marker-to-active span | Piecewise map sets genomic breakpoint hazard `r(x)*L(G)` |
| State | One arrangement label per lineage | Segment ancestry plus one background label per active lineage and a marker coordinate | One dated local tree |
| Origin boundary | Uses recorded origin; if absent, compatibility fallback `root_age + root_extension` | `structured_origin_age=None` means `_structured()` remains true at all times | Not applicable |

The legacy implementation is an effective structured coalescent, not a
physical recombination simulator: `simulate_genealogy_reference` applies the
switch clock directly. The ARG implementation attempts to produce that clock
from marker-mediated crossover events. The neutral SMC′ uses a guide tree and
does not retain the ancestral segment graph. Thus the ARG is not automatically
ground truth for the legacy model, and neither is an automatic oracle for a
new SMC′ projection.

The Moran implementation uses `M=2Ne`, state-changing count rates
`q_plus=(M/2)*p_sel*(1-p)` and
`q_minus=(M/2)*(1-p_sel)*p`, with
`p_sel=p(1+s)/(1+sp)`. It draws branch histories and initializes daughters
from the parent terminal frequency. Wright–Fisher and diffusion histories are
also available through `simulate_population_history`. These histories can be
generated once and shared across loci, but the simulator must retain their
population-specific event boundaries and split initialization when evaluating
backward rates.

### 1.1 Does the current ARG produce the stated one-locus switching rate?

For a lineage in state `h`, `MSRCTwoBackground.recombination_hazard` sums the
two parental-state measures over its active-material hull, extended to the
stored anchor when applicable. At a position where the crossover changes the
focal locus's marker association, the cross-parent component is
`r0(x) * p_g(t) * eta(h,g,x)`. Integrating only crossover points separating
focal coordinate `u` from anchor `z` gives

```
q_h_to_g(u,t) = p_g(t) * integral_[min(u,z), max(u,z)]
                         r0(x) * eta(h,g,x) dx,  g != h.
```

Therefore, **the intended intensity is produced under the implementation's
marker-present, one-active-focal-block assumptions**: a shared rate map, the
configured anchor, and a lineage segment state spanning the anchor and focal
ancestry. It is an effective event intensity, not by itself a complete
description of the resulting joint genealogy. At positions inside active
material the code instead splits the segment and assigns parental backgrounds
to the two sides. Events outside active neutral material are recorded as
`background_transfer`; the whole lineage's background label is updated, and
the segment tree is not split.

This distinction matters: the same integrated intensity can be obtained by
different hidden marker genealogies and phase histories. The current event
rules do not retain an explicit marker genealogy or a full segment-by-segment
ancestral phase state after every background transfer. Consequently the
intensity audit is not a proof that conditional local-tree transitions, or
joint multi-position marginals, match the intended recombination process.

### 1.2 Existing boundary conventions that need a decision

The inherited finite-copy floor gives a positive rate even for an explicitly
supplied lineage in an absent background. For example, at `p=0`, a state-1
pair receives rate
`1/(2Ne * max(0,1/(2Ne))) = 1` per generation. The stated biological rule
“lineages in different backgrounds cannot coalesce” is respected, but a
lineage in background 1 can coalesce even though that background has zero
frequency. This is a regularization convention, not the literal limit of the
boxed rate. It must be retained as a named approximation or replaced with a
boundary rule before scientific claims. A Bernoulli terminal draw avoids
sampling the absent state at fixation/loss, but explicit/conditioned states
and ancestral transitions can still encounter this case.

There is also a concrete origin default mismatch: the legacy simulator uses
`root.age + root_extension` when no origin is attached to the history, whereas
`MSRCTwoBackground._structured()` treats `None` as structured at every time.
Matched marginal validation must pass an explicit shared origin or reconcile
these defaults. At speciation the existing history samplers initialize
daughter populations from their parent frequency (binomial for finite-copy
Moran/WF variants, with an optional split policy); the backward lineage
population changes at the species boundary. A linked simulator must use the
same realized branch-specific history and boundary values in all reference
methods.

## 2. The desired single-locus marginal process

Within population `B`, the requested coalescence rate for a compatible pair
in background `h` is

```
c_(B,h)(t) = 1 / (2 Ne_B(t) p_(B,h)(t)),
p_(B,1)=p_B,  p_(B,0)=1-p_B.
```

The legacy implementation uses the finite-copy floor described above. At the
origin age it drops the background partition and uses the ordinary MSC. The
legacy backward switching generator is

```
Q_B(t) = [[-rho*p, rho*p], [rho*(1-p), -rho*(1-p)]].
```

For an unstructured genetic exchange process with a marker at `z`, the
candidate physical-recombination construction for a focal locus `u` gives

```
q_(0->1)(u,t) = p_B(t) * integral_between(u,z) r0(x) eta(0,1,x) dx,
q_(1->0)(u,t) = (1-p_B(t)) * integral_between(u,z) r0(x) eta(1,0,x) dx.
```

When the integral is the same scalar `rho`, this agrees with the legacy
switching rates. When the map, suppression, or marker distance varies, each
focal coordinate has its own effective rate. This agreement concerns one
focal lineage's switching intensity only; it does not establish equality of
genealogy distributions or of a chromosome-wide joint process.

At `p=1/2`, the switching generator has stationary distribution `(1/2,1/2)`.
At variable `p(t)`, the realized history is an exogenous shared path in this
task, so conditional genealogy draws must not independently regenerate it.
At fixation/loss, the exact two-background coalescent requires a boundary
convention for absent states; the current floor is not a neutral mathematical
consequence of the displayed rate.

## 3. Why the proposed local-tree state is insufficient

Let the proposed SMC state at coordinate `x` be

```
S_x = (G_x, x, species-population labels, branch background labels, p_B(t), z).
```

It does not contain which ancestral lineage at the anchor `z` is the marker's
ancestor, the ancestral segment/phase links between `x` and `z`, or the joint
ancestral recombination history that created those links. Those variables
determine which background parent a physical crossover assigns, what material
is on each side, and which guide branches are available to the detached focal
lineage.

### Counterexample: identical local state, different next transition law

Take two focal positions `u` and `v` on the same side of `z`, with identical
maps and the same realized constant `p=1/2`. Consider two valid ancestral
configurations that project to the same current focal local tree `G_u`, the
same focal branch-background labels, the same species populations, and the
same `x=u`:

1. In configuration A, the focal ancestry and anchor ancestry still occupy
   one ancestral lineage over the interval separating them.
2. In configuration B, an earlier recombination has separated the anchor
   ancestry; the focal local tree is unchanged, but the anchor is on a
   different ancestral lineage and has a different re-coalescence history.

Both configurations have the same proposed `S_u`, and the supplied one-locus
rate formula gives the same instantaneous `q_(h->g)` for both. But following
a crossover, the chance that the focal lineage re-coalesces with the marker
ancestry before the next species/frequency boundary differs: configuration A
has the anchor ancestry on the shared lineage, while configuration B has a
separate marker lineage with its own coalescence opportunities and background
history. The resulting probability that `G_v` shares a topology/time with
`G_u` therefore differs. A Markov transition kernel conditioned only on
`S_u` would have to assign one law to both cases, so it cannot reproduce both
ancestral configurations. This is a direct failure of the Markov sufficiency
criterion, independent of the number of simulations or the choice of
re-coalescence rate formula.

Equivalently, single-site transition rates do not identify the joint process:
many couplings of the per-site generators have the same one-site marginals
but different cross-site genealogical autocorrelation and tract-length
distributions. The requested acceptance criterion requires both, so selecting
one coupling without deriving it from the segment/marker ancestry would be
ad hoc.

## 4. State required for a well-defined extension

An exact Markov ancestral state sufficient for this model must at least retain
the active ancestral recombination graph over the focal interval and marker:

- all ancestral lineages and species-tree populations;
- each lineage's ancestral material intervals, including the marker;
- background/arrangement state for each marker-linked ancestral haplotype;
- the linkage/phase relation between the marker and each active interval;
- the shared realized frequency history and its population-specific event
  boundaries;
- recombination and coalescence nodes needed to determine future carriers.

For a general map and many scanned positions, the state must retain enough
segment ancestry to update all focal/marker associations. This approaches the
full MSRC-ARG state. A reduced SMC′ state may be possible, but it requires a
specified projection and a derived conditional distribution over omitted
marker ancestry given the retained state. The current papers and repository
implementation do not provide that conditional kernel. Keeping only a marker
gene tree is not automatically sufficient if the focal and marker ancestry
share hidden recombination/phase links; this needs a separate proof.

The sequential scan would need to condition on that augmented state, sample
physical recombination locations from the joint rate measure, update the
ancestral material/background assignment on both sides, and then apply
species-compatible, background-compatible re-coalescence over every change
boundary. SMC′ can omit some ancestral graph detail only after its conditional
transition is explicitly derived. Treating an arbitrary SPR or coloring a
local tree is not such a derivation.

## 5. Relationship among the three models

- **Legacy single-locus MSRC:** conditional structured-coalescent model with
  direct switching rates `rho p` and `rho(1-p)`, plus background-dependent
  coalescence. Its switching clock is effective and does not encode the
  physical breakpoint shared by neighboring loci.
- **Current MSRC-ARG:** adds physical breakpoint positions, segment ancestry,
  and a marker-mediated parental-background rule. Under the matching anchor,
  map, history, and origin assumptions, its single-focal switch intensity
  reduces to the integrated formula above. It is an explicit effective ARG
  model, but its background/phase abstractions still require direct marginal
  validation and should not automatically be treated as ground truth.
- **Neutral SMC′:** retains one guide genealogy and uses cut/re-coalescence
  transitions with population constraints. It has no background or marker
  state. Its successful neutral marginal validation says nothing by itself
  about a structured transition.
- **Proposed MSRC-SMC′:** not defined yet. A valid version must specify what
  ancestral marker/phase state is retained or the conditional law used to
  integrate it out, and show that its one-site marginal agrees with the
  explicitly selected target model before making spatial claims.

## 6. Phase A exit criteria for future implementation

Before coding a background SMC′ kernel, settle these mathematical choices:

1. Is the target one-locus process the legacy effective switch model, the
   current marker-mediated ARG, or a new forward haplotype model? State which
   process defines the marginal truth.
2. Specify the marker inheritance and parental haplotype model, including
   whether recombination between marker and focal material can change phase
   without changing focal genealogy.
3. Prove a sufficient reduced state (or retain the full marked segment ARG)
   and derive its conditional transition law. Show that the law is consistent
   across two or more focal positions, not just one-site switching rates.
4. Resolve the absent-background coalescence convention at `p=0,1` and align
   origin-age defaults across implementations.
5. Establish matched single-locus reference settings: same frequency path,
   species split history, sampled background states, map-integrated `rho`,
   origin, and coalescence boundary rule. Only then compare quartet
   probabilities, `D=q2-q1`, and coalescence-time distributions.

Until these are resolved, a two-background local-tree SMC′ implementation
would be under-specified. No implementation, spatial-correlation claim,
performance benchmark, or MSRC-SMC′ validation result is reported here.
