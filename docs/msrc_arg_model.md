# Two-background MSRC ancestral recombination graph

This document specifies the first structured-ARG extension. It is an explicit
effective model, not a complete inversion heterokaryotype model: it has two
background states, one background-defining marker, no gene conversion, and no
multihaplotype forward process.

## State

Each active ARG lineage has a species-tree population `b`, a background state
`h in {0,1}`, and the neutral ancestral segments already used by the Hudson
engine. One realized frequency history `p_b(t)` is shared by every segment in
the chromosome. Terminal states are sampled once from the history, or supplied
explicitly for reproducibility.

The rearrangement-origin age is read from the history. Younger times use the
structured process; at and older than that age, background labels are discarded
and the ordinary neutral coalescent is used. This boundary is not inferred from
the end of a numerical history unless the configuration explicitly makes them
the same.

## Coalescence

Within a structured population, two lineages can coalesce only when their
background states agree. For background 1 frequency `p` and background 0
frequency `1-p`, the finite-copy convention inherited from
`structured_coalescent.py` is

\[
c_{b,h}(t)=\frac{1}{2N_b\max(p_{b,h}(t),1/(2N_b))}.
\]

After the origin boundary, the rate is `1/(2N_b)` for every same-population
pair. Disjoint neutral segments remain eligible for common-ancestor events;
only overlapping material receives a new local ancestry node.

## Physical recombination and background inheritance

The baseline map is `r0(x)`. For a lineage in state `h`, a candidate parental
background `g` has local rate

\[
r_0(x)p_{b,g}(t)\eta(h,g,x),
\qquad
\eta(h,g,x)=\begin{cases}1&h=g\\\kappa_{cross}(x)&h\ne g.\end{cases}
\]

The engine integrates this same measure to select an event and samples both
position and `g` from it. Within-background events use the neutral Hudson
split. Cross-background events use the background-defining marker `z`: when a
cut is left of `z`, the left parent receives `g` and the right parent retains
`h`; when a cut is right of `z`, the left parent receives `h` and the right
parent receives `g`. This is a marker-inheritance convention, not a claim
about the detailed viability of inversion heterozygotes.

The eligible neutral crossover span is the existing Hudson span of unresolved
material. If the marker is stored but lies outside that span, the expanded
marker-to-active interval also permits a `background_transfer` event. Such an
event changes the background association of the active neutral lineage,
retains its neutral segment/node structure, and does not create artificial
genealogy edges. It is the explicit replacement for the old independent
switching event; no additional `kappa_s` process is added.

For a focal position `u`, separating it from marker `z` gives the target
single-locus rates

\[
q_{0\to1}(u,t)=p_b(t)\int_{\min(u,z)}^{\max(u,z)}r_0(x)\eta(0,1,x)dx,
\]
\[
q_{1\to0}(u,t)=[1-p_b(t)]\int_{\min(u,z)}^{\max(u,z)}r_0(x)\eta(1,0,x)dx.
\]

With symmetric suppression these reduce to the documented
`rho_eff(u) p` and `rho_eff(u)(1-p)` equations. The equality applies only when
the marker-transfer geometry, random parental pairing, and focal-locus
conditioning match these assumptions.

## Boundaries and limitations

Species boundaries and frequency-history boundaries are deterministic event
boundaries. At species boundaries, lineages move to the parent population and
retain their states. At frequency boundaries, hazards are recalculated. The
origin boundary removes the biological background partition. The model does
not implement multiple rearrangements, gene conversion, crossover viability,
phase, or a full forward haplotype process. These are explicit extension
points rather than hidden approximations.
