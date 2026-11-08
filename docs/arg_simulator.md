# Neutral ARG simulator

The experimental `mode: arg` path implements a continuous-time backward
Hudson process over a species tree. Samples begin with one full-chromosome
segment. An active lineage stores sorted, non-overlapping half-open segments,
the current node ancestral to each segment, its species-tree population, and a
background state. The neutral state is deliberately a one-state implementation
of the background interface.

For a population branch with effective size `N`, every pair in that population
has coalescence rate `1/(2N)`, including pairs whose material is disjoint. A
coalescence creates a parent node and edges for all material in both lineages.
Each lineage's recombination hazard is

`integral over ancestral material r(x) dx`.

A sampled crossover splits the material at its physical position into two
parental lineages. Species-tree boundaries move lineages to their parent
population; coalescence is never allowed across populations. The root branch
extension ensures that all material resolves.

`tskit` is used only after the ancestral process: node and edge tables are
constructed from the event history, sorted, and loaded as a tree sequence.
`arg_events.csv` retains species-boundary, coalescence, and physical
recombination events, including recombinations that do not change a marginal
tree. `local_gene_trees.nwk` and its CSV companion enumerate marginal tree
intervals.

## Scope and validation

Implemented: neutral multi-species ARGs, arbitrary ultrametric `SpeciesTree`
topologies, branch-specific effective sizes, constant and piecewise-constant
recombination maps, reproducible NumPy seeds, event logging, tskit export, and
local-tree iteration. The included unit tests check segment boundaries,
nonnegative rates, termination, recombination, and tskit loading.

This first version has not yet included a completed msprime statistical
benchmark suite; `msprime` is intentionally not a runtime dependency. The
diagnostic script reports adjacent topology agreement and shared pairwise MRCA
times while treating intervals from one chromosome as linked observations, not
independent replicates. A future benchmark should compare replicate summaries
and confidence intervals against msprime under matching species-tree and
recombination settings.

The implementation records the full ARG through unary recombination-parent
nodes and event metadata. A tskit tree sequence represents marginal ancestry,
not every event annotation, which is why the separate event log is retained.

## MSRC extension plan

The `backgrounds.BackgroundModel` interface separates species populations from
genetic backgrounds and from physical segments. A two-background model can
replace the neutral coalescence method with
`1/(2*N_b*p_bh(t))`, return background-specific states on recombination, and
add position-dependent transition rules. Such transitions should be attached
to physical crossover events (or to an explicitly defined background process),
not added as independent neutral switching events. Segment-associated state
objects allow multiple rearrangements without changing the ARG segment
machinery.

Known limitations are the lack of a background-dependent implementation, no
demography changes beyond the existing species-tree branch sizes, and no
large-scale performance optimization. The engine prioritizes transparent event
rates and valid ancestry tables over matching msprime's optimized runtime.
