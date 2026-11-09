# Neutral ARG simulator

The experimental `mode: arg` path implements a continuous-time backward
Hudson process over a species tree. Samples begin with one full-chromosome
segment. An active lineage stores sorted, non-overlapping half-open segments,
the current node ancestral to each segment, its species-tree population, and a
background state. The neutral state is deliberately a one-state implementation
of the background interface.

For a population branch with effective size `N`, every eligible pair in that
population has coalescence rate `1/(2N)`, including pairs with disjoint active
material. At a coalescence, only intervals carried by both lineages acquire a
new local parent node. Disjoint intervals retain their existing local roots
while the carrier lineages are merged. This is important: a common ancestor of
different genomic regions must not create a fictitious unary branch above a
completed marginal genealogy.

For a lineage with active material spanning `[L, R]`, the ancestry-changing
recombination hazard is

`integral_L^R r(x) dx`.

This includes gaps between disjoint active ancestral segments when material
exists on both sides of the gap, and excludes stored sequence outside the
leftmost/rightmost unresolved active bases. The exact same span is used for
the hazard and for map-weighted breakpoint sampling. A sampled cut is drawn
only from positive-rate pieces of this span and is accepted only when both
parental material sets are non-empty. Physical
crossovers that do not change the active ancestry are not inserted into the
event process; the event log distinguishes the ancestry-changing events that
are simulated from marginal-tree changes that happen to be observable.

After each coalescence, the engine counts active carriers over all segment
boundaries. An interval with one remaining carrier is marked resolved and is
removed from future hazards, but its nodes and edges remain in the output.
Lineages with no active material are retired. The process terminates when all
genomic intervals are resolved. Species-tree boundaries move active lineages
to their parent population; coalescence is never allowed across populations.

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

`msprime` is intentionally not a production dependency. The optional
`validation/arg/validate_neutral.py` program compares independent replicate
summaries against msprime for two samples, a four-population species tree,
heterogeneous recombination, and 4/8/16-taxon total branch lengths. The
diagnostic script samples fixed genomic positions and reports binned quartet
agreement, Pearson correlation of pairwise MRCA times, and Pearson correlation
of internal branch lengths. Positions from one chromosome are not treated as
independent replicate confidence intervals.

The implementation retains unary recombination-parent nodes in the full
tree-sequence representation when they are needed to preserve event ancestry;
the gene-tree exporter calls `simplify()` so ASTRAL receives ordinary marginal
trees. A tskit tree sequence represents marginal ancestry, not every event
annotation, which is why the separate event log is retained.

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
