# Neutral multispecies SMC′ simulator

`msrcsim.smc` is an isolated linked-locus simulator. It leaves the independent
MSC simulator and `msrcsim.arg` (the full ancestral recombination graph
benchmark) intact.

## Model and conventions

Species-tree branch lengths, genealogy node ages, and population sizes use
generations and diploid effective sizes. Each compatible pair coalesces at
`1/(2 Ne)` per generation. The first local tree is drawn by the existing
ordinary MSC simulator. The current tree is represented as a dated rooted
Newick tree; tips are one sampled lineage per species at time zero.

At genomic position `x`, each edge of the local genealogy contributes its
length in generations to `L_G`. The next breakpoint is drawn by inverting the
integrated map hazard for an Exp(1) variate, equivalent to
`Delta_x ~ Exp(r * L_G)` for a constant rate. Piecewise map boundaries are
integrated exactly, including zero-rate intervals. At the breakpoint, an edge
is selected proportional to its eligible length, and a time uniformly on that
edge is selected as the cut.

The cut subtree is detached. Its lineage is advanced backward through species
tree epochs and re-coalesces with guide-tree lineages only when both occupy the
same species-tree population. Pair hazards are `1/(2 Ne)` per compatible guide
lineage. The detached edge remains an eligible target, implementing the SMC′
invisible recombination outcome. Re-coalescence to that edge leaves the local
tree unchanged. Other outcomes are classified as no tree change, branch-time
change, or topology change. A recombination event is recorded even when its
tree is unchanged; adjacent equal trees are retained as distinct genomic
intervals.

This is an SMC′ approximation, not the exact ARG. The ARG permits a detached
lineage to coalesce with any ancestral lineage carrying the segment; SMC′
restricts its targets to lineages in the current marginal genealogy, with the
same-edge event retained. McKenzie & Eaton (2026) derive event probabilities
for precisely this SMC′ process embedded in the MSC and account for barriers
between populations. He, Scornavacca & Chan (2026) instead describe a
conditional guide-tree sampler at fixed genomic distances in which uncoalesced
locus lineages may coalesce with each other. That algorithm motivates linked
gene-tree simulation but is not the variable-distance SMC′ transition used
here.

## Recombination map

Map rates are per base per generation. `breakpoints` partitions the
half-open chromosome `[0,L)` and `rates` supplies one rate per interval. A
region `[a,b)` with reduced or zero rate is a region of suppressed
recombination. It changes spatial persistence through the breakpoint hazard;
it does not by itself change the marginal quartet distribution.

## Extension interface and limits

The package boundary is intentionally separate from tree representation and
recombination-map handling so a later augmented lineage state can hold
background histories. The intended extension points are `p(t)`, a two-state
background, background-specific coalescence, backward background switching,
selection in the frequency process, and shared demographic trajectories.
None is active in this neutral implementation. Genomic recombination is a
spatial event; background switching is an evolutionary-time event. No
independent switching clock is introduced here. A future MSRC transition
must derive both events from a joint ancestral process, not equate them by
name or add an uncoupled rate.

Known limitations: one sample per species, ultrametric species trees, no
migration or serial sampling, no recombination on terminal sequence outside
the sampled marginal genealogy, and no full-MSRC integration. The initial MSC
tree is produced by the project's existing implementation and serialized to
Newick before transitions. The edge-based transition is experimental and its
broader validation **fails**: several higher-recombination cases differ from
independent MSC marginal root-age and other genealogy summaries, including
after multiple-test correction. Do not treat the transition as a validated
MSC-SMC′ simulator or extend it to MSRC-SMC until this is resolved. See
[`msrc_smc_validation.md`](msrc_smc_validation.md) for the full audit, tests,
and reproduction commands.

## CLI and output

Run `msrc-sim-smc --config config.yaml --output smc_output`. Inputs use the
normal project `species_tree` block and a `smc` block with chromosome length,
seed, and map breakpoints/rates. Outputs are `local_genealogies.csv`,
`local_gene_trees.nwk`, `recombination_events.csv`, `genealogy_changes.csv`,
and `metadata.json`.

## References and code

- He, Scornavacca & Chan (2026), [The effect of gene tree dependence on summary methods for species tree inference](https://peercommunityjournal.org/articles/10.24072/pcjournal.694/); simulation code: [ASTI-GTD](https://github.com/WAHE1/ASTI-GTD/tree/main/code).
- McKenzie & Eaton (2026), [Estimating Waiting Distances between Genealogy Changes under a Multi-Species Extension of the Sequentially Markov Coalescent](https://academic.oup.com/sysbio/article/75/2/366/8250006).
- Their published calculations are implemented in [ipcoal](https://github.com/eaton-lab/ipcoal); reproducible validation notebooks are in [waiting-distance-code](https://github.com/eaton-lab/waiting-distance-code). These implement event probability and waiting-distance calculations, rather than a standalone chromosome simulator.
- He et al.'s code release is [ASTI-GTD v1.0.0](https://github.com/WAHE1/ASTI-GTD/tree/v1.0.0), archived at [Zenodo](https://zenodo.org/records/18930409). The shell could not resolve GitHub to inspect a local checkout during this implementation.
