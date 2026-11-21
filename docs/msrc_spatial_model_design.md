# Framework investigation for spatial MSRC

## Recommendation

Use a trajectory-conditioned, two-background ancestral-material process
(Hudson-style structured ARG), with a versioned **background-reservoir
inheritance rule**, as the reference spatial MSRC model. Build and validate
this reference before deriving any SMC approximation. Store local genealogies
as a tree sequence. This is an effective coalescent conditional on a supplied
frequency trajectory, not an exact finite-population reproductive pedigree.

The new prototype is `msrcsim.spatial_reservoir`. It reuses the existing ARG
scheduler and segment operations through a new opt-in background adapter;
existing simulators remain unchanged. The complete rates, proof of the
one-site generator, restrictions, and sampler are in
[the specification](msrc_spatial_recommended_model.md). The implementation and
validation stages are in [the plan](msrc_spatial_implementation_plan.md).

## Scientific requirements and the key distinction

The primary marginal target is the legacy MSRC conditional structured
coalescent, with a shared realized frequency history, explicit origin time,
same-background coalescence and marker-distance-dependent switching. The
spatial target is an explicitly stated coupling of these marginals through
shared ancestral chromosomes and crossover points. Preserving each marginal
does not uniquely identify that coupling or establish biological realism.

A local-tree process is Markov along genomic position only under additional
assumptions. A set of ancestral chromosomes with material intervals is Markov
in evolutionary time. We can use the second representation without deriving
the first. Keeping active ancestry is computationally more expensive, but
provides a generator whose projection to a focal locus is directly checkable.

The earlier audits overreached when they treated missing hidden-marker
conditional distributions as a barrier to any useful implementation. Their
examples demonstrate nonuniqueness of exact projected transitions; they do
not rule out choosing a complete, effective generative model. In particular,
the marker-associated target-count example assumes a particular hidden
coalescence construction without establishing that its two configurations
are in the same retained-state fiber. It is an illustration, not a proof of
impossibility. The recommendation here instead specifies an explicit
reservoir assumption and proves its one-site generator.

## Literature and source audit

Hudson and Kaplan (1988), *The Coalescent Process in Models With Selection and
Recombination*, develop genealogies of linked neutral loci conditional on
background structure at a selected site. This establishes the conceptual
foundation for joint ancestry over loci rather than independent recoloring.
[Paper](https://doi.org/10.1093/genetics/120.3.831).

Ewing and Hermisson (2010), *MSMS*, implement a forward frequency-trajectory
stage followed by a structured ARG conditional on that trajectory. Their
methods use inverse-frequency coalescence, competing event clocks, geographic
structure, and piecewise constant trajectories. This is the closest published
architectural precedent; our suppression and origin conventions remain
explicit model choices, and our adapter is not a port or independent
validation of MSMS. [Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC2916717/).

Steinrücken, Paul and Song (2013), *A sequentially Markov conditional
sampling distribution for structured populations with migration and
recombination*, derive a structured conditional sampling distribution using
an approximation to ancestral absorption histories. This shows that useful
structured sequential models exist, but their trunk/absorption approximation
is not a turnkey generator for arbitrary MSRC sample genealogies.
[Derivation](https://pmc.ncbi.nlm.nih.gov/articles/PMC3532580/).

McKenzie and Eaton (2026) analyze waiting distances under a multispecies SMC.
Our corrected neutral SMC′ is the relevant baseline. Its successful finite
neutral marginal checks do not imply that adding state-dependent switching
and cut hazards preserves a structured target.
[Paper](https://academic.oup.com/sysbio/article/75/2/366/8250006).
He, Scornavacca and Chan (2026) provide a fixed-separation conditional
guide-tree sampler; it is another coupling, with its own conditional
assumptions. [Paper](https://peercommunityjournal.org/articles/10.24072/pcjournal.694/).
The publisher pages for these two papers were inaccessible in this research
session (403/redirect error); their prior repository audits were read, and no
new full-method inspection is claimed here.

Deng, Song and Nielsen (2021), *The distribution of waiting distances in ancestral
recombination graphs*, and Wilton et al. (2015), *The SMC′ is a highly accurate
approximation to the ancestral recombination graph*, motivate explicit
no-change, tree-change and topology-change records and comparison of
multi-position laws rather than judging an SMC by its single-site agreement.
[Waiting distances](https://pmc.ncbi.nlm.nih.gov/articles/PMC8556732/),
[SMC′ accuracy](https://pmc.ncbi.nlm.nih.gov/articles/PMC4423375/).

A further caveat comes from Eriksson et al., *An Accurate Model for Genetic
Hitchhiking*: standard background coalescent rates can differ from finite
Moran genealogies near small background sizes. This is a reason to name our
model an effective conditional coalescent, and validate founder/boundary
behavior separately. [Analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC2206092/).

Baumdicker et al. (2022), *Efficient ancestry and mutation simulation with
msprime 1.0*, supports the choice of tree-sequence output and supplies neutral
comparison methods. Selected-history ancestry is not simply obtained by
setting two static population sizes in msprime.
[Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9176297/),
[simulation documentation](https://tskit.dev/msprime/docs/stable/ancestry.html).

## Candidate comparison

| Candidate | Generator/state and inheritance | Marginal target | Spatial properties and cost | Decision |
|---|---|---|---|---|
| Trajectory-conditioned reservoir ARG | Ancestral intervals, population and background per carrier; crossovers split material, marker-distal parent draws background from shared reservoir | Exact legacy one-site generator on admissible histories with matched rates and inactive floor; proof available | Shared carrier ancestry creates a consistent multi-locus process; active lineage/interval count can grow rapidly | Primary reference and working prototype |
| Explicit finite Moran/WF pedigree with meiosis and marker | Entire parental haplotypes or forward tree sequence; background inherited at actual marker | Generally differs from legacy continuous-time coalescent at finite copy number, selection and founder boundaries | Most mechanistic phase and linkage; cost scales with population size × generations plus recording | Later biological benchmark, especially at fixation/founding |
| Structured SMC′ with branchwise background histories | Guide genealogy plus migration/background paths; prune and re-coalesce in matching states | Requires detailed-balance/projection derivation for position-dependent switching and background-dependent cut rate; cannot assume preservation | Attractive eventual scaling; can distort branch times and marker associations | Derive only after reference is validated |
| Conditional sampling/trunk HMM | Fixed trunk plus absorption population, time and migration path; structured CSD transition | Targets conditional samples under its own approximation, not automatically joint MSRC genealogies | Efficient HMM/SMC structure; trunk and multi-sample order can bias genealogy law | Useful algorithmic precedent; not selected as simulator |
| Marginal-preserving refresh or MH coupling | At each coordinate target MSRC density; accept/reject a tree/path proposal | Possible to preserve a fixed target by detailed balance, or independent draws | Physical interpretation of distance and tract laws needs separate calibration; position-varying targets require transport/weights | Diagnostic coupling, insufficient as biological reference |
| Particle approximation of hidden marker ancestry | Particles over ancestral segment state conditioned on guide histories | Can converge to a specified reference if weights/proposals are correct | Particle and reference costs; degeneracy for long chromosomes | Consider only if a sequential implementation is later necessary |

The choice favors a provable marginal projection and explicit joint law over
an unmeasured expectation that SMC is faster. A complete exact genomic
conditional law can be defined by disintegration of the ARG distribution,
but that does not give an efficient sampler. We need not solve that
conditional-sampling problem to simulate chromosomes backward in time.

## Existing implementations: what actually differs

All three existing models use generations, diploid effective sizes, and
species boundaries. Legacy MSRC and current MSRC-ARG share the finite-copy
coalescence floor. Legacy uses a direct effective switch clock; the current
ARG uses a position-and-parental-state event measure. Neutral SMC′ omits
backgrounds and retains a dated guide tree.

The current MSRC-ARG extends its recombination domain to the marker only when
`_stored_anchor(segments)` is true. Once that condition fails, a focal lineage
can lose its marker-to-focal switch rate. For example, marker `z=0`, material
`[80,90)`, constant map `r`, and focal `u=85` have desired cross rate
`85*r*k*p_g`. The old model integrates only over `[80,90)`, allowing a switch
at the focal point only over `[80,85)`: rate `5*r*k*p_g`. This is a specific
generator discrepancy, not an expected consequence of ordinary ARG
recombination or finite chromosome length. The new adapter always retains
the marker-to-carrier exchange domain; it never assumes a stored marker is
required for the unobserved reservoir to exist.

At zero frequencies the old floor gives a positive rate for an impossible
occupied background. Missing origin also differs: legacy falls back to root
extension; current ARG remains structured indefinitely. The prototype requires
an explicit origin, uses exact positive-frequency rates, and raises on an
occupied absent background. It does not silently force transitions or reject
and redraw complete histories.

## What remains uncertain

The generator establishes the intended effective single-locus law; it does
not establish the correct finite-population selected-marker law. The
reservoir assumption discards dependence on unobserved marker ancestry and
parental phase. The spatial process is coherent by construction but its
biological tract lengths need a separate pedigree/MSMS comparison.
Frequency jumps into an absent occupied background need a founder/bridge
model. The explicit legacy “erase labels at origin” convention is usable,
but is distinct from forcing all copies of a unique newly arisen selected
allele through one founder. Both choices affect root ages and quartet laws.
These are named model choices to resolve, rather than missing SMC machinery.
