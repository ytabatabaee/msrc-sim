# Project Status

The simulator now has two supported modes:

- Legacy quartet mode: the frozen four-taxon MSRC validation path and existing
  command-line tools remain available for reproducible manuscript experiments.
- Arbitrary-tree mode: `msrc-sim-arbitrary` runs Wright-Fisher or Moran
  structured-coalescent core over an arbitrary rooted species tree and writes
  n-tip gene trees plus optional quartet summaries.

Current local extension:

- Completed an optional speciation-associated daughter-split kernel with
  reciprocal log-odds parameter `eta`, shared by Wright-Fisher and Moran.
- Added arbitrary-tree, spatial, frozen-history, provenance, tests, docs, and
  one example config for the opt-in daughter-split feature.
- Scientific default remains legacy unbiased initialization; configs without
  `speciation_split` and explicit `eta=0` retain the historical simulator.
- New future scientific axis: eta != 0 daughter assortment at species splits.

Current validation:

- Full test suite: `163 passed, 6 warnings` on 2026-09-16 after the corrected
  Moran selection-scaling rerun.
- Legacy regression manifest:
  `tests/regression/legacy_quartet_baseline.json`.
- Arbitrary examples:
  `examples/output/10taxa` and `examples/output/50taxa`.
- Validation report:
  `results/arbitrary_tree_validation.md`.
- Moran validation report:
  `results/moran_validation.md`.

Known limitations before empirical use:

- The 10- and 50-taxon examples are implementation checks, not biological
  validation studies.
- Wright-Fisher and Moran are alternative stochastic generators of frequency
  histories; process-dependent quantities such as persistence and terminal
  pattern probabilities are not claimed to be invariant.
- Quartet-specific analysis, exact pattern probabilities, and model-comparison
  utilities remain quartet layers by design.
- Large stochastic validation and performance tuning should be run separately
  for empirical-scale analyses.

## Prevalence Parameter-Sensitivity Grid

New forward-process prevalence study added under
`results/prevalence_grid/` and `figures/prevalence_grid_*`. This is a
parameter-sensitivity experiment only; it does not modify empirical analyses,
quartet-distortion experiments, threshold experiments, pseudo-replication
experiments, or correction experiments.

Rationale:

- The biological/theoretical question is whether a structural polymorphism can
  persist long enough to generate phylogenetically relevant discordant 2:2
  sorting outside one tuned parameter setting.
- The study separates root-branch population-genetic persistence from terminal
  2:2 sorting and emphasizes unconditional discordant 2:2 probability as the
  more relevant compound outcome.

Parameter grid:

- Primary Wright-Fisher grid: Ne = 20, 50, 100; p0 = 0.05, 0.10, 0.25, 0.50;
  scaled branch duration tau = T/(2Ne) = 0.10, 0.30, 0.60, 1.00; scaled
  genic directional selection gamma = 2Ne s = -5, -2, 0, +2, +5; 5,000
  trajectories per cell.
- Corrected matched Moran robustness grid: Ne = 20, 100; p0 = 0.10, 0.50;
  tau = 0.30, 0.60, 1.00; shared diffusion-scale `gamma_target` = -5, -2,
  0, +2, +5; 5,000 trajectories per cell.
- The old Moran robustness grid is preserved as
  `moran_matched_grid_old_unscaled_selection.tsv`. It used the same nominal
  `gamma = 2Ne s` as Wright-Fisher. Those raw Moran simulations are valid for
  their stated raw `s`, but they under-matched the WF diffusion-scale
  directional-selection strength by approximately a factor of two.
- Supplementary raw-s grid: raw s = -0.01, 0, +0.01 across Ne, labeled
  `supplement_raw_s`, because fixed raw s changes gamma with Ne.

Simulator conventions:

- Both processes use M = 2Ne chromosome copies.
- Wright-Fisher uses Binomial(M, p_sel) updates with
  p_sel = p(1+s)/(1+s p).
- Moran uses the same genic reproductive probability with CTMC rates
  q+ = (M/2) p_sel(1-p), q- = (M/2)(1-p_sel)p.
- The implemented Moran drift is
  dE[p]/dt = (q+ - q-)/M = s p(1-p)/[2(1+s p)], which is approximately
  (s/2)p(1-p) under weak selection. Therefore corrected WF/Moran comparisons
  use `s_WF = gamma_target/(2Ne)` and `s_Moran = gamma_target/Ne`.
- This is a diffusion-scale matching rule, not an exact finite-s equality of
  the WF and Moran generators.
- Daughter branches are initialized by binomial sampling from parental
  terminal frequency, matching the existing prevalence machinery.

Current focal setting:

- Legacy manuscript-facing WF/Moran prevalence setting: Ne = 20, p0 = 0.5,
  internal tau = 0.6, root extension = 0.5, s = 0, 2,500 attempts.
- Legacy prevalence counted persistence as both daughter branches A and B
  segregating. The new grid's primary P_persist is root-branch terminal
  0 < p(T) < 1, so `focal_reproduction.tsv` reports the legacy comparison
  separately.

Results:

- New focal primary Wright-Fisher cell, Ne = 20, p0 = 0.5, tau = 0.6,
  gamma = 0: P_persist = 0.769 and unconditional P_discordant_2to2 = 0.032.
- Largest primary-grid unconditional P_discordant_2to2: 0.171 at Ne = 50,
  p0 = 0.5, tau = 0.1, gamma = 0.
- Weakest primary-grid discordant outcome: 0.000 in several low-frequency,
  long/selected cells, including Ne = 20, p0 = 0.05, tau = 0.6, gamma = -5.
- At Ne = 50 and neutrality, P_discordant_2to2 declines from 0.171 to 0.009
  as tau increases from 0.1 to 1.0 for p0 = 0.5; lower p0 values are much
  less productive.
- Neutral P_persist curves approximately collapse across Ne after scaling by
  2Ne; the largest p0 = 0.5 within-tau spread across Ne is 0.033.
- Moran and Wright-Fisher are close under neutrality for matched cells, but
  strong directional selection required process-specific raw-s matching.
  After correction, selected-cell mean absolute WF/Moran differences dropped
  from 0.151 to 0.023 for P_persist and from 0.011 to 0.001 for unconditional
  P_discordant_2to2 in the cells available for old-vs-corrected comparison.
  Remaining differences are process dependence, not a known parameterization
  bug.
- Exact small-population validation of loss/persistence/fixation agrees with
  Monte Carlo within ordinary sampling error; corrected selected/neutral Moran
  checks have maximum absolute error 0.0049 and all exact P_loss/P_persist/P_fix
  values inside Wilson intervals.

Interpretation and manuscript implications:

- The mechanism is not confined to the single legacy focal setting, but
  phylogenetically relevant discordant 2:2 sorting is much narrower than mere
  polymorphism persistence.
- Intermediate/high initial frequencies and short scaled branch durations
  produce the strongest discordant 2:2 prevalence.
- Moderate positive selection can rescue low-frequency A1 from immediate loss
  at short times, but long/strong directional selection increases absorption
  into loss or fixation.
- Do not interpret Ne = 20, 50, 100 as biologically realistic population
  sizes. They are computational scaling experiments. Organism-specific
  prevalence requires calibrated Ne, selection, and divergence-time estimates.
- The LaTeX manuscript has not been modified for this correction yet.
- Empirical analyses are unchanged.
- The completed full Wright-Fisher grid is unchanged and was not rerun.

Outputs:

- `results/prevalence_grid/prevalence_grid.tsv`
- `results/prevalence_grid/prevalence_summary.md`
- `results/prevalence_grid/selection_scaling_derivation.md`
- `results/prevalence_grid/moran_matched_grid_old_unscaled_selection.tsv`
- `results/prevalence_grid/corrected_moran_matched_grid.tsv`
- `results/prevalence_grid/wf_moran_corrected_comparison.tsv`
- `results/prevalence_grid/corrected_moran_exact_validation.tsv`
- `results/prevalence_grid/exact_validation.tsv`
- `results/prevalence_grid/focal_reproduction.tsv`
- `results/prevalence_grid/calibration_template.tsv`
- `figures/prevalence_grid_main.pdf`
- `figures/prevalence_grid_supplement.pdf`
- `figures/prevalence_selection_main.pdf`
- `figures/prevalence_selection_supplement.pdf`

Next step:

- Inspect the selection-focused figures, then update the manuscript text to
  quote root-terminal persistence and legacy daughter-persistence with clearly
  separated definitions.

## Daughter-Sorting Eta Sensitivity

Completed:

- The optional speciation daughter-sorting kernel was already implemented and
  validated before this experiment.
- Added the first focused eta sensitivity experiment under
  `experiments/daughter_sorting_eta/` and `results/daughter_sorting_eta/`.
- Froze the baseline prospectively: symmetric `((1,2)A,(3,4)B)ROOT`
  Wright-Fisher setting with Ne = 20, p0 = 0.5, internal tau = 0.6,
  root extension = 0.5, neutral selection, and descendant-split sorting only.
- Main run used 5,000 forward histories per eta for eta = 0, 0.5, 1, 2, 4.
- Sign-symmetry check used eta = +2 and eta = -2.
- Secondary genealogy check used up to 50 persistent discordant histories per
  main eta and 300 gene genealogies per retained history.
- Optional controls completed: root-only sorting with 2,000 histories per eta
  for eta = 0, 1, 2, 4; matched neutral Moran robustness with 2,000 histories
  per eta for eta = 0, 2, 4.

Main results:

- P(persistence) stayed near 0.30 across descendant eta values:
  0.2980, 0.2946, 0.3010, 0.3008, 0.2950.
- P(discordant strict 2:2) increased with eta:
  0.0352, 0.0400, 0.0602, 0.1238, 0.2442.
- P(discordant strict 2:2 | persistence) increased strongly:
  0.1181, 0.1358, 0.2000, 0.4116, 0.8278.
- D without P count was zero in all main and sign-symmetry cells.
- Mean realized daughter-frequency divergence increased from 0.0710 at eta 0
  to 0.9041 at eta 4.
- P(opposite enrichment) increased from 0.1752 at eta 0 to 0.9991 at eta 4.
- Sign symmetry behaved as expected: eta +2 favored pattern 1010
  (`P=0.1192`) while eta -2 favored 0101 (`P=0.1130`); aggregate discordant
  strict 2:2 was 0.1238 versus 0.1178.
- In the genealogy check, the arrangement-associated topology was dominant in
  all retained persistent discordant histories. Mean arrangement-minus-species
  quartet margin increased from 0.7385 at eta 0 to 0.8345 at eta 4.
- A later review found that this first daughter-sorting genealogy check used
  the forward-experiment cross-arrangement fraction 0.1, not the historical
  WF/Moran robustness genealogy setting 0.002. Those first margins therefore
  are useful as an internal check but are not directly comparable to the
  historical approximately 0.994 robustness margins.
- Revised matched-genealogy outputs under `results/daughter_sorting_eta/revised/`
  use the exact historical backward settings recovered from
  `experiments/wf_moran_robustness/config.yaml` and
  `results/wf_moran_msrc_robustness/run_metadata.json`:
  recombination rate 0.02, effective cross-arrangement fraction 0.002,
  50 histories, and 300 genealogies per history. The eta=0 reproduction gave
  mean margin 0.993533, matching the historical Wright-Fisher target.
- The revised matched-genealogy margins for eta 0, 0.5, 1, 2, and 4 were
  0.9930, 0.9940, 0.9953, 0.9945, and 0.9962, with the arrangement-associated
  topology dominant in 50/50 histories for every eta.
- Root-only sorting mainly reduced descendant persistence in this baseline
  rather than producing the same descendant-discordant response: P(persistence)
  fell from 0.3055 at eta 0 to 0.0025 at eta 4.
- Matched neutral Moran robustness was qualitatively similar for eta:
  P(discordant | persistence) was 0.1190, 0.3632, and 0.8127 for eta 0, 2,
  and 4.

Decisions:

- The baseline was frozen before viewing eta results.
- Descendant-split sorting is the primary analysis; root-only sorting is a
  location-control result.
- Large crossed eta x Ne/selection/kappa/branch-duration grids are deferred.
- Previous eta=0 manuscript results remain unchanged.
- The structural-theorem check still held exactly in the preserved and revised
  outputs: no discordant strict 2:2 replicate occurred without persistence.

Next steps:

- Decide whether the descendant-split eta result belongs in the main MSRC
  manuscript figure/table or in a focused supplement.
- If included, use the revised matched-genealogy figures/tables rather than the
  original genealogy check for comparison to historical robustness margins.
- If included, run a later crossed eta x persistence-parameter grid to show
  robustness across Ne, p0, branch duration, selection, and kappa.
# Phase: neutral multispecies SMC′ prototype (2026-11-16)

Completed:

- Added isolated `msrcsim.smc` implementation, `msrc-sim-smc` CLI, piecewise
  recombination maps, local-tree/event/metadata exports, four-taxon config,
  demo generator, and `docs/msrc_smc.md`.
- Initial genealogy uses the existing MSC simulator. Recombination positions
  invert the integrated `r(x) * L_G` hazard, including rate-map boundaries.
- Transitions cut a length-weighted genealogy edge, sample a uniform time on
  the edge, and re-coalesce against guide lineages only in a shared species
  population. The detached edge and ancestral guide-tree stem remain
  re-coalescence targets, retaining SMC′ no-change outcomes.
- Classified recombinations as no genealogy change, branch-time change, or
  topology change. Consecutive equal trees are preserved as separate spans.
- Added explicit future-state documentation for p(t), two backgrounds,
  background-specific coalescence, background switching, selection, and shared
  demographic histories. No MSRC background switching is implemented.
- `msrc-arg` and production simulation paths were not modified.

Validation run:

- Full repository suite `.venv/bin/python -m pytest -q -p no:cacheprovider`:
  passed. One existing warning reports a regression threshold outside the
  sampled validation grid.
- `.venv/bin/python -m pytest -q -p no:cacheprovider tests/test_msrc_smc.py`:
  11 passed. Coverage includes zero rate, seed reproducibility, constant-map
  equivalence, map boundary integration, suppressed intervals, turnover and
  topology autocorrelation under increased recombination, and tree validity.
- The marginal-distribution check compares genomic-span-weighted SMC′ local
  trees to 600 independent MSC draws on a four-taxon tree. Topology frequencies
  and mean root age pass the stated Monte Carlo tolerance. The weighting is by
  interval width, corresponding to sampling a random chromosome position.
- Demo command `MPLCONFIGDIR=/private/tmp/msrc-mpl .venv/bin/python
  scripts/demo_msrc_smc.py` completed for constant rates 5e-9, 1e-8, 2e-8
  and a map suppressed to 1e-9 over [50,000, 150,000). Recombination counts
  were 90, 206, 444, and 148; tree-change counts were 52, 132, 264, and 80.
  Outputs are in `smc_demo/`; figure: `smc_demo/smc_map_comparison.png` and PDF.
- Attempts to retrieve the associated GitHub source with `git ls-remote`
  failed because shell network DNS was unavailable. The paper and code release
  metadata were reviewed through the paper's linked GitHub/Zenodo references;
  source code itself was not locally audited.

Mathematical decisions and limitations:

- The selected process is the multispecies SMC′ approximation described and
  analyzed by McKenzie & Eaton (2026), not exact ARG. He et al. (2026) use a
  fixed-separation conditional guide-tree process which permits the detached
  locus lineages to coalesce with each other; that is not the transition
  implemented here.
- Times are generations and all compatible pair rates are `1/(2Ne)`.
  Recombination map rates are per base per generation; the spatial hazard is
  `r(x) * total_genealogy_branch_length`.
- The SMC′ marginal test is a finite seeded check, not a proof. The transition
  is experimental pending broader parameterized and ARG comparisons.
- McKenzie & Eaton's code was reviewed through the public ipcoal and
  waiting-distance-code repository pages: these contain analytical event
  probability/waiting-distance functions and reproducible notebooks, rather
  than a direct chromosome transition simulator. He et al.'s code release was
  identified as ASTI-GTD v1.0.0 and its Zenodo archive, but could not be cloned
  for local source inspection because shell DNS could not resolve GitHub.
- Current sampling supports one contemporaneous lineage per species, with an
  ultrametric rooted species tree. Migration, serial samples, full ARG, and
  all MSRC background behavior are out of scope.
- No claim is made that lowering recombination changes marginal quartet
  frequencies; it changes spatial dependence and tract lengths.

Phase 2:

1. Audit and formalize the SMC′ transition against the McKenzie–Eaton
   supplementary derivations and an independent implementation (msprime/full
   ARG in ultrametric neutral cases).
2. Expand the MSC marginal checks over topology, internal coalescence times,
   species-tree shapes, and branch-specific Ne.
3. Design the augmented lineage state and joint backward process for genetic
   backgrounds, p(t), selection, switching, and genomic recombination before
   implementing any switching clock.
4. Only then add an MSRC-SMC mode and compare it with `msrc-arg`; retain
   `msrc-arg` as the validation benchmark.

## Neutral SMC′ Validation Follow-up (2026-11-19)

Corrected and completed:

- Fixed root-age bias in `msrcsim/smc/engine.py`: the transition interval
  partition now includes the original cut edge's parent time, when its
  re-coalescence target expires. A stale target set had allowed impossible
  invisible events past that time.
- Added the two-tip analytic target-expiry and conditional re-coalescence-time
  tests, including 12,000 fixed-guide transitions; added four- and eight-taxon
  high-recombination marginal comparisons and a regression test for absent
  quartet topologies being counted as zero.
- Fixed the validation summarizer's absent-quartet handling. Previously, it
  omitted chromosomes with zero count for a topology and biased its frequency
  summaries upward. All corrected results below use the corrected summary.
- Compared against independent MSC draws, msprime SMC′, msprime Hudson ARG,
  and `msrc-arg` over seven settings. Added integrated-hazard, spatial
  autocorrelation, suppressed-region, and runtime/memory checks. No changes
  were made to `msrc-arg` or production simulation code.
- Saved compact audit tables and figures at
  `validation/msrc_smc_corrected/`; raw chromosome/event tables were kept in
  the temporary run directory and are not committed.

Validation results:

- Corrected full repository suite: pending final rerun for this commit.
- Marginal audit: 40 linked chromosomes and 1,000 independent MSC reference
  trees in each of seven 2-, 3-, 4-, and 8-taxon cases. Of 439 comparisons per
  method, Holm-significant counts were 0 for this SMC′ implementation,
  msprime SMC′, msprime Hudson ARG, and `msrc-arg`. Nominal 95% intervals
  excluded zero 16, 18, 2, and 26 times respectively, consistent with
  multiple-comparison noise.
- Four-taxon long-branch root age: independent MSC mean 1008.0 generations;
  corrected SMC′ 1026.1, chromosome-level 95% difference CI [-2.4, 38.6].
  Eight-taxon long-branch means were 142516.4 and 142619.9, with CI
  [-2228.1, 2435.1].
- Integrated genomic hazards: 23,821 samples, mean 0.9974 for expected 1,
  exponential KS p=0.731. Spatial topology agreement and coalescence-time
  correlation decline with genomic distance and increase in recombination,
  tracking msprime references.
- Suppressed-region experiment: inside recombination density increased from
  zero at κ=0 through 0.000003, 0.000073, and 0.000362 to 0.000827
  events/bp at κ=0.01, 0.1, 0.5, and 1. Marginal quartet differences were
  not significant in this finite sample; lowering the map rate is treated as
  a spatial-dependence effect, not a marginal-frequency change.
- Diagnostic 4–32-taxon runtime: corrected SMC′ median approximately 0.004 s
  at 4 taxa and 0.86 s at 32 taxa; `msrc-arg` approximately 0.006 s and
  0.09 s. Current SMC′ implementation is slower for larger trees in this
  benchmark. Maximum traced Python allocation was about 0.66 MiB.
- Demo event totals: 90, 206, 444, and 148 recombinations for three constant
  rates and one map with a suppressed region; genealogy changes 52, 132, 264,
  and 80.

Mathematical decision and limits:

- This is the guide-genealogy multispecies SMC′ approximation with same-edge
  invisible events and compatibility constrained by species populations. It
  is not the exact ARG and is not the fixed-distance conditional process of
  He, Scornavacca & Chan. Detailed transition and two-tip CDF derivation are
  in `docs/msrc_smc_validation.md`.
- The kernel has conditional empirical support across the audited neutral
  grid; this is not a proof for all trees or parameters. The authors'
  `ipcoal` waiting-distance implementation was not run, so validation is
  conditional pending that direct comparison and expanded tree/Ne settings.
- The implementation supports one contemporaneous sample per species and
  ultrametric species trees. No migration, serial sampling, full ancestral
  background state, or two-background MSRC behavior is implemented.
- Genomic recombination and backward background switching remain separate
  process dimensions. No independent switching clock is introduced.

Next steps:

1. Compare event waiting distances and conditional outcomes directly to the
   paper-associated `ipcoal` calculations and expand species-tree/Ne cases.
2. Profile and optimize transition reconstruction at larger taxon counts.
3. Derive an augmented ancestral state and joint process for backgrounds,
   p(t), selection, switching, recombination, and shared demography before
   starting the neutral-kernel-based MSRC-SMC extension.
4. Keep `msrc-arg` unchanged as a validation benchmark.

## Two-background MSRC-SMC′ Phase A Gate (2026-11-20)

Status: **FAIL / implementation blocked by an under-specified kernel**. No
two-background SMC′ code was added. This follows the task's explicit Phase A
gate: do not implement if the retained state cannot define a coherent
transition.

Completed audit:

- Added `docs/msrc_smc_mathematical_model.md`, reviewing the neutral SMC′,
  legacy single-locus structured coalescent, current MSRC-ARG background
  model, and shared Moran/Wright–Fisher history conventions.
- Verified the current marker-mediated ARG's intended one-focal switching
  intensity under its marker-present assumptions: integrate
  `r0(x) * p_g(t) * eta(h,g,x)` between focal position and anchor. This
  intensity alone does not specify the conditional joint ancestry of focal
  and marker positions.
- Documented a counterexample: distinct marker/phase ancestry configurations
  can have the same local tree, focal background, population, map, and shared
  frequency history, but different conditional local-tree transitions at a
  neighboring position. Thus the proposed local-tree-plus-color state is not
  Markov sufficient for the requested spatial coupling.
- Identified two existing matching hazards: absent-origin handling differs
  (legacy falls back to root age plus extension; ARG remains structured when
  origin is missing), and the shared finite-copy floor allows positive
  coalescence rates for an explicitly supplied lineage in a zero-frequency
  background.

No tests were run because this commit changes only the mathematical audit and
project status. No claims of MSRC-SMC′ implementation, spatial dependence,
marginal equivalence, or performance are made. `msrc-arg`, neutral SMC′, and
production simulation code remain unchanged.

Required next steps before implementation:

1. Choose the marginal target: legacy effective switching, current
   marker-mediated ARG, or a specified forward haplotype process.
2. Define marker inheritance and ancestral phase, and derive either a
   sufficient reduced state plus its conditional kernel or retain the full
   marked segment ARG.
3. Resolve coalescence at background fixation/loss and align origin-age
   defaults before matched marginal testing.
4. Only after Phase A passes, implement limiting-case tests and compare
   topology probabilities, `D=q2-q1`, coalescence-time distributions, and
   chromosome dependence against appropriately matched references.

## Marker-Aware Approximation Assessment (2026-11-21)

The follow-up task explicitly permits a reduced approximation, so candidate
states were assessed rather than requiring exact ARG equivalence. The
minimal binary marker-linkage flag fails: for the same focal tree and flag,
one versus two compatible marker-associated guide lineages give no-recoal
probabilities 0.1353 versus 0.0183 over a `2Ne` interval at `p_h=0.5`.
Consequently it does not define a unique transition law.

A marker gene tree and branch-association map resolves the present target
count but does not define the conditional ancestry of a third position from
the current two-locus state. The provided equations define only per-locus
switch rates, not the hidden multi-position conditional distribution. A
particle approximation would need a specified proposal, weights, and error
analysis; those are not defined in the repository. See
`docs/msrc_smc_marker_aware_validation.md` for the worked calculation and
candidate-state assessment.

This commit adds analysis only. No simulator, tests, or benchmark were added;
neutral SMC′, legacy single-locus MSRC, and MSRC-ARG remain unchanged. Overall
readiness remains **FAIL** until a reproducible hidden-state conditional law
is defined. The next step is to specify and calibrate that law before coding.

## Spatial MSRC Framework Decision and Reference Prototype (2026-11-21)

**Recommendation: trajectory-conditioned structured ancestral material with
a background-reservoir inheritance rule, before pursuing another SMC kernel.**
Status is **CONDITIONAL PASS for the fixed-interior-history reference
prototype**. This supersedes the preceding blanket implementation blocker:
failure of an exact local-tree Markov state does not prevent a coherent
evolutionary-time ancestral-material process. The prior target-count example
is illustrative and does not establish a general impossibility theorem.

Mathematical decisions:

- Primary marginal truth is the legacy conditional MSRC generator, with
  an explicit origin age and positive occupied background sizes. Under the
  reservoir rule, each candidate crossover assigns the marker-distal parent
  a background from the shared realized population-frequency path.
- Always integrate marker-to-carrier exchange, even when no sampled segment
  retains the marker. The current ARG's `_stored_anchor` condition removes
  this measure in some states; a constructed focal-85/marker-0/material-80:90
  case gives rate 5 rather than 85 map-distance units.
- A generator-projection argument yields exactly the legacy one-site switch
  rates and same-background coalescence rates on admissible histories where
  the finite-copy floor is inactive. The full active carrier partition also
  defines consistent multi-position laws without a genomic SMC state.
- Marker phase is approximated by an exchangeable unobserved reservoir.
  Origin erasure preserves the legacy convention but is distinct from a
  unique physical founder. Occupied zero-frequency backgrounds raise an
  explicit error; no forced switch, silent history rejection, or floor is
  introduced into the new mode.

Implementation:

- Added isolated `msrcsim/spatial_reservoir/` API with an opt-in background
  adapter, reusing the unchanged existing ARG scheduler/segment machinery.
  No production simulator, neutral SMC′ transition, or legacy MSRC behavior
  was modified. Existing Moran/WF/frozen history objects can be supplied,
  but arbitrary generated paths hitting incompatible boundaries are not
  yet a validated application domain.
- Added the requested design, complete specification, and implementation
  plan at `docs/msrc_spatial_{model_design,recommended_model,implementation_plan}.md`.
  Added `docs/msrc_spatial_prototype_validation.md` and seeded audit scripts,
  compact machine-readable tables, profiles, and figures under
  `validation/spatial_reservoir/`.

Validation performed:

- Focused regression run: **39 passed in 41.10 seconds**, covering the new
  prototype, existing MSRC-ARG, corrected neutral SMC′, and its marginal tests.
  The full repository suite was not rerun for this isolated extension.
- Primary audit: 18 two-/four-taxon cases, each with 1,000 independent
  chromosome replicates per method, p=0.2/0.5/0.8, k=0/0.1/1, and matched
  histories, terminal states, map-integrated switching, and origin. No
  Holm-significant deviations in 63 reservoir-vs-legacy, 63 old-ARG-vs-legacy,
  or 9 reservoir-vs-independent-two-tip-CTMC comparisons.
- Reservoir quartet point differences were at most 0.012. Of 63 primary
  difference intervals, 61 lay within preset practical tolerances (0.05 for
  quartet frequencies, 10% for time means). The two unresolved time intervals
  in four-taxon p=0.2/k=1 were refined to 10,000 replicates; all five refined
  comparisons fit the tolerances. Root mean difference 1.13 generations,
  95% interval [-0.27,2.54]. The refinement includes the original 1,000 draws.
- Direct marked-event audit: 216,000 proposals across p/background/k settings;
  maximum standardized switching-rate discrepancy 2.381. Constructed absent
  marker domains, origin transition, complete cross suppression, nonzero
  same-background recombination, zero maps, seed reproducibility, and
  fixation/invalid-map errors are covered by tests.
- Spatial diagnostic at 50-bp separation and p=0.5: root correlations 0.570,
  0.537, and 0.363 for k=0,0.1,1, with chromosome-bootstrap intervals. These
  establish dependence in the specified model, not correct biological phase
  dynamics or a multi-distance spatial validation against an independent ARG.
- Diagnostic benchmark, 4/8/16/32 taxa and two short lengths: reservoir
  median times 0.0044/0.0084/0.0222/0.1012 seconds with allocation tracing;
  comparable old ARG medians 0.0045/0.0079/0.0235/0.1058. Neutral SMC′ medians
  0.0032/0.0170/0.0254/2.4023 but generates different event workloads. No
  general speed advantage is claimed. Profiling identifies active-material
  retirement, rate evaluation, and interval operations as main costs.

Remaining limits and next steps:

1. Define founder/boundary bridges before general Moran/WF paths crossing
   fixation/loss; keep path failures explicit to avoid conditioning bias.
2. Expand to internal times, D, dwell/switch counts, time-varying shared
   histories, branch-specific sizes, and an independent selected-ancestry or
   finite-pedigree implementation. Do not treat the old ARG as ground truth.
3. Validate multi-distance correlations, tree identity, tract lengths, and
   map-boundary effects, using independent chromosomes and restriction tests.
4. Optimize only against the now-specified generator and then derive any SMC
   approximation against this reference, measuring both marginal and spatial
   errors. Genome-scale performance is not established.

## Structured Multispecies SMC′ Prototype (2026-11-21)

Scientific status: **CONDITIONAL PASS for an explicitly frozen-reference
window; FAIL as a general whole-chromosome MSRC replacement**. A working
sequential tree/path implementation is now available. It is not another ARG.
The known spatial and off-reference errors below are retained in the report.

Mathematical model:

- Retain one dated local tree, background jump histories on every branch,
  a root-stem history to the origin, and one shared frequency-history object.
  An inherited-state prune/re-coalescence proposal uses the MSRC switching
  and compatible-pair rates, guide population/background constraints, and
  the original-edge target only until its old parent time.
- Candidate spatial rate is the integrated branch-time measure
  `r0(x) * [p_h(t)+k(x)*p_other(t)]`, reverting to r0 after origin.
  Map boundaries integrate the hazard; they do not cause tree changes.
- The target marginal is held at a declared reference position u*: rho is
  the map-integrated marker-to-u* cross-exchange rate. A forward/reverse
  path-density calculation gives a continuous-position MH balance equation;
  total-tree cut-length normalizers must not be inserted into that equation.
- For the inherited-state proposal, the ratio cancels analytically. Default
  fast mode uses the cancellation; `audit_density: true` computes all terms.
  This is a conditional background-path closure and a frozen-window
  approximation, not a literal census of physical ARG crossovers.
- `strict_v1` requires explicit origin and rejects occupied absent or
  sub-one-copy backgrounds. Origin comparison uses the legacy 1e-10 numerical
  convention, including the Moran nextafter origin record. No legacy floor
  behavior or existing production implementation was changed.

Software:

- Added `msrcsim/smc/structured.py`, a runnable module CLI, constant/frozen/
  Moran/WF shared-history inputs, explicit or Bernoulli terminal states,
  date/path outputs, and independent proposal/genealogy/topology/background
  flags. Neutral mode delegates directly to corrected neutral SMC′ with
  identical RNG semantics.
- Added `examples/msrc_smc_structured.yaml`, saved window realization and
  figure under `validation/msrc_smc_structured/example/`, transition tests,
  audited marginal/spatial/history/background/benchmark programs, and compact
  result tables under `validation/msrc_smc_structured/results/`.
- Added full model and quantitative validation documents at
  `docs/msrc_smc_structured_{model,validation}.md`.

Validation completed:

- Final focused new/neutral/reservoir/original-ARG regression suites:
  **56 passed in 44.92 seconds**. Full repository suite was not rerun.
- Primary grid: 30 two-/four-taxon settings × 1,000 independent chromosomes
  per method, p=0.2/0.5/0.8 and k=0/0.01/0.1/0.5/1. Matched coordinate-50
  marginal comparisons cover quartets, D, root/pair/internal times. No
  Holm-significant discrepancies among 150 structured-SMC versus legacy,
  150 reservoir versus legacy, and 15 independent analytic two-lineage tests.
- SMC quartets differed at most 0.011 and root means at most 3.8% in point
  estimates. Initial practical intervals passed 147/150 SMC comparisons.
  The three unresolved D cases were expanded to 3,000 replicates; all 21 SMC
  comparisons in those expanded cases fit the original tolerance. Refined
  runs include the first 1,000, so they are not independent replications.
- Background/dwell audit: 90 between-method comparisons and 45 integrated-
  exposure rate tests, no Holm-significant differences. Only 34/45 SMC
  background practical intervals fit 10% tolerance; remaining rare-count/
  dwell precision is unresolved, not declared equivalent.
- Actual shared-history checks: 500 SMC and legacy draws each for a realized
  Moran path (154 records) and WF path (29), origin age 8, history seed 93.
  No nominal p<0.05 in 14 comparisons. These are admitted conditional paths,
  not validation of every generated trajectory or founder boundary.
- Spatial comparisons include genealogy/topology agreement, rooted RF,
  pair/root correlation, topology agreement above its independent-marginal
  baseline, completed/censored waits and tract lengths. A five-k boundary
  experiment adds 500 independent pairs per setting across both edges of a
  suppressed region. Uncertainty resamples complete chromosomes.

Scientific discrepancies and readiness limits:

- Broad uniform-map p=0.5/k=1: at 90-bp separation pair-time correlation is
  0.497 for SMC versus 0.151 for reservoir ARG. Excess persistence is
  substantial and reported, not calibrated away.
- Freezing rho fails the intended position-specific marginals across a broad
  chromosome. At position 5, q2 is 0.842 versus reservoir 0.977, difference
  -0.135, 95% CI [-0.159,-0.111]; at position 95, 0.856 versus 0.751,
  difference +0.105, CI [0.070,0.140]. This is why whole-chromosome readiness
  is FAIL despite reference-coordinate marginal support.
- Within a strongly suppressed region whose marker-path integral is nearly
  constant, the approximation is more plausible. At 350-bp separation in
  the boundary experiment, correlation errors were -0.054 for k=0 and
  -0.011 for k=0.01, but intervals remain wide; no spatial equivalence claimed.

Performance:

- Profiling identified repeated species-LCA traversal and full density scans.
  Cached clade populations and target density preserved a 20-seed local-tree/
  background output hash; the analytic fast path is also directly tested
  against audited seeded states.
- Diagnostic 32-taxon/1,000-bp runtime went from 84.2 s in uncached audit mode
  to 2.13 s in cached audit mode and 1.19 s in analytic fast mode. Reservoir
  ARG was about 0.097 s in the final small benchmark. No overall speed
  advantage or genome-scale claim is supported. Marks and stochastic workload
  counts differ between algorithms; one-seed timings are descriptive.

Next steps:

1. Establish usable window-size/rho-variation and spatial-error limits for
   nearly constant marker-path integrals, especially complete/strong cross
   suppression with ongoing within-background recombination.
2. Derive position-varying marginal transport and a more accurate conditional
   parental-background closure; local invariance at each x alone is insufficient.
3. Specify genuine fixation/loss founder bridges and test broader shared
   Moran/WF paths. Keep boundary failures explicit and avoid silent retries.
4. Optimize cut-history integration only against the stated kernel; benchmark
   replicated longer chromosomes and matched observables before speed claims.
