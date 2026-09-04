# MSRC Simulator

`msrc-sim` simulates quartet gene-tree distributions under the Multi-Species
Rearrangement Coalescent (MSRC) model. It is designed for experiments where a
chromosomal rearrangement, such as an inversion, arises in a population,
evolves forward in time with a Wright-Fisher process, and then affects
backward-time genealogies through arrangement-dependent coalescence and
recombination.

The simulator currently focuses on four sampled taxa. It can be used to:

- simulate a single mechanistic MSRC history and its gene trees;
- compute empirical and exact quartet probabilities for a fixed structured
  interval;
- run independent evolutionary replicates and estimate prevalence statistics;
- condition replicates on persistence or terminal arrangement patterns;
- run multidimensional parameter grids for prevalence analyses;
- run spatially ordered MSRC and pulse-hybridization comparator simulations.

The independent unit in prevalence analyses is an evolutionary replicate, not a
locus. Each accepted replicate draws one rearrangement frequency history, and
the loci within that replicate estimate the quartet distribution conditional on
that history.

## Installation

Clone the repository and install it in editable mode:

```bash
git clone https://github.com/ytabatabaee/msrc-sim.git
cd msrc-sim
pip install -e .
```

For development and tests, install the test dependencies:

```bash
pip install -e ".[test]"
pytest
```

The package requires Python 3.9 or later, NumPy, SciPy, and PyYAML.

## Commands

The package installs these command-line programs. Simulation commands read YAML
configuration files.

```bash
msrc-sim --config <config.yaml>
```

Runs one simulation. The config `mode` can be `mechanistic` or `conditional`.

```bash
msrc-sim-replicates --config <config.yaml>
```

Runs independent evolutionary replicates from a `replicate_experiment`
configuration.

```bash
msrc-sim-grid --config <config.yaml>
```

Runs a parameter grid from a `parameter_grid` configuration. Each grid cell is
run as a replicate experiment.

```bash
msrc-sim-compare --input <replicate_summary.csv> --output <comparison.csv>
msrc-sim-freeze-history --config <config.yaml> --output <history.yaml>
msrc-sim-replay-history --history <history.yaml> --num-loci <n> --output <dir>
msrc-sim-plot --input <replicate_summary.csv> --output <figures>
msrc-sim-plot-history --run-dir <run_output> --output <history.pdf>
```

Compare quartet vectors, freeze and replay realized rearrangement histories,
create automated replicate/model-comparison figures, and render a static
Wright-Fisher frequency-history tree. Mechanistic simulations can also render
that history figure automatically with `output.make_history_plot: true`.

```bash
msrc-sim-spatial --config examples/spatial_inversion.yaml
msrc-sim-plot-spatial --input spatial_output --format png
msrc-sim-spatial-summarize --input loci.csv --breakpoints 25000000,65000000 --window-loci 50 --step-loci 10 --output summary_dir
msrc-sim-hybridization --config examples/pulse_hybridization.yaml
msrc-sim-plot-spatial-compare --msrc-dir spatial_output --hyb-dir hybridization_output --output comparison/spatial_compare.png
msrc-sim-find-matched-hybridization --msrc-dir spatial_output --hyb-grid-dir hybridization_grid --output matched_hybridization.json
msrc-sim-match-hybridization --msrc-dir spatial_output --major-topology '12|34' --introgressed-topology '13|24' --output matched_hybridization.yaml
msrc-sim-spatial-distinguishability --msrc-dir spatial_output --matched-hybridization-yaml matched_hybridization.yaml --output-dir distinguishability
msrc-sim-spatial-identifiability --msrc-dir spatial_output --output-dir identifiability --replicates 100
```

Run the v0.7.0 spatially ordered locus prototype, replot a spatial output
directory, summarize an external ordered table with `position` and `topology`
columns, run the pulse-hybridization spatial comparator, or make a matched
side-by-side MSRC versus hybridization spatial comparison figure. When the
scientific question requires comparable bag-of-genes outcomes, use the matcher
to select a precomputed hybridization run with a marginal quartet vector close
to the MSRC run before plotting.

Equivalent script wrappers are provided in `scripts/`:

```bash
python scripts/simulate_msrc.py --config examples/mechanistic_balanced.yaml
python scripts/simulate_replicates.py --config examples/replicates_unconditional.yaml
python scripts/simulate_parameter_grid.py --config examples/parameter_grid.yaml
```

## Simulation Modes

### Mechanistic MSRC Simulation

Mechanistic mode simulates a rearrangement history forward through a dated
quartet species tree and then simulates locus genealogies backward through that
realized history.

```bash
msrc-sim --config examples/mechanistic_balanced.yaml
```

The main configuration sections are:

- `mode`: set to `mechanistic`;
- `seed`: random seed;
- `num_loci`: number of loci to simulate;
- `species_tree`: ultrametric four-taxon Newick tree, root extension, and
  effective population sizes;
- `rearrangement`: rearrangement type, origin branch, origin time, initial copy
  count, and selection coefficient;
- `recombination`: baseline recombination rate and effective cross-arrangement
  fraction;
- `sampling`: currently one sample per species;
- `output`: output directory and which intermediate files to record.

Example:

```yaml
mode: mechanistic
seed: 12345
num_loci: 200
species_tree:
  newick: "((1:100,2:100)A:50,(3:100,4:100)B:50)ROOT;"
  time_units: generations
  root_extension: 500
  default_effective_population_size: 1000
rearrangement:
  id: inv_1
  type: inversion
  origin_branch: ROOT
  origin_time_from_branch_start: 120
  initial_copy_count: 20
  selection:
    model: genic
    coefficient: 0.0
recombination:
  baseline_rate: 0.01
  effective_cross_arrangement_fraction: 0.05
sampling:
  samples_per_species: 1
output:
  directory: balanced_output
  record_frequency_history: true
  record_sampled_arrangements: true
  record_gene_trees: true
  record_backward_events: true
  make_history_plot: true
  history_plot:
    filename: wright_fisher_history.pdf
    glyphs_per_row: 12
    max_rows_per_branch: 30
```

Species trees must be ultrametric and must have exactly four sampled taxa.
Internal node names are used as branch identifiers, so named internal nodes such
as `A`, `B`, and `ROOT` are recommended.

### Conditional Quartet Simulation

Conditional mode simulates quartet outcomes for a fixed four-lineage
arrangement configuration over one structured interval. It also computes the
exact conditional quartet matrix `H_m(t)` for the same model parameters.

```bash
msrc-sim --config examples/conditional_quartet.yaml
```

Example:

```yaml
mode: conditional
seed: 12345
num_loci: 100000
structured_interval:
  duration: 1.0
  configuration: "1010"
  migration:
    m01: 0.05
    m10: 0.05
  coalescence:
    lambda0: 1.0
    lambda1: 1.0
output:
  directory: conditional_output
```

The `configuration` string gives the arrangement state of the four lineages.
For example, `1010` means lineages 1 and 3 carry state `1`, while lineages 2
and 4 carry state `0`.

### Replicate Experiments

Replicate experiments simulate many independent rearrangement histories. For
each accepted history, the simulator samples terminal arrangements and then
simulates a fixed number of loci to estimate the replicate's quartet
distribution.

```bash
msrc-sim-replicates --config examples/replicates_unconditional.yaml
```

The `experiment` section controls the number of accepted replicates, loci per
replicate, and prevalence statistics:

```yaml
experiment:
  replicates: 20
  loci_per_replicate: 200
  asymmetry_threshold: 0.10
  persistence_target_branches: [ROOT]
```

The simulator records both attempted and accepted histories. This is important
for conditioned experiments because rare conditioning events should remain
visible in the reported acceptance rate.

### Conditioning

Replicate experiments can be unconditioned or conditioned on features of the
forward rearrangement history.

Unconditional simulation:

```yaml
conditioning:
  mode: none
```

Require the rearrangement to remain segregating at specified branch ends:

```yaml
conditioning:
  mode: persistent_polymorphism
  require_segregating_at: [ROOT]
  max_attempts: 10000
```

Require one of a set of sampled terminal arrangement patterns:

```yaml
conditioning:
  mode: terminal_pattern
  accepted_patterns: ["1010", "0101"]
  max_attempts: 100000
```

Require both persistence and terminal pattern conditions:

```yaml
conditioning:
  mode: persistent_and_pattern
  require_segregating_at: [ROOT]
  accepted_patterns: ["1010", "0101"]
  max_attempts: 100000
```

### Parameter Grids

Parameter-grid mode runs a replicate experiment for every combination of values
in `parameter_grid`. Grid keys are dotted paths into the base YAML
configuration.

```bash
msrc-sim-grid --config examples/parameter_grid.yaml
```

Example:

```yaml
parameter_grid:
  species_tree.default_effective_population_size: [25, 50]
  rearrangement.initial_copy_count: [10, 20]
  recombination.effective_cross_arrangement_fraction: [0.01, 0.1]
```

Each grid cell is written to its own `cell_####` directory, and the grid-level
summary is written to `parameter_grid_summary.csv`.

### Spatially Ordered Loci

Spatial mode simulates one realized forward rearrangement history, samples one
terminal arrangement pattern, places loci along a chromosome, and uses a local
model at each ordered locus. Loci inside the configured rearrangement interval
use the mechanistic MSRC genealogy simulator with the shared realized
Wright-Fisher history and sampled terminal arrangements. Loci outside the
interval use an ordinary unstructured MSC genealogy on the same species tree.

```bash
msrc-sim-spatial --config examples/spatial_inversion.yaml
```

Spatial mode is designed for comparing a bag-of-genes quartet vector with a
chromosome-wide quartet-support track. If genomic positions are discarded, the
analysis retains only the total counts of `12|34`, `13|24`, and `14|23`.
With ordered positions retained, the output can show whether quartet support is
localized around a physical rearrangement interval and its known structural
breakpoints.

This is not a linked-locus ARG simulator. In v0.7.0, loci are conditionally
independent given the realized rearrangement history and local region
parameters. The simulator does not generate topology autocorrelation or
tract-length information from a multi-locus genealogy.

Spatial mode supports frozen-history replay:

```bash
msrc-sim-spatial --config spatial.yaml --history frozen_history.yaml
```

This allows different genomic or recombination configurations to be compared
while holding the same forward rearrangement history and sampled terminal
arrangements fixed.

### Pulse-Hybridization Spatial Comparator

`msrc-sim-hybridization` implements a spatial pulse-hybridization comparator.
A hybridization event introduces a fraction `gamma` of ancestry from an
alternative parental history. Recombination during the subsequent `h`
generations breaks that ancestry into tracts along the chromosome. Local
quartet topologies are sampled from the MSC distribution associated with the
local ancestry state. The model therefore captures spatially correlated
introgressed ancestry but is not a full ARG or full multispecies-network-
coalescent simulator.

The causal interpretation is explicit:

- hybridization = origin of introgressed ancestry;
- recombination = fragmentation of that ancestry into tracts.

For each genomic position, the latent ancestry state is `0` for the major/native
parental history and `1` for the introgressed parental history. Conditional on
that state, the local quartet topology is sampled independently from
`msc_probabilities(topology, internal_branch_length)`. The ancestry mosaic is
linked along the chromosome, but this release does not simulate a full ARG or
full linked multispecies-network coalescent. Conditional on the ancestry state
at a locus, the local quartet topology is sampled independently from the
corresponding MSC distribution.

The genome-wide bag-of-genes vector is written as `expected_marginal_q`:

```text
(1 - gamma) * q_major + gamma * q_introgressed
```

This is the same two-tree mixture formula used by `msrcsim.model_fitting`.
The ordered output keeps the spatial information needed to ask whether MSRC and
hybridization can produce similar genome-wide quartet-frequency vectors while
producing different chromosome-wide profiles.

### MSRC vs Hybridization Spatial Comparison

`msrc-sim-plot-spatial-compare` creates matched spatial comparison plots for an
MSRC spatial run and a pulse-hybridization run. The figure is designed to show
how two models may have similar genome-wide averaged quartet frequencies while
differing in the spatial organization of local quartet support along a
chromosome. In the MSRC model, local signal can be associated with rearranged
structural intervals; in the hybridization model, local signal can follow
recombined ancestry tracts.

This is a visualization and reporting tool for studying identifiability. A
bag-of-genes comparison alone may be non-identifiable because it discards
genomic order, while spatially ordered loci may provide additional
identifiability through the organization of `q1(x)`, `q2(x)`, and `q3(x)`.
The command does not claim that MSRC and hybridization are always spatially
distinguishable.

Example:

```bash
msrc-sim-spatial --config examples/spatial_inversion.yaml
msrc-sim-hybridization --config examples/pulse_hybridization.yaml

msrc-sim-plot-spatial-compare \
  --msrc-dir spatial_output \
  --hyb-dir hybridization_output \
  --output comparison/spatial_compare.png
```

If window files are absent, windows can be recomputed from locus-level outputs:

```bash
msrc-sim-plot-spatial-compare \
  --msrc-dir spatial_output \
  --hyb-dir hybridization_output \
  --window-size-loci 50 \
  --step-loci 10 \
  --output comparison/spatial_compare.pdf
```

Use `--formats png,pdf,svg` to write multiple figure formats from the same
comparison. The command also writes `spatial_compare_summary.json`, recording
the compared marginal quartet vectors and their L1 difference. By default it
also writes a compact `spatial_compare_overlay.<format>` figure that overlays
the model-specific local fraction and one topology-support track.

For calibrated comparisons, first run a grid or collection of hybridization
simulations, then select the run whose marginal quartet vector is closest to
the MSRC run:

```bash
msrc-sim-find-matched-hybridization \
  --msrc-dir spatial_output \
  --hyb-grid-dir hybridization_grid_results \
  --metric l1 \
  --output matched_hybridization.json

msrc-sim-plot-spatial-compare \
  --msrc-dir spatial_output \
  --matched-hybridization-json matched_hybridization.json \
  --output comparison/spatial_compare.png
```

The matcher scans precomputed hybridization run directories, ranks candidates
by distance between marginal quartet vectors, and records the selected
`best_hyb_dir`. It does not simulate new hybridization parameters; it selects
from the runs already present under `--hyb-grid-dir`.

For exact expected bag-of-genes calibration, use
`msrc-sim-match-hybridization`. This fits `gamma`, the major-parent MSC branch
length, and the introgressed-parent MSC branch length so the expected two-tree
mixture vector

```text
(1 - gamma) q_MSC(T_major, t_major) + gamma q_MSC(T_intro, t_intro)
```

matches the MSRC marginal quartet vector as closely as possible:

```bash
msrc-sim-match-hybridization \
  --msrc-dir spatial_output \
  --major-topology '12|34' \
  --introgressed-topology '13|24' \
  --output matched_hybridization.yaml
```

The YAML records the target vector, fitted vector, fitted `gamma`, `t_major`,
`t_introgressed`, and L1/Euclidean errors. When the target vector lies on the
chosen two-tree mixture surface, the expected-vector fit should be essentially
exact; finite simulated chromosomes can still have sampling error in their
observed quartet frequencies.

The fitted expected mixture can be used in an age/recombination experiment that
holds bag-of-genes quartet frequencies fixed while changing ancestry-tract
structure:

```bash
msrc-sim-spatial-distinguishability \
  --msrc-dir spatial_output \
  --matched-hybridization-yaml matched_hybridization.yaml \
  --output-dir distinguishability \
  --h-values 10,25,50,100,250,500,1000 \
  --r-multipliers 0.25,0.5,1,2 \
  --baseline-r 1e-8 \
  --replicates 100
```

This writes one feature row per model/parameter/replicate to
`spatial_distinguishability.csv`, including marginal quartet frequencies,
tract summaries, topology change-point counts, same-topology probabilities at
configured lags, longest focal-topology interval, focal support concentrated in
the model-specific interval or tract, and MSRC breakpoint-distance summaries.
It also trains a simple reproducible logistic classifier from spatial features
only and writes accuracy/ROC-AUC summaries plus a four-panel figure. The
analysis is intended to reveal regimes where spatial organization separates
the models and regimes where long hybridization tracts can make them difficult
to distinguish; it is not a universal identifiability claim.

`msrc-sim-spatial-identifiability` runs the paper-oriented matched experiment
using the closed-form T1/T2 mixture. For a target MSRC quartet vector
`q = (q1, q2, q3)` and selected major/introgressed topologies, it computes the
mathematically feasible `gamma` interval, samples interior gamma values across
that interval, and sets:

```text
d_major = q_major - q_shared
d_intro = q_intro - q_shared
t_major = -log(1 - d_major / (1 - gamma))
t_intro = -log(1 - d_intro / gamma)
```

The main spatial grid is `eta = E[L_intro] / L_R`, where `L_R` is the
rearranged interval length. For each `(gamma, eta)`, the command solves
`h*r = 1 / (eta * L_R * (1 - gamma))`, simulates matched hybridization
chromosomes, generates independent MSRC resampled chromosomes on the same locus
positions, and trains two reproducible classifiers:

- a bag-of-genes baseline using only matched `q1`, `q2`, and `q3`;
- a spatial classifier using topology-run, autocorrelation, off-arm, focal-run,
  and breakpoint-distance features.

Hybridization-specific latent quantities such as ancestry-tract counts and
realized introgressed fraction are written as diagnostics but excluded from the
main spatial classifier.

Example:

```bash
msrc-sim-spatial-identifiability \
  --msrc-dir spatial_output_q2_dominant \
  --major-topology '12|34' \
  --introgressed-topology '13|24' \
  --num-gamma 10 \
  --eta-values 0.02,0.05,0.1,0.25,0.5,1,2,5,10 \
  --replicates 100 \
  --output-dir identifiability_q2_dominant
```

The result is designed to show that identical expected bag-of-genes quartet
frequencies can be non-identifying, while spatial organization can restore
information about the mechanism in some parameter regimes.
The command rejects very small replicate counts by default because ROC-AUC and
accuracy estimates are not interpretable with only a few chromosomes per
class. Use at least the default `--replicates 100` for analysis. The
`--allow-small-sample` flag exists only for debugging and automated smoke
tests.

Example:

```yaml
mode: pulse_hybridization
seed: 12345
chromosome:
  length_bp: 100000000
  num_loci: 5000
  locus_positions:
    mode: evenly_spaced
hybridization:
  gamma: 0.25
  generations_since_pulse: 200
  major:
    topology: "12|34"
    internal_branch_length: 0.5
  introgressed:
    topology: "13|24"
    internal_branch_length: 1.5
  donor: "3"
  recipient: "2"
recombination:
  rate_per_bp_per_generation: 1.0e-8
windows:
  loci_per_window: 50
  step_loci: 10
```

The `donor` and `recipient` fields are recorded as biological metadata. In this
release, the explicitly supplied parental topology and internal branch length
fields are authoritative.

## Outputs

### Mechanistic Outputs

Mechanistic runs write files to `output.directory`. Depending on the output
flags, the directory can contain:

- `config.resolved.yaml`: YAML configuration after defaults are applied;
- `frequency_history.csv`: forward Wright-Fisher frequency path on each branch;
- `sampled_arrangements.csv`: sampled terminal arrangement state for each taxon;
- `wright_fisher_history.png` or `.pdf`: optional static visualization created
  by `msrc-sim-plot-history`;
- `true_gene_trees.nwk`: simulated true gene trees in Newick format;
- `coalescence_times.csv`: coalescence times for each locus;
- `coalescence_events.csv`: coalescence-event records;
- `genealogy_events.csv`: optional full backward-event log for selected loci;
- `summary.json`: topology counts, topology frequencies, taxa, and sampled
  arrangements.

### Conditional Outputs

Conditional runs write:

- `config.resolved.yaml`: YAML configuration after defaults are applied;
- `quartet_probabilities.csv`: topology counts, empirical probabilities, exact
  probabilities, and absolute errors;
- `summary.json`: the same conditional summary in JSON format.

The three quartet topologies are reported as `12|34`, `13|24`, and `14|23`.

### Wright-Fisher History Visualization

The Wright-Fisher history plot uses `frequency_history.csv` to render the
realized rearrangement-frequency trajectory directly on the species tree.
Colored arrows summarize the proportion of ancestral and rearranged chromosomes
at selected generations. The visualization does not resimulate the process and
does not treat displayed arrows as individual chromosome copies.

Set `output.make_history_plot: true` in a mechanistic configuration to write
`wright_fisher_history.png` automatically, or set
`output.history_plot.filename` to choose a different output file such as
`wright_fisher_history.pdf`.

Existing runs can be replotted with explicit display controls:

```bash
msrc-sim-plot-history \
  --run-dir balanced_output \
  --output balanced_output/wright_fisher_history.png \
  --glyphs-per-row 14 \
  --max-rows-per-branch 35 \
  --tip-order 1,2,3,4
```

`frequency_history.csv` is the only required data file. When present,
`sampled_arrangements.csv` is used for terminal taxon labels; otherwise the
plot labels terminal population frequencies without implying sampled
chromosome states.

For a compact example where all four present-day populations remain
polymorphic:

```bash
msrc-sim --config examples/mechanistic_balanced_polymorphic_tips.yaml
```

That example enables `output.make_history_plot: true`, so the run writes
`polymorphic_tips_output/wright_fisher_history.png` alongside
`frequency_history.csv`.

For a higher-variance example with mixed terminal fates:

```bash
msrc-sim --config examples/mechanistic_balanced_mixed_fates.yaml
```

This run writes `mixed_fates_output/wright_fisher_history_mixed_fates.png` and
shows a combination of loss, low-frequency polymorphism, and A1 fixation at
the present-day tips.

### Replicate Outputs

Replicate experiments write:

- `replicate_summary.csv`: one row per attempted history, including rejected
  attempts;
- `simplex_points.csv`: accepted replicate quartet probabilities and simplex
  coordinates;
- `prevalence_summary.json`: acceptance rate, terminal pattern counts,
  prevalence estimates, and 95% Wilson intervals;
- `config.resolved.yaml`: the experiment configuration.

The prevalence summary includes estimates for:

- persistent polymorphism on target branches;
- 2:2 terminal arrangement patterns;
- asymmetric quartet distributions;
- discordant-topology dominance.

### Grid Outputs

Parameter grids write:

- `parameter_grid_summary.csv`: one row per parameter combination;
- `cell_####/replicate_summary.csv`: replicate-level records for each cell;
- `cell_####/simplex_points.csv`: accepted quartet-simplex points for each cell;
- `cell_####/prevalence_summary.json`: prevalence summary for each cell;
- `cell_####/config.resolved.yaml`: resolved cell configuration.

### Spatial Outputs

Spatial runs write:

- `spatial_loci.csv`: one row per ordered locus with `locus_id`, `position`,
  region, local model, topology, terminal pattern, switch count, and
  coalescence-time summary;
- `spatial_gene_trees.tsv`: optional position-aware gene trees with
  `locus_id`, `position`, and `newick`;
- `spatial_windows.csv`: sliding-window topology counts, local quartet
  concordance factors, nearest-MSC-arm distance, off-arm statistic, and
  fraction of loci inside the rearrangement interval;
- `spatial_summary.json`: whole-chromosome bag-of-genes vector, inside/outside
  vectors, spatial contrasts, breakpoint-aligned jump summaries, strongest
  adjacent-window jumps, and metadata stating that linked loci are not modeled;
- `spatial_quartet_profile.png`: chromosome-wide quartet-support track when
  plotting dependencies are available;
- `spatial_bag_inside_outside.png`: comparison of overall, inside, and outside
  quartet vectors when plotting dependencies are available.

### Hybridization Outputs

Pulse-hybridization runs write:

- `hybridization_tracts.csv`: ancestry tracts with start, end, length, state,
  and `major`/`introgressed` labels;
- `hybridization_loci.csv`: ordered loci with `position_bp`, `model`,
  ancestry state, tract ID, sampled topology, topology label, and conditional
  MSC probabilities;
- `hybridization_windows.csv`: moving-window `q1`, `q2`, `q3`,
  introgressed fraction, dominant topology, nearest-MSC-arm distance, and
  off-arm difference;
- `hybridization_spatial_autocorrelation.csv`: same-topology and same-ancestry
  probabilities at configured genomic lags;
- `hybridization_summary.json`: expected and observed marginal quartet
  vectors, realized ancestry fractions, tract-length summaries, topology
  counts/frequencies, parental parameters, and model-limitation metadata;
- `hybridization_spatial_profile.png`: ancestry mosaic, moving quartet support,
  and moving introgressed fraction when plotting is enabled.

Spatial comparison runs write:

- `spatial_compare.png`, `.pdf`, or `.svg`: a two-column MSRC versus
  hybridization figure with latent interval structure, local quartet support,
  local model-specific fraction, and bag-of-genes summaries;
- `spatial_compare_summary.json`: compared marginal quartet vectors, their L1
  difference, run sizes, chromosome length, figure paths, and a warning when
  the marginal vectors are not closely matched;
- `spatial_compare_overlay.<format>`: optional direct overlay of local feature
  fraction and a focal quartet-support track.
- `matched_hybridization.json`: when written by
  `msrc-sim-find-matched-hybridization`, the selected `best_hyb_dir`, MSRC and
  hybridization marginal quartet vectors, the matching distance, ranked
  candidate runs, and a warning flag when no close match was found.
- `matched_hybridization.yaml`: when written by
  `msrc-sim-match-hybridization`, exact expected-mixture fit parameters
  `gamma`, `t_major`, `t_introgressed`, target and fitted quartet vectors, and
  L1/Euclidean fit errors;
- `spatial_distinguishability.csv`: one row per MSRC bootstrap or
  hybridization chromosome replicate with matched marginal-q fields and spatial
  features;
- `spatial_distinguishability_metrics.csv`: train/test classification accuracy
  and ROC-AUC by hybridization age and recombination rate;
- `spatial_distinguishability.png`: example matched profiles, latent tracks,
  accuracy versus age, and accuracy heatmap over `h x r`.
- `spatial_identifiability_replicates.csv`: one row per matched MSRC or
  hybridization chromosome replicate across the closed-form `gamma x eta`
  grid;
- `spatial_identifiability_summary.csv`: bag-of-genes and spatial classifier
  accuracy/ROC-AUC with cross-validation confidence intervals for each
  `(gamma, eta)` cell;
- `identifiability_heatmap.pdf`: spatial ROC-AUC heatmap over `gamma x eta`;
- `auc_vs_eta.pdf`: spatial ROC-AUC versus `eta` for the gamma grid;
- `bag_vs_spatial_auc.pdf`: bag-of-genes baseline versus spatial classifier;
- `example_profiles.pdf`: matched MSRC/HYB example profiles for small,
  intermediate, and large `eta`.

## Configuration Reference

Common fields:

- `mode`: one of `mechanistic`, `conditional`, `replicate_experiment`,
  `parameter_grid`, `spatial`, or `pulse_hybridization`;
- `seed`: random seed, defaulting to `1` for single-run modes;
- `num_loci`: number of loci for single-run modes;
- `output.directory`: output directory.

Species-tree fields:

- `species_tree.newick`: ultrametric four-taxon Newick tree;
- `species_tree.root_extension`: length of the population above the root;
- `species_tree.default_effective_population_size`: default diploid effective
  population size used by branches;
- `species_tree.branch_parameters.<branch>.effective_population_size`: optional
  branch-specific effective population size.

Rearrangement fields:

- `rearrangement.id`: identifier written to outputs;
- `rearrangement.type`: rearrangement type, such as `inversion`;
- `rearrangement.origin_branch`: branch where the rearrangement originates;
- `rearrangement.origin_time_from_branch_start`: forward-time origin location on
  the origin branch;
- `rearrangement.initial_copy_count`: initial number of rearranged chromosomes;
- `rearrangement.selection.coefficient`: genic selection coefficient.

Recombination fields:

- `recombination.baseline_rate`: baseline switching/recombination rate;
- `recombination.effective_cross_arrangement_fraction`: fraction of the
  baseline rate retained across arrangements.

Output flags for mechanistic runs:

- `record_resolved_config`;
- `record_frequency_history`;
- `record_sampled_arrangements`;
- `record_gene_trees`;
- `record_coalescence_times`;
- `record_backward_events`;
- `event_log_loci.first_n`.
- `make_history_plot`: render `frequency_history.csv` as a static species-tree
  history figure after a mechanistic run, default `false`;
- `history_plot.filename`, `.format`, `.dpi`, `.glyphs_per_row`,
  `.max_rows_per_branch`, `.tip_order`, `.width_mode`,
  `.show_frequency_trace`, and `.title`: optional display controls for the
  automatic history figure.

Spatial fields:

- `genome.length`: chromosome length in integer coordinate units;
- `genome.loci.count`: number of ordered loci;
- `genome.loci.placement`: `evenly_spaced` or `uniform_random`;
- `genome.rearrangement_interval.start` and `.end`: one physical interval,
  requiring `0 <= start < end <= genome.length`;
- `genome.inside_model.type`: `msrc`;
- `genome.inside_model.effective_cross_arrangement_fraction`: local
  cross-arrangement fraction inside the interval;
- `genome.outside_model.type`: `msc`;
- `spatial_summary.window_loci`, `.step_loci`, and
  `.breakpoint_bandwidth_loci`: sliding-window and breakpoint summary sizes.

Pulse-hybridization fields:

- `chromosome.length_bp`: chromosome length in base pairs;
- `chromosome.num_loci`: number of ordered loci;
- `chromosome.locus_positions.mode`: `evenly_spaced`, `random_uniform`, or
  `file`;
- `hybridization.gamma`: initial introgressed ancestry fraction created by the
  pulse;
- `hybridization.generations_since_pulse`: generations of recombination after
  the pulse;
- `hybridization.major` and `.introgressed`: authoritative parental topology
  labels (`12|34`, `13|24`, `14|23`) and internal branch lengths;
- `hybridization.donor` and `.recipient`: metadata only in this release;
- `recombination.rate_per_bp_per_generation`: recombination rate per bp per
  generation for ancestry-tract breakpoints;
- `windows.loci_per_window` and `.step_loci`: moving-window sizes;
- `spatial_statistics.lags_bp`: genomic lags for same-topology and
  same-ancestry probabilities.

## Examples

The `examples/` directory contains ready-to-run configurations:

- `mechanistic_balanced.yaml`: mechanistic simulation on a balanced quartet tree;
- `mechanistic_balanced_mixed_fates.yaml`: short balanced-tree example with a
  mix of terminal loss, polymorphism, and A1 fixation;
- `mechanistic_balanced_polymorphic_tips.yaml`: short neutral balanced-tree
  example designed to keep the terminal populations segregating;
- `mechanistic_unbalanced.yaml`: mechanistic simulation on an unbalanced quartet
  tree;
- `conditional_quartet.yaml`: fixed-configuration conditional quartet
  simulation;
- `replicates_unconditional.yaml`: unconditioned prevalence experiment;
- `replicates_conditioned.yaml`: terminal-pattern-conditioned prevalence
  experiment;
- `parameter_grid.yaml`: multidimensional parameter grid;
- `spatial_inversion.yaml`: ordered-locus spatial inversion prototype;
- `pulse_hybridization.yaml`: spatial pulse-hybridization comparator.

## Development

Run the test suite with:

```bash
pytest
```

Release notes are maintained in [CHANGELOG.md](CHANGELOG.md).
