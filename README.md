# msrc-sim

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
- run multidimensional parameter grids for prevalence analyses.

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
```

Run the v0.7.0 spatially ordered locus prototype, replot a spatial output
directory, or summarize an external ordered table with `position` and
`topology` columns.

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

## Configuration Reference

Common fields:

- `mode`: one of `mechanistic`, `conditional`, `replicate_experiment`,
  `parameter_grid`, or `spatial`;
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

## Examples

The `examples/` directory contains ready-to-run configurations:

- `mechanistic_balanced.yaml`: mechanistic simulation on a balanced quartet tree;
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
- `spatial_inversion.yaml`: ordered-locus spatial inversion prototype.

## Development

Run the test suite with:

```bash
pytest
```

## Version 0.5.0: off-arm and quartet-model comparison

Version 0.5.0 makes backward genealogy simulation piecewise exact with respect
to the generation-by-generation Wright–Fisher frequency path. A proposed
Gillespie event can no longer cross a frequency-change boundary while retaining
outdated rates.

Replicate outputs now include topology counts, distance to the nearest MSC arm,
an off-arm contrast and confidence interval, maximum-likelihood MSC fits, and a
quartet-level two-tree introgression-mixture fit.

Run model comparison on an existing replicate table with:

```bash
msrc-sim-compare \
  --input replicate_output/replicate_summary.csv \
  --output model_comparison.csv
```

The comparison distinguishes a strong alternative-tree signal from a genuine
off-arm signal. For a vector whose best MSC topology is `T2`, the off-arm
contrast compares the two minor probabilities, `q1 - q3`. The two-tree mixture
is a quartet-level model-comparison device; its fitted parameters should not be
interpreted as uniquely identifiable demographic estimates from one quartet.


## v0.5.1 model-comparison interpretation

The comparison output deliberately separates three questions:

- `network_representable`: can the two-tree mixture reproduce the observed quartet vector geometrically?
- `off_arm_supported`: does the empirical vector significantly violate the nearest single-tree MSC arm at the 0.05 level?
- `network_aic_preferred` / `network_strongly_preferred`: is the network likelihood worth its additional parameters (`delta AIC < 0` / `< -4`)?

`best_network_boundary_warning` is true when the fitted gamma is close to 0 or 1, a branch length is close to zero, or a branch length reaches the optimization ceiling. Such a fit can still be geometrically valid, but its parameters should not be described as a well-interior introgression estimate. Terminal patterns are always written as four-character strings such as `0101`.


## v0.6.0: frozen-history validation and automated figures

Freeze one accepted evolutionary history:

```bash
msrc-sim-freeze-history --config examples/replicates_conditioned.yaml --output frozen_history.yaml
```

Replay genealogy simulations on exactly that history:

```bash
msrc-sim-replay-history --history frozen_history.yaml --num-loci 100000 --seed 7 --output replay_100k
```

The replay output includes `replay_summary.json`, `branch_history_summary.csv`, and `true_gene_trees.nwk`.
This supports locus-count convergence experiments without resimulating the Wright–Fisher trajectory.

Create automated figures from a replicate or comparison table:

```bash
msrc-sim-plot --input replicate_output/replicate_summary.csv --output figures --format png
```

The command creates quartet-simplex, off-arm-distance, terminal-pattern-prevalence, and MSC-versus-network AIC figures. `msrc-sim-compare` now also adds Benjamini–Hochberg adjusted off-arm q-values.

## v0.7.0: spatial profiles

Version 0.7.0 adds ordered genomic loci and rearrangement-interval quartet
profiles. The primary spatial figure is analogous to empirical chromosome-wide
quartet-support plots: local `q1(x)`, `q2(x)`, and `q3(x)` tracks are plotted
against genomic position, with the rearrangement interval shaded and structural
breakpoints marked.

The scientific target is spatial identifiability. Genome-averaged quartet
counts can be non-identifying when two mechanisms produce the same average
quartet vector. Ordered profiles can contain additional information through the
alignment of local quartet support with rearrangement breakpoints, orientation,
recombination suppression, and arrangement-state partition. Localization alone
is not claimed to uniquely identify MSRC, because introgression/network models
can also generate spatial ancestry patterns.
