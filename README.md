# msrc-sim

`msrc-sim` simulates quartet gene-tree distributions under the Multi-Species
Rearrangement Coalescent (MSRC) model. It is designed for experiments where a
chromosomal rearrangement, such as an inversion, arises in a species tree,
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
- compare quartet probability vectors with MSC arms and two-tree mixture fits.

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

The package installs four command-line programs. The simulation commands read
YAML configuration files.

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
msrc-sim-compare --input <replicate_summary.csv> --output <model_comparison.csv>
```

Fits MSC and two-tree quartet-mixture models to accepted quartet-count rows
from `replicate_summary.csv` or an equivalent table.

Equivalent script wrappers are provided in `scripts/`:

```bash
python scripts/simulate_msrc.py --config examples/mechanistic_balanced.yaml
python scripts/simulate_replicates.py --config examples/replicates_unconditional.yaml
python scripts/simulate_parameter_grid.py --config examples/parameter_grid.yaml
python scripts/compare_quartet_models.py --input replicate_output/replicate_summary.csv --output model_comparison.csv
```

## Simulation Modes

### Mechanistic MSRC Simulation

Mechanistic mode simulates a rearrangement history forward through a dated
quartet species tree and then simulates locus genealogies backward through that
realized history. Backward genealogy simulation is piecewise exact with respect
to the generation-by-generation Wright-Fisher frequency path, so a proposed
Gillespie event cannot cross a frequency-change boundary while retaining
outdated rates.

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

### Quartet Model Comparison

The `msrc-sim-compare` command fits model summaries to an existing table of
quartet counts or probabilities. It accepts `replicate_summary.csv` from a
replicate experiment, or any CSV with either `n1`, `n2`, and `n3` topology
counts or `num_loci` plus `q1`, `q2`, and `q3` probabilities.

```bash
msrc-sim-compare \
  --input replicate_output/replicate_summary.csv \
  --output model_comparison.csv
```

The comparison reports the nearest MSC arm, an off-arm contrast with standard
error, z-score, p-value, and 95% confidence interval, maximum-likelihood MSC
fits, and a quartet-level two-tree mixture fit. It is intended to distinguish a
strong alternative-tree signal from a genuine off-arm signal.

For a vector whose best MSC topology is `T2`, the off-arm contrast compares the
two minor probabilities, `q1 - q3`. The two-tree mixture is a quartet-level
model-comparison device; its fitted parameters should not be interpreted as
uniquely identifiable demographic estimates from one quartet.

## Outputs

### Mechanistic Outputs

Mechanistic runs write files to `output.directory`. Depending on the output
flags, the directory can contain:

- `config.resolved.yaml`: YAML configuration after defaults are applied;
- `frequency_history.csv`: forward Wright-Fisher frequency path on each branch;
- `sampled_arrangements.csv`: sampled terminal arrangement state for each taxon;
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

### Replicate Outputs

Replicate experiments write:

- `replicate_summary.csv`: one row per attempted history, including rejected
  attempts;
- `simplex_points.csv`: accepted replicate quartet probabilities and simplex
  coordinates, with selected MSC-distance, off-arm, and model-classification
  fields;
- `prevalence_summary.json`: acceptance rate, terminal pattern counts,
  prevalence estimates, and 95% Wilson intervals;
- `config.resolved.yaml`: the experiment configuration.

Accepted rows in `replicate_summary.csv` include topology counts (`n1`, `n2`,
`n3`), quartet probabilities (`q1`, `q2`, `q3`), distance to the nearest MSC
arm, off-arm statistics, best MSC and two-tree mixture fits, and a model
classification.

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

### Model-Comparison Outputs

`msrc-sim-compare` writes one CSV row per accepted input row. The output keeps
the input columns and appends nearest-MSC-arm fields, off-arm confidence
statistics, best MSC fit fields, best two-tree mixture fit fields,
`delta_aic_network_vs_msc`, `network_loglik_gain`, and
`model_classification`.

## Configuration Reference

Common fields:

- `mode`: one of `mechanistic`, `conditional`, `replicate_experiment`, or
  `parameter_grid`;
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

## Examples

The `examples/` directory contains ready-to-run configurations:

- `mechanistic_balanced.yaml`: mechanistic simulation on a balanced quartet tree;
- `mechanistic_unbalanced.yaml`: mechanistic simulation on an unbalanced quartet
  tree;
- `conditional_quartet.yaml`: fixed-configuration conditional quartet
  simulation;
- `replicates_unconditional.yaml`: unconditioned prevalence experiment;
- `replicates_conditioned.yaml`: terminal-pattern-conditioned prevalence
  experiment;
- `parameter_grid.yaml`: multidimensional parameter grid.
