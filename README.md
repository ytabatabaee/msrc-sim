# msrc-sim

`msrc-sim` is a simulator for the **Multi-Species Recombination-partitioned Coalescent (MSRC)** model. MSRC models the effect of **recombination suppression** on gene tree distributions across a species tree. When a genomic region carries partially isolated ancestry backgrounds — due to chromosomal inversions, supergenes, or any other mechanism that suppresses recombination — lineages within that region remain associated with their background rather than freely recombining. This produces spatially coherent, biased genealogical signal that can mislead standard species tree methods. `msrc-sim` simulates this process by running a forward-time model of background-frequency evolution (Wright-Fisher or Moran) and then generating gene trees through a backward-time structured coalescent conditional on that history.

## Requirements

- Python ≥ 3.9
- NumPy, SciPy, PyYAML

## Installation

```bash
git clone https://github.com/ytabatabaee/msrc-sim.git
cd msrc-sim
pip install -e .
```

To run tests:

```bash
pip install -e ".[test]"
pytest
```

## Commands

`msrc-sim` installs several command-line programs. All simulation commands read YAML configuration files; the `examples/` directory contains ready-to-run configurations for each mode.

| Command | Description |
|---|---|
| `msrc-sim --config <config.yaml>` | Run a single mechanistic or conditional simulation |
| `msrc-sim-replicates --config <config.yaml>` | Run independent evolutionary replicates |
| `msrc-sim-grid --config <config.yaml>` | Run a multidimensional parameter grid |
| `msrc-sim-arbitrary --config <config.yaml>` | Simulate gene trees on an arbitrary n-taxon species tree |
| `msrc-sim-spatial --config <config.yaml>` | Simulate spatially ordered loci along a chromosome |
| `msrc-sim-hybridization --config <config.yaml>` | Run a pulse-hybridization spatial comparator |
| `msrc-sim-pattern-probabilities --config <config.yaml>` | Compare exact quartet probabilities (theory vs. simulation) |
| `msrc-sim-plot-history --run-dir <dir>` | Plot the realized structural-frequency history |
| `msrc-sim-compare --input <csv> --output <csv>` | Compare quartet vectors across runs |

## Simulation Modes

### Mechanistic simulation

Simulates a single rearrangement/recombination-suppression history forward through a dated species tree, then simulates gene trees backward through that realized history.

```bash
msrc-sim --config examples/mechanistic_balanced.yaml
```

Key configuration sections:

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

output:
  directory: my_output
  record_gene_trees: true
  record_frequency_history: true
  make_history_plot: true
```

The species tree must be ultrametric. Internal node names (e.g., `A`, `B`, `ROOT`) are used as branch identifiers and are recommended. The `effective_cross_arrangement_fraction` (κ) controls how much recombination suppression is applied: κ = 1 means no suppression; κ ≪ 1 means strong suppression and long correlated blocks.

**Outputs** (depending on flags): `frequency_history.csv`, `sampled_arrangements.csv`, `true_gene_trees.nwk`, `coalescence_times.csv`, `summary.json`, `wright_fisher_history.png`.

### Replicate experiments

Runs many independent evolutionary replicates, each drawing its own rearrangement-frequency history and simulating loci conditional on it. Useful for estimating prevalence of persistence, asymmetric quartet distributions, and discordant-topology dominance.

```bash
msrc-sim-replicates --config examples/replicates_unconditional.yaml
```

Replicates can be unconditioned, or conditioned on features of the forward history:

```yaml
# Unconditioned
conditioning:
  mode: none

# Require the rearrangement to remain polymorphic at the root
conditioning:
  mode: persistent_polymorphism
  require_segregating_at: [ROOT]
  max_attempts: 10000

# Require a specific terminal arrangement pattern
conditioning:
  mode: terminal_pattern
  accepted_patterns: ["1010", "0101"]
  max_attempts: 100000
```

**Outputs**: `replicate_summary.csv`, `simplex_points.csv`, `prevalence_summary.json`.

### Parameter grids

Sweeps a replicate experiment over all combinations of specified parameters.

```bash
msrc-sim-grid --config examples/parameter_grid.yaml
```

```yaml
parameter_grid:
  species_tree.default_effective_population_size: [25, 50, 100]
  rearrangement.initial_copy_count: [10, 20]
  recombination.effective_cross_arrangement_fraction: [0.01, 0.05, 0.1]
```

Each grid cell is written to its own `cell_####/` directory. Grid-level results are in `parameter_grid_summary.csv`.

### Arbitrary-tree mode

Simulates gene trees on a rooted species tree with any number of taxa. Accepts either a YAML config or a direct Newick file.

```bash
msrc-sim-arbitrary --config examples/arbitrary_10taxa.yaml
msrc-sim-arbitrary --species-tree my_tree.nwk --output my_output --num-loci 100 --origin-branch ROOT
```

**Outputs**: `species_tree.nwk`, `local_gene_trees.nwk`, `arrangement_history.tsv`, `terminal_states.tsv`, `quartet_summaries.tsv`, `simulation_summary.json`.

### Conditional quartet simulation

Simulates quartet outcomes for a fixed four-lineage arrangement configuration over a single structured interval, and computes the exact conditional quartet matrix for comparison.

```bash
msrc-sim --config examples/conditional_quartet.yaml
```

```yaml
mode: conditional
seed: 12345
num_loci: 100000
structured_interval:
  duration: 1.0
  configuration: "1010"   # lineages 1 and 3 carry state 1, lineages 2 and 4 carry state 0
  migration:
    m01: 0.05
    m10: 0.05
  coalescence:
    lambda0: 1.0
    lambda1: 1.0
output:
  directory: conditional_output
```

**Outputs**: `quartet_probabilities.csv` (empirical and exact probabilities with errors), `summary.json`.

### Spatially ordered loci

Simulates ordered loci along a chromosome. Loci inside the recombination-suppressed interval use the MSRC genealogy generator; loci outside use the standard MSC. Useful for studying whether quartet support shifts are spatially localized around the suppressed region.

```bash
msrc-sim-spatial --config examples/spatial_inversion.yaml
```

**Outputs**: `spatial_loci.csv`, `spatial_windows.csv`, `spatial_summary.json`, `spatial_gene_trees.tsv`, `spatial_quartet_profile.png`.

### Pulse-hybridization comparator

Simulates a spatial pulse-hybridization model for comparison with MSRC. Since MSRC suppression and introgression can produce similar genome-wide quartet frequencies, the spatial comparison helps study whether the two models are distinguishable from their chromosomal profiles.

```bash
msrc-sim-hybridization --config examples/pulse_hybridization.yaml
```

To compare MSRC and hybridization runs side by side:

```bash
msrc-sim-plot-spatial-compare \
  --msrc-dir spatial_output \
  --hyb-dir hybridization_output \
  --output comparison/spatial_compare.png
```

## Forward Population Process

The default forward process is Wright-Fisher. Moran dynamics are also available:

```yaml
population_process:
  model: wright_fisher  # default

population_process:
  model: moran
```

Both models converge to the same diffusion limit under weak selection and are useful for robustness checks. The `effective_cross_arrangement_fraction` (κ) is the key parameter governing recombination suppression: it multiplies the baseline recombination rate inside the suppressed interval, reducing the rate at which lineages switch between ancestry backgrounds.

## Examples

The `examples/` directory contains ready-to-run configs:

- `mechanistic_balanced.yaml` — balanced quartet tree, mechanistic simulation
- `mechanistic_unbalanced.yaml` — unbalanced quartet tree
- `replicates_unconditional.yaml` — unconditioned prevalence experiment
- `replicates_conditioned.yaml` — terminal-pattern-conditioned experiment
- `parameter_grid.yaml` — multidimensional parameter sweep
- `conditional_quartet.yaml` — fixed-configuration conditional simulation
- `spatial_inversion.yaml` — ordered-locus spatial prototype
- `pulse_hybridization.yaml` — spatial hybridization comparator
- `arbitrary_10taxa.yaml`, `arbitrary_50taxa.yaml` — arbitrary-tree examples

## Citation

If you use `msrc-sim`, please cite:

> Yasamin Tabatabaee, Siavash Mirarab, and Jian Ma. *Suppressed recombination generates introgression-like signals and biases phylogenomic inference.* 2026.
