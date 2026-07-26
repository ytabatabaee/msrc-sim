# msrc-sim 0.4.0

This release adds the prevalence-analysis milestone for the two-state Multi-Species Rearrangement Coalescent (MSRC).

## Existing functionality

- exact conditional quartet matrix `H_m(t)`;
- balanced or unbalanced dated-Newick quartet species trees;
- a single-origin Wright–Fisher inversion process;
- backward arrangement-structured genealogy simulation;
- coalescence and optional full backward-event logs.

## New functionality

- independent evolutionary replicates, each with a new inversion-frequency history;
- unconditional, persistent-polymorphism, terminal-pattern, and combined conditioning;
- explicit attempted/accepted replicate counts and conditioning acceptance rates;
- prevalence estimates with 95% Wilson intervals;
- raw quartet-simplex points;
- reproducible multidimensional parameter grids.

The independent unit for prevalence is an **evolutionary replicate**, not a locus. Multiple loci within a replicate estimate the quartet distribution conditional on one realized inversion history.

## Install and test

```bash
pip install -e ".[test]"
pytest
```

## Single mechanistic run

```bash
msrc-sim --config examples/mechanistic_balanced.yaml
```

## Unconditional replicate experiment

```bash
msrc-sim-replicates --config examples/replicates_unconditional.yaml
```

## Pattern-conditioned experiment

```bash
msrc-sim-replicates --config examples/replicates_conditioned.yaml
```

## Parameter grid

```bash
msrc-sim-grid --config examples/parameter_grid.yaml
```

## Replicate outputs

- `replicate_summary.csv`: one row per attempted history, including rejected attempts;
- `simplex_points.csv`: accepted quartet probability vectors and simplex coordinates;
- `prevalence_summary.json`: persistence, 2:2 sorting, asymmetry, and discordant-dominance estimates;
- `config.resolved.yaml`.

## Conditioning

```yaml
conditioning:
  mode: none
```

```yaml
conditioning:
  mode: persistent_polymorphism
  require_segregating_at: [ROOT]
  max_attempts: 10000
```

```yaml
conditioning:
  mode: terminal_pattern
  accepted_patterns: ["1010", "0101"]
  max_attempts: 100000
```

The simulator always reports the number of attempted and accepted histories, so conditioning does not hide rarity.


## Backward-compatible single-run CLI

The `msrc-sim` command supports both standalone modes:

```bash
msrc-sim --config examples/conditional_quartet.yaml
msrc-sim --config examples/mechanistic_balanced.yaml
```

Use `mode: conditional` for a fixed four-lineage arrangement configuration and `mode: mechanistic` for a forward Wright--Fisher history followed by backward structured coalescence.
