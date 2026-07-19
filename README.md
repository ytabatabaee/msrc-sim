# msrc-sim — Milestone 2

This release contains two simulation modes for the two-state
**Multi-Species Rearrangement Coalescent (MSRC)**.

1. **Conditional quartet mode** reproduces the exact 16-state CTMC matrix
   \(H_m(t)\) and simulates quartet topologies by a Gillespie process.
2. **Mechanistic quartet mode** simulates one neutral inversion forward through
   a balanced quartet species tree using Wright–Fisher drift, samples terminal
   arrangement states, and then simulates local genealogies backward through
   the same species tree under arrangement-dependent switching and coalescence.

The second mode implements the map

\[
\Theta_R \longrightarrow \mathcal P_R
\longrightarrow Z
\longrightarrow G_1,\ldots,G_L,
\]

where \(\mathcal P_R\) is one realized inversion-frequency history shared by
all loci in the simulated inversion block.

## Install and test

```bash
pip install -e ".[test]"
pytest
```

## Conditional example

```bash
msrc-sim --config examples/conditional_quartet.yaml
```

## Mechanistic example

```bash
msrc-sim --config examples/mechanistic_quartet.yaml
```

The mechanistic output directory contains the resolved configuration, the
frequency history, sampled terminal arrangements, true gene trees, quartet
counts, and run metadata.

## Scope

This milestone supports a balanced quartet species tree with integer branch
lengths in generations. It is intentionally not yet a full linked-genome or
ARG simulator. Loci share one inversion-frequency history and one sampled
terminal arrangement pattern, but genealogies are conditionally independent.
