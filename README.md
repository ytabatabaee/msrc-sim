# msrc-sim

Milestone 1 implementation of the conditional quartet simulator for the two-state rearrangement-structured coalescent (MSRC).

## Model

Four lineages coexist in one ancestral population for a structured interval of duration `t`.

- Each lineage switches arrangement background: `0 -> 1` at rate `m01`, and `1 -> 0` at rate `m10`.
- Two lineages coalesce only when they occupy the same arrangement.
- Same-arrangement pairwise coalescence rates are `lambda0` and `lambda1`.
- If no first coalescence occurs before time `t`, the unresolved quartet enters an exchangeable ancestor and each topology receives probability `1/3`.

Rows of the exact matrix `H_m(t)` correspond to initial configurations `0000,...,1111`. Columns are `12|34`, `13|24`, and `14|23`.

## Install and test

```bash
pip install -e ".[test]"
pytest
```

## Run

```bash
msrc-sim --config examples/conditional_quartet.yaml
```

or

```bash
python scripts/simulate_msrc.py --config examples/conditional_quartet.yaml
```

For symmetric rates, the exact asymmetry from `1010` is

\[
P(13|24\mid1010)-P(14|23\mid1010)
=
\frac{\lambda}{\lambda+2m}\left(1-e^{-2(\lambda+2m)t}\right).
\]
