# msrc-sim 0.3.0

This release adds a general dated-Newick quartet engine (balanced or unbalanced), branch-specific effective population sizes, a single-origin forward Wright–Fisher inversion process, backward MSRC genealogy simulation, full optional event logging, explicit coalescence times, and expanded frequency-history metadata.

## Install

```bash
pip install -e ".[test]"
pytest
```

## Run

```bash
msrc-sim --config examples/mechanistic_balanced.yaml
msrc-sim --config examples/mechanistic_unbalanced.yaml
```

## Main outputs

- `config.resolved.yaml`
- `frequency_history.csv` with status, A0/A1 counts, frequencies, origin/start/end indicators, selection coefficient, and population size
- `sampled_arrangements.csv`
- `true_gene_trees.nwk`
- `coalescence_times.csv`
- `coalescence_events.csv`
- optional `genealogy_events.csv`
- `summary.json`

The branch identifier is the child node label for each Newick edge. The root label identifies the ancestral root-extension population. `origin_time_from_branch_start` is measured forward from the older end of the specified branch.
