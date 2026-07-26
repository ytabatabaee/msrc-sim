from __future__ import annotations
import argparse
import json
import csv
from pathlib import Path
import numpy as np
import yaml

from .config import load_config
from .species_tree import SpeciesTree
from .rearrangement import Rearrangement
from .wright_fisher import simulate_frequency_history
from .structured_coalescent import simulate_genealogy
from .io import write_outputs
from .conditional import simulate_conditional
from .analytic import TOPOLOGY_NAMES


def _run_conditional(c):
    result = simulate_conditional(c)
    out = Path(c['output']['directory'])
    out.mkdir(parents=True, exist_ok=True)
    if c['output'].get('record_resolved_config', True):
        with open(out / 'config.resolved.yaml', 'w') as f:
            yaml.safe_dump(c, f, sort_keys=False)
    with open(out / 'quartet_probabilities.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['topology', 'count', 'empirical_probability', 'exact_probability', 'absolute_error'])
        for i, topology in enumerate(TOPOLOGY_NAMES):
            w.writerow([topology, int(result.counts[i]), float(result.empirical_probabilities[i]), float(result.exact_probabilities[i]), float(result.absolute_error[i])])
    summary = {
        'version': '0.4.1',
        'mode': 'conditional',
        'configuration': ''.join(map(str, result.configuration)),
        'num_loci': int(c['num_loci']),
        'counts': result.counts.tolist(),
        'empirical_probabilities': result.empirical_probabilities.tolist(),
        'exact_probabilities': result.exact_probabilities.tolist(),
        'absolute_error': result.absolute_error.tolist(),
    }
    with open(out / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Configuration: {summary['configuration']}")
    print(f"Empirical: {summary['empirical_probabilities']}")
    print(f"Exact: {summary['exact_probabilities']}")
    print(f"Wrote {out}")


def _run_mechanistic(c):
    rng = np.random.default_rng(c['seed'])
    st = c['species_tree']
    tree = SpeciesTree(st['newick'], st['default_effective_population_size'], st['root_extension'], st.get('branch_parameters'))
    rr = c['rearrangement']
    sel = rr.get('selection', {}).get('coefficient', rr.get('selection_coefficient', 0.0))
    rearr = Rearrangement(rr.get('id', 'inv_1'), rr['type'], rr['origin_branch'], int(rr['origin_time_from_branch_start']), int(rr.get('initial_copy_count', 1)), float(sel))
    hist = simulate_frequency_history(tree, rearr, rng)
    sampled = {t: int(rng.random() < hist.terminal_frequency(t)) for t in tree.taxa}
    rec = c['recombination']
    base = rec.get('baseline_rate', rec.get('rate'))
    frac = rec.get('effective_cross_arrangement_fraction', rec.get('suppression_factor'))
    logopt = c['output']['event_log_loci']
    first_n = logopt.get('first_n', c['num_loci']) if isinstance(logopt, dict) else c['num_loci']
    results = [simulate_genealogy(i, tree, hist, sampled, float(base), float(frac), rng, c['output']['record_backward_events'] and i < first_n) for i in range(c['num_loci'])]
    out = write_outputs(c, tree, rearr, hist, sampled, results)
    print(f"Wrote {out}")


def main():
    ap = argparse.ArgumentParser(description='Run one MSRC simulation')
    ap.add_argument('--config', required=True)
    a = ap.parse_args()
    c = load_config(a.config)
    if c['mode'] == 'conditional':
        _run_conditional(c)
    else:
        _run_mechanistic(c)


if __name__ == '__main__':
    main()
