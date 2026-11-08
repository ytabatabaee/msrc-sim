"""Small reproducible neutral ARG checks against msprime.

msprime is imported only by this validation program; it is not used by
``msrcsim.arg``. The reported standard errors are replicate-level, so intervals
from one chromosome are not incorrectly treated as independent replicates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import msprime
import tskit

from msrcsim.arg.engine import ARGConfig, simulate_arg
from msrcsim.arg.output.tree_sequence import build_tree_sequence
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap
from msrcsim.species_tree import SpeciesTree


def _topology(tree):
    splits = []
    for node in tree.nodes():
        leaves = tuple(sorted(tree.samples(node)))
        if len(leaves) == 2:
            other = tuple(sorted(set(range(4)) - set(leaves)))
            splits.append(min(leaves, other))
    return min(splits)


def _own_two_taxon(seed, recombination_rate):
    tree = SpeciesTree("(A:0,B:0)ROOT;", 100, 5000)
    result = simulate_arg(tree, ARGConfig(100_000, PiecewiseConstantRecombinationMap.constant(100_000, recombination_rate)), np.random.default_rng(seed))
    ts = build_tree_sequence(result)
    values = [t.time(t.mrca(0, 1)) for t in ts.trees()]
    return float(np.mean(values)), ts.num_trees


def _own_two_series(seed, rec_map):
    tree = SpeciesTree("(A:0,B:0)ROOT;", 100, 5000)
    result = simulate_arg(tree, ARGConfig(100_000, rec_map), np.random.default_rng(seed))
    ts = build_tree_sequence(result)
    return np.array([ts.at(float(x)).time(ts.at(float(x)).mrca(0, 1)) for x in (25_000, 50_000, 75_000)]), ts.num_trees


def _msprime_two_taxon(seed, recombination_rate):
    # With ploidy=1 msprime uses rate 1/N; N=200 therefore matches this
    # simulator's diploid neutral rate 1/(2 * 100).
    ts = msprime.sim_ancestry(2, sequence_length=100_000, recombination_rate=recombination_rate,
                              population_size=200, ploidy=1, random_seed=seed)
    values = [t.time(t.mrca(0, 1)) for t in ts.trees()]
    return float(np.mean(values)), ts.num_trees


def _msprime_two_series(seed, rec_map):
    ts = msprime.sim_ancestry(2, sequence_length=100_000, recombination_rate=rec_map,
                              population_size=200, ploidy=1, random_seed=seed)
    return np.array([ts.at(float(x)).time(ts.at(float(x)).mrca(0, 1)) for x in (25_000, 50_000, 75_000)]), ts.num_trees


def _species_tree_replicate(seed):
    tree = SpeciesTree("((A:10,B:10)AB:50,(C:10,D:10)CD:50)ROOT;", 100, 500)
    result = simulate_arg(tree, ARGConfig(1_000, PiecewiseConstantRecombinationMap.constant(1_000, 0)), np.random.default_rng(seed))
    return _topology(build_tree_sequence(result).first())


def _msprime_species_replicate(seed):
    demography = msprime.Demography()
    for name in ("A", "B", "C", "D", "AB", "CD", "ROOT"):
        demography.add_population(name=name, initial_size=200)
    demography.add_population_split(time=10, derived=["A", "B"], ancestral="AB")
    demography.add_population_split(time=10, derived=["C", "D"], ancestral="CD")
    demography.add_population_split(time=60, derived=["AB", "CD"], ancestral="ROOT")
    ts = msprime.sim_ancestry({x: 1 for x in "ABCD"}, demography=demography, ploidy=1, sequence_length=1_000, random_seed=seed)
    return _topology(ts.first())


def _total_branch_replicate(seed, n, reference=False):
    if reference:
        ts = msprime.sim_ancestry(n, sequence_length=1_000, population_size=200, ploidy=1, random_seed=seed)
    else:
        taxa = [f"T{i}" for i in range(n)]
        tree = SpeciesTree("(" + ",".join(f"{x}:0" for x in taxa) + ")ROOT;", 100, 5000)
        result = simulate_arg(tree, ARGConfig(1_000, PiecewiseConstantRecombinationMap.constant(1_000, 0)), np.random.default_rng(seed))
        ts = build_tree_sequence(result)
    return ts.first().total_branch_length


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261108)
    parser.add_argument("--output", type=Path, default=Path("validation/arg/neutral_validation.json"))
    args = parser.parse_args()
    own_times, ref_times, own_trees, ref_trees = [], [], [], []
    own_series, ref_series = [], []
    own_map_trees, ref_map_trees = [], []
    own_branch, ref_branch = {}, {}
    map_own = PiecewiseConstantRecombinationMap(100_000, (50_000,), (1e-8, 1e-9))
    map_ref = msprime.RateMap(position=[0, 50_000, 100_000], rate=[1e-8, 1e-9])
    for i in range(args.replicates):
        seed = args.seed + i
        a, b = _own_two_taxon(seed, 1e-8)
        c, d = _msprime_two_taxon(seed, 1e-8)
        own_times.append(a); ref_times.append(c); own_trees.append(b); ref_trees.append(d)
        a, b = _own_two_series(seed, PiecewiseConstantRecombinationMap.constant(100_000, 1e-8))
        c, d = _msprime_two_series(seed, 1e-8)
        own_series.append(a); ref_series.append(c)
        a, b = _own_two_series(seed, map_own)
        c, d = _msprime_two_series(seed, map_ref)
        own_map_trees.append(b); ref_map_trees.append(d)
    own_top = [_species_tree_replicate(args.seed + i) for i in range(args.replicates)]
    ref_top = [_msprime_species_replicate(args.seed + i) for i in range(args.replicates)]
    for n in (4, 8, 16):
        own_branch[n] = [_total_branch_replicate(args.seed + i, n) for i in range(args.replicates)]
        ref_branch[n] = [_total_branch_replicate(args.seed + i, n, True) for i in range(args.replicates)]
    # The two terminal ancestral branches each have length 50. Here ``t`` in
    # the quartet formula is their combined no-coalescence length: 100/(2N).
    expected = 1 - 2 / 3 * np.exp(-100 / 200)
    summary = {
        "replicates": args.replicates,
        "two_sample": {"own_mean_mrca": float(np.mean(own_times)), "msprime_mean_mrca": float(np.mean(ref_times)),
                        "theoretical_mean": 200.0, "own_mean_trees": float(np.mean(own_trees)), "msprime_mean_trees": float(np.mean(ref_trees))},
        "two_locus": {"own_corr_25k_75k": float(np.corrcoef(np.asarray(own_series).T)[0, 2]),
                       "msprime_corr_25k_75k": float(np.corrcoef(np.asarray(ref_series).T)[0, 2]),
                       "own_map_mean_trees": float(np.mean(own_map_trees)), "msprime_map_mean_trees": float(np.mean(ref_map_trees))},
        "species_tree": {"expected_concordant": float(expected),
                          "own_concordant": float(np.mean([x == (0, 1) for x in own_top])),
                          "msprime_concordant": float(np.mean([x == (0, 1) for x in ref_top])),
                          "own_total_branch_means": {str(n): float(np.mean(own_branch[n])) for n in own_branch},
                          "msprime_total_branch_means": {str(n): float(np.mean(ref_branch[n])) for n in ref_branch}},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
