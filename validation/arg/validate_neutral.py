"""Independent statistical validation of the neutral ARG against msprime."""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import msprime
import tskit

from msrcsim.arg.engine import ARGConfig, simulate_arg
from msrcsim.arg.output.tree_sequence import build_tree_sequence
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap
from msrcsim.species_tree import SpeciesTree

LENGTH = 100_000
N_E = 100
RATE = 1e-7
POSITIONS = (1_000, 2_000, 5_000, 10_000, 25_000, 50_000)


def own_ts(seed, n=2, rate=RATE, rec_map=None, tree=None):
    if tree is None:
        taxa = [f"T{i}" for i in range(n)]
        tree = SpeciesTree("(" + ",".join(f"{x}:0" for x in taxa) + ")ROOT;", N_E, 5_000)
    rec_map = rec_map or PiecewiseConstantRecombinationMap.constant(LENGTH, rate)
    # Marginal statistics must not count unary recombination-provenance
    # breakpoints. The production result remains unsimplified; validation
    # simplifies only this analysis view.
    return build_tree_sequence(simulate_arg(tree, ARGConfig(LENGTH, rec_map), np.random.default_rng(seed))).simplify()


def ref_ts(seed, n=2, rate=RATE, heterogeneous=False):
    kwargs = dict(samples=n, ploidy=1, sequence_length=LENGTH,
                  population_size=2 * N_E, random_seed=seed)
    kwargs["recombination_rate"] = (
        msprime.RateMap(position=[0, 50_000, LENGTH], rate=[1e-7, 1e-9])
        if heterogeneous else rate
    )
    return msprime.sim_ancestry(**kwargs).simplify()


def mrca_time(ts, x, a=0, b=1):
    tree = ts.at(float(x)); node = tree.mrca(a, b)
    if node == tskit.NULL:
        raise AssertionError("marginal tree is not rooted")
    return float(tree.time(node))


def quartet_topology(tree):
    """Return a root-independent unrooted quartet split."""
    samples = frozenset(tree.samples())
    if len(samples) != 4:
        raise ValueError("requires four samples")
    splits = []
    for u in tree.nodes():
        if u == tree.root:
            continue
        side = frozenset(tree.samples(u))
        if len(side) == 2:
            other = samples - side
            splits.append(tuple(sorted((tuple(sorted(side)), tuple(sorted(other))))))
    if not splits:
        raise AssertionError("unresolved quartet")
    return min(splits)


def quartet_replicate(seed, internal_length, reference):
    if reference:
        demography = msprime.Demography()
        for name in ("A", "B", "C", "D", "AB", "CD", "ROOT"):
            demography.add_population(name=name, initial_size=2 * N_E)
        # msprime requires sampled populations to be active at time zero;
        # epsilon separates sampling from the first instantaneous split.
        epsilon = 1e-6
        demography.add_population_split(time=epsilon, derived=["A", "B"], ancestral="AB")
        demography.add_population_split(time=epsilon, derived=["C", "D"], ancestral="CD")
        demography.add_population_split(time=internal_length + epsilon, derived=["AB", "CD"], ancestral="ROOT")
        ts = msprime.sim_ancestry({x: 1 for x in "ABCD"}, demography=demography,
                                  ploidy=1, sequence_length=1_000, random_seed=seed)
    else:
        tree = SpeciesTree(f"((A:0,B:0)AB:{internal_length},(C:0,D:0)CD:{internal_length})ROOT;",
                           N_E, 5_000)
        result = simulate_arg(tree, ARGConfig(1_000, PiecewiseConstantRecombinationMap.constant(1_000, 0)),
                              np.random.default_rng(seed))
        ts = build_tree_sequence(result)
    return quartet_topology(ts.first())


def mean_ci(values):
    x = np.asarray(values, dtype=float); mean = float(np.mean(x))
    se = float(np.std(x, ddof=1) / np.sqrt(len(x))) if len(x) > 1 else 0.0
    return {"mean": mean, "variance": float(np.var(x, ddof=1)) if len(x) > 1 else 0.0,
            "se": se, "ci95": [mean - 1.96 * se, mean + 1.96 * se]}


def difference_ci(a, b):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    estimate = float(np.mean(a) - np.mean(b))
    se = float(np.sqrt(np.var(a, ddof=1) / len(a) + np.var(b, ddof=1) / len(b)))
    return {"estimate": estimate, "se": se, "ci95": [estimate - 1.96 * se, estimate + 1.96 * se]}


def proportion_ci(values):
    x = np.asarray(values, dtype=float); p = float(np.mean(x)); se = np.sqrt(p * (1 - p) / len(x))
    return {"estimate": p, "se": float(se), "ci95": [max(0.0, p - 1.96 * se), min(1.0, p + 1.96 * se)]}


def correlation_ci(x, y):
    r = float(np.corrcoef(x, y)[0, 1]); z = np.arctanh(np.clip(r, -.999999, .999999))
    se = 1 / np.sqrt(max(1, len(x) - 3))
    return {"estimate": r, "ci95": [float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))]}


def paired_correlation_difference(own, ref, seed, bootstrap=400):
    own = np.asarray(own); ref = np.asarray(ref)
    observed = float(np.corrcoef(own[:, 0], own[:, 1])[0, 1] - np.corrcoef(ref[:, 0], ref[:, 1])[0, 1])
    rng = np.random.default_rng(seed); values = []
    for _ in range(bootstrap):
        own_ix = rng.integers(0, len(own), len(own)); ref_ix = rng.integers(0, len(ref), len(ref))
        values.append(np.corrcoef(own[own_ix, 0], own[own_ix, 1])[0, 1] -
                      np.corrcoef(ref[ref_ix, 0], ref[ref_ix, 1])[0, 1])
    return {"estimate": observed, "se": float(np.std(values, ddof=1)),
            "ci95": [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]}


def summarize_two_sample(seeds):
    own = [mrca_time(own_ts(s), 50_000) for s in seeds]
    ref = [mrca_time(ref_ts(s), 50_000) for s in seeds]
    return {"own": mean_ci(own), "msprime": mean_ci(ref), "difference": difference_ci(own, ref),
            "theoretical": {"mean": 2 * N_E, "variance": (2 * N_E) ** 2}}


def summarize_two_locus(seeds):
    own_pairs = {d: [] for d in POSITIONS[1:]}; ref_pairs = {d: [] for d in POSITIONS[1:]}
    own_counts = []; ref_counts = []; own_lengths = []; ref_lengths = []
    own_same = {d: [] for d in POSITIONS[1:]}; ref_same = {d: [] for d in POSITIONS[1:]}
    for s in seeds:
        a = own_ts(s); b = ref_ts(s); a0 = mrca_time(a, POSITIONS[0]); b0 = mrca_time(b, POSITIONS[0])
        for d in POSITIONS[1:]:
            at = mrca_time(a, POSITIONS[0] + d); bt = mrca_time(b, POSITIONS[0] + d)
            own_pairs[d].append((a0, at)); ref_pairs[d].append((b0, bt))
            own_same[d].append(int(a0 == at)); ref_same[d].append(int(b0 == bt))
        own_counts.append(a.num_trees); ref_counts.append(b.num_trees)
        own_lengths.extend(t.interval.right - t.interval.left for t in a.trees())
        ref_lengths.extend(t.interval.right - t.interval.left for t in b.trees())
    distances = []
    for d in POSITIONS[1:]:
        own = np.asarray(own_pairs[d]); ref = np.asarray(ref_pairs[d])
        distances.append({"distance_bp": d, "own_correlation": correlation_ci(own[:, 0], own[:, 1]),
                          "msprime_correlation": correlation_ci(ref[:, 0], ref[:, 1]),
                          "correlation_difference": paired_correlation_difference(own, ref, d),
                          "own_covariance": float(np.cov(own[:, 0], own[:, 1], ddof=1)[0, 1]),
                          "msprime_covariance": float(np.cov(ref[:, 0], ref[:, 1], ddof=1)[0, 1]),
                          "own_identical": proportion_ci(own_same[d]), "msprime_identical": proportion_ci(ref_same[d]),
                          "identical_difference": difference_ci(own_same[d], ref_same[d])})
    return {"distances": distances, "own_marginal_tree_count": mean_ci(own_counts),
            "msprime_marginal_tree_count": mean_ci(ref_counts), "own_tree_interval_lengths": mean_ci(own_lengths),
            "msprime_tree_interval_lengths": mean_ci(ref_lengths)}


def summarize_quartets(seeds):
    out = {}
    for label, internal in {"high_ILS": 10.0, "moderate_ILS": 100.0, "low_ILS": 600.0}.items():
        own = [quartet_replicate(s, internal, False) for s in seeds]
        ref = [quartet_replicate(s, internal, True) for s in seeds]
        topologies = sorted(set(own + ref)); concordant = next(t for t in topologies if (0, 1) in t)
        # This is a symmetric four-species tree: both AB and CD have an
        # ancestral branch of this length before the root. The probability
        # that neither pair coalesces is exp(-2*t), not exp(-t).
        t = internal / (2 * N_E); expected = {concordant: 1 - 2 / 3 * np.exp(-2 * t)}
        expected.update({k: 1 / 3 * np.exp(-2 * t) for k in topologies if k != concordant})
        out[label] = {"internal_length_generations": internal, "expected": {str(k): float(v) for k, v in expected.items()},
                      "own": {str(k): proportion_ci([x == k for x in own]) for k in topologies},
                      "msprime": {str(k): proportion_ci([x == k for x in ref]) for k in topologies},
                      "difference": {str(k): difference_ci([x == k for x in own], [x == k for x in ref]) for k in topologies}}
    return out


def summarize_larger_samples(seeds):
    out = {}
    for n in (8, 16):
        own_branch = []; ref_branch = []; own_counts = []; ref_counts = []; own_pairs = []; ref_pairs = []
        for s in seeds:
            a = own_ts(s, n); b = ref_ts(s, n)
            own_branch.append(a.first().total_branch_length); ref_branch.append(b.first().total_branch_length)
            own_counts.append(a.num_trees); ref_counts.append(b.num_trees)
            own_pairs.append((mrca_time(a, POSITIONS[0], 0, 1), mrca_time(a, POSITIONS[-1], 0, 1)))
            ref_pairs.append((mrca_time(b, POSITIONS[0], 0, 1), mrca_time(b, POSITIONS[-1], 0, 1)))
        out[str(n)] = {"own_total_branch_length": mean_ci(own_branch), "msprime_total_branch_length": mean_ci(ref_branch),
                       "total_branch_difference": difference_ci(own_branch, ref_branch), "own_marginal_trees": mean_ci(own_counts),
                       "msprime_marginal_trees": mean_ci(ref_counts), "correlation_difference": paired_correlation_difference(own_pairs, ref_pairs, n)}
    return out


def summarize_heterogeneous(seeds):
    own_map = PiecewiseConstantRecombinationMap(LENGTH, (50_000,), (1e-7, 1e-9))
    own_counts = []; ref_counts = []; own_times = []; ref_times = []
    for s in seeds:
        a = own_ts(s, rec_map=own_map); b = ref_ts(s, heterogeneous=True)
        own_counts.append(a.num_trees); ref_counts.append(b.num_trees)
        own_times.append(mrca_time(a, 25_000)); ref_times.append(mrca_time(b, 25_000))
    return {"own_marginal_trees": mean_ci(own_counts), "msprime_marginal_trees": mean_ci(ref_counts),
            "marginal_tree_difference": difference_ci(own_counts, ref_counts), "own_mrca": mean_ci(own_times),
            "msprime_mrca": mean_ci(ref_times), "mrca_difference": difference_ci(own_times, ref_times)}


def write_figures(summary, directory):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    d = summary["two_locus"]["distances"]; x = [r["distance_bp"] for r in d]
    plt.figure(figsize=(7, 4))
    for key, label in (("own_correlation", "msrcsim"), ("msprime_correlation", "msprime")):
        y = [r[key]["estimate"] for r in d]; lo = [r[key]["ci95"][0] for r in d]; hi = [r[key]["ci95"][1] for r in d]
        plt.plot(x, y, marker="o", label=label); plt.fill_between(x, lo, hi, alpha=.2)
    plt.xlabel("Separation (bp)"); plt.ylabel("MRCA-time correlation"); plt.legend(); plt.tight_layout()
    plt.savefig(directory / "two_locus_correlation.png", dpi=160); plt.close()
    labels = list(summary["quartets"]); fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), sharey=True)
    for ax, label in zip(axes, labels):
        q = summary["quartets"][label]; keys = list(q["own"]); xx = np.arange(len(keys))
        ax.bar(xx - .18, [q["own"][k]["estimate"] for k in keys], .36, label="msrcsim")
        ax.bar(xx + .18, [q["msprime"][k]["estimate"] for k in keys], .36, label="msprime")
        ax.set_title(label); ax.set_xticks(xx); ax.set_xticklabels([str(i + 1) for i in range(len(keys))])
    axes[0].set_ylabel("Topology frequency"); axes[-1].legend(); fig.tight_layout()
    fig.savefig(directory / "quartet_topology_frequencies.png", dpi=160); plt.close(fig)
    return True


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261108)
    parser.add_argument("--output", type=Path, default=Path("validation/arg/results/neutral_validation_1000.json"))
    parser.add_argument("--figures", type=Path, default=Path("validation/arg/results/figures"))
    args = parser.parse_args(argv)
    if args.replicates < 1000:
        raise SystemExit("primary validation requires at least 1000 independent replicates")
    seeds = [args.seed + i for i in range(args.replicates)]; started = time.perf_counter()
    summary = {"software": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
                             "msprime": msprime.__version__, "tskit": tskit.__version__}, "seed": args.seed,
               "replicates": args.replicates, "N_e": N_E, "length_bp": LENGTH, "recombination_rate": RATE}
    summary["two_sample"] = summarize_two_sample(seeds)
    summary["two_locus"] = summarize_two_locus(seeds)
    summary["quartets"] = summarize_quartets(seeds)
    summary["larger_samples"] = summarize_larger_samples(seeds)
    summary["heterogeneous_map"] = summarize_heterogeneous(seeds)
    summary["elapsed_seconds"] = time.perf_counter() - started
    args.output.parent.mkdir(parents=True, exist_ok=True); args.figures.mkdir(parents=True, exist_ok=True)
    summary["figures_generated"] = write_figures(summary, args.figures)
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
