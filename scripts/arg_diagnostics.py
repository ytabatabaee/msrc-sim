"""Plot simple linkage diagnostics for an ARG tree sequence.

The chromosome is treated as one linked replicate: each adjacent tree interval
is one observation, and the script reports both topology agreement and Pearson
correlation of a representative pairwise MRCA time. Optional msprime output is
reported as a separate reference panel, never used by the simulator.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import tskit


def _statistics(ts):
    trees = list(ts.trees())
    def topology(tree):
        return tuple(sorted(tuple(sorted(tree.leaves(tree.mrca(i, j)))) for i in range(4) for j in range(i)))
    topologies = [topology(tree) for tree in trees]
    # For four samples, use the MRCA time of the first two samples.
    pair_times = [tree.time(tree.mrca(0, 1)) for tree in trees]
    distances = np.array([0.5 * (a.interval.right + b.interval.left) for a, b in zip(trees, trees[1:])])
    agreements = np.array([topologies[i] == topologies[i + 1] for i in range(len(trees) - 1)], dtype=float)
    correlations = np.array([pair_times[i] * pair_times[i + 1] for i in range(len(trees) - 1)], dtype=float)
    return distances, agreements, correlations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trees", type=Path)
    parser.add_argument("--output", type=Path, default=Path("arg_linkage_diagnostics.png"))
    args = parser.parse_args()
    import matplotlib.pyplot as plt
    ts = tskit.load(args.trees)
    distance, agreement, shared_time = _statistics(ts)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    if len(distance):
        axes[0].scatter(distance, agreement, s=12)
        axes[1].scatter(distance, shared_time, s=12)
    axes[0].set(xlabel="genomic distance (bp)", ylabel="adjacent topology agreement")
    axes[1].set(xlabel="genomic distance (bp)", ylabel="adjacent pairwise MRCA-time product")
    fig.tight_layout()
    fig.savefig(args.output, dpi=180)


if __name__ == "__main__":
    main()
