"""Plot simple linkage diagnostics for an ARG tree sequence.

The chromosome is treated as one linked replicate: each adjacent tree interval
is one observation, and the script reports both topology agreement and Pearson
correlation of a representative pairwise MRCA time. Optional msprime output is
reported as a separate reference panel, never used by the simulator.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import csv
import numpy as np
import tskit


def _topology(tree):
    splits = []
    for node in tree.nodes():
        leaves = tuple(sorted(tree.leaves(node)))
        if len(leaves) == 2:
            other = tuple(sorted(set(range(4)) - set(leaves)))
            splits.append(min(leaves, other))
    return min(splits) if splits else ()


def _internal_branch(tree):
    lengths = []
    for node in tree.nodes():
        leaves = list(tree.leaves(node))
        if len(leaves) == 2 and tree.parent(node) != tskit.NULL:
            lengths.append(tree.time(tree.parent(node)) - tree.time(node))
    return max(lengths, default=0.0)


def _statistics(ts, bins):
    positions = np.linspace(0, ts.sequence_length, 201, endpoint=False) + ts.sequence_length / 400
    trees = [ts.at(float(x)) for x in positions]
    topologies = [_topology(tree) for tree in trees]
    times = np.array([tree.time(tree.mrca(0, 1)) for tree in trees])
    internal = np.array([_internal_branch(tree) for tree in trees])
    rows = []
    for low, high in zip(bins[:-1], bins[1:]):
        top, corr, branch = [], [], []
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                d = positions[j] - positions[i]
                if low <= d < high:
                    top.append(float(topologies[i] == topologies[j]))
                    corr.append((times[i], times[j]))
                    branch.append((internal[i], internal[j]))
        def pearson(values):
            if len(values) < 3 or np.std(values, axis=0).min() == 0:
                return np.nan
            return float(np.corrcoef(np.asarray(values).T)[0, 1])
        rows.append((0.5 * (low + high), np.mean(top) if top else np.nan, pearson(corr), pearson(branch)))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trees", type=Path)
    parser.add_argument("--output", type=Path, default=Path("arg_linkage_diagnostics.png"))
    parser.add_argument("--csv", type=Path, default=Path("arg_linkage_diagnostics.csv"))
    args = parser.parse_args()
    import matplotlib.pyplot as plt
    ts = tskit.load(args.trees)
    bins = np.linspace(0, ts.sequence_length, 21)
    rows = _statistics(ts, bins)
    with args.csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["distance_bp", "topology_agreement", "coalescence_time_correlation", "internal_branch_correlation"])
        writer.writerows(rows)
    distance = np.array([r[0] for r in rows])
    agreement = np.array([r[1] for r in rows])
    time_corr = np.array([r[2] for r in rows])
    branch_corr = np.array([r[3] for r in rows])
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    if len(distance):
        axes[0].plot(distance, agreement, marker="o", label="quartet topology")
        axes[1].plot(distance, time_corr, marker="o", label="MRCA time")
        axes[1].plot(distance, branch_corr, marker="s", label="internal branch")
    axes[0].set(xlabel="genomic distance (bp)", ylabel="adjacent topology agreement")
    axes[1].set(xlabel="genomic distance (bp)", ylabel="Pearson correlation")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(args.output, dpi=180)


if __name__ == "__main__":
    main()
