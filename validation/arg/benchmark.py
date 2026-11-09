"""Runtime and memory benchmark for the repository ARG and msprime reference."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import time
import tracemalloc
import numpy as np
import msprime

from msrcsim.arg.engine import ARGConfig, simulate_arg
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap
from msrcsim.species_tree import SpeciesTree


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=Path("validation/arg/results/benchmark.csv"))
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--seed", type=int, default=20261108)
    args = ap.parse_args(argv)
    rows = []
    for n in (4, 8, 16):
        taxa = [f"T{i}" for i in range(n)]
        tree = SpeciesTree("(" + ",".join(f"{x}:0" for x in taxa) + ")ROOT;", 100, 5000)
        for length, rate in ((10_000, 1e-8), (100_000, 1e-8), (100_000, 1e-7)):
            for replicate in range(args.replicates):
                seed = args.seed + 1000 * n + replicate
                tracemalloc.start()
                start = time.perf_counter()
                result = simulate_arg(tree, ARGConfig(length, PiecewiseConstantRecombinationMap.constant(length, rate)), np.random.default_rng(seed))
                simulation_time = time.perf_counter() - start
                sim_peak = tracemalloc.get_traced_memory()[1]
                start = time.perf_counter()
                from msrcsim.arg.output.tree_sequence import build_tree_sequence
                ts_own = build_tree_sequence(result)
                output_time = time.perf_counter() - start
                _, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                rows.append({"engine": "msrcsim", "taxa": n, "length_bp": length, "rate": rate,
                             "replicate": replicate, "simulation_seconds": simulation_time,
                             "tree_sequence_seconds": output_time, "runtime_seconds": simulation_time + output_time,
                             "peak_memory_bytes": peak, "events": len(result.events),
                             "recombination_events": sum(e.event_type == "recombination" for e in result.events),
                             "coalescence_events": sum(e.event_type == "coalescence" for e in result.events),
                             "marginal_trees": ts_own.num_trees, "peak_lineages": result.peak_lineages,
                             "peak_segments": result.peak_segments})
                start = time.perf_counter()
                ts = msprime.sim_ancestry(samples=n, ploidy=1, sequence_length=length,
                                          recombination_rate=rate, population_size=200,
                                          random_seed=seed)
                reference_time = time.perf_counter() - start
                rows.append({"engine": "msprime", "taxa": n, "length_bp": length, "rate": rate,
                             "replicate": replicate, "simulation_seconds": reference_time,
                             "tree_sequence_seconds": 0.0, "runtime_seconds": reference_time,
                             "peak_memory_bytes": "", "events": "", "recombination_events": "",
                             "coalescence_events": "", "marginal_trees": ts.num_trees,
                             "peak_lineages": "", "peak_segments": ""})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    main()
