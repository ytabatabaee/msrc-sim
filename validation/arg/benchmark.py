"""Small benchmark for the repository ARG and optional msprime reference."""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=Path("validation/arg/benchmark.csv"))
    args = ap.parse_args()
    rows = []
    for n in (4, 8, 16):
        taxa = [f"T{i}" for i in range(n)]
        tree = SpeciesTree("(" + ",".join(f"{x}:0" for x in taxa) + ")ROOT;", 100, 5000)
        for length, rate in ((10_000, 1e-8), (100_000, 1e-8), (100_000, 1e-7)):
            tracemalloc.start()
            start = time.perf_counter()
            result = simulate_arg(tree, ARGConfig(length, PiecewiseConstantRecombinationMap.constant(length, rate)), np.random.default_rng(42))
            runtime = time.perf_counter() - start
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            rows.append({"engine": "msrcsim", "taxa": n, "length_bp": length, "rate": rate,
                         "runtime_seconds": runtime, "peak_memory_bytes": peak,
                         "events": len(result.events), "peak_lineages": result.peak_lineages})
            start = time.perf_counter()
            ts = msprime.sim_ancestry(n, sequence_length=length, recombination_rate=rate,
                                      population_size=100, random_seed=42)
            rows.append({"engine": "msprime", "taxa": n, "length_bp": length, "rate": rate,
                         "runtime_seconds": time.perf_counter() - start, "peak_memory_bytes": "",
                         "events": ts.num_trees, "peak_lineages": ""})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    main()
