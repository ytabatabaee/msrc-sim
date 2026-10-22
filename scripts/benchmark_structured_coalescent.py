#!/usr/bin/env python3
"""Benchmark the reference and grouped structured-coalescent implementations."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from msrcsim import SpeciesTree  # noqa: E402
from msrcsim.population_process import simulate_population_history  # noqa: E402
from msrcsim.rearrangement import Rearrangement  # noqa: E402
from msrcsim.structured_coalescent import (  # noqa: E402
    simulate_genealogy,
    simulate_genealogy_reference,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--rho", type=float, default=0.02)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    species = config["species_tree"]
    tree = SpeciesTree(species["newick"], species["default_effective_population_size"],
                       species["root_extension"], ultrametric_tolerance=1e-6)
    rr = config["rearrangement"]
    rearrangement = Rearrangement(rr["id"], rr["type"], rr["origin_branch"],
                                  int(rr["origin_time_from_branch_start"]),
                                  int(rr["initial_copy_count"]), 0.0)
    history = simulate_population_history(tree, rearrangement, np.random.default_rng(1), config)
    sampled = {taxon: int(history.terminal_frequency(taxon) > 0.0) for taxon in tree.taxa}
    rows = []
    for name, function in (("reference", simulate_genealogy_reference), ("optimized", simulate_genealogy)):
        started = time.perf_counter()
        events = 0
        coalescences = 0
        for replicate in range(args.replicates):
            result = function(replicate, tree, history, sampled, args.rho, 1.0,
                              np.random.default_rng(1000 + replicate), record_events=True)
            events += len(result.events)
            coalescences += len(result.coalescences)
        elapsed = time.perf_counter() - started
        rows.append({"implementation": name, "replicates": args.replicates,
                     "seconds": elapsed, "seconds_per_genealogy": elapsed / args.replicates,
                     "mean_events": events / args.replicates,
                     "mean_coalescences": coalescences / args.replicates})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"config": str(args.config), "rho": args.rho,
                                       "rows": rows}, indent=2) + "\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
