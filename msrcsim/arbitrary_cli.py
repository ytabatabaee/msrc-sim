from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from .arbitrary import load_arbitrary_config, run_arbitrary_simulation


def main() -> None:
    parser = argparse.ArgumentParser(description="Run an arbitrary-tree MSRC simulation")
    parser.add_argument("--config", help="YAML arbitrary-tree simulation config")
    parser.add_argument("--species-tree", help="Newick file for an arbitrary rooted species tree")
    parser.add_argument("--output", default="arbitrary_tree_output")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--num-loci", type=int, default=20)
    parser.add_argument("--default-ne", type=int, default=100)
    parser.add_argument("--root-extension", type=float, default=200)
    parser.add_argument("--origin-branch", default="ROOT")
    parser.add_argument("--origin-time-from-branch-start", type=int, default=10)
    parser.add_argument("--initial-copy-count", type=int, default=20)
    parser.add_argument("--population-process", choices=["wright_fisher", "moran"], default=None)
    args = parser.parse_args()

    if args.config:
        config = load_arbitrary_config(args.config)
    else:
        if not args.species_tree:
            parser.error("--config or --species-tree is required")
        config = {
            "mode": "arbitrary_tree",
            "seed": args.seed,
            "num_loci": args.num_loci,
            "species_tree": {
                "path": args.species_tree,
                "root_extension": args.root_extension,
                "default_effective_population_size": args.default_ne,
            },
            "rearrangement": {
                "id": "inv_1",
                "type": "inversion",
                "origin_branch": args.origin_branch,
                "origin_time_from_branch_start": args.origin_time_from_branch_start,
                "initial_copy_count": args.initial_copy_count,
                "selection": {"model": "genic", "coefficient": 0.0},
            },
            "recombination": {
                "baseline_rate": 0.01,
                "effective_cross_arrangement_fraction": 0.05,
            },
            "output": {"directory": args.output},
            "population_process": {"model": args.population_process or "wright_fisher"},
        }
    if args.population_process:
        config.setdefault("population_process", {})
        config["population_process"]["model"] = args.population_process
    config.setdefault("output", {})
    config["output"].setdefault("directory", args.output)
    out = run_arbitrary_simulation(config)
    print(f"Wrote {Path(out)}")


if __name__ == "__main__":
    main()
