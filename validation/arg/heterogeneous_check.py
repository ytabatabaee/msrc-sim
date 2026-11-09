"""Supplementary high-replicate heterogeneous-map marginal check."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from validation.arg.validate_neutral import summarize_heterogeneous


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=31_000_000)
    parser.add_argument("--output", type=Path, default=Path("validation/arg/results/heterogeneous_3000.json"))
    args = parser.parse_args(argv)
    summary = {"replicates": args.replicates, "seed": args.seed,
               "heterogeneous_map": summarize_heterogeneous(range(args.seed, args.seed + args.replicates))}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
