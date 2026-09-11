from __future__ import annotations

import argparse

from .config import load_config
from .pattern_probabilities import run_pattern_probability_analysis


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute exact arrangement-pattern probabilities and compare to simulation"
    )
    parser.add_argument("--config", required=True, help="Mechanistic msrc-sim YAML config")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument(
        "--monte-carlo-replicates",
        type=int,
        default=20000,
        help="Number of Wright-Fisher histories to simulate for comparison",
    )
    parser.add_argument("--seed", type=int, help="Monte Carlo seed; defaults to config seed")
    parser.add_argument(
        "--max-states",
        type=int,
        default=2001,
        help="Maximum allowed WF count states, equal to 2Ne+1",
    )
    args = parser.parse_args()
    config = load_config(args.config)
    if config["mode"] != "mechanistic":
        raise ValueError("Pattern probabilities require a mechanistic configuration")
    result = run_pattern_probability_analysis(
        config,
        args.output,
        monte_carlo_replicates=args.monte_carlo_replicates,
        seed=args.seed,
        max_states=args.max_states,
    )
    theory = result["theory"]
    print(f"Wrote {result['output']}")
    print(f"Exact DP runtime: {theory.runtime_seconds:.3f}s")
    print(f"Maximum WF state count: {theory.max_state_count}")
    print(f"Maximum absolute discrepancy: {result['max_abs_error']:.6g}")


if __name__ == "__main__":
    main()
