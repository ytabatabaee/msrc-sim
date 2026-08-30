from __future__ import annotations

import argparse

import yaml

from .hybridization import simulate_hybridization
from .hybridization_plot import plot_hybridization_output


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a spatial pulse-hybridization ancestry-tract simulation")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    out = simulate_hybridization(config)
    if config.get("output", {}).get("make_plots", True):
        try:
            plot_hybridization_output(out, fmt="png")
        except Exception as exc:
            print(f"Plot generation failed: {exc}")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
