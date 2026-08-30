from __future__ import annotations

import argparse

import yaml

from .spatial import simulate_spatial
from .spatial_plotting import plot_spatial_output


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a spatially ordered MSRC/MSC chromosome simulation")
    parser.add_argument("--config", required=True)
    parser.add_argument("--history", default=None, help="Override config history.frozen_history")
    args = parser.parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    if args.history:
        config.setdefault("history", {})["frozen_history"] = args.history
    out = simulate_spatial(config)
    if config.get("output", {}).get("make_plots", True):
        try:
            plot_spatial_output(out, fmt="png")
        except Exception as exc:
            print(f"Plot generation failed: {exc}")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
