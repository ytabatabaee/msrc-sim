from __future__ import annotations
import argparse
import yaml
from .experiments import simulate_parameter_grid


def main() -> None:
    parser = argparse.ArgumentParser(description="Run an MSRC prevalence parameter grid")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    rows, out = simulate_parameter_grid(config)
    print(f"Completed {len(rows)} parameter cells")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
