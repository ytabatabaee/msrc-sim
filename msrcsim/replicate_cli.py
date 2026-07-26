from __future__ import annotations
import argparse
import yaml
from .experiments import simulate_replicates, write_replicate_outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Run independent MSRC evolutionary replicates")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    records, summary = simulate_replicates(config)
    out = write_replicate_outputs(config, records, summary)
    print(f"Wrote {out}")
    print(f"Accepted {summary['accepted_replicates']} of {summary['attempted_replicates']} attempts")


if __name__ == "__main__":
    main()
