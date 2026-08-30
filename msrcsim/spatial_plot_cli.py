from __future__ import annotations

import argparse

from .spatial_plotting import plot_spatial_output


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot spatial quartet profile outputs")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default=None)
    parser.add_argument("--format", choices=["png", "pdf", "svg"], default="png")
    args = parser.parse_args()
    figures = plot_spatial_output(args.input, args.output, args.format)
    print(f"Wrote {', '.join(str(p) for p in figures)}")


if __name__ == "__main__":
    main()

