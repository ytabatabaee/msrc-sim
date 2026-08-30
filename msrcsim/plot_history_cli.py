from __future__ import annotations

from pathlib import Path
import argparse
import json

from .frequency_history_plot import (
    A0_COLOR,
    A1_COLOR,
    plot_frequency_history_tree,
    read_frequency_history,
    read_sampled_arrangements,
)


def _infer_format(output: Path, explicit: str | None) -> str:
    if explicit:
        return explicit
    suffix = output.suffix.lower().lstrip(".")
    if suffix in {"png", "pdf", "svg"}:
        return suffix
    return "png"


def _manifest_path(output: Path) -> Path:
    return output.with_name(f"{output.stem}.figure.json")


def _resolve_inputs(args: argparse.Namespace) -> tuple[Path, Path | None, Path | None, Path]:
    run_dir = Path(args.run_dir) if args.run_dir else None
    if run_dir:
        frequency_history = Path(args.frequency_history) if args.frequency_history else run_dir / "frequency_history.csv"
        config = Path(args.config) if args.config else run_dir / "config.resolved.yaml"
        sampled = Path(args.sampled_arrangements) if args.sampled_arrangements else run_dir / "sampled_arrangements.csv"
        output = Path(args.output) if args.output else run_dir / "wright_fisher_history.png"
    else:
        if not args.frequency_history:
            raise SystemExit("--frequency-history is required unless --run-dir is supplied")
        frequency_history = Path(args.frequency_history)
        config = Path(args.config) if args.config else None
        sampled = Path(args.sampled_arrangements) if args.sampled_arrangements else None
        output = Path(args.output) if args.output else Path("wright_fisher_history.png")
    return frequency_history, config if config and config.exists() else None, sampled if sampled and sampled.exists() else None, output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot a Wright-Fisher rearrangement-frequency history on the species tree")
    parser.add_argument("--run-dir", help="Run directory containing frequency_history.csv and optional metadata files")
    parser.add_argument("--frequency-history", help="Path to frequency_history.csv")
    parser.add_argument("--config", help="Optional config.resolved.yaml path")
    parser.add_argument("--sampled-arrangements", help="Optional sampled_arrangements.csv path")
    parser.add_argument("--output", help="Output figure path")
    parser.add_argument("--format", choices=["png", "pdf", "svg"], help="Figure format; inferred from --output suffix by default")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--glyphs-per-row", type=int, default=12)
    parser.add_argument("--max-rows-per-branch", type=int, default=30)
    parser.add_argument("--all-generations", action="store_true", help="Disable generation downsampling")
    parser.add_argument("--width-mode", choices=["constant", "ne"], default="constant")
    parser.add_argument("--tip-order", help="Comma-separated terminal branch/taxon order, for example 1,2,3,4")
    parser.add_argument("--hide-internal-labels", action="store_true")
    parser.add_argument("--show-frequency-trace", action="store_true")
    parser.add_argument("--title")
    parser.add_argument("--no-manifest", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    frequency_history, config, sampled_path, output = _resolve_inputs(args)
    if not frequency_history.exists():
        raise SystemExit(f"frequency history not found: {frequency_history}")
    fmt = _infer_format(output, args.format)
    if output.suffix.lower().lstrip(".") not in {"png", "pdf", "svg"}:
        output = output.with_suffix(f".{fmt}")
    output.parent.mkdir(parents=True, exist_ok=True)

    records = read_frequency_history(frequency_history)
    sampled = read_sampled_arrangements(sampled_path)
    tip_order = [piece.strip() for piece in args.tip_order.split(",") if piece.strip()] if args.tip_order else None
    max_rows = None if args.all_generations else args.max_rows_per_branch
    fig, _ax = plot_frequency_history_tree(
        records,
        sampled_arrangements=sampled,
        tip_order=tip_order,
        glyphs_per_row=args.glyphs_per_row,
        max_rows_per_branch=max_rows,
        width_mode=args.width_mode,
        show_internal_labels=not args.hide_internal_labels,
        show_frequency_trace=args.show_frequency_trace,
        title=args.title,
    )
    fig.savefig(output, format=fmt, dpi=args.dpi, bbox_inches="tight")
    metadata = getattr(fig, "_msrc_history_metadata", {})
    import matplotlib.pyplot as plt

    plt.close(fig)

    if not args.no_manifest:
        origin = next((r for r in records if r.get("is_origin")), None)
        manifest = {
            "frequency_history": str(frequency_history),
            "config": str(config) if config else None,
            "sampled_arrangements": str(sampled_path) if sampled_path else None,
            "figure": str(output),
            "num_records_total": len(records),
            "num_records_displayed": metadata.get("num_records_displayed"),
            "glyphs_per_row": args.glyphs_per_row,
            "max_rows_per_branch": max_rows,
            "width_mode": args.width_mode,
            "tip_order": metadata.get("tip_order", tip_order),
            "origin_branch": str(origin["branch_id"]) if origin else None,
            "origin_age": float(origin["absolute_age"]) if origin else None,
            "arrangement_colors": {"A0": A0_COLOR, "A1": A1_COLOR},
        }
        _manifest_path(output).write_text(json.dumps(manifest, indent=2))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
