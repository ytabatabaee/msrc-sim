from __future__ import annotations

from pathlib import Path
from typing import Iterable
import argparse
import json
import sys

from .spatial_compare import marginal_l1_distance
from .spatial_compare_io import load_spatial_model_run
from .spatial_compare_plot import make_spatial_compare_figure, make_spatial_overlay_figure


SUPPORTED_FORMATS = {"png", "pdf", "svg"}


def _parse_formats(output: Path, requested: str | None) -> list[str]:
    if requested:
        formats = [piece.strip().lower().lstrip(".") for piece in requested.split(",") if piece.strip()]
    else:
        suffix = output.suffix.lower().lstrip(".")
        formats = [suffix or "png"]
    bad = [fmt for fmt in formats if fmt not in SUPPORTED_FORMATS]
    if bad:
        raise ValueError(f"unsupported output format(s): {', '.join(bad)}")
    return list(dict.fromkeys(formats))


def _outputs_for_formats(output: Path, formats: Iterable[str]) -> list[Path]:
    base = output.with_suffix("") if output.suffix else output
    return [base.with_suffix(f".{fmt}") for fmt in formats]


def _summary_path(output: Path) -> Path:
    base = output.with_suffix("") if output.suffix else output
    return base.parent / f"{base.name}_summary.json"


def _overlay_path(output: Path, fmt: str) -> Path:
    base = output.with_suffix("") if output.suffix else output
    return base.parent / f"{base.name}_overlay.{fmt}"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot matched spatial comparison panels for MSRC and pulse-hybridization runs")
    parser.add_argument("--msrc-dir", required=True, help="Directory containing MSRC spatial outputs")
    parser.add_argument("--hyb-dir", help="Directory containing pulse-hybridization outputs")
    parser.add_argument("--matched-hybridization-json", help="JSON written by msrc-sim-find-matched-hybridization; supplies best_hyb_dir when --hyb-dir is omitted")
    parser.add_argument("--output", required=True, help="Output figure path, ending in .png, .pdf, or .svg")
    parser.add_argument("--formats", help="Comma-separated extra/requested formats: png,pdf,svg. Defaults to the --output suffix.")
    parser.add_argument("--window-size-loci", type=int, default=50, help="Window size when windows must be recomputed from loci")
    parser.add_argument("--step-loci", type=int, default=10, help="Window step when windows must be recomputed from loci")
    parser.add_argument("--title", help="Optional figure title")
    parser.add_argument("--marginal-match-threshold", type=float, default=0.05, help="L1 threshold used for mismatch warning")
    parser.add_argument("--overlay", action=argparse.BooleanOptionalAction, default=True, help="Also write direct overlay figure")
    parser.add_argument("--focal-topology", type=int, choices=[1, 2, 3], help="Topology index for overlay support track, using 1-based q numbering")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    output = Path(args.output)
    try:
        formats = _parse_formats(output, args.formats)
        hyb_dir = args.hyb_dir
        if hyb_dir is None and args.matched_hybridization_json:
            with Path(args.matched_hybridization_json).open() as handle:
                matched = json.load(handle)
            hyb_dir = matched.get("best_hyb_dir")
        if hyb_dir is None:
            raise ValueError("--hyb-dir is required unless --matched-hybridization-json supplies best_hyb_dir")
        msrc_run = load_spatial_model_run(args.msrc_dir, "msrc", window_size_loci=args.window_size_loci, step_loci=args.step_loci)
        hyb_run = load_spatial_model_run(hyb_dir, "hybridization", window_size_loci=args.window_size_loci, step_loci=args.step_loci)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig = make_spatial_compare_figure(msrc_run, hyb_run, title=args.title)
        figure_paths = _outputs_for_formats(output, formats)
        for path in figure_paths:
            fig.savefig(path, dpi=300, bbox_inches="tight")
        import matplotlib.pyplot as plt

        plt.close(fig)
        overlay_paths: list[Path] = []
        if args.overlay:
            overlay_fig = make_spatial_overlay_figure(msrc_run, hyb_run, focal_topology=(args.focal_topology - 1) if args.focal_topology else None)
            for fmt in formats:
                path = _overlay_path(output, fmt)
                overlay_fig.savefig(path, dpi=300, bbox_inches="tight")
                overlay_paths.append(path)
            plt.close(overlay_fig)
        l1 = marginal_l1_distance(msrc_run.marginal_q, hyb_run.marginal_q)
        summary = {
            "msrc_dir": str(msrc_run.run_dir),
            "hyb_dir": str(hyb_run.run_dir),
            "chromosome_length_bp": max(msrc_run.chromosome_length_bp, hyb_run.chromosome_length_bp),
            "msrc_marginal_q": [float(x) for x in msrc_run.marginal_q],
            "hyb_marginal_q": [float(x) for x in hyb_run.marginal_q],
            "marginal_q_l1_difference": float(l1),
            "marginal_match_threshold": float(args.marginal_match_threshold),
            "warning": "Spatial comparison generated, but marginal quartet vectors are not closely matched." if l1 > args.marginal_match_threshold else None,
            "msrc_num_windows": len(msrc_run.windows),
            "hyb_num_windows": len(hyb_run.windows),
            "msrc_num_loci": len(msrc_run.loci),
            "hyb_num_loci": len(hyb_run.loci),
            "has_overlay": bool(overlay_paths),
            "figures": [str(path) for path in figure_paths],
            "overlay_figures": [str(path) for path in overlay_paths],
        }
        summary_file = _summary_path(output)
        summary_file.write_text(json.dumps(summary, indent=2))
        if summary["warning"]:
            print(summary["warning"], file=sys.stderr)
        print(f"Wrote {', '.join(str(path) for path in figure_paths)}")
        print(f"Wrote {summary_file}")
        return 0
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
