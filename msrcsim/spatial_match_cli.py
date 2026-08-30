from __future__ import annotations

from pathlib import Path
from typing import Any
import argparse
import json
import math

from .spatial_compare_io import load_spatial_model_run


def _distance(a: tuple[float, float, float], b: tuple[float, float, float], metric: str) -> float:
    diffs = [float(x) - float(y) for x, y in zip(a, b)]
    if metric == "l1":
        return float(sum(abs(x) for x in diffs))
    if metric == "l2":
        return float(math.sqrt(sum(x * x for x in diffs)))
    if metric == "linf":
        return float(max(abs(x) for x in diffs))
    raise ValueError("metric must be l1, l2, or linf")


def _candidate_dirs(grid_dir: Path) -> list[Path]:
    dirs: list[Path] = []
    if (grid_dir / "hybridization_summary.json").exists() or (grid_dir / "hyb_summary.json").exists():
        dirs.append(grid_dir)
    for path in sorted(grid_dir.rglob("*")):
        if path.is_dir() and ((path / "hybridization_summary.json").exists() or (path / "hyb_summary.json").exists()):
            dirs.append(path)
    return list(dict.fromkeys(dirs))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Find the precomputed hybridization run whose marginal quartet vector best matches an MSRC spatial run")
    parser.add_argument("--msrc-dir", required=True, help="Directory containing MSRC spatial outputs")
    parser.add_argument("--hyb-grid-dir", required=True, help="Directory containing one or more hybridization run directories")
    parser.add_argument("--metric", choices=["l1", "l2", "linf"], default="l1", help="Distance metric for quartet-vector matching")
    parser.add_argument("--output", required=True, help="Output JSON path")
    parser.add_argument("--top-n", type=int, default=10, help="Number of ranked candidates to include in the JSON")
    parser.add_argument("--window-size-loci", type=int, default=50, help="Window size if candidate windows must be recomputed")
    parser.add_argument("--step-loci", type=int, default=10, help="Window step if candidate windows must be recomputed")
    parser.add_argument("--match-threshold", type=float, default=0.05, help="Distance threshold used for the matched warning flag")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        msrc_run = load_spatial_model_run(args.msrc_dir, "msrc", window_size_loci=args.window_size_loci, step_loci=args.step_loci)
        candidates: list[dict[str, Any]] = []
        for path in _candidate_dirs(Path(args.hyb_grid_dir)):
            try:
                run = load_spatial_model_run(path, "hybridization", window_size_loci=args.window_size_loci, step_loci=args.step_loci)
            except Exception as exc:
                candidates.append({"hyb_dir": str(path), "error": str(exc)})
                continue
            dist = _distance(msrc_run.marginal_q, run.marginal_q, args.metric)
            candidates.append({
                "hyb_dir": str(path),
                "hyb_marginal_q": [float(x) for x in run.marginal_q],
                "distance": float(dist),
                "num_loci": len(run.loci),
                "num_windows": len(run.windows),
            })
        valid = [c for c in candidates if "distance" in c]
        if not valid:
            raise FileNotFoundError("No loadable hybridization runs found under --hyb-grid-dir")
        ranked = sorted(valid, key=lambda row: float(row["distance"]))
        best = ranked[0]
        result = {
            "msrc_dir": str(msrc_run.run_dir),
            "hyb_grid_dir": str(Path(args.hyb_grid_dir)),
            "metric": args.metric,
            "msrc_marginal_q": [float(x) for x in msrc_run.marginal_q],
            "best_hyb_dir": best["hyb_dir"],
            "best_hyb_marginal_q": best["hyb_marginal_q"],
            "best_distance": best["distance"],
            "match_threshold": float(args.match_threshold),
            "is_close_match": float(best["distance"]) <= float(args.match_threshold),
            "warning": None if float(best["distance"]) <= float(args.match_threshold) else "Best hybridization run is not closely matched to the MSRC marginal quartet vector.",
            "ranked_candidates": ranked[:max(1, int(args.top_n))],
            "num_candidates_scanned": len(valid),
            "num_candidates_with_errors": len(candidates) - len(valid),
        }
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2))
        if result["warning"]:
            print(result["warning"])
        print(f"Best hybridization run: {best['hyb_dir']} ({args.metric}={float(best['distance']):.6g})")
        print(f"Wrote {output}")
        return 0
    except Exception as exc:
        print(str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
