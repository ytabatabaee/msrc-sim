from __future__ import annotations

from pathlib import Path
from typing import Any
import argparse
import csv
import json
import math

import numpy as np
import yaml
from scipy.optimize import minimize

from .analytic import TOPOLOGY_NAMES
from .hybridization import simulate_hybridization
from .spatial_compare import SpatialModelRun
from .spatial_compare_io import compute_windows_from_loci, load_spatial_model_run
from .spatial_compare_plot import QUARTET_COLORS, HYB_INTRO, HYB_MAJOR, MSRC_FEATURE, configure_matplotlib_cache
from .spatial_features import extract_spatial_features


def _parse_floats(text: str) -> list[float]:
    return [float(piece.strip()) for piece in text.split(",") if piece.strip()]


def _load_match(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict) or "hybridization" not in data:
        raise ValueError("matched parameters YAML must contain a hybridization section")
    return data


def _topology_counts(loci: list[dict[str, Any]]) -> np.ndarray:
    return np.bincount(np.asarray([int(row["topology_index"]) for row in loci], dtype=int), minlength=3)[:3]


def _bootstrap_msrc_run(template: SpatialModelRun, rng: np.random.Generator, window_size_loci: int, step_loci: int, replicate: int) -> SpatialModelRun:
    inside = [row for row in template.loci if bool(row.get("is_inside_rearranged_interval", False))]
    outside = [row for row in template.loci if not bool(row.get("is_inside_rearranged_interval", False))]
    q_inside = _topology_counts(inside) / max(1, len(inside)) if inside else np.asarray(template.marginal_q, dtype=float)
    q_outside = _topology_counts(outside) / max(1, len(outside)) if outside else np.asarray(template.marginal_q, dtype=float)
    loci = []
    for row in template.loci:
        new_row = dict(row)
        q = q_inside if bool(row.get("is_inside_rearranged_interval", False)) else q_outside
        top = int(rng.choice(3, p=q / q.sum()))
        new_row["topology_index"] = top
        new_row["topology_label"] = TOPOLOGY_NAMES[top]
        loci.append(new_row)
    windows = compute_windows_from_loci(loci, "msrc", window_size_loci, step_loci)
    counts = _topology_counts(loci)
    marginal = counts / counts.sum()
    summary = dict(template.summary)
    summary.update({"overall_q1": float(marginal[0]), "overall_q2": float(marginal[1]), "overall_q3": float(marginal[2]), "bootstrap_replicate": replicate})
    return SpatialModelRun("MSRC", template.chromosome_length_bp, loci, windows, summary, (float(marginal[0]), float(marginal[1]), float(marginal[2])), template.feature_intervals, template.run_dir)


def _hybridization_config(match: dict[str, Any], h: float, r: float, seed: int, output_dir: Path, num_loci: int, window_size_loci: int, step_loci: int, lags: list[float]) -> dict[str, Any]:
    hyb = match["hybridization"]
    template = match.get("simulation_template", {})
    chromosome = dict(template.get("chromosome", {}))
    chromosome.setdefault("length_bp", 100000000.0)
    chromosome["num_loci"] = int(num_loci)
    chromosome.setdefault("locus_positions", {"mode": "evenly_spaced"})
    return {
        "mode": "pulse_hybridization",
        "seed": int(seed),
        "chromosome": chromosome,
        "hybridization": {
            "gamma": float(hyb["gamma"]),
            "generations_since_pulse": float(h),
            "major": {"topology": int(hyb["major_topology"]), "internal_branch_length": float(hyb["t_major"])},
            "introgressed": {"topology": int(hyb["introgressed_topology"]), "internal_branch_length": float(hyb["t_introgressed"])},
            "donor": hyb.get("donor", "donor"),
            "recipient": hyb.get("recipient", "recipient"),
        },
        "recombination": {"rate_per_bp_per_generation": float(r)},
        "windows": {"loci_per_window": int(window_size_loci), "step_loci": int(step_loci)},
        "spatial_statistics": {"lags_bp": [float(x) for x in lags]},
        "output": {
            "directory": str(output_dir),
            "record_tracts": True,
            "record_loci": True,
            "record_windows": True,
            "make_plots": False,
        },
    }


def _auc(y_true: np.ndarray, scores: np.ndarray) -> float:
    pos = scores[y_true == 1]
    neg = scores[y_true == 0]
    if len(pos) == 0 or len(neg) == 0:
        return math.nan
    wins = 0.0
    for p in pos:
        wins += float(np.sum(p > neg)) + 0.5 * float(np.sum(p == neg))
    return float(wins / (len(pos) * len(neg)))


def _sigmoid(values: np.ndarray) -> np.ndarray:
    out = np.empty_like(values, dtype=float)
    positive = values >= 0.0
    out[positive] = 1.0 / (1.0 + np.exp(-values[positive]))
    exp_values = np.exp(values[~positive])
    out[~positive] = exp_values / (1.0 + exp_values)
    return out


def _ci95(values: np.ndarray) -> tuple[float, float]:
    if len(values) == 0:
        return (math.nan, math.nan)
    if len(values) == 1:
        value = float(values[0])
        return (value, value)
    mean = float(np.mean(values))
    half_width = float(1.96 * np.std(values, ddof=1) / math.sqrt(len(values)))
    return (max(0.0, mean - half_width), min(1.0, mean + half_width))


def _fit_logistic(train_x: np.ndarray, train_y: np.ndarray) -> np.ndarray:
    x = np.column_stack([np.ones(train_x.shape[0]), train_x])

    def objective(beta: np.ndarray) -> float:
        z = x @ beta
        return float(np.sum(np.logaddexp(0.0, z) - train_y * z) + 0.5 * 1e-4 * np.dot(beta[1:], beta[1:]))

    result = minimize(objective, np.zeros(x.shape[1]), method="BFGS")
    return np.asarray(result.x, dtype=float)


def _classifier_metrics(rows: list[dict[str, Any]], feature_cols: list[str], seed: int) -> list[dict[str, Any]]:
    rng = np.random.default_rng(seed)
    out = []
    combos = sorted({(float(row["h"]), float(row["r"])) for row in rows})
    for h, r in combos:
        subset = [row for row in rows if float(row["h"]) == h and float(row["r"]) == r]
        if len({row["model"] for row in subset}) < 2:
            continue
        x = np.asarray([[float(row[col]) if row[col] not in ("", None) and np.isfinite(float(row[col])) else 0.0 for col in feature_cols] for row in subset], dtype=float)
        y = np.asarray([1 if row["model"] == "hybridization" else 0 for row in subset], dtype=int)
        indices = np.arange(len(y))
        train_idx = []
        test_idx = []
        for label in (0, 1):
            label_idx = indices[y == label]
            rng.shuffle(label_idx)
            split = max(1, int(round(0.7 * len(label_idx))))
            if split >= len(label_idx):
                split = len(label_idx) - 1
            train_idx.extend(label_idx[:split])
            test_idx.extend(label_idx[split:])
        train_idx = np.asarray(train_idx, dtype=int)
        test_idx = np.asarray(test_idx, dtype=int)
        mean = x[train_idx].mean(axis=0)
        sd = x[train_idx].std(axis=0)
        sd[sd == 0.0] = 1.0
        train_x = (x[train_idx] - mean) / sd
        test_x = (x[test_idx] - mean) / sd
        beta = _fit_logistic(train_x, y[train_idx])
        logits = np.column_stack([np.ones(test_x.shape[0]), test_x]) @ beta
        scores = _sigmoid(logits)
        pred = (scores >= 0.5).astype(int)
        accuracy = np.asarray([float(np.mean(pred == y[test_idx]))], dtype=float)
        auc = np.asarray([_auc(y[test_idx], scores)], dtype=float)
        acc_low, acc_high = _ci95(accuracy)
        auc_low, auc_high = _ci95(auc[np.isfinite(auc)])
        out.append({"h": h, "r": r, "accuracy": float(accuracy[0]), "accuracy_ci95_low": acc_low, "accuracy_ci95_high": acc_high, "roc_auc": float(auc[0]), "roc_auc_ci95_low": auc_low, "roc_auc_ci95_high": auc_high, "n_train": int(len(train_idx)), "n_test": int(len(test_idx))})
    return out


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    preferred = ["model", "h", "r", "replicate", "marginal_q1", "marginal_q2", "marginal_q3"]
    keys = sorted({key for row in rows for key in row})
    fieldnames = [key for key in preferred if key in keys] + [key for key in keys if key not in preferred]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _plot_experiment(path: Path, msrc_run: SpatialModelRun, hyb_run: SpatialModelRun, metrics: list[dict[str, Any]]) -> None:
    configure_matplotlib_cache()
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), constrained_layout=True)
    x_m = [float(row["center_bp"]) for row in msrc_run.windows]
    x_h = [float(row["center_bp"]) for row in hyb_run.windows]
    axes[0, 0].plot(x_m, [float(row["q1"]) for row in msrc_run.windows], color=QUARTET_COLORS[0], linewidth=1.2)
    axes[0, 0].plot(x_m, [float(row["q2"]) for row in msrc_run.windows], color=QUARTET_COLORS[1], linewidth=1.2)
    axes[0, 0].plot(x_h, [float(row["q1"]) for row in hyb_run.windows], color=QUARTET_COLORS[0], linestyle="--", linewidth=1.2)
    axes[0, 0].plot(x_h, [float(row["q2"]) for row in hyb_run.windows], color=QUARTET_COLORS[1], linestyle="--", linewidth=1.2)
    axes[0, 0].set_title("A. matched spatial profiles")
    axes[0, 0].set_ylabel("Quartet support")
    axes[0, 1].broken_barh([(0, msrc_run.chromosome_length_bp)], (0.65, 0.18), facecolors="#e6e6e6")
    for interval in msrc_run.feature_intervals:
        axes[0, 1].broken_barh([(float(interval["start_bp"]), float(interval["end_bp"]) - float(interval["start_bp"]))], (0.62, 0.24), facecolors=MSRC_FEATURE)
    for interval in hyb_run.feature_intervals:
        color = HYB_INTRO if interval.get("state") == 1 or str(interval.get("label", "")) == "introgressed" else HYB_MAJOR
        axes[0, 1].broken_barh([(float(interval["start_bp"]), float(interval["end_bp"]) - float(interval["start_bp"]))], (0.22, 0.24), facecolors=color)
    axes[0, 1].set_yticks([0.74, 0.34], ["MSRC", "HYB"])
    axes[0, 1].set_title("B. structural and ancestry tracks")
    by_h: dict[float, list[float]] = {}
    for row in metrics:
        by_h.setdefault(float(row["h"]), []).append(float(row["accuracy"]))
    hs = sorted(by_h)
    axes[1, 0].plot(hs, [float(np.mean(by_h[h])) for h in hs], marker="o", color="0.15")
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_ylim(0, 1.02)
    axes[1, 0].set_title("C. accuracy vs hybridization age")
    axes[1, 0].set_xlabel("h")
    axes[1, 0].set_ylabel("Test accuracy")
    hs = sorted({float(row["h"]) for row in metrics})
    rs = sorted({float(row["r"]) for row in metrics})
    heat = np.full((len(rs), len(hs)), np.nan)
    for row in metrics:
        heat[rs.index(float(row["r"])), hs.index(float(row["h"]))] = float(row["accuracy"])
    image = axes[1, 1].imshow(heat, vmin=0, vmax=1, aspect="auto", origin="lower", cmap="viridis")
    axes[1, 1].set_xticks(range(len(hs)), [f"{h:g}" for h in hs], rotation=45)
    axes[1, 1].set_yticks(range(len(rs)), [f"{r:.1e}" for r in rs])
    axes[1, 1].set_title("D. accuracy over h x r")
    axes[1, 1].set_xlabel("h")
    axes[1, 1].set_ylabel("r")
    fig.colorbar(image, ax=axes[1, 1], label="Accuracy")
    fig.suptitle("Bag-of-genes quartet frequencies are deliberately matched; separation uses spatial organization only")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run a matched MSRC-vs-hybridization spatial distinguishability experiment")
    parser.add_argument("--msrc-dir", required=True)
    parser.add_argument("--matched-hybridization-yaml", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--h-values", default="10,25,50,100,250,500,1000")
    parser.add_argument("--r-multipliers", default="0.25,0.5,1,2")
    parser.add_argument("--baseline-r", type=float, default=1e-8)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--num-loci", type=int, default=5000)
    parser.add_argument("--window-size-loci", type=int, default=50)
    parser.add_argument("--step-loci", type=int, default=10)
    parser.add_argument("--lags-bp", default="10000,100000,1000000")
    parser.add_argument("--seed", type=int, default=1)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        out = Path(args.output_dir)
        if int(args.replicates) < 2:
            raise ValueError("--replicates must be at least 2 for train/test classification")
        match = _load_match(Path(args.matched_hybridization_yaml))
        h_values = _parse_floats(args.h_values)
        r_values = [float(args.baseline_r) * x for x in _parse_floats(args.r_multipliers)]
        lags = _parse_floats(args.lags_bp)
        rng = np.random.default_rng(args.seed)
        msrc_template = load_spatial_model_run(args.msrc_dir, "msrc", window_size_loci=args.window_size_loci, step_loci=args.step_loci)
        focal = int(match["hybridization"]["introgressed_topology"])
        rows: list[dict[str, Any]] = []
        example_hyb: SpatialModelRun | None = None
        for h in h_values:
            for r in r_values:
                for rep in range(int(args.replicates)):
                    msrc_run = _bootstrap_msrc_run(msrc_template, rng, args.window_size_loci, args.step_loci, rep)
                    msrc_features = extract_spatial_features(msrc_run, focal_topology=focal, lags_bp=lags)
                    rows.append({"model": "msrc", "h": h, "r": r, "replicate": rep, **msrc_features})
                    run_dir = out / "hybridization_runs" / f"h_{h:g}_r_{r:.3e}" / f"rep_{rep:04d}"
                    cfg = _hybridization_config(match, h, r, int(rng.integers(1, 2**31 - 1)), run_dir, args.num_loci, args.window_size_loci, args.step_loci, lags)
                    simulate_hybridization(cfg)
                    hyb_run = load_spatial_model_run(run_dir, "hybridization", window_size_loci=args.window_size_loci, step_loci=args.step_loci, marginal_q_source="expected")
                    hyb_run.summary["expected_mean_introgressed_tract_length_bp"] = None if h <= 0 or r <= 0 or float(match["hybridization"]["gamma"]) >= 1.0 else 1.0 / (h * r * max(1e-15, 1.0 - float(match["hybridization"]["gamma"])))
                    if example_hyb is None:
                        example_hyb = hyb_run
                    hyb_features = extract_spatial_features(hyb_run, focal_topology=focal, lags_bp=lags)
                    rows.append({"model": "hybridization", "h": h, "r": r, "replicate": rep, **hyb_features})
        feature_cols = [
            key for key in rows[0]
            if key not in {"model", "h", "r", "replicate", "marginal_q1", "marginal_q2", "marginal_q3"}
            and isinstance(rows[0][key], (int, float))
        ]
        metrics = _classifier_metrics(rows, feature_cols, args.seed)
        _write_csv(out / "spatial_distinguishability.csv", rows)
        _write_csv(out / "spatial_distinguishability_metrics.csv", metrics)
        assert example_hyb is not None
        _plot_experiment(out / "spatial_distinguishability.png", msrc_template, example_hyb, metrics)
        (out / "spatial_distinguishability_summary.json").write_text(json.dumps({
            "msrc_dir": str(msrc_template.run_dir),
            "matched_hybridization_yaml": str(args.matched_hybridization_yaml),
            "target_q": match.get("target_q"),
            "fitted_q": match.get("fitted_q"),
            "matched_l1_error": match.get("l1_error"),
            "feature_columns": feature_cols,
            "num_rows": len(rows),
            "num_metrics": len(metrics),
            "note": "Bag-of-genes quartet frequencies are deliberately matched in expectation; distinguishability is estimated from spatial features only.",
        }, indent=2))
        print(f"Wrote {out / 'spatial_distinguishability.csv'}")
        print(f"Wrote {out / 'spatial_distinguishability_metrics.csv'}")
        print(f"Wrote {out / 'spatial_distinguishability.png'}")
        return 0
    except Exception as exc:
        print(str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
