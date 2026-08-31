from __future__ import annotations

from pathlib import Path
from typing import Any
import argparse
import csv
import json
import tempfile

import numpy as np

from .hybridization import simulate_hybridization
from .hybridization import topology_index
from .spatial_compare import SpatialModelRun
from .spatial_compare_io import load_spatial_model_run
from .spatial_compare_plot import configure_matplotlib_cache, make_spatial_compare_figure
from .spatial_identifiability import (
    BAG_FEATURE_COLUMNS,
    eta_to_hr,
    exact_match_for_gamma,
    gamma_grid,
    identifiability_features,
    rearranged_interval_length,
    resample_msrc_chromosome,
    spatial_classifier_columns,
    stratified_cv_metrics,
)


def _parse_floats(text: str) -> list[float]:
    return [float(piece.strip()) for piece in text.split(",") if piece.strip()]


def _write_locus_positions(path: Path, msrc_run: SpatialModelRun) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["locus_id", "position_bp"])
        writer.writeheader()
        for row in msrc_run.loci:
            writer.writerow({"locus_id": row["locus_id"], "position_bp": row["position_bp"]})


def _hybridization_config(
    *,
    run_dir: Path,
    positions_file: Path,
    chromosome_length: float,
    num_loci: int,
    match: Any,
    h: float,
    r: float,
    window_size_loci: int,
    step_loci: int,
    lags_bp: list[float],
    seed: int,
) -> dict[str, Any]:
    return {
        "mode": "pulse_hybridization",
        "seed": int(seed),
        "chromosome": {
            "length_bp": float(chromosome_length),
            "num_loci": int(num_loci),
            "locus_positions": {"mode": "file", "path": str(positions_file), "column": "position_bp"},
        },
        "hybridization": {
            "gamma": float(match.gamma),
            "generations_since_pulse": float(h),
            "major": {"topology": int(match.major_topology), "internal_branch_length": float(match.t_major)},
            "introgressed": {"topology": int(match.introgressed_topology), "internal_branch_length": float(match.t_introgressed)},
            "donor": "donor",
            "recipient": "recipient",
        },
        "recombination": {"rate_per_bp_per_generation": float(r)},
        "windows": {"loci_per_window": int(window_size_loci), "step_loci": int(step_loci)},
        "spatial_statistics": {"lags_bp": [float(x) for x in lags_bp]},
        "output": {
            "directory": str(run_dir),
            "record_tracts": True,
            "record_loci": True,
            "record_windows": True,
            "make_plots": False,
        },
    }


def _write_csv(path: Path, rows: list[dict[str, Any]], preferred: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({key for row in rows for key in row})
    fieldnames = [key for key in preferred if key in keys] + [key for key in keys if key not in preferred]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _plot_heatmap(path: Path, summary_rows: list[dict[str, Any]]) -> None:
    configure_matplotlib_cache()
    import matplotlib.pyplot as plt

    gammas = sorted({float(row["gamma"]) for row in summary_rows})
    etas = sorted({float(row["eta"]) for row in summary_rows})
    heat = np.full((len(gammas), len(etas)), np.nan)
    for row in summary_rows:
        heat[gammas.index(float(row["gamma"])), etas.index(float(row["eta"]))] = float(row["spatial_roc_auc"])
    fig, ax = plt.subplots(figsize=(8.2, 5.4), constrained_layout=True)
    image = ax.imshow(heat, origin="lower", aspect="auto", vmin=0.0, vmax=1.0, cmap="viridis")
    ax.set_xticks(range(len(etas)), [f"{eta:g}" for eta in etas], rotation=45)
    ax.set_yticks(range(len(gammas)), [f"{gamma:.3f}" for gamma in gammas])
    ax.set_xlabel("eta = E[introgressed tract length] / rearranged interval length")
    ax.set_ylabel("gamma")
    ax.set_title("Spatial ROC-AUC for exactly matched expected quartet vectors")
    ax.axhline(-10, color="white", label="chance = 0.5")
    fig.colorbar(image, ax=ax, label="Spatial ROC-AUC")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _plot_auc_vs_eta(path: Path, summary_rows: list[dict[str, Any]]) -> None:
    configure_matplotlib_cache()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.4, 5.0), constrained_layout=True)
    gammas = sorted({float(row["gamma"]) for row in summary_rows})
    for gamma in gammas:
        rows = sorted([row for row in summary_rows if float(row["gamma"]) == gamma], key=lambda row: float(row["eta"]))
        ax.plot([float(row["eta"]) for row in rows], [float(row["spatial_roc_auc"]) for row in rows], marker="o", linewidth=1.2, label=f"gamma={gamma:.3f}")
    ax.axhline(0.5, color="0.4", linestyle="--", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("eta")
    ax.set_ylabel("Spatial ROC-AUC")
    ax.set_title("Spatial distinguishability as tract scale changes")
    ax.legend(frameon=False, fontsize=7, ncols=2)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _plot_bag_vs_spatial(path: Path, summary_rows: list[dict[str, Any]]) -> None:
    configure_matplotlib_cache()
    import matplotlib.pyplot as plt

    etas = sorted({float(row["eta"]) for row in summary_rows})
    bag = []
    spatial = []
    for eta in etas:
        rows = [row for row in summary_rows if float(row["eta"]) == eta]
        bag.append(float(np.mean([float(row["bag_roc_auc"]) for row in rows])))
        spatial.append(float(np.mean([float(row["spatial_roc_auc"]) for row in rows])))
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    ax.plot(etas, bag, marker="o", color="0.45", label="bag-of-genes q")
    ax.plot(etas, spatial, marker="o", color="#3b6fb6", label="spatial features")
    ax.axhline(0.5, color="0.4", linestyle="--", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("eta")
    ax.set_ylabel("ROC-AUC")
    ax.set_title("Matched bag-of-genes baseline versus spatial classifier")
    ax.legend(frameon=False)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _plot_examples(path: Path, examples: dict[str, tuple[SpatialModelRun, SpatialModelRun]]) -> None:
    configure_matplotlib_cache()
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt

    with PdfPages(path) as pdf:
        for label in ("eta << 1", "eta ~ 1", "eta >> 1"):
            if label not in examples:
                continue
            msrc_run, hyb_run = examples[label]
            fig = make_spatial_compare_figure(msrc_run, hyb_run, title=f"Example profiles: {label}; expected bag-of-genes q is matched")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Quantify spatial identifiability for exactly matched MSRC and pulse-hybridization quartet vectors")
    parser.add_argument("--msrc-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--major-topology", default="12|34")
    parser.add_argument("--introgressed-topology", default="13|24")
    parser.add_argument("--num-gamma", type=int, default=10)
    parser.add_argument("--eta-values", default="0.02,0.05,0.1,0.25,0.5,1,2,5,10")
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--min-replicates-for-analysis", type=int, default=20)
    parser.add_argument("--allow-small-sample", action="store_true", help="Allow debug/smoke runs below the analysis replicate minimum")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--window-size-loci", type=int, default=50)
    parser.add_argument("--step-loci", type=int, default=10)
    parser.add_argument("--lags-bp", default="10000,100000,1000000")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--h", type=float, default=1.0, help="Hybridization age used with r=hr/h; only the product h*r controls tract lengths")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.replicates < 2:
            raise ValueError("--replicates must be at least 2")
        if args.replicates < args.min_replicates_for_analysis and not args.allow_small_sample:
            raise ValueError(
                f"--replicates must be at least {args.min_replicates_for_analysis} for interpretable analysis; "
                "use --allow-small-sample only for debug/smoke tests"
            )
        if args.h <= 0.0:
            raise ValueError("--h must be positive")
        output = Path(args.output_dir)
        output.mkdir(parents=True, exist_ok=True)
        msrc_template = load_spatial_model_run(args.msrc_dir, "msrc", window_size_loci=args.window_size_loci, step_loci=args.step_loci)
        if not msrc_template.loci:
            raise ValueError("MSRC locus-level output is required for spatial-identifiability replicates")
        major = topology_index(args.major_topology)
        intro = topology_index(args.introgressed_topology)
        target_q = msrc_template.marginal_q
        gammas = gamma_grid(target_q, args.num_gamma, major, intro)
        etas = _parse_floats(args.eta_values)
        lags = _parse_floats(args.lags_bp)
        l_r = rearranged_interval_length(msrc_template)
        positions_file = output / "locus_positions.csv"
        _write_locus_positions(positions_file, msrc_template)
        rng = np.random.default_rng(args.seed)
        rows: list[dict[str, Any]] = []
        example_candidates: dict[float, tuple[SpatialModelRun, SpatialModelRun]] = {}
        with tempfile.TemporaryDirectory(prefix="msrcsim-identifiability-") as tmp:
            tmp_root = Path(tmp)
            for gamma in gammas:
                match = exact_match_for_gamma(target_q, gamma, major, intro)
                for eta in etas:
                    hr = eta_to_hr(eta, l_r, gamma)
                    r = hr / float(args.h)
                    for replicate in range(args.replicates):
                        msrc_run = resample_msrc_chromosome(msrc_template, rng, args.window_size_loci, args.step_loci)
                        msrc_features = identifiability_features(msrc_run, msrc_template, focal_topology=intro, lags_bp=lags)
                        msrc_features.update({"marginal_q1": float(target_q[0]), "marginal_q2": float(target_q[1]), "marginal_q3": float(target_q[2])})
                        common = {
                            "gamma": float(gamma),
                            "eta": float(eta),
                            "hr": float(hr),
                            "h": float(args.h),
                            "r": float(r),
                            "replicate": int(replicate),
                            "target_q1": float(target_q[0]),
                            "target_q2": float(target_q[1]),
                            "target_q3": float(target_q[2]),
                            "t_major": float(match.t_major),
                            "t_introgressed": float(match.t_introgressed),
                            "expected_introgressed_tract_length_bp": float(eta * l_r),
                        }
                        rows.append({"model": "msrc", **common, **msrc_features})
                        run_dir = tmp_root / f"g_{gamma:.6f}_e_{eta:.6f}_r_{replicate}"
                        cfg = _hybridization_config(
                            run_dir=run_dir,
                            positions_file=positions_file,
                            chromosome_length=msrc_template.chromosome_length_bp,
                            num_loci=len(msrc_template.loci),
                            match=match,
                            h=float(args.h),
                            r=r,
                            window_size_loci=args.window_size_loci,
                            step_loci=args.step_loci,
                            lags_bp=lags,
                            seed=int(rng.integers(1, 2**31 - 1)),
                        )
                        simulate_hybridization(cfg)
                        hyb_run = load_spatial_model_run(run_dir, "hybridization", window_size_loci=args.window_size_loci, step_loci=args.step_loci, marginal_q_source="expected")
                        hyb_run.summary["expected_mean_introgressed_tract_length_bp"] = eta * l_r
                        hyb_features = identifiability_features(hyb_run, msrc_template, focal_topology=intro, lags_bp=lags)
                        hyb_features.update({"marginal_q1": float(target_q[0]), "marginal_q2": float(target_q[1]), "marginal_q3": float(target_q[2])})
                        rows.append({"model": "hybridization", **common, **hyb_features})
                        if replicate == 0:
                            example_candidates.setdefault(float(eta), (msrc_run, hyb_run))
        spatial_cols = spatial_classifier_columns(rows)
        summary_rows = []
        for gamma in gammas:
            for eta in etas:
                subset = [row for row in rows if float(row["gamma"]) == float(gamma) and float(row["eta"]) == float(eta)]
                bag = stratified_cv_metrics(subset, BAG_FEATURE_COLUMNS, folds=args.folds, seed=args.seed)
                spatial = stratified_cv_metrics(subset, spatial_cols, folds=args.folds, seed=args.seed)
                summary_rows.append({
                    "gamma": float(gamma),
                    "eta": float(eta),
                    "hr": eta_to_hr(eta, l_r, gamma),
                    "h": float(args.h),
                    "r": eta_to_hr(eta, l_r, gamma) / float(args.h),
                    "replicates_per_model": int(args.replicates),
                    "small_sample_warning": (
                        f"Only {args.replicates} replicates per model; AUC/accuracy are not scientifically interpretable."
                        if args.replicates < args.min_replicates_for_analysis else ""
                    ),
                    "bag_accuracy": bag["accuracy"],
                    "bag_accuracy_ci95_low": bag["accuracy_ci95_low"],
                    "bag_accuracy_ci95_high": bag["accuracy_ci95_high"],
                    "bag_roc_auc": bag["roc_auc"],
                    "bag_roc_auc_ci95_low": bag["roc_auc_ci95_low"],
                    "bag_roc_auc_ci95_high": bag["roc_auc_ci95_high"],
                    "spatial_accuracy": spatial["accuracy"],
                    "spatial_accuracy_ci95_low": spatial["accuracy_ci95_low"],
                    "spatial_accuracy_ci95_high": spatial["accuracy_ci95_high"],
                    "spatial_roc_auc": spatial["roc_auc"],
                    "spatial_roc_auc_ci95_low": spatial["roc_auc_ci95_low"],
                    "spatial_roc_auc_ci95_high": spatial["roc_auc_ci95_high"],
                })
        _write_csv(
            output / "spatial_identifiability_replicates.csv",
            rows,
            ["model", "gamma", "eta", "hr", "h", "r", "replicate", "marginal_q1", "marginal_q2", "marginal_q3"],
        )
        _write_csv(output / "spatial_identifiability_summary.csv", summary_rows, ["gamma", "eta", "hr", "h", "r"])
        _plot_heatmap(output / "identifiability_heatmap.pdf", summary_rows)
        _plot_auc_vs_eta(output / "auc_vs_eta.pdf", summary_rows)
        _plot_bag_vs_spatial(output / "bag_vs_spatial_auc.pdf", summary_rows)
        selected_examples: dict[str, tuple[SpatialModelRun, SpatialModelRun]] = {}
        eta_sorted = sorted(example_candidates)
        if eta_sorted:
            selected_examples["eta << 1"] = example_candidates[min(eta_sorted)]
            selected_examples["eta ~ 1"] = example_candidates[min(eta_sorted, key=lambda x: abs(np.log(x)))]
            selected_examples["eta >> 1"] = example_candidates[max(eta_sorted)]
            _plot_examples(output / "example_profiles.pdf", selected_examples)
        metadata = {
            "msrc_dir": str(msrc_template.run_dir),
            "target_q": [float(x) for x in target_q],
            "rearranged_interval_length": l_r,
            "major_topology": major,
            "introgressed_topology": intro,
            "gamma_values": gammas,
            "eta_values": etas,
            "replicates_per_model": args.replicates,
            "min_replicates_for_analysis": args.min_replicates_for_analysis,
            "small_sample_warning": (
                f"Only {args.replicates} replicates per model per cell; outputs are for debugging, not inference."
                if args.replicates < args.min_replicates_for_analysis else None
            ),
            "bag_feature_columns": BAG_FEATURE_COLUMNS,
            "spatial_feature_columns": spatial_cols,
            "excluded_from_spatial_classifier": ["model", "gamma", "eta", "hr", "h", "r", "replicate", "t_major", "t_introgressed", "expected_introgressed_tract_length_bp", "marginal_q*", "observed_marginal_q*", "num_ancestry_tracts", "realized_introgressed_fraction", "mean_introgressed_tract_length"],
            "msrc_replicate_mode": "conditional_resampling_from_fixed_run_inside_outside_topology_distributions",
            "interpretation": "Expected bag-of-genes quartet frequencies are deliberately matched. Spatial distinguishability is parameter-regime dependent and comes only from spatial organization features.",
        }
        (output / "spatial_identifiability_metadata.json").write_text(json.dumps(metadata, indent=2))
        print(f"Wrote {output / 'spatial_identifiability_replicates.csv'}")
        print(f"Wrote {output / 'spatial_identifiability_summary.csv'}")
        print(f"Wrote {output / 'identifiability_heatmap.pdf'}")
        print(f"Wrote {output / 'auc_vs_eta.pdf'}")
        print(f"Wrote {output / 'bag_vs_spatial_auc.pdf'}")
        print(f"Wrote {output / 'example_profiles.pdf'}")
        return 0
    except Exception as exc:
        print(str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
