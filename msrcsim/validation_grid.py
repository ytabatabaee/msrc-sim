from __future__ import annotations

import argparse
import csv
import warnings
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np

from .model_fitting import msc_probabilities
from .robustness import (
    binomial_confidence_interval,
    dominant_quartet_threshold,
    infer_with_strategy,
    interpolate_first_crossing,
    msrc_probabilities_from_beta,
    support_fraction,
)
from .robustness_cli import STRATEGIES, _linkage_diagnostics, _rows_for_fraction


THRESHOLD_REPLICATE_FIELDS = [
    "seed", "replicate_id", "tau", "beta", "rearrangement_fraction",
    "theory_threshold", "support_t1", "support_t2", "support_t3",
    "support_margin_t1_minus_t2", "inferred_topology_index", "recovered_true_t1",
    "num_windows", "num_genealogy_blocks", "kappa", "mode",
]
THRESHOLD_SUMMARY_FIELDS = [
    "tau", "beta", "theory_threshold", "empirical_support_threshold",
    "empirical_recovery50_threshold", "abs_error", "n_replicates",
]
KAPPA_REPLICATE_FIELDS = [
    "seed", "replicate_id", "kappa", "rearrangement_fraction", "region",
    "num_windows", "num_genealogy_blocks", "mean_block_length_bp",
    "median_block_length_bp", "breakpoint_density_per_bp",
    "observed_inside_outside_rate_ratio", "target_kappa",
    "absolute_calibration_error", "relative_calibration_error",
    "fraction_rearrangement_intervals_with_1_block",
    "fraction_rearrangement_intervals_with_2_blocks",
    "fraction_rearrangement_intervals_with_3plus_blocks",
]
KAPPA_SUMMARY_FIELDS = [
    "kappa", "rearrangement_fraction", "n_replicates",
    "inside_breakpoint_density", "outside_breakpoint_density", "observed_ratio",
    "target_kappa", "absolute_calibration_error", "relative_calibration_error",
    "mean_block_length_inside", "median_block_length_inside",
    "mean_block_length_outside", "median_block_length_outside",
    "fraction_rearrangement_intervals_with_1_block",
    "fraction_rearrangement_intervals_with_2_blocks",
    "fraction_rearrangement_intervals_with_3plus_blocks",
]
CORRECTION_SUMMARY_FIELDS = [
    "tau", "beta", "kappa", "rearrangement_fraction", "strategy",
    "n_replicates", "recovery_probability", "ci_low", "ci_high",
    "regime",
]


def _parse_floats(text: str | None, default: Iterable[float]) -> list[float]:
    if text is None or text == "":
        return [float(x) for x in default]
    return [float(x) for x in text.split(",") if x.strip()]


def default_rearrangement_fractions() -> list[float]:
    return [round(i * 0.025, 3) for i in range(0, 29)]


def _seed(base_seed: int, *parts: int) -> int:
    return int(np.random.SeedSequence([int(base_seed), *map(int, parts)]).generate_state(1)[0])


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def _cell_rows(
    *,
    tau: float,
    beta: float,
    fractions: list[float],
    replicates: int,
    base_seed: int,
    cell_index: int,
    chrom: str,
    chrom_length: float,
    windows: int,
    block_windows: int,
    kappa: float,
    mode: str,
) -> list[dict[str, Any]]:
    q_msc = msc_probabilities(0, tau)
    q_msrc = msrc_probabilities_from_beta(beta)
    theory = dominant_quartet_threshold(tau, beta)
    out: list[dict[str, Any]] = []
    for replicate_id in range(replicates):
        for fraction_index, fraction in enumerate(fractions):
            seed = _seed(base_seed, 10, cell_index, replicate_id, fraction_index)
            rows = _rows_for_fraction(
                fraction,
                replicate_id=replicate_id,
                chrom=chrom,
                length=chrom_length,
                windows=windows,
                block_windows=block_windows,
                msrc_block_windows=1,
                kappa=kappa,
                rng=np.random.default_rng(seed),
                q_msc=q_msc,
                q_msrc=q_msrc,
                soft_probability_mode="oracle",
                soft_sensitivity=1.0,
                soft_specificity=1.0,
                soft_noise_sd=0.0,
            )
            result = infer_with_strategy(rows, "all_windows")
            s1 = support_fraction(result, 0)
            s2 = support_fraction(result, 1)
            s3 = support_fraction(result, 2)
            out.append({
                "seed": seed,
                "replicate_id": replicate_id,
                "tau": float(tau),
                "beta": float(beta),
                "rearrangement_fraction": float(fraction),
                "theory_threshold": theory,
                "support_t1": s1,
                "support_t2": s2,
                "support_t3": s3,
                "support_margin_t1_minus_t2": s1 - s2,
                "inferred_topology_index": int(result.inferred_topology_index),
                "recovered_true_t1": result.inferred_topology_index == 0,
                "num_windows": len(rows),
                "num_genealogy_blocks": len({int(row["block_id"]) for row in rows}),
                "kappa": float(kappa),
                "mode": mode,
            })
    return out


def summarize_threshold_grid(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for tau, beta in sorted({(float(row["tau"]), float(row["beta"])) for row in rows}):
        chunk = [row for row in rows if float(row["tau"]) == tau and float(row["beta"]) == beta]
        fractions = sorted({float(row["rearrangement_fraction"]) for row in chunk})
        margins = []
        recoveries = []
        for fraction in fractions:
            sub = [row for row in chunk if float(row["rearrangement_fraction"]) == fraction]
            margins.append(float(np.mean([float(row["support_margin_t1_minus_t2"]) for row in sub])))
            recoveries.append(float(np.mean([bool(row["recovered_true_t1"]) for row in sub])))
        support_threshold = interpolate_first_crossing(fractions, margins, 0.0)
        recovery_threshold = interpolate_first_crossing(fractions, recoveries, 0.5)
        theory = dominant_quartet_threshold(tau, beta)
        if support_threshold is None:
            warnings.warn(f"no support crossing inside sampled grid for tau={tau}, beta={beta}", RuntimeWarning)
        if recovery_threshold is None:
            warnings.warn(f"no recovery=0.5 crossing inside sampled grid for tau={tau}, beta={beta}", RuntimeWarning)
        out.append({
            "tau": tau,
            "beta": beta,
            "theory_threshold": theory,
            "empirical_support_threshold": "" if support_threshold is None else support_threshold,
            "empirical_recovery50_threshold": "" if recovery_threshold is None else recovery_threshold,
            "abs_error": "" if support_threshold is None else abs(support_threshold - theory),
            "n_replicates": len({int(row["replicate_id"]) for row in chunk}),
        })
    return out


def run_threshold_grid(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    taus = _parse_floats(args.tau_values, [0.02, 0.05, 0.1, 0.25, 0.5, 1.0])
    betas = _parse_floats(args.beta_values, [0.1, 0.25, 0.5, 0.75, 1.0])
    fractions = _parse_floats(args.rearrangement_fractions, default_rearrangement_fractions())
    rows: list[dict[str, Any]] = []
    cell_index = 0
    for tau in taus:
        for beta in betas:
            rows.extend(_cell_rows(
                tau=tau,
                beta=beta,
                fractions=fractions,
                replicates=int(args.threshold_replicates),
                base_seed=int(args.seed),
                cell_index=cell_index,
                chrom=args.chrom,
                chrom_length=float(args.chrom_length),
                windows=int(args.windows),
                block_windows=int(args.block_windows),
                kappa=float(args.kappa),
                mode="threshold_grid",
            ))
            cell_index += 1
    return rows, summarize_threshold_grid(rows)


def run_kappa_grid(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    kappas = _parse_floats(args.kappa_values, [1.0, 0.5, 0.25, 0.1, 0.05])
    fractions = _parse_floats(args.kappa_rearrangement_fractions, [0.1, 0.25, 0.5])
    q_msc = msc_probabilities(0, float(args.kappa_tau))
    q_msrc = msrc_probabilities_from_beta(float(args.kappa_beta))
    replicate_rows: list[dict[str, Any]] = []
    cell_index = 0
    for kappa in kappas:
        for fraction in fractions:
            for replicate_id in range(int(args.kappa_replicates)):
                seed = _seed(int(args.seed), 20, cell_index, replicate_id)
                rows = _rows_for_fraction(
                    fraction,
                    replicate_id=replicate_id,
                    chrom=args.chrom,
                    length=float(args.chrom_length),
                    windows=int(args.windows),
                    block_windows=int(args.block_windows),
                    msrc_block_windows=1,
                    kappa=float(kappa),
                    rng=np.random.default_rng(seed),
                    q_msc=q_msc,
                    q_msrc=q_msrc,
                    soft_probability_mode="oracle",
                    soft_sensitivity=1.0,
                    soft_specificity=1.0,
                    soft_noise_sd=0.0,
                )
                diags = _linkage_diagnostics(rows, fraction=fraction, replicate_id=replicate_id)
                ratio = next((float(row["observed_inside_outside_rate_ratio"]) for row in diags if np.isfinite(float(row["observed_inside_outside_rate_ratio"]))), float("nan"))
                abs_error = abs(ratio - kappa) if np.isfinite(ratio) else float("nan")
                rel_error = abs_error / kappa if np.isfinite(abs_error) and kappa > 0 else float("nan")
                for row in diags:
                    item = dict(row)
                    item["seed"] = seed
                    item["kappa"] = float(kappa)
                    item["target_kappa"] = float(kappa)
                    item["absolute_calibration_error"] = abs_error
                    item["relative_calibration_error"] = rel_error
                    replicate_rows.append(item)
            cell_index += 1
    return replicate_rows, summarize_kappa_grid(replicate_rows)


def _finite_mean(values: list[float]) -> float:
    finite = [value for value in values if np.isfinite(value)]
    return float(np.mean(finite)) if finite else float("nan")


def summarize_kappa_grid(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    keys = sorted({(float(row["kappa"]), float(row["rearrangement_fraction"])) for row in rows})
    for kappa, fraction in keys:
        chunk = [row for row in rows if float(row["kappa"]) == kappa and float(row["rearrangement_fraction"]) == fraction]
        inside = [row for row in chunk if row["region"] == "inside"]
        outside = [row for row in chunk if row["region"] == "outside"]
        inside_density = _finite_mean([float(row["breakpoint_density_per_bp"]) for row in inside])
        outside_density = _finite_mean([float(row["breakpoint_density_per_bp"]) for row in outside])
        ratio = inside_density / outside_density if np.isfinite(inside_density) and np.isfinite(outside_density) and outside_density > 0 else float("nan")
        abs_error = abs(ratio - kappa) if np.isfinite(ratio) else float("nan")
        rel_error = abs_error / kappa if np.isfinite(abs_error) and kappa > 0 else float("nan")
        out.append({
            "kappa": kappa,
            "rearrangement_fraction": fraction,
            "n_replicates": len({int(row["replicate_id"]) for row in chunk}),
            "inside_breakpoint_density": inside_density,
            "outside_breakpoint_density": outside_density,
            "observed_ratio": ratio,
            "target_kappa": kappa,
            "absolute_calibration_error": abs_error,
            "relative_calibration_error": rel_error,
            "mean_block_length_inside": _finite_mean([float(row["mean_block_length_bp"]) for row in inside]),
            "median_block_length_inside": _finite_mean([float(row["median_block_length_bp"]) for row in inside]),
            "mean_block_length_outside": _finite_mean([float(row["mean_block_length_bp"]) for row in outside]),
            "median_block_length_outside": _finite_mean([float(row["median_block_length_bp"]) for row in outside]),
            "fraction_rearrangement_intervals_with_1_block": _finite_mean([float(row["fraction_rearrangement_intervals_with_1_block"]) for row in inside]),
            "fraction_rearrangement_intervals_with_2_blocks": _finite_mean([float(row["fraction_rearrangement_intervals_with_2_blocks"]) for row in inside]),
            "fraction_rearrangement_intervals_with_3plus_blocks": _finite_mean([float(row["fraction_rearrangement_intervals_with_3plus_blocks"]) for row in inside]),
        })
    return out


def _representative_cells(args: argparse.Namespace) -> list[tuple[str, float, float, float]]:
    if args.correction_cells:
        cells = []
        for i, text in enumerate(args.correction_cells.split(";")):
            tau, beta, kappa = (float(x) for x in text.split(","))
            cells.append((f"cell_{i}", tau, beta, kappa))
        return cells
    return [
        ("easy_long_tau", 1.0, 0.25, 0.25),
        ("hard_short_tau", 0.05, 0.75, 0.25),
        ("linked_hard", 0.1, 0.5, 0.05),
    ]


def run_correction_grid(args: argparse.Namespace) -> list[dict[str, Any]]:
    fractions = _parse_floats(args.correction_rearrangement_fractions, [0.0, 0.1, 0.25, 0.4, 0.55, 0.7])
    out: list[dict[str, Any]] = []
    for cell_index, (regime, tau, beta, kappa) in enumerate(_representative_cells(args)):
        q_msc = msc_probabilities(0, tau)
        q_msrc = msrc_probabilities_from_beta(beta)
        for fraction_index, fraction in enumerate(fractions):
            successes = {strategy: 0 for strategy in STRATEGIES}
            for replicate_id in range(int(args.correction_replicates)):
                seed = _seed(int(args.seed), 30, cell_index, fraction_index, replicate_id)
                rows = _rows_for_fraction(
                    fraction,
                    replicate_id=replicate_id,
                    chrom=args.chrom,
                    length=float(args.chrom_length),
                    windows=int(args.windows),
                    block_windows=int(args.block_windows),
                    msrc_block_windows=1,
                    kappa=float(kappa),
                    rng=np.random.default_rng(seed),
                    q_msc=q_msc,
                    q_msrc=q_msrc,
                    soft_probability_mode=args.soft_probability_mode,
                    soft_sensitivity=float(args.soft_sensitivity),
                    soft_specificity=float(args.soft_specificity),
                    soft_noise_sd=float(args.soft_noise_sd),
                )
                for strategy in STRATEGIES:
                    successes[strategy] += infer_with_strategy(rows, strategy).inferred_topology_index == 0
            for strategy in STRATEGIES:
                low, high = binomial_confidence_interval(successes[strategy], int(args.correction_replicates))
                out.append({
                    "tau": tau,
                    "beta": beta,
                    "kappa": kappa,
                    "rearrangement_fraction": fraction,
                    "strategy": strategy,
                    "n_replicates": int(args.correction_replicates),
                    "recovery_probability": successes[strategy] / int(args.correction_replicates),
                    "ci_low": low,
                    "ci_high": high,
                    "regime": regime,
                })
    return out


def _plot_thresholds(summary: list[dict[str, Any]], path: Path) -> None:
    xs = [float(row["theory_threshold"]) for row in summary if row["empirical_support_threshold"] != ""]
    ys = [float(row["empirical_support_threshold"]) for row in summary if row["empirical_support_threshold"] != ""]
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    ax.scatter(xs, ys, s=24, color="#2c7fb8")
    lim = [0.0, min(1.0, max(xs + ys + [0.7]) + 0.05)]
    ax.plot(lim, lim, color="black", linewidth=1.0)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Analytic threshold")
    ax.set_ylabel("Empirical support-crossing threshold")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_threshold_heatmap(summary: list[dict[str, Any]], path: Path) -> None:
    taus = sorted({float(row["tau"]) for row in summary})
    betas = sorted({float(row["beta"]) for row in summary})
    data = np.full((len(taus), len(betas)), np.nan)
    for row in summary:
        if row["empirical_support_threshold"] != "":
            i = taus.index(float(row["tau"]))
            j = betas.index(float(row["beta"]))
            data[i, j] = float(row["empirical_support_threshold"]) - float(row["theory_threshold"])
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    im = ax.imshow(data, origin="lower", aspect="auto", cmap="coolwarm")
    ax.set_xticks(range(len(betas)), [str(x) for x in betas])
    ax.set_yticks(range(len(taus)), [str(x) for x in taus])
    ax.set_xlabel("beta")
    ax.set_ylabel("tau")
    fig.colorbar(im, ax=ax, label="Empirical - theoretical threshold")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_recovery_phase(rows: list[dict[str, Any]], path: Path) -> None:
    betas = sorted({float(row["beta"]) for row in rows})
    taus = sorted({float(row["tau"]) for row in rows})
    fractions = sorted({float(row["rearrangement_fraction"]) for row in rows})
    fig, axes = plt.subplots(1, len(betas), figsize=(3.0 * len(betas), 4.6), squeeze=False)
    for ax, beta in zip(axes[0], betas):
        data = np.full((len(taus), len(fractions)), np.nan)
        for i, tau in enumerate(taus):
            for j, fraction in enumerate(fractions):
                chunk = [row for row in rows if float(row["beta"]) == beta and float(row["tau"]) == tau and float(row["rearrangement_fraction"]) == fraction]
                if chunk:
                    data[i, j] = np.mean([bool(row["recovered_true_t1"]) for row in chunk])
        im = ax.imshow(data, origin="lower", aspect="auto", vmin=0.0, vmax=1.0, cmap="viridis")
        ax.set_title(f"beta={beta:g}")
        ax.set_xticks(range(len(fractions))[::max(1, len(fractions) // 5)], [f"{fractions[k]:.2f}" for k in range(len(fractions))[::max(1, len(fractions) // 5)]], rotation=45)
        ax.set_yticks(range(len(taus)), [str(x) for x in taus])
        ax.set_xlabel("f_R")
        ax.set_ylabel("tau")
    fig.colorbar(im, ax=axes.ravel().tolist(), label="Naive P(T1)")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _plot_kappa(summary: list[dict[str, Any]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    fractions = sorted({float(row["rearrangement_fraction"]) for row in summary})
    for fraction in fractions:
        chunk = [row for row in summary if float(row["rearrangement_fraction"]) == fraction]
        xs = [float(row["target_kappa"]) for row in chunk]
        ys = [float(row["observed_ratio"]) for row in chunk]
        ax.plot(xs, ys, marker="o", linewidth=1.4, label=f"f_R={fraction:g}")
    ax.plot([0.0, 1.0], [0.0, 1.0], color="black", linewidth=1.0)
    ax.set_xlabel("Target kappa")
    ax.set_ylabel("Observed inside/outside breakpoint-rate ratio")
    ax.set_xlim(0.0, 1.02)
    ax.set_ylim(0.0, 1.02)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_corrections(summary: list[dict[str, Any]], path: Path) -> None:
    regimes = sorted({row["regime"] for row in summary})
    colors = {
        "all_windows": "#d95f02",
        "oracle_filter": "#1b9e77",
        "genealogy_block_collapse": "#7570b3",
        "rearrangement_interval_collapse": "#e7298a",
        "soft_weight": "#66a61e",
    }
    fig, axes = plt.subplots(1, len(regimes), figsize=(4.8 * len(regimes), 4.4), squeeze=False)
    for ax, regime in zip(axes[0], regimes):
        for strategy in STRATEGIES:
            chunk = [row for row in summary if row["regime"] == regime and row["strategy"] == strategy]
            xs = [float(row["rearrangement_fraction"]) for row in chunk]
            ys = np.asarray([float(row["recovery_probability"]) for row in chunk])
            low = np.asarray([float(row["ci_low"]) for row in chunk])
            high = np.asarray([float(row["ci_high"]) for row in chunk])
            ax.plot(xs, ys, marker="o", markersize=2.5, linewidth=1.4, color=colors[strategy], label=strategy)
            ax.fill_between(xs, low, high, color=colors[strategy], alpha=0.12, linewidth=0)
        ax.set_title(regime)
        ax.set_xlabel("f_R")
        ax.set_ylabel("P(T1)")
        ax.set_ylim(0.0, 1.0)
    axes[0][-1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def run_validation_grid(args: argparse.Namespace) -> Path:
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    threshold_rows, threshold_summary = run_threshold_grid(args)
    kappa_rows, kappa_summary = run_kappa_grid(args)
    correction_summary = run_correction_grid(args)
    _write_csv(out / "threshold_grid_replicates.csv", THRESHOLD_REPLICATE_FIELDS, threshold_rows)
    _write_csv(out / "threshold_grid_summary.csv", THRESHOLD_SUMMARY_FIELDS, threshold_summary)
    _write_csv(out / "kappa_calibration_replicates.csv", KAPPA_REPLICATE_FIELDS, kappa_rows)
    _write_csv(out / "kappa_calibration_summary.csv", KAPPA_SUMMARY_FIELDS, kappa_summary)
    _write_csv(out / "correction_grid_summary.csv", CORRECTION_SUMMARY_FIELDS, correction_summary)
    _plot_thresholds(threshold_summary, out / "theory_vs_empirical_threshold.pdf")
    _plot_threshold_heatmap(threshold_summary, out / "threshold_error_heatmap.pdf")
    _plot_recovery_phase(threshold_rows, out / "recovery_phase_diagram.pdf")
    _plot_kappa(kappa_summary, out / "kappa_calibration.pdf")
    _plot_corrections(correction_summary, out / "correction_methods_grid.pdf")
    return out
