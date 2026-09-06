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
THRESHOLD_REGRESSION_FIELDS = [
    "tau", "beta", "theory_threshold",
    "regression_threshold", "regression_ci_low", "regression_ci_high",
    "interpolation_threshold",
    "intercept", "slope", "r2", "abs_error", "n_replicates",
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
CONSISTENCY_SUMMARY_FIELDS = [
    "tau", "beta", "theory_threshold", "regime", "epsilon", "n_blocks",
    "strategy", "n_replicates", "recovery_probability", "ci_low", "ci_high",
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


def _is_inside(value: float | None, xs: Iterable[float]) -> bool:
    if value is None or not np.isfinite(value):
        return False
    vals = [float(x) for x in xs]
    return min(vals) <= float(value) <= max(vals)


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


def fit_margin_regression(xs: Iterable[float], ys: Iterable[float]) -> dict[str, float | None]:
    x = np.asarray([float(value) for value in xs], dtype=float)
    y = np.asarray([float(value) for value in ys], dtype=float)
    if x.size != y.size or x.size < 2:
        raise ValueError("regression threshold requires at least two paired points")
    slope, intercept = np.polyfit(x, y, 1)
    fitted = intercept + slope * x
    ss_res = float(np.sum((y - fitted) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    r2 = 1.0 if ss_tot == 0.0 and ss_res == 0.0 else 0.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    threshold = None if slope == 0.0 else float(-intercept / slope)
    return {
        "intercept": float(intercept),
        "slope": float(slope),
        "r2": float(r2),
        "threshold": threshold,
    }


def _bootstrap_regression_threshold(
    chunk: list[dict[str, Any]],
    fractions: list[float],
    *,
    seed: int,
    n_bootstrap: int,
) -> tuple[float | None, float | None]:
    if n_bootstrap <= 0:
        return None, None
    rng = np.random.default_rng(seed)
    by_fraction = {
        fraction: [float(row["support_margin_t1_minus_t2"]) for row in chunk if float(row["rearrangement_fraction"]) == fraction]
        for fraction in fractions
    }
    draws = []
    for _ in range(n_bootstrap):
        means = []
        for fraction in fractions:
            values = np.asarray(by_fraction[fraction], dtype=float)
            if values.size == 0:
                means.append(float("nan"))
            else:
                means.append(float(np.mean(rng.choice(values, size=values.size, replace=True))))
        if any(not np.isfinite(value) for value in means):
            continue
        fit = fit_margin_regression(fractions, means)
        threshold = fit["threshold"]
        if _is_inside(threshold, fractions):
            draws.append(float(threshold))
    if not draws:
        return None, None
    low, high = np.percentile(np.asarray(draws, dtype=float), [2.5, 97.5])
    return float(low), float(high)


def summarize_threshold_regression(
    rows: list[dict[str, Any]],
    *,
    bootstrap_replicates: int = 500,
    seed: int = 1,
) -> list[dict[str, Any]]:
    out = []
    for cell_index, (tau, beta) in enumerate(sorted({(float(row["tau"]), float(row["beta"])) for row in rows})):
        chunk = [row for row in rows if float(row["tau"]) == tau and float(row["beta"]) == beta]
        fractions = sorted({float(row["rearrangement_fraction"]) for row in chunk})
        margins = []
        for fraction in fractions:
            sub = [row for row in chunk if float(row["rearrangement_fraction"]) == fraction]
            margins.append(float(np.mean([float(row["support_margin_t1_minus_t2"]) for row in sub])))
        interpolation_threshold = interpolate_first_crossing(fractions, margins, 0.0)
        fit = fit_margin_regression(fractions, margins)
        raw_threshold = fit["threshold"]
        regression_threshold = float(raw_threshold) if _is_inside(raw_threshold, fractions) else None
        ci_low, ci_high = _bootstrap_regression_threshold(
            chunk,
            fractions,
            seed=_seed(seed, 40, cell_index),
            n_bootstrap=bootstrap_replicates,
        )
        theory = dominant_quartet_threshold(tau, beta)
        if regression_threshold is None:
            warnings.warn(f"regression threshold outside sampled grid for tau={tau}, beta={beta}", RuntimeWarning)
        out.append({
            "tau": tau,
            "beta": beta,
            "theory_threshold": theory,
            "regression_threshold": "" if regression_threshold is None else regression_threshold,
            "regression_ci_low": "" if ci_low is None else max(0.0, min(1.0, ci_low)),
            "regression_ci_high": "" if ci_high is None else max(0.0, min(1.0, ci_high)),
            "interpolation_threshold": "" if interpolation_threshold is None else interpolation_threshold,
            "intercept": float(fit["intercept"]),
            "slope": float(fit["slope"]),
            "r2": float(fit["r2"]),
            "abs_error": "" if regression_threshold is None else abs(regression_threshold - theory),
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


def _independent_block_rows(
    *,
    tau: float,
    beta: float,
    epsilon: float,
    n_blocks: int,
    rng: np.random.Generator,
    soft_probability_mode: str,
    soft_sensitivity: float,
    soft_specificity: float,
    soft_noise_sd: float,
) -> list[dict[str, Any]]:
    q_msc = msc_probabilities(0, tau)
    q_msrc = msrc_probabilities_from_beta(beta)
    n_rearranged = int(round(float(epsilon) * int(n_blocks)))
    rows = []
    for block_id in range(int(n_blocks)):
        is_rearranged = block_id < n_rearranged
        topology = int(rng.choice(3, p=q_msrc if is_rearranged else q_msc))
        if soft_probability_mode == "oracle":
            p_msrc = 1.0 if is_rearranged else 0.0
        elif soft_probability_mode == "noisy":
            mean = soft_sensitivity if is_rearranged else 1.0 - soft_specificity
            p_msrc = float(np.clip(mean + rng.normal(0.0, soft_noise_sd), 0.0, 1.0))
        else:
            raise ValueError("soft probability mode must be oracle or noisy")
        rows.append({
            "block_id": block_id,
            "topology_index": topology,
            "topology": str(topology),
            "is_rearranged": is_rearranged,
            "rearrangement_id": "independent_msrc_interval" if is_rearranged else "",
            "msrc_probability": p_msrc,
            "weight": 1.0 - p_msrc,
        })
    rng.shuffle(rows)
    return rows


def run_consistency_benchmark(args: argparse.Namespace) -> list[dict[str, Any]]:
    n_blocks_values = [int(x) for x in _parse_floats(args.consistency_n_blocks, [10, 25, 50, 100, 250, 500, 1000])]
    tau = float(args.consistency_tau)
    beta = float(args.consistency_beta)
    threshold = dominant_quartet_threshold(tau, beta)
    regimes = [
        ("pure_msc", 0.0),
        ("below_threshold", max(0.0, min(1.0, 0.75 * threshold))),
        ("above_threshold", max(0.0, min(1.0, 1.25 * threshold))),
    ]
    out = []
    for regime_index, (regime, epsilon) in enumerate(regimes):
        for n_index, n_blocks in enumerate(n_blocks_values):
            successes = {strategy: 0 for strategy in STRATEGIES}
            for replicate_id in range(int(args.consistency_replicates)):
                seed = _seed(int(args.seed), 50, regime_index, n_index, replicate_id)
                rows = _independent_block_rows(
                    tau=tau,
                    beta=beta,
                    epsilon=epsilon,
                    n_blocks=n_blocks,
                    rng=np.random.default_rng(seed),
                    soft_probability_mode=args.consistency_soft_probability_mode,
                    soft_sensitivity=float(args.consistency_soft_sensitivity),
                    soft_specificity=float(args.consistency_soft_specificity),
                    soft_noise_sd=float(args.consistency_soft_noise_sd),
                )
                for strategy in STRATEGIES:
                    successes[strategy] += infer_with_strategy(rows, strategy).inferred_topology_index == 0
            for strategy in STRATEGIES:
                low, high = binomial_confidence_interval(successes[strategy], int(args.consistency_replicates))
                out.append({
                    "tau": tau,
                    "beta": beta,
                    "theory_threshold": threshold,
                    "regime": regime,
                    "epsilon": epsilon,
                    "n_blocks": n_blocks,
                    "strategy": strategy,
                    "n_replicates": int(args.consistency_replicates),
                    "recovery_probability": successes[strategy] / int(args.consistency_replicates),
                    "ci_low": low,
                    "ci_high": high,
                })
    return out


def _plot_thresholds(summary: list[dict[str, Any]], path: Path) -> None:
    threshold_field = "regression_threshold" if summary and "regression_threshold" in summary[0] else "empirical_support_threshold"
    xs = [float(row["theory_threshold"]) for row in summary if row[threshold_field] != ""]
    ys = [float(row[threshold_field]) for row in summary if row[threshold_field] != ""]
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    if threshold_field == "regression_threshold":
        yerr_low = []
        yerr_high = []
        for row in summary:
            if row[threshold_field] == "":
                continue
            y = float(row[threshold_field])
            low = float(row["regression_ci_low"]) if row["regression_ci_low"] != "" else y
            high = float(row["regression_ci_high"]) if row["regression_ci_high"] != "" else y
            yerr_low.append(max(0.0, y - low))
            yerr_high.append(max(0.0, high - y))
        ax.errorbar(xs, ys, yerr=[yerr_low, yerr_high], fmt="o", markersize=4, color="#2c7fb8", ecolor="#7bccc4", elinewidth=1.0, capsize=2)
    else:
        ax.scatter(xs, ys, s=24, color="#2c7fb8")
    lim = [0.0, min(1.0, max(xs + ys + [0.7]) + 0.05)]
    ax.plot(lim, lim, color="black", linewidth=1.0)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Analytic threshold")
    ax.set_ylabel("Empirical regression threshold" if threshold_field == "regression_threshold" else "Empirical support-crossing threshold")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_threshold_heatmap(summary: list[dict[str, Any]], path: Path) -> None:
    threshold_field = "regression_threshold" if summary and "regression_threshold" in summary[0] else "empirical_support_threshold"
    taus = sorted({float(row["tau"]) for row in summary})
    betas = sorted({float(row["beta"]) for row in summary})
    data = np.full((len(taus), len(betas)), np.nan)
    for row in summary:
        if row[threshold_field] != "":
            i = taus.index(float(row["tau"]))
            j = betas.index(float(row["beta"]))
            data[i, j] = float(row[threshold_field]) - float(row["theory_threshold"])
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


def _plot_consistency(summary: list[dict[str, Any]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    colors = {"pure_msc": "#1b9e77", "below_threshold": "#7570b3", "above_threshold": "#d95f02"}
    for regime in ("pure_msc", "below_threshold", "above_threshold"):
        chunk = [row for row in summary if row["strategy"] == "all_windows" and row["regime"] == regime]
        xs = [int(row["n_blocks"]) for row in chunk]
        ys = np.asarray([float(row["recovery_probability"]) for row in chunk])
        low = np.asarray([float(row["ci_low"]) for row in chunk])
        high = np.asarray([float(row["ci_high"]) for row in chunk])
        ax.plot(xs, ys, marker="o", linewidth=1.5, color=colors[regime], label=regime)
        ax.fill_between(xs, low, high, color=colors[regime], alpha=0.12, linewidth=0)
    ax.set_xscale("log")
    ax.set_xlabel("n_blocks")
    ax.set_ylabel("P(inferred quartet = T1)")
    ax.set_ylim(0.0, 1.0)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_consistency_corrections(summary: list[dict[str, Any]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    colors = {
        "all_windows": "#d95f02",
        "oracle_filter": "#1b9e77",
        "genealogy_block_collapse": "#7570b3",
        "rearrangement_interval_collapse": "#e7298a",
        "soft_weight": "#66a61e",
    }
    for strategy in STRATEGIES:
        chunk = [row for row in summary if row["regime"] == "above_threshold" and row["strategy"] == strategy]
        xs = [int(row["n_blocks"]) for row in chunk]
        ys = np.asarray([float(row["recovery_probability"]) for row in chunk])
        low = np.asarray([float(row["ci_low"]) for row in chunk])
        high = np.asarray([float(row["ci_high"]) for row in chunk])
        ax.plot(xs, ys, marker="o", linewidth=1.5, color=colors[strategy], label=strategy)
        ax.fill_between(xs, low, high, color=colors[strategy], alpha=0.12, linewidth=0)
    ax.set_xscale("log")
    ax.set_xlabel("n_blocks")
    ax.set_ylabel("P(inferred quartet = T1)")
    ax.set_ylim(0.0, 1.0)
    ax.legend(fontsize=8)
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
    regression_summary = summarize_threshold_regression(
        threshold_rows,
        bootstrap_replicates=int(args.threshold_bootstrap_replicates),
        seed=int(args.seed),
    )
    kappa_rows, kappa_summary = run_kappa_grid(args)
    correction_summary = run_correction_grid(args)
    consistency_summary = run_consistency_benchmark(args)
    _write_csv(out / "threshold_grid_replicates.csv", THRESHOLD_REPLICATE_FIELDS, threshold_rows)
    _write_csv(out / "threshold_grid_summary.csv", THRESHOLD_SUMMARY_FIELDS, threshold_summary)
    _write_csv(out / "threshold_regression_summary.csv", THRESHOLD_REGRESSION_FIELDS, regression_summary)
    _write_csv(out / "kappa_calibration_replicates.csv", KAPPA_REPLICATE_FIELDS, kappa_rows)
    _write_csv(out / "kappa_calibration_summary.csv", KAPPA_SUMMARY_FIELDS, kappa_summary)
    _write_csv(out / "correction_grid_summary.csv", CORRECTION_SUMMARY_FIELDS, correction_summary)
    _write_csv(out / "consistency_summary.csv", CONSISTENCY_SUMMARY_FIELDS, consistency_summary)
    _plot_thresholds(regression_summary, out / "theory_vs_empirical_threshold.pdf")
    _plot_threshold_heatmap(regression_summary, out / "threshold_error_heatmap.pdf")
    _plot_recovery_phase(threshold_rows, out / "recovery_phase_diagram.pdf")
    _plot_kappa(kappa_summary, out / "kappa_calibration.pdf")
    _plot_corrections(correction_summary, out / "correction_methods_grid.pdf")
    _plot_consistency(consistency_summary, out / "consistency_vs_nblocks.pdf")
    _plot_consistency_corrections(consistency_summary, out / "consistency_corrections.pdf")
    return out
