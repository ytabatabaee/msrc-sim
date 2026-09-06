from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from .analytic import TOPOLOGY_NAMES
from .model_fitting import msc_probabilities
from .robustness import (
    binomial_confidence_interval,
    infer_with_strategy,
    support_fraction,
    theoretical_flip_threshold,
)


STRATEGIES = [
    "all_windows",
    "oracle_filter",
    "genealogy_block_collapse",
    "rearrangement_interval_collapse",
    "soft_weight",
]

ROBUSTNESS_FIELDS = [
    "rearrangement_fraction", "replicate_id", "mode", "strategy",
    "inferred_topology", "inferred_topology_index", "recovered_true_t1",
    "support_t1", "support_t2", "support_t3", "total_weight",
    "num_windows", "num_genealogy_blocks", "num_rearrangement_intervals",
    "naive_flip", "theoretical_flip_threshold",
]
RECOVERY_FIELDS = [
    "rearrangement_fraction", "mode", "strategy", "replicates",
    "recovery_probability", "ci_low", "ci_high",
]
GENEALOGY_FIELDS = [
    "rearrangement_fraction", "replicate_id", "window_id", "chrom", "start",
    "end", "midpoint", "block_id", "topology", "topology_index",
    "is_rearranged", "rearrangement_id", "msrc_probability", "weight",
]


def _parse_floats(text: str) -> list[float]:
    return [float(x) for x in text.split(",") if x.strip()]


def _default_fractions() -> list[float]:
    return [i / 100.0 for i in range(0, 61, 5)]


def _sample_probabilities(
    is_rearranged: bool,
    *,
    mode: str,
    sensitivity: float,
    specificity: float,
    noise_sd: float,
    rng: np.random.Generator,
) -> float:
    oracle = 1.0 if is_rearranged else 0.0
    if mode == "oracle":
        return oracle
    if mode != "noisy":
        raise ValueError("soft probability mode must be oracle or noisy")
    mean = sensitivity if is_rearranged else 1.0 - specificity
    return float(np.clip(mean + rng.normal(0.0, noise_sd), 0.0, 1.0))


def _newick_for_topology(topology: int) -> str:
    if topology == 0:
        return "((1,2),(3,4));"
    if topology == 1:
        return "((1,3),(2,4));"
    return "((1,4),(2,3));"


def _make_window_skeleton(chrom: str, length: float, windows: int) -> list[dict[str, Any]]:
    window_size = length / windows
    rows = []
    for window_id in range(windows):
        start = window_id * window_size
        end = min(length, (window_id + 1) * window_size)
        rows.append({
            "window_id": int(window_id),
            "chrom": chrom,
            "start": float(start),
            "end": float(end),
            "midpoint": float((start + end) / 2.0),
        })
    return rows


def _assign_background(
    rows: list[dict[str, Any]],
    *,
    q_msc: np.ndarray,
    block_windows: int,
    rng: np.random.Generator,
) -> list[dict[str, Any]]:
    out = []
    current_key = None
    block_topology = 0
    for row in rows:
        key = int(row["window_id"]) // max(1, block_windows)
        if key != current_key:
            current_key = key
            block_topology = int(rng.choice(3, p=q_msc))
        item = dict(row)
        item["base_block_key"] = key
        item["base_topology_index"] = block_topology
        out.append(item)
    return out


def _rows_for_fraction(
    fraction: float,
    *,
    replicate_id: int,
    chrom: str,
    length: float,
    windows: int,
    block_windows: int,
    msrc_block_windows: int,
    rng: np.random.Generator,
    q_msc: np.ndarray,
    q_msrc: np.ndarray,
    soft_probability_mode: str,
    soft_sensitivity: float,
    soft_specificity: float,
    soft_noise_sd: float,
    baseline_rows: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    skeleton = baseline_rows if baseline_rows is not None else _assign_background(
        _make_window_skeleton(chrom, length, windows),
        q_msc=q_msc,
        block_windows=block_windows,
        rng=rng,
    )
    rearranged_windows = int(round(windows * fraction))
    start_window = max(0, (windows - rearranged_windows) // 2)
    end_window = start_window + rearranged_windows
    rows: list[dict[str, Any]] = []
    current_key: tuple[str, int] | None = None
    block_id = -1
    block_topology = 0
    for base in skeleton:
        window_id = int(base["window_id"])
        is_rearranged = start_window <= window_id < end_window
        if is_rearranged:
            local_index = window_id - start_window
            key = ("rearrangement", local_index // max(1, msrc_block_windows))
            if key != current_key:
                block_id += 1
                current_key = key
                block_topology = int(rng.choice(3, p=q_msrc))
            rearrangement_id = "central_msrc_interval"
        else:
            key = ("background", int(base["base_block_key"]))
            if key != current_key:
                block_id += 1
                current_key = key
            block_topology = int(base["base_topology_index"])
            rearrangement_id = ""
        p_msrc = _sample_probabilities(
            is_rearranged,
            mode=soft_probability_mode,
            sensitivity=soft_sensitivity,
            specificity=soft_specificity,
            noise_sd=soft_noise_sd,
            rng=rng,
        )
        item = {
            "rearrangement_fraction": float(fraction),
            "replicate_id": int(replicate_id),
            "window_id": window_id,
            "chrom": str(base["chrom"]),
            "start": float(base["start"]),
            "end": float(base["end"]),
            "midpoint": float(base["midpoint"]),
            "block_id": int(block_id),
            "topology": TOPOLOGY_NAMES[block_topology],
            "topology_index": int(block_topology),
            "is_rearranged": bool(is_rearranged),
            "rearrangement_id": rearrangement_id,
            "msrc_probability": float(p_msrc),
            "weight": float(1.0 - p_msrc),
        }
        rows.append(item)
    return rows


def _summarize_replicate(
    rows: list[dict[str, Any]],
    *,
    fraction: float,
    replicate_id: int,
    mode: str,
    threshold: float | None,
) -> list[dict[str, Any]]:
    naive = infer_with_strategy(rows, "all_windows")
    naive_flip = support_fraction(naive, 1) > support_fraction(naive, 0)
    num_intervals = len({row["rearrangement_id"] for row in rows if row["rearrangement_id"]})
    out = []
    for strategy in STRATEGIES:
        result = infer_with_strategy(rows, strategy)
        out.append({
            "rearrangement_fraction": float(fraction),
            "replicate_id": int(replicate_id),
            "mode": mode,
            "strategy": strategy,
            "inferred_topology": result.inferred_topology,
            "inferred_topology_index": int(result.inferred_topology_index),
            "recovered_true_t1": result.inferred_topology_index == 0,
            "support_t1": support_fraction(result, 0),
            "support_t2": support_fraction(result, 1),
            "support_t3": support_fraction(result, 2),
            "total_weight": result.total_weight,
            "num_windows": len(rows),
            "num_genealogy_blocks": len({int(row["block_id"]) for row in rows}),
            "num_rearrangement_intervals": num_intervals,
            "naive_flip": bool(naive_flip),
            "theoretical_flip_threshold": "" if threshold is None else float(threshold),
        })
    return out


def _aggregate_recovery(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    keys = sorted({(float(row["rearrangement_fraction"]), row["mode"], row["strategy"]) for row in rows})
    for fraction, mode, strategy in keys:
        chunk = [row for row in rows if float(row["rearrangement_fraction"]) == fraction and row["mode"] == mode and row["strategy"] == strategy]
        successes = sum(bool(row["recovered_true_t1"]) for row in chunk)
        low, high = binomial_confidence_interval(successes, len(chunk))
        out.append({
            "rearrangement_fraction": fraction,
            "mode": mode,
            "strategy": strategy,
            "replicates": len(chunk),
            "recovery_probability": successes / len(chunk) if chunk else float("nan"),
            "ci_low": low,
            "ci_high": high,
        })
    return out


def run_benchmark(args: argparse.Namespace) -> Path:
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(int(args.seed))
    fractions = _parse_floats(args.rearrangement_fractions) if args.rearrangement_fractions else _default_fractions()
    q_msc = msc_probabilities(0, float(args.t1_branch))
    q_msrc = np.asarray([(1.0 - float(args.t2_probability)) / 2.0, float(args.t2_probability), (1.0 - float(args.t2_probability)) / 2.0], dtype=float)
    threshold = theoretical_flip_threshold(q_msc[0], q_msc[1], q_msrc[0], q_msrc[1])
    all_genealogies: list[dict[str, Any]] = []
    result_rows: list[dict[str, Any]] = []
    for replicate_id in range(int(args.replicates)):
        baseline_rows = None
        if args.mode == "paired":
            baseline_rows = _assign_background(
                _make_window_skeleton(args.chrom, float(args.chrom_length), int(args.windows)),
                q_msc=q_msc,
                block_windows=int(args.block_windows),
                rng=rng,
            )
        for fraction in fractions:
            rows = _rows_for_fraction(
                fraction,
                replicate_id=replicate_id,
                chrom=args.chrom,
                length=float(args.chrom_length),
                windows=int(args.windows),
                block_windows=int(args.block_windows),
                msrc_block_windows=int(args.msrc_block_windows),
                rng=rng,
                q_msc=q_msc,
                q_msrc=q_msrc,
                soft_probability_mode=args.soft_probability_mode,
                soft_sensitivity=float(args.soft_sensitivity),
                soft_specificity=float(args.soft_specificity),
                soft_noise_sd=float(args.soft_noise_sd),
                baseline_rows=baseline_rows,
            )
            all_genealogies.extend(rows)
            result_rows.extend(_summarize_replicate(rows, fraction=fraction, replicate_id=replicate_id, mode=args.mode, threshold=threshold))
    recovery_rows = _aggregate_recovery(result_rows)
    _write_outputs(out, all_genealogies, result_rows, recovery_rows)
    _plot_support(result_rows, out / "quartet_support_vs_rearrangement_fraction.pdf", threshold)
    _plot_recovery(recovery_rows, out / "species_tree_recovery_vs_rearrangement_fraction.pdf", threshold)
    return out


def _write_outputs(
    out: Path,
    genealogies: list[dict[str, Any]],
    result_rows: list[dict[str, Any]],
    recovery_rows: list[dict[str, Any]],
) -> None:
    with (out / "spatial_genealogies.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GENEALOGY_FIELDS)
        writer.writeheader()
        writer.writerows({k: row[k] for k in GENEALOGY_FIELDS} for row in genealogies)
    with (out / "species_tree_robustness.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROBUSTNESS_FIELDS)
        writer.writeheader()
        writer.writerows(result_rows)
    with (out / "species_tree_recovery.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RECOVERY_FIELDS)
        writer.writeheader()
        writer.writerows(recovery_rows)
    with (out / "gene_trees.nwk").open("w") as handle:
        for row in genealogies:
            handle.write(_newick_for_topology(int(row["topology_index"])) + "\n")


def _mean_by_fraction(rows: list[dict[str, Any]], strategy: str, field: str) -> tuple[list[float], list[float]]:
    fractions = sorted({float(row["rearrangement_fraction"]) for row in rows})
    xs = []
    ys = []
    for fraction in fractions:
        chunk = [float(row[field]) for row in rows if float(row["rearrangement_fraction"]) == fraction and row["strategy"] == strategy]
        xs.append(fraction)
        ys.append(float(np.mean(chunk)))
    return xs, ys


def _plot_support(rows: list[dict[str, Any]], path: Path, threshold: float | None) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    colors = {"support_t1": "#1b9e77", "support_t2": "#d95f02", "support_t3": "#7570b3"}
    styles = {
        "all_windows": "-",
        "oracle_filter": "--",
        "genealogy_block_collapse": "-.",
        "rearrangement_interval_collapse": (0, (3, 1, 1, 1)),
        "soft_weight": ":",
    }
    for strategy, linestyle in styles.items():
        for field, color in colors.items():
            xs, ys = _mean_by_fraction(rows, strategy, field)
            ax.plot(xs, ys, linestyle=linestyle, color=color, linewidth=1.2, label=f"{strategy} {field[-2:].upper()}")
    if threshold is not None:
        ax.axvline(threshold, color="black", linewidth=1.0)
    ax.set_xlabel("Rearrangement fraction")
    ax.set_ylabel("Mean weighted quartet support")
    ax.set_ylim(0.0, 1.0)
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_recovery(rows: list[dict[str, Any]], path: Path, threshold: float | None) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    colors = {
        "all_windows": "#d95f02",
        "oracle_filter": "#1b9e77",
        "genealogy_block_collapse": "#7570b3",
        "rearrangement_interval_collapse": "#e7298a",
        "soft_weight": "#66a61e",
    }
    for strategy in STRATEGIES:
        chunk = [row for row in rows if row["strategy"] == strategy]
        xs = [float(row["rearrangement_fraction"]) for row in chunk]
        ys = np.asarray([float(row["recovery_probability"]) for row in chunk])
        low = np.asarray([float(row["ci_low"]) for row in chunk])
        high = np.asarray([float(row["ci_high"]) for row in chunk])
        ax.plot(xs, ys, marker="o", markersize=2.5, linewidth=1.5, color=colors[strategy], label=strategy)
        ax.fill_between(xs, low, high, color=colors[strategy], alpha=0.12, linewidth=0)
    if threshold is not None:
        ax.axvline(threshold, color="black", linewidth=1.0)
    ax.set_xlabel("Rearrangement fraction")
    ax.set_ylabel("P(inferred quartet = T1)")
    ax.set_ylim(0.0, 1.0)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the v0.8.1 linked-spatial species-tree robustness benchmark")
    parser.add_argument("--output-dir", default="robustness_output")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--mode", choices=["independent", "paired"], default="independent")
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--chrom-length", type=float, default=1_000_000.0)
    parser.add_argument("--windows", type=int, default=400)
    parser.add_argument("--block-windows", type=int, default=20, help="Background MSC genealogy block size in windows")
    parser.add_argument("--msrc-block-windows", type=int, default=5, help="MSRC genealogy block size in windows inside the rearrangement")
    parser.add_argument("--rearrangement-fractions", default=",".join(f"{x:.2f}" for x in _default_fractions()))
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--t1-branch", type=float, default=1.2)
    parser.add_argument("--t2-probability", type=float, default=0.95)
    parser.add_argument("--soft-probability-mode", choices=["oracle", "noisy"], default="noisy")
    parser.add_argument("--soft-sensitivity", type=float, default=0.85)
    parser.add_argument("--soft-specificity", type=float, default=0.90)
    parser.add_argument("--soft-noise-sd", type=float, default=0.05)
    args = parser.parse_args()
    out = run_benchmark(args)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
