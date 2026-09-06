from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from .analytic import TOPOLOGY_NAMES
from .genomic import GenomicInterval
from .linked_spatial import generate_genealogy_breakpoints
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
    "end", "midpoint", "block_id", "block_start", "block_end", "topology", "topology_index",
    "is_rearranged", "rearrangement_id", "msrc_probability", "weight",
]
LINKAGE_DIAGNOSTIC_FIELDS = [
    "rearrangement_fraction", "replicate_id", "region", "num_windows",
    "num_genealogy_blocks", "mean_block_length_bp", "median_block_length_bp",
    "breakpoint_density_per_bp", "observed_inside_outside_rate_ratio", "expected_kappa",
    "fraction_rearrangement_intervals_with_1_block",
    "fraction_rearrangement_intervals_with_2_blocks",
    "fraction_rearrangement_intervals_with_3plus_blocks",
]
LINKAGE_SUMMARY_FIELDS = [
    "rearrangement_fraction", "region", "replicates", "mean_num_genealogy_blocks",
    "mean_block_length_bp", "median_block_length_bp", "breakpoint_density_per_bp",
    "observed_inside_outside_rate_ratio", "expected_kappa",
    "fraction_rearrangement_intervals_with_1_block",
    "fraction_rearrangement_intervals_with_2_blocks",
    "fraction_rearrangement_intervals_with_3plus_blocks",
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
    if not rows:
        return []
    length = max(float(row["end"]) for row in rows)
    chrom = str(rows[0]["chrom"])
    window_size = length / len(rows)
    background_rate = 1.0 / (max(1, block_windows) * window_size)
    blocks = _sample_background_blocks(chrom=chrom, length=length, background_rate=background_rate, q_msc=q_msc, rng=rng)
    out = []
    seen_blocks: set[int] = set()
    for row in rows:
        block = _block_at(blocks, float(row["midpoint"]))
        block_idx = blocks.index(block)
        item = dict(row)
        item["base_block_key"] = block_idx
        item["base_block_start"] = float(block["start"])
        item["base_block_end"] = float(block["end"])
        item["base_block_first"] = block_idx not in seen_blocks
        item["base_topology_index"] = int(block["topology_index"])
        seen_blocks.add(block_idx)
        out.append(item)
    return out


def _sample_blocks(
    *,
    chrom: str,
    length: float,
    intervals: list[GenomicInterval],
    background_rate: float,
    kappa: float,
    q_msc: np.ndarray,
    q_msrc: np.ndarray,
    rng: np.random.Generator,
) -> list[dict[str, Any]]:
    breakpoints = generate_genealogy_breakpoints(length, intervals, background_rate, kappa, rng)
    blocks = []
    for start, end in zip(breakpoints, breakpoints[1:]):
        midpoint = (start + end) / 2.0
        interval = next((item for item in intervals if item.contains(midpoint)), None)
        is_rearranged = interval is not None
        topology = int(rng.choice(3, p=q_msrc if is_rearranged else q_msc))
        blocks.append({
            "chrom": chrom,
            "start": float(start),
            "end": float(end),
            "is_rearranged": bool(is_rearranged),
            "rearrangement_id": interval.interval_id if interval else "",
            "topology_index": topology,
        })
    return blocks


def _sample_background_blocks(
    *,
    chrom: str,
    length: float,
    background_rate: float,
    q_msc: np.ndarray,
    rng: np.random.Generator,
) -> list[dict[str, Any]]:
    breakpoints = generate_genealogy_breakpoints(length, [], background_rate, 1.0, rng)
    blocks = []
    for start, end in zip(breakpoints, breakpoints[1:]):
        blocks.append({
            "chrom": chrom,
            "start": float(start),
            "end": float(end),
            "is_rearranged": False,
            "rearrangement_id": "",
            "topology_index": int(rng.choice(3, p=q_msc)),
        })
    return blocks


def _block_at(blocks: list[dict[str, Any]], position: float) -> dict[str, Any]:
    starts = np.asarray([float(block["start"]) for block in blocks], dtype=float)
    idx = int(np.searchsorted(starts, position, side="right") - 1)
    return blocks[min(max(idx, 0), len(blocks) - 1)]


def _split_baseline_blocks(
    baseline_blocks: list[dict[str, Any]],
    intervals: list[GenomicInterval],
    q_msrc: np.ndarray,
    rng: np.random.Generator,
    background_rate: float,
    kappa: float,
    length: float,
) -> list[dict[str, Any]]:
    if not intervals:
        return [dict(block) for block in baseline_blocks]
    points = {0.0, float(length)}
    for block in baseline_blocks:
        for point in (float(block["start"]), float(block["end"])):
            if not any(interval.start < point < interval.end for interval in intervals):
                points.add(point)
    for interval in intervals:
        points.add(interval.start)
        points.add(interval.end)
        inside_bps = generate_genealogy_breakpoints(interval.end - interval.start, [], background_rate * kappa, 1.0, rng)
        for bp in inside_bps[1:-1]:
            points.add(interval.start + bp)
    ordered = sorted(points)
    out = []
    for start, end in zip(ordered, ordered[1:]):
        midpoint = (start + end) / 2.0
        interval = next((item for item in intervals if item.contains(midpoint)), None)
        if interval is None:
            base = _block_at(baseline_blocks, midpoint)
            topology = int(base["topology_index"])
            rearrangement_id = ""
            is_rearranged = False
        else:
            topology = int(rng.choice(3, p=q_msrc))
            rearrangement_id = interval.interval_id
            is_rearranged = True
        out.append({
            "chrom": baseline_blocks[0]["chrom"],
            "start": float(start),
            "end": float(end),
            "is_rearranged": is_rearranged,
            "rearrangement_id": rearrangement_id,
            "topology_index": topology,
        })
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
    kappa: float = 0.25,
    baseline_rows: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    skeleton = baseline_rows if baseline_rows is not None else _make_window_skeleton(chrom, length, windows)
    rearranged_windows = int(round(windows * fraction))
    start_window = max(0, (windows - rearranged_windows) // 2)
    end_window = start_window + rearranged_windows
    window_size = length / windows
    background_rate = 1.0 / (max(1, block_windows) * window_size)
    interval_start = start_window * window_size
    interval_end = end_window * window_size
    intervals = []
    if interval_end > interval_start:
        intervals.append(GenomicInterval(chrom, interval_start, interval_end, "central_msrc_interval"))
    if baseline_rows is None:
        blocks = _sample_blocks(
            chrom=chrom,
            length=length,
            intervals=intervals,
            background_rate=background_rate,
            kappa=kappa,
            q_msc=q_msc,
            q_msrc=q_msrc,
            rng=rng,
        )
    elif baseline_rows and "base_block_start" in baseline_rows[0]:
        baseline_blocks = [
            {
                "chrom": chrom,
                "start": float(row["base_block_start"]),
                "end": float(row["base_block_end"]),
                "topology_index": int(row["base_topology_index"]),
            }
            for row in baseline_rows
            if bool(row.get("base_block_first", False))
        ]
        blocks = _split_baseline_blocks(baseline_blocks, intervals, q_msrc, rng, background_rate, kappa, length)
    else:
        blocks = _sample_blocks(
            chrom=chrom,
            length=length,
            intervals=intervals,
            background_rate=background_rate,
            kappa=kappa,
            q_msc=q_msc,
            q_msrc=q_msrc,
            rng=rng,
        )
    rows: list[dict[str, Any]] = []
    used_blocks: dict[int, int] = {}
    for base in skeleton:
        window_id = int(base["window_id"])
        block_idx = blocks.index(_block_at(blocks, float(base["midpoint"])))
        if block_idx not in used_blocks:
            used_blocks[block_idx] = len(used_blocks)
        block = blocks[block_idx]
        block_id = used_blocks[block_idx]
        is_rearranged = bool(block["is_rearranged"])
        block_topology = int(block["topology_index"])
        rearrangement_id = str(block["rearrangement_id"])
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
            "block_start": float(block["start"]),
            "block_end": float(block["end"]),
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


def _linkage_diagnostics(rows: list[dict[str, Any]], *, fraction: float, replicate_id: int) -> list[dict[str, Any]]:
    out = []
    by_region: dict[str, dict[str, Any]] = {}
    ordered_rows = sorted(rows, key=lambda row: int(row["window_id"]))
    for region, rearranged in (("inside", True), ("outside", False)):
        chunk = [row for row in rows if bool(row["is_rearranged"]) is rearranged]
        block_ids = sorted({int(row["block_id"]) for row in chunk})
        if chunk:
            length = sum(float(row["end"]) - float(row["start"]) for row in chunk)
            blocks = len(block_ids)
            breakpoints = 0
            run_blocks: set[int] = set()
            for row in ordered_rows:
                if bool(row["is_rearranged"]) is rearranged:
                    run_blocks.add(int(row["block_id"]))
                elif run_blocks:
                    breakpoints += max(0, len(run_blocks) - 1)
                    run_blocks = set()
            if run_blocks:
                breakpoints += max(0, len(run_blocks) - 1)
            block_lengths = [
                max(float(row["block_end"]) - float(row["block_start"]) for row in chunk if int(row["block_id"]) == block_id)
                for block_id in block_ids
            ]
            mean_length = float(np.mean(block_lengths)) if block_lengths else float("nan")
            median_length = float(np.median(block_lengths)) if block_lengths else float("nan")
            density = breakpoints / length if length > 0.0 else float("nan")
        else:
            blocks = 0
            mean_length = float("nan")
            median_length = float("nan")
            density = float("nan")
        by_region[region] = {
            "rearrangement_fraction": float(fraction),
            "replicate_id": int(replicate_id),
            "region": region,
            "num_windows": int(len(chunk)),
            "num_genealogy_blocks": int(blocks),
            "mean_block_length_bp": float(mean_length),
            "median_block_length_bp": float(median_length),
            "breakpoint_density_per_bp": float(density),
        }
    inside_density = float(by_region["inside"]["breakpoint_density_per_bp"])
    outside_density = float(by_region["outside"]["breakpoint_density_per_bp"])
    ratio = inside_density / outside_density if np.isfinite(inside_density) and np.isfinite(outside_density) and outside_density > 0.0 else float("nan")
    interval_counts: dict[str, set[int]] = {}
    for row in rows:
        if row["rearrangement_id"]:
            interval_counts.setdefault(str(row["rearrangement_id"]), set()).add(int(row["block_id"]))
    counts = [len(value) for value in interval_counts.values()]
    denom = len(counts)
    fractions = {
        "fraction_rearrangement_intervals_with_1_block": sum(value == 1 for value in counts) / denom if denom else float("nan"),
        "fraction_rearrangement_intervals_with_2_blocks": sum(value == 2 for value in counts) / denom if denom else float("nan"),
        "fraction_rearrangement_intervals_with_3plus_blocks": sum(value >= 3 for value in counts) / denom if denom else float("nan"),
    }
    expected = ""
    for row in by_region.values():
        row["observed_inside_outside_rate_ratio"] = ratio
        row["expected_kappa"] = expected
        row.update(fractions)
        out.append(row)
    return out


def _aggregate_linkage_diagnostics(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    def mean_or_nan(values: list[float]) -> float:
        finite = [value for value in values if np.isfinite(value)]
        return float(np.mean(finite)) if finite else float("nan")

    out = []
    keys = sorted({(float(row["rearrangement_fraction"]), row["region"]) for row in rows})
    for fraction, region in keys:
        chunk = [row for row in rows if float(row["rearrangement_fraction"]) == fraction and row["region"] == region]
        finite_ratio = [float(row["observed_inside_outside_rate_ratio"]) for row in chunk if row["observed_inside_outside_rate_ratio"] != "" and np.isfinite(float(row["observed_inside_outside_rate_ratio"]))]
        expected = next((row["expected_kappa"] for row in chunk if row["expected_kappa"] != ""), "")
        out.append({
            "rearrangement_fraction": fraction,
            "region": region,
            "replicates": len(chunk),
            "mean_num_genealogy_blocks": float(np.mean([float(row["num_genealogy_blocks"]) for row in chunk])),
            "mean_block_length_bp": mean_or_nan([float(row["mean_block_length_bp"]) for row in chunk]),
            "median_block_length_bp": mean_or_nan([float(row["median_block_length_bp"]) for row in chunk]),
            "breakpoint_density_per_bp": mean_or_nan([float(row["breakpoint_density_per_bp"]) for row in chunk]),
            "observed_inside_outside_rate_ratio": float(np.mean(finite_ratio)) if finite_ratio else float("nan"),
            "expected_kappa": expected,
            "fraction_rearrangement_intervals_with_1_block": mean_or_nan([float(row["fraction_rearrangement_intervals_with_1_block"]) for row in chunk]),
            "fraction_rearrangement_intervals_with_2_blocks": mean_or_nan([float(row["fraction_rearrangement_intervals_with_2_blocks"]) for row in chunk]),
            "fraction_rearrangement_intervals_with_3plus_blocks": mean_or_nan([float(row["fraction_rearrangement_intervals_with_3plus_blocks"]) for row in chunk]),
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
    diagnostic_rows: list[dict[str, Any]] = []
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
                kappa=float(args.kappa),
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
            replicate_diagnostics = _linkage_diagnostics(rows, fraction=fraction, replicate_id=replicate_id)
            for row in replicate_diagnostics:
                row["expected_kappa"] = float(args.kappa)
            diagnostic_rows.extend(replicate_diagnostics)
    recovery_rows = _aggregate_recovery(result_rows)
    diagnostic_summary_rows = _aggregate_linkage_diagnostics(diagnostic_rows)
    _write_outputs(out, all_genealogies, result_rows, recovery_rows, diagnostic_rows, diagnostic_summary_rows)
    _plot_support(result_rows, out / "quartet_support_vs_rearrangement_fraction.pdf", threshold)
    _plot_recovery(recovery_rows, out / "species_tree_recovery_vs_rearrangement_fraction.pdf", threshold)
    return out


def _write_outputs(
    out: Path,
    genealogies: list[dict[str, Any]],
    result_rows: list[dict[str, Any]],
    recovery_rows: list[dict[str, Any]],
    diagnostic_rows: list[dict[str, Any]],
    diagnostic_summary_rows: list[dict[str, Any]],
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
    with (out / "spatial_linkage_diagnostics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LINKAGE_DIAGNOSTIC_FIELDS)
        writer.writeheader()
        writer.writerows(diagnostic_rows)
    with (out / "spatial_linkage_diagnostic_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LINKAGE_SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(diagnostic_summary_rows)
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
    parser = argparse.ArgumentParser(description="Run the linked-spatial species-tree robustness benchmark")
    parser.add_argument("--output-dir", default="robustness_output")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--mode", choices=["independent", "paired"], default="independent")
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--chrom-length", type=float, default=1_000_000.0)
    parser.add_argument("--windows", type=int, default=400)
    parser.add_argument("--block-windows", type=int, default=20, help="Background MSC genealogy block size in windows")
    parser.add_argument("--msrc-block-windows", type=int, default=5, help="MSRC genealogy block size in windows inside the rearrangement")
    parser.add_argument("--kappa", type=float, default=0.25, help="Multiplier on the rearranged genealogy breakpoint rate")
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
