from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from .analytic import TOPOLOGY_NAMES
from .model_fitting import msc_probabilities
from .robustness import infer_with_strategy, support_fraction


ROBUSTNESS_FIELDS = [
    "rearrangement_fraction", "strategy", "inferred_topology", "inferred_topology_index",
    "support_t1", "support_t2", "support_t3", "total_weight", "naive_flip",
]
GENEALOGY_FIELDS = [
    "rearrangement_fraction", "window_id", "chrom", "start", "end", "midpoint",
    "block_id", "topology", "topology_index", "is_rearranged", "rearrangement_id",
    "msrc_probability", "weight",
]


def _parse_floats(text: str) -> list[float]:
    return [float(x) for x in text.split(",") if x.strip()]


def _rows_for_fraction(
    fraction: float,
    *,
    chrom: str,
    length: float,
    windows: int,
    block_windows: int,
    rng: np.random.Generator,
    t1_branch: float,
    t2_probability: float,
) -> list[dict[str, Any]]:
    q_msc = msc_probabilities(0, t1_branch)
    q_msrc = np.asarray([(1.0 - t2_probability) / 2.0, t2_probability, (1.0 - t2_probability) / 2.0], dtype=float)
    window_size = length / windows
    rearranged_windows = int(round(windows * fraction))
    start_window = max(0, (windows - rearranged_windows) // 2)
    end_window = start_window + rearranged_windows
    rows: list[dict[str, Any]] = []
    current_block_key: tuple[bool, int] | None = None
    block_id = -1
    for window_id in range(windows):
        is_rearranged = start_window <= window_id < end_window
        local_index = (window_id - start_window) if is_rearranged else window_id
        block_key = (is_rearranged, local_index // max(1, block_windows))
        if block_key != current_block_key:
            block_id += 1
            current_block_key = block_key
            q = q_msrc if is_rearranged else q_msc
            block_topology = int(rng.choice(3, p=q))
        start = window_id * window_size
        end = min(length, (window_id + 1) * window_size)
        p_msrc = 1.0 if is_rearranged else 0.0
        rows.append({
            "rearrangement_fraction": float(fraction),
            "window_id": int(window_id),
            "chrom": chrom,
            "start": float(start),
            "end": float(end),
            "midpoint": float((start + end) / 2.0),
            "block_id": int(block_id),
            "topology": TOPOLOGY_NAMES[block_topology],
            "topology_index": int(block_topology),
            "is_rearranged": bool(is_rearranged),
            "rearrangement_id": "msrc_block" if is_rearranged else "",
            "msrc_probability": p_msrc,
            "weight": 1.0 - p_msrc,
        })
    return rows


def run_benchmark(args: argparse.Namespace) -> Path:
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(int(args.seed))
    fractions = _parse_floats(args.rearrangement_fractions)
    all_genealogies: list[dict[str, Any]] = []
    result_rows: list[dict[str, Any]] = []
    strategies = ["all_windows", "oracle_filter", "block_collapse", "soft_weight"]
    for fraction in fractions:
        rows = _rows_for_fraction(
            fraction,
            chrom=args.chrom,
            length=float(args.chrom_length),
            windows=int(args.windows),
            block_windows=int(args.block_windows),
            rng=rng,
            t1_branch=float(args.t1_branch),
            t2_probability=float(args.t2_probability),
        )
        all_genealogies.extend(rows)
        naive = infer_with_strategy(rows, "all_windows")
        naive_flip = support_fraction(naive, 1) > support_fraction(naive, 0)
        for strategy in strategies:
            result = infer_with_strategy(rows, strategy)
            result_rows.append({
                "rearrangement_fraction": float(fraction),
                "strategy": strategy,
                "inferred_topology": result.inferred_topology,
                "inferred_topology_index": int(result.inferred_topology_index),
                "support_t1": support_fraction(result, 0),
                "support_t2": support_fraction(result, 1),
                "support_t3": support_fraction(result, 2),
                "total_weight": result.total_weight,
                "naive_flip": bool(naive_flip),
            })
    with (out / "spatial_genealogies.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GENEALOGY_FIELDS)
        writer.writeheader()
        writer.writerows({k: row[k] for k in GENEALOGY_FIELDS} for row in all_genealogies)
    with (out / "species_tree_robustness.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROBUSTNESS_FIELDS)
        writer.writeheader()
        writer.writerows(result_rows)
    with (out / "gene_trees.nwk").open("w") as handle:
        for row in all_genealogies:
            top = int(row["topology_index"])
            if top == 0:
                handle.write("((1,2),(3,4));\n")
            elif top == 1:
                handle.write("((1,3),(2,4));\n")
            else:
                handle.write("((1,4),(2,3));\n")
    _plot(result_rows, out / "species_tree_robustness.png")
    return out


def _plot(rows: list[dict[str, Any]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    styles = {
        "all_windows": "-",
        "oracle_filter": "--",
        "block_collapse": "-.",
        "soft_weight": ":",
    }
    for strategy, linestyle in styles.items():
        chunk = [row for row in rows if row["strategy"] == strategy]
        xs = [float(row["rearrangement_fraction"]) for row in chunk]
        ax.plot(xs, [float(row["support_t1"]) for row in chunk], linestyle=linestyle, color="#1b9e77", label=f"{strategy} T1")
        ax.plot(xs, [float(row["support_t2"]) for row in chunk], linestyle=linestyle, color="#d95f02", label=f"{strategy} T2")
    flips = [row for row in rows if row["strategy"] == "all_windows" and row["naive_flip"]]
    if flips:
        ax.axvline(float(flips[0]["rearrangement_fraction"]), color="black", linewidth=1.0)
    ax.set_xlabel("Rearrangement fraction")
    ax.set_ylabel("Quartet support")
    ax.set_ylim(0.0, 1.0)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the v0.8.0 linked-spatial species-tree robustness benchmark")
    parser.add_argument("--output-dir", default="robustness_output")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--chrom-length", type=float, default=1_000_000.0)
    parser.add_argument("--windows", type=int, default=400)
    parser.add_argument("--block-windows", type=int, default=20)
    parser.add_argument("--rearrangement-fractions", default="0,0.1,0.2,0.3,0.4,0.5,0.6")
    parser.add_argument("--t1-branch", type=float, default=1.2)
    parser.add_argument("--t2-probability", type=float, default=0.95)
    args = parser.parse_args()
    out = run_benchmark(args)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
