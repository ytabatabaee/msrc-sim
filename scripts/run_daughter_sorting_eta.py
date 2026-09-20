from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping

import matplotlib.pyplot as plt
import numpy as np
import yaml

from msrcsim.conditioning import branch_end_status, terminal_pattern
from msrcsim.moran import TIME_CONVENTION
from msrcsim.population_process import simulate_population_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.species_tree import SpeciesTree
from msrcsim.speciation_split import daughter_target_frequency
from msrcsim.statistics import wilson_interval
from msrcsim.structured_coalescent import simulate_genealogy


DISCORDANT_PATTERNS = {"0101", "1010", "0110", "1001"}
ALL_PATTERNS = [format(i, "04b") for i in range(16)]
SPECIES_TREE_TOPOLOGY = 0
ARRANGEMENT_TOPOLOGY_BY_PATTERN = {
    "0101": 1,
    "1010": 1,
    "0110": 2,
    "1001": 2,
}


def _git_commit() -> str | None:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return None


def _git_status() -> str:
    try:
        return subprocess.check_output(["git", "status", "--short"], text=True, stderr=subprocess.DEVNULL)
    except Exception as exc:
        return f"unavailable: {exc}"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_csv(path: Path, rows: list[Mapping[str, Any]], fields: list[str] | None = None) -> None:
    if not rows:
        path.write_text("")
        return
    fields = fields or list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fields} for row in rows)


def _tree_and_rearrangement(config: Mapping[str, Any]) -> tuple[SpeciesTree, Rearrangement]:
    base = config["baseline"]
    tree = SpeciesTree(
        base["species_tree_newick"],
        int(base["default_effective_population_size"]),
        float(base["root_extension_generations"]),
    )
    rr = base["rearrangement"]
    rearrangement = Rearrangement(
        str(rr["id"]),
        str(rr["type"]),
        str(rr["origin_branch"]),
        int(rr["origin_time_from_branch_start"]),
        int(rr["initial_copy_count"]),
        float(rr["selection_coefficient"]),
    )
    return tree, rearrangement


def _base_sim_config(config: Mapping[str, Any], *, process: str = "wright_fisher") -> dict[str, Any]:
    base = config["baseline"]
    rr = base["rearrangement"]
    return {
        "species_tree": {
            "newick": base["species_tree_newick"],
            "root_extension": base["root_extension_generations"],
            "default_effective_population_size": base["default_effective_population_size"],
        },
        "rearrangement": {
            "id": rr["id"],
            "type": rr["type"],
            "origin_branch": rr["origin_branch"],
            "origin_time_from_branch_start": rr["origin_time_from_branch_start"],
            "initial_copy_count": rr["initial_copy_count"],
            "selection": {"model": "genic", "coefficient": rr["selection_coefficient"]},
        },
        "recombination": dict(base["recombination"]),
        "population_process": {"model": process},
    }


def _split_policy_config(config: Mapping[str, Any], eta: float, mode: str) -> dict[str, Any]:
    if mode == "descendant":
        policy = yaml.safe_load(yaml.safe_dump(config["eta_design"]["primary_split_policy"]))
        policy["by_split"]["A"]["eta"] = float(eta)
        policy["by_split"]["B"]["eta"] = float(eta)
        return policy
    if mode == "root_only":
        policy = yaml.safe_load(yaml.safe_dump(config["eta_design"]["root_control_split_policy"]))
        policy["by_split"]["ROOT"]["eta"] = float(eta)
        return policy
    raise ValueError("split policy mode must be descendant or root_only")


def _simulate_one(
    tree: SpeciesTree,
    rearrangement: Rearrangement,
    base_config: Mapping[str, Any],
    eta: float,
    rng: np.random.Generator,
    *,
    split_mode: str,
) -> tuple[Any, dict[str, int], str, bool, bool, bool]:
    sim_config = dict(base_config)
    sim_config["speciation_split"] = _split_policy_config(base_config["_experiment_config"], eta, split_mode)
    history = simulate_population_history(tree, rearrangement, rng, sim_config)
    sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
    pattern = terminal_pattern(sampled, tree.taxa)
    persistent = branch_end_status(history, "A") == "segregating" and branch_end_status(history, "B") == "segregating"
    strict = pattern.count("1") == 2
    discordant = pattern in DISCORDANT_PATTERNS
    return history, sampled, pattern, persistent, strict, discordant


def _event_diagnostics(history: Any, split_ids: set[str]) -> tuple[list[dict[str, Any]], dict[str, float]]:
    by_parent: dict[str, dict[str, Any]] = {}
    rows = []
    for event in getattr(history, "split_events", []):
        row = asdict(event) if not isinstance(event, Mapping) else dict(event)
        if row["parent_branch_id"] in split_ids and row["parent_state"] == "segregating":
            by_parent.setdefault(str(row["parent_branch_id"]), {})[str(row["child_side_or_sign"])] = row
    for parent, sides in by_parent.items():
        if "positive" not in sides or "negative" not in sides:
            continue
        pos = sides["positive"]
        neg = sides["negative"]
        pos_freq = float(pos["realized_child_A1_count"]) / float(pos["child_total"])
        neg_freq = float(neg["realized_child_A1_count"]) / float(neg["child_total"])
        parent_freq = float(pos["parent_frequency_A1"])
        rows.append({
            "split": parent,
            "parent_frequency_A1": parent_freq,
            "target_frequency_positive": float(pos["target_frequency_A1"]),
            "target_frequency_negative": float(neg["target_frequency_A1"]),
            "realized_positive_frequency": pos_freq,
            "realized_negative_frequency": neg_freq,
            "realized_daughter_difference": abs(pos_freq - neg_freq),
            "opposite_enrichment": bool(pos_freq > parent_freq and neg_freq < parent_freq),
        })
    if rows:
        aggregate = {
            "mean_realized_daughter_difference": float(np.mean([r["realized_daughter_difference"] for r in rows])),
            "opposite_enrichment_p": float(np.mean([r["opposite_enrichment"] for r in rows])),
            "n_polymorphic_split_events": int(len(rows)),
        }
    else:
        aggregate = {
            "mean_realized_daughter_difference": float("nan"),
            "opposite_enrichment_p": float("nan"),
            "n_polymorphic_split_events": 0,
        }
    return rows, aggregate


def _prop_fields(prefix: str, successes: int, total: int) -> dict[str, Any]:
    lo, hi = wilson_interval(int(successes), int(total))
    return {
        f"{prefix}_n": int(successes),
        f"{prefix}_p": float(successes) / float(total) if total else float("nan"),
        f"{prefix}_ci_low": lo,
        f"{prefix}_ci_high": hi,
    }


def _summarize_eta(
    *,
    eta: float,
    n_histories: int,
    pattern_counts: Counter,
    persistence_n: int,
    strict_n: int,
    discordant_n: int,
    persistent_discordant_n: int,
    d_without_p_count: int,
    diagnostic: Mapping[str, float],
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "eta": float(eta),
        "odds_ratio_exp_eta": math.exp(float(eta)),
        "n_histories": int(n_histories),
    }
    row.update(_prop_fields("persistence", persistence_n, n_histories))
    row.update(_prop_fields("strict_2_2", strict_n, n_histories))
    d_fields = _prop_fields("discordant_2_2", discordant_n, n_histories)
    row.update(d_fields)
    lo, hi = wilson_interval(persistent_discordant_n, persistence_n)
    row.update({
        "persistent_discordant_n": int(persistent_discordant_n),
        "p_discordant_given_persistent": persistent_discordant_n / persistence_n if persistence_n else float("nan"),
        "conditional_ci_low": lo,
        "conditional_ci_high": hi,
    })
    lo, hi = wilson_interval(discordant_n, strict_n)
    row.update({
        "p_discordant_given_2_2": discordant_n / strict_n if strict_n else float("nan"),
        "discordant_given_2_2_ci_low": lo,
        "discordant_given_2_2_ci_high": hi,
        "d_without_p_count": int(d_without_p_count),
        "mean_realized_daughter_difference": diagnostic["mean_realized_daughter_difference"],
        "opposite_enrichment_p": diagnostic["opposite_enrichment_p"],
        "n_polymorphic_split_events": diagnostic["n_polymorphic_split_events"],
    })
    for pattern in ALL_PATTERNS:
        row[f"pattern_{pattern}_n"] = int(pattern_counts.get(pattern, 0))
        row[f"pattern_{pattern}_p"] = int(pattern_counts.get(pattern, 0)) / n_histories if n_histories else float("nan")
    return row


def run_forward_sweep(
    config: Mapping[str, Any],
    output: Path,
    *,
    etas: Iterable[float],
    histories_per_eta: int,
    process: str,
    split_mode: str,
    master_ss: np.random.SeedSequence,
    retain_genealogy_histories: bool = False,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[float, list[tuple[Any, dict[str, int], str]]]]:
    tree, rearrangement = _tree_and_rearrangement(config)
    base_config = _base_sim_config(config, process=process)
    base_config["_experiment_config"] = config
    summary_rows = []
    replicate_rows = []
    diagnostic_rows = []
    retained: dict[float, list[tuple[Any, dict[str, int], str]]] = defaultdict(list)
    cell_streams = master_ss.spawn(len(list(etas)))
    for eta, cell_ss in zip(list(etas), cell_streams):
        rng = np.random.default_rng(cell_ss)
        pattern_counts: Counter = Counter()
        persistence_n = strict_n = discordant_n = persistent_discordant_n = d_without_p_count = 0
        split_diag_values = []
        cell_seed = int(cell_ss.generate_state(1)[0])
        for replicate_id in range(int(histories_per_eta)):
            history, sampled, pattern, persistent, strict, discordant = _simulate_one(
                tree, rearrangement, base_config, float(eta), rng, split_mode=split_mode
            )
            diagnostic_splits = {"ROOT"} if split_mode == "root_only" else {"A", "B"}
            event_rows, aggregate_diag = _event_diagnostics(history, diagnostic_splits)
            for row in event_rows:
                row.update({"eta": float(eta), "replicate_id": replicate_id, "process": process, "split_mode": split_mode})
                diagnostic_rows.append(row)
            if aggregate_diag["n_polymorphic_split_events"]:
                split_diag_values.append(aggregate_diag)
            pattern_counts[pattern] += 1
            persistence_n += int(persistent)
            strict_n += int(strict)
            discordant_n += int(discordant)
            persistent_discordant_n += int(persistent and discordant)
            d_without_p_count += int(discordant and not persistent)
            replicate_rows.append({
                "process": process,
                "split_mode": split_mode,
                "eta": float(eta),
                "eta_cell_seed": cell_seed,
                "replicate_id": replicate_id,
                "persistent": persistent,
                "strict_2_2": strict,
                "discordant_2_2": discordant,
                "terminal_pattern": pattern,
                "sample_1": sampled["1"],
                "sample_2": sampled["2"],
                "sample_3": sampled["3"],
                "sample_4": sampled["4"],
                "A_end_status": branch_end_status(history, "A"),
                "B_end_status": branch_end_status(history, "B"),
            })
            if retain_genealogy_histories and persistent and discordant:
                limit = int(config["experiment"]["genealogy_histories_per_eta"])
                if len(retained[float(eta)]) < limit:
                    retained[float(eta)].append((history, sampled, pattern))
        if split_diag_values:
            diag = {
                "mean_realized_daughter_difference": float(np.mean([x["mean_realized_daughter_difference"] for x in split_diag_values])),
                "opposite_enrichment_p": float(np.mean([x["opposite_enrichment_p"] for x in split_diag_values])),
                "n_polymorphic_split_events": int(sum(x["n_polymorphic_split_events"] for x in split_diag_values)),
            }
        else:
            diag = {"mean_realized_daughter_difference": float("nan"), "opposite_enrichment_p": float("nan"), "n_polymorphic_split_events": 0}
        row = _summarize_eta(
            eta=float(eta),
            n_histories=int(histories_per_eta),
            pattern_counts=pattern_counts,
            persistence_n=persistence_n,
            strict_n=strict_n,
            discordant_n=discordant_n,
            persistent_discordant_n=persistent_discordant_n,
            d_without_p_count=d_without_p_count,
            diagnostic=diag,
        )
        row.update({"process": process, "split_mode": split_mode, "eta_cell_seed": cell_seed})
        summary_rows.append(row)
    return summary_rows, replicate_rows, diagnostic_rows, retained


def run_genealogy_check(
    config: Mapping[str, Any],
    retained: Mapping[float, list[tuple[Any, dict[str, int], str]]],
    master_ss: np.random.SeedSequence,
) -> list[dict[str, Any]]:
    tree, _rearrangement = _tree_and_rearrangement(config)
    rec = config["baseline"]["recombination"]
    base_rate = float(rec["baseline_rate"])
    fraction = float(rec["effective_cross_arrangement_fraction"])
    rows = []
    eta_streams = master_ss.spawn(len(retained))
    for (eta, histories), eta_ss in zip(sorted(retained.items()), eta_streams):
        history_streams = eta_ss.spawn(len(histories))
        for history_id, ((history, sampled, pattern), hist_ss) in enumerate(zip(histories, history_streams)):
            rng = np.random.default_rng(hist_ss)
            counts = [0, 0, 0]
            loci = int(config["experiment"]["genealogy_loci_per_history"])
            for locus_id in range(loci):
                result = simulate_genealogy(locus_id, tree, history, sampled, base_rate, fraction, rng, record_events=False)
                counts[int(result.topology_index)] += 1
            q = [count / loci for count in counts]
            arrangement_topology = ARRANGEMENT_TOPOLOGY_BY_PATTERN[pattern]
            rows.append({
                "eta": float(eta),
                "history_id": history_id,
                "terminal_pattern": pattern,
                "num_loci": loci,
                "n1": counts[0],
                "n2": counts[1],
                "n3": counts[2],
                "q1": q[0],
                "q2": q[1],
                "q3": q[2],
                "species_tree_topology": SPECIES_TREE_TOPOLOGY,
                "arrangement_associated_topology": arrangement_topology,
                "dominant_topology": int(np.argmax(q)),
                "arrangement_associated_dominant": int(np.argmax(q)) == arrangement_topology,
                "arrangement_minus_species_margin": q[arrangement_topology] - q[SPECIES_TREE_TOPOLOGY],
            })
    return rows


def theoretical_targets(etas: Iterable[float]) -> list[dict[str, Any]]:
    rows = []
    for eta in etas:
        for p in (0.25, 0.5, 0.75):
            rows.append({
                "eta": float(eta),
                "parent_frequency_A1": p,
                "target_frequency_positive": daughter_target_frequency(p, float(eta)),
                "target_frequency_negative": daughter_target_frequency(p, -float(eta)),
            })
    return rows


def _plot_main(summary: list[Mapping[str, Any]], out: Path) -> None:
    rows = [r for r in summary if r["process"] == "wright_fisher" and r["split_mode"] == "descendant"]
    etas = np.array([float(r["eta"]) for r in rows])
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.2), constrained_layout=True)
    panels = [
        ("persistence_p", "persistence_ci_low", "persistence_ci_high", "P(persistence)"),
        ("discordant_2_2_p", "discordant_2_2_ci_low", "discordant_2_2_ci_high", "P(discordant strict 2:2)"),
        ("p_discordant_given_persistent", "conditional_ci_low", "conditional_ci_high", "P(discordant | persistence)"),
        ("mean_realized_daughter_difference", None, None, "Mean daughter-frequency divergence"),
    ]
    for ax, (value, low, high, ylabel) in zip(axes.flat, panels):
        y = np.array([float(r[value]) for r in rows])
        if low:
            yerr = np.vstack([y - np.array([float(r[low]) for r in rows]), np.array([float(r[high]) for r in rows]) - y])
            ax.errorbar(etas, y, yerr=yerr, marker="o", color="#1f77b4", linewidth=1.6, capsize=3)
        else:
            ax.plot(etas, y, marker="o", color="#1f77b4", linewidth=1.6)
        ax.set_xlabel("eta")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25, linewidth=0.6)
    fig.suptitle("Speciation daughter sorting at descendant splits")
    fig.savefig(out / "daughter_sorting_eta_main.pdf", bbox_inches="tight")
    fig.savefig(out / "daughter_sorting_eta_main.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_patterns(pattern_rows: list[Mapping[str, Any]], out: Path) -> None:
    rows = [r for r in pattern_rows if r["process"] == "wright_fisher" and r["split_mode"] in {"descendant", "sign_symmetry"}]
    by_eta = {float(r["eta"]): r for r in rows}
    etas = sorted(by_eta)
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    for pattern, color in (("1010", "#1f77b4"), ("0101", "#d62728"), ("discordant_total", "#2ca02c")):
        if pattern == "discordant_total":
            y = [float(by_eta[e]["discordant_2_2_p"]) for e in etas]
            label = "discordant strict 2:2"
        else:
            y = [float(by_eta[e][f"pattern_{pattern}_p"]) for e in etas]
            label = pattern
        ax.plot(etas, y, marker="o", linewidth=1.8, label=label, color=color)
    ax.set_xlabel("eta")
    ax.set_ylabel("Terminal pattern frequency")
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.legend(frameon=False)
    ax.set_title("Labeled terminal patterns under daughter sorting")
    fig.savefig(out / "daughter_sorting_eta_patterns.pdf", bbox_inches="tight")
    fig.savefig(out / "daughter_sorting_eta_patterns.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_genealogy(genealogy_summary: list[Mapping[str, Any]], out: Path) -> None:
    if not genealogy_summary:
        return
    rows = sorted(genealogy_summary, key=lambda r: float(r["eta"]))
    etas = [float(r["eta"]) for r in rows]
    margins = [float(r["mean_arrangement_minus_species_margin"]) for r in rows]
    dominant = [float(r["arrangement_associated_dominant_fraction"]) for r in rows]
    fig, ax1 = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    ax1.plot(etas, margins, marker="o", color="#1f77b4", label="arrangement - species margin")
    ax1.set_xlabel("eta")
    ax1.set_ylabel("Mean quartet margin")
    ax1.grid(alpha=0.25, linewidth=0.6)
    ax2 = ax1.twinx()
    ax2.plot(etas, dominant, marker="s", color="#d62728", label="arrangement dominant fraction")
    ax2.set_ylabel("Dominant fraction")
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines], frameon=False, loc="best")
    ax1.set_title("Genealogy check for persistent discordant histories")
    fig.savefig(out / "daughter_sorting_eta_genealogy.pdf", bbox_inches="tight")
    fig.savefig(out / "daughter_sorting_eta_genealogy.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _genealogy_summary(rows: list[Mapping[str, Any]]) -> list[dict[str, Any]]:
    out = []
    by_eta: dict[float, list[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        by_eta[float(row["eta"])].append(row)
    for eta, chunk in sorted(by_eta.items()):
        out.append({
            "eta": eta,
            "num_histories": len(chunk),
            "arrangement_associated_dominant_n": sum(bool(r["arrangement_associated_dominant"]) for r in chunk),
            "arrangement_associated_dominant_fraction": sum(bool(r["arrangement_associated_dominant"]) for r in chunk) / len(chunk) if chunk else float("nan"),
            "mean_arrangement_minus_species_margin": float(np.mean([float(r["arrangement_minus_species_margin"]) for r in chunk])) if chunk else float("nan"),
            "mean_q_species": float(np.mean([float(r["q1"]) for r in chunk])) if chunk else float("nan"),
        })
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the focused daughter-sorting eta experiment.")
    parser.add_argument("--config", default="experiments/daughter_sorting_eta/config.yaml")
    parser.add_argument("--output", default="results/daughter_sorting_eta")
    args = parser.parse_args()

    start = perf_counter()
    config_path = Path(args.config)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    config = yaml.safe_load(config_path.read_text())
    master_seed = int(config["experiment"]["master_seed"])
    master = np.random.SeedSequence(master_seed)
    streams = master.spawn(6)

    main_etas = [float(x) for x in config["eta_design"]["main_eta_values"]]
    sign_etas = [-2.0]
    root_etas = [float(x) for x in config["eta_design"]["root_control_eta_values"]]
    moran_etas = [float(x) for x in config["eta_design"]["moran_robustness_eta_values"]]

    summary, replicates, split_diag, retained = run_forward_sweep(
        config,
        output,
        etas=main_etas,
        histories_per_eta=int(config["experiment"]["main_histories_per_eta"]),
        process="wright_fisher",
        split_mode="descendant",
        master_ss=streams[0],
        retain_genealogy_histories=True,
    )
    sign_summary, sign_reps, sign_diag, _sign_retained = run_forward_sweep(
        config,
        output,
        etas=sign_etas,
        histories_per_eta=int(config["experiment"]["sign_symmetry_histories_per_eta"]),
        process="wright_fisher",
        split_mode="descendant",
        master_ss=streams[1],
        retain_genealogy_histories=False,
    )
    for row in sign_summary:
        row["split_mode"] = "sign_symmetry"
    for row in sign_reps:
        row["split_mode"] = "sign_symmetry"
    for row in sign_diag:
        row["split_mode"] = "sign_symmetry"

    if any(int(row["d_without_p_count"]) != 0 for row in summary + sign_summary):
        bad = [row for row in summary + sign_summary if int(row["d_without_p_count"]) != 0]
        raise RuntimeError(f"Discordant strict 2:2 without persistence observed: {bad}")

    root_summary, root_reps, root_diag, _ = run_forward_sweep(
        config,
        output,
        etas=root_etas,
        histories_per_eta=int(config["experiment"]["root_control_histories_per_eta"]),
        process="wright_fisher",
        split_mode="root_only",
        master_ss=streams[2],
        retain_genealogy_histories=False,
    )
    moran_summary, moran_reps, moran_diag, _ = run_forward_sweep(
        config,
        output,
        etas=moran_etas,
        histories_per_eta=int(config["experiment"]["moran_robustness_histories_per_eta"]),
        process="moran",
        split_mode="descendant",
        master_ss=streams[3],
        retain_genealogy_histories=False,
    )

    genealogy_rows = run_genealogy_check(config, retained, streams[4])
    genealogy_summary = _genealogy_summary(genealogy_rows)

    all_summary = summary + sign_summary + root_summary + moran_summary
    all_replicates = replicates + sign_reps + root_reps + moran_reps
    all_diag = split_diag + sign_diag + root_diag + moran_diag

    terminal_rows = []
    for row in all_summary:
        for pattern in ALL_PATTERNS:
            terminal_rows.append({
                "process": row["process"],
                "split_mode": row["split_mode"],
                "eta": row["eta"],
                "pattern": pattern,
                "count": row[f"pattern_{pattern}_n"],
                "frequency": row[f"pattern_{pattern}_p"],
            })
    sign_rows = [row for row in all_summary if row["process"] == "wright_fisher" and row["split_mode"] in {"descendant", "sign_symmetry"} and float(row["eta"]) in {-2.0, 2.0}]
    target_rows = theoretical_targets(sorted(set(main_etas + [-2.0])))

    frozen_config = output / "frozen_config.yaml"
    frozen_config.write_text(yaml.safe_dump(config, sort_keys=False))
    _write_csv(output / "eta_summary.csv", all_summary)
    _write_csv(output / "eta_replicates.csv", all_replicates)
    _write_csv(output / "eta_terminal_patterns.csv", terminal_rows)
    _write_csv(output / "eta_split_diagnostics.csv", all_diag)
    _write_csv(output / "eta_sign_symmetry.csv", sign_rows)
    _write_csv(output / "eta_theoretical_targets.csv", target_rows)
    _write_csv(output / "eta_genealogy_check.csv", genealogy_rows)
    _write_csv(output / "eta_genealogy_summary.csv", genealogy_summary)
    _write_csv(output / "root_sorting_control.csv", root_summary)
    _write_csv(output / "moran_robustness.csv", moran_summary)

    _plot_main(summary, output)
    _plot_patterns(summary + sign_summary, output)
    _plot_genealogy(genealogy_summary, output)

    runtime = perf_counter() - start
    files_for_checksums = [
        "frozen_config.yaml",
        "eta_summary.csv",
        "eta_replicates.csv",
        "eta_terminal_patterns.csv",
        "eta_split_diagnostics.csv",
        "eta_sign_symmetry.csv",
        "eta_theoretical_targets.csv",
        "eta_genealogy_check.csv",
        "eta_genealogy_summary.csv",
        "root_sorting_control.csv",
        "moran_robustness.csv",
        "daughter_sorting_eta_main.pdf",
        "daughter_sorting_eta_main.png",
        "daughter_sorting_eta_patterns.pdf",
        "daughter_sorting_eta_patterns.png",
        "daughter_sorting_eta_genealogy.pdf",
        "daughter_sorting_eta_genealogy.png",
    ]
    provenance = {
        "command": f"python scripts/run_daughter_sorting_eta.py --config {config_path} --output {output}",
        "git_commit": _git_commit(),
        "dirty_working_tree_status": _git_status(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "matplotlib_version": plt.matplotlib.__version__,
        "master_seed": master_seed,
        "runtime_seconds": runtime,
        "moran_time_convention": TIME_CONVENTION,
        "config": config,
        "checksums_sha256": {name: _sha256(output / name) for name in files_for_checksums if (output / name).exists()},
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True))

    readme = output / "README.md"
    readme.write_text(
        "# Daughter-Sorting Eta Experiment\n\n"
        "Run command:\n\n"
        f"```bash\npython scripts/run_daughter_sorting_eta.py --config {config_path} --output {output}\n```\n\n"
        "Primary outputs are `eta_summary.csv`, `eta_terminal_patterns.csv`, "
        "`eta_split_diagnostics.csv`, `eta_genealogy_summary.csv`, and the "
        "`daughter_sorting_eta_*.pdf/png` figures. The primary analysis applies "
        "eta only at descendant splits A and B with positive children 1 and 3.\n"
    )
    print(f"Wrote {output}")
    print(f"Runtime seconds: {runtime:.2f}")


if __name__ == "__main__":
    main()
