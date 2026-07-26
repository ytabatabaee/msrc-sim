from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Mapping, Iterable
import copy
import csv
import itertools
import json

import numpy as np
import yaml

from .conditioning import evaluate_conditioning, terminal_pattern, branch_end_status
from .rearrangement import Rearrangement
from .simplex import barycentric_to_cartesian
from .species_tree import SpeciesTree
from .statistics import summarize_replicates
from .structured_coalescent import simulate_genealogy
from .wright_fisher import simulate_frequency_history


@dataclass(frozen=True)
class ReplicateRecord:
    replicate_id: int
    attempt_id: int
    seed: int
    accepted: bool
    conditioning_reason: str
    origin_branch: str
    origin_time_from_branch_start: int
    initial_copy_count: int
    selection_coefficient: float
    terminal_pattern: str
    root_end_status: str
    persisted_to_target: bool
    is_2_2_pattern: bool
    num_loci: int
    q1: float
    q2: float
    q3: float
    q2_minus_q3: float
    absolute_asymmetry: float
    dominant_topology: int
    simplex_x: float
    simplex_y: float
    mean_switches_per_locus: float
    mean_coalescence_time: float


def _tree_from_config(config: Mapping[str, Any]) -> SpeciesTree:
    st = config["species_tree"]
    return SpeciesTree(
        st["newick"],
        st["default_effective_population_size"],
        st["root_extension"],
        st.get("branch_parameters"),
    )


def _rearrangement_from_config(config: Mapping[str, Any]) -> Rearrangement:
    rr = config["rearrangement"]
    sel = rr.get("selection", {}).get("coefficient", rr.get("selection_coefficient", 0.0))
    return Rearrangement(
        rr.get("id", "inv_1"),
        rr["type"],
        rr["origin_branch"],
        int(rr["origin_time_from_branch_start"]),
        int(rr.get("initial_copy_count", 1)),
        float(sel),
    )


def _recombination(config: Mapping[str, Any]) -> tuple[float, float]:
    rec = config["recombination"]
    base = rec.get("baseline_rate", rec.get("rate"))
    frac = rec.get("effective_cross_arrangement_fraction", rec.get("suppression_factor"))
    if base is None or frac is None:
        raise ValueError("Recombination requires baseline_rate and effective_cross_arrangement_fraction")
    return float(base), float(frac)


def _is_two_two(pattern: str) -> bool:
    return pattern.count("1") == 2


def simulate_replicates(config: Mapping[str, Any]) -> tuple[list[ReplicateRecord], dict[str, Any]]:
    """Simulate independent inversion histories and conditional locus genealogies."""
    tree = _tree_from_config(config)
    rearrangement = _rearrangement_from_config(config)
    base_rate, fraction = _recombination(config)

    exp = config.get("experiment", {})
    requested = int(exp.get("replicates", config.get("replicates", 100)))
    loci_per_replicate = int(exp.get("loci_per_replicate", config.get("num_loci", 1000)))
    threshold = float(exp.get("asymmetry_threshold", 0.05))
    conditioning = config.get("conditioning", {"mode": "none"})
    max_attempts = int(conditioning.get("max_attempts", max(requested, requested * 100)))
    target_branches = exp.get("persistence_target_branches", conditioning.get("require_segregating_at", []))
    if isinstance(target_branches, str):
        target_branches = [target_branches]

    master = np.random.SeedSequence(int(config.get("seed", 1)))
    attempt_sequences = master.spawn(max_attempts)
    records: list[ReplicateRecord] = []
    accepted_count = 0

    for attempt_id, ss in enumerate(attempt_sequences):
        rng = np.random.default_rng(ss)
        history = simulate_frequency_history(tree, rearrangement, rng)
        sampled = {t: int(rng.random() < history.terminal_frequency(t)) for t in tree.taxa}
        cond = evaluate_conditioning(conditioning, history, sampled, tree.taxa)

        if not cond.accepted:
            records.append(
                ReplicateRecord(
                    replicate_id=-1,
                    attempt_id=attempt_id,
                    seed=int(ss.generate_state(1)[0]),
                    accepted=False,
                    conditioning_reason=cond.reason,
                    origin_branch=rearrangement.origin_branch,
                    origin_time_from_branch_start=rearrangement.origin_time_from_branch_start,
                    initial_copy_count=rearrangement.initial_copy_count,
                    selection_coefficient=rearrangement.selection_coefficient,
                    terminal_pattern=terminal_pattern(sampled, tree.taxa),
                    root_end_status=branch_end_status(history, tree.root.name),
                    persisted_to_target=False,
                    is_2_2_pattern=_is_two_two(terminal_pattern(sampled, tree.taxa)),
                    num_loci=0,
                    q1=float("nan"), q2=float("nan"), q3=float("nan"),
                    q2_minus_q3=float("nan"), absolute_asymmetry=float("nan"),
                    dominant_topology=-1, simplex_x=float("nan"), simplex_y=float("nan"),
                    mean_switches_per_locus=float("nan"), mean_coalescence_time=float("nan"),
                )
            )
            continue

        results = [
            simulate_genealogy(
                i, tree, history, sampled, base_rate, fraction, rng, record_events=True
            )
            for i in range(loci_per_replicate)
        ]
        counts = np.bincount([r.topology_index for r in results], minlength=3)
        q = counts / loci_per_replicate
        x, y = barycentric_to_cartesian(float(q[0]), float(q[1]), float(q[2]))
        switches = [sum(e.event_type == "switch" for e in r.events) for r in results]
        coal_times = [t for r in results for t in r.coalescence_times]
        persisted = bool(target_branches) and all(branch_end_status(history, b) == "segregating" for b in target_branches)
        pattern = terminal_pattern(sampled, tree.taxa)
        rep_id = accepted_count
        accepted_count += 1
        records.append(
            ReplicateRecord(
                replicate_id=rep_id,
                attempt_id=attempt_id,
                seed=int(ss.generate_state(1)[0]),
                accepted=True,
                conditioning_reason=cond.reason,
                origin_branch=rearrangement.origin_branch,
                origin_time_from_branch_start=rearrangement.origin_time_from_branch_start,
                initial_copy_count=rearrangement.initial_copy_count,
                selection_coefficient=rearrangement.selection_coefficient,
                terminal_pattern=pattern,
                root_end_status=branch_end_status(history, tree.root.name),
                persisted_to_target=persisted,
                is_2_2_pattern=_is_two_two(pattern),
                num_loci=loci_per_replicate,
                q1=float(q[0]), q2=float(q[1]), q3=float(q[2]),
                q2_minus_q3=float(q[1] - q[2]),
                absolute_asymmetry=float(abs(q[1] - q[2])),
                dominant_topology=int(np.argmax(q)),
                simplex_x=x, simplex_y=y,
                mean_switches_per_locus=float(np.mean(switches)),
                mean_coalescence_time=float(np.mean(coal_times)),
            )
        )
        if accepted_count >= requested:
            break

    if accepted_count < requested:
        raise RuntimeError(
            f"Conditioning accepted only {accepted_count}/{requested} replicates "
            f"after {max_attempts} attempts"
        )

    summary = summarize_replicates([asdict(r) for r in records], threshold)
    summary["requested_replicates"] = requested
    summary["loci_per_replicate"] = loci_per_replicate
    return records, summary


def write_replicate_outputs(config: Mapping[str, Any], records: Iterable[ReplicateRecord], summary: Mapping[str, Any]) -> Path:
    out = Path(config.get("output", {}).get("directory", "replicate_output"))
    out.mkdir(parents=True, exist_ok=True)
    rows = [asdict(r) for r in records]
    with (out / "replicate_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader(); writer.writerows(rows)
    accepted = [r for r in rows if r["accepted"]]
    with (out / "simplex_points.csv").open("w", newline="") as handle:
        fields = ["replicate_id", "q1", "q2", "q3", "simplex_x", "simplex_y", "terminal_pattern", "q2_minus_q3", "dominant_topology"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows({k: r[k] for k in fields} for r in accepted)
    with (out / "prevalence_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    with (out / "config.resolved.yaml").open("w") as handle:
        yaml.safe_dump(dict(config), handle, sort_keys=False)
    return out


def _set_dotted(config: dict[str, Any], path: str, value: Any) -> None:
    parts = path.split(".")
    target = config
    for part in parts[:-1]:
        target = target.setdefault(part, {})
    target[parts[-1]] = value


def parameter_combinations(grid: Mapping[str, list[Any]]) -> Iterable[dict[str, Any]]:
    keys = list(grid)
    for values in itertools.product(*(grid[k] for k in keys)):
        yield dict(zip(keys, values))


def simulate_parameter_grid(config: Mapping[str, Any]) -> tuple[list[dict[str, Any]], Path]:
    grid = config.get("parameter_grid")
    if not grid:
        raise ValueError("parameter_grid configuration is required")
    base = copy.deepcopy(dict(config))
    base.pop("parameter_grid", None)
    output_root = Path(config.get("output", {}).get("directory", "grid_output"))
    output_root.mkdir(parents=True, exist_ok=True)
    aggregate: list[dict[str, Any]] = []

    for cell_id, params in enumerate(parameter_combinations(grid)):
        cell = copy.deepcopy(base)
        for path, value in params.items():
            _set_dotted(cell, path, value)
        cell_dir = output_root / f"cell_{cell_id:04d}"
        cell.setdefault("output", {})["directory"] = str(cell_dir)
        records, summary = simulate_replicates(cell)
        write_replicate_outputs(cell, records, summary)
        row: dict[str, Any] = {"cell_id": cell_id, **params}
        row.update({
            "attempted_replicates": summary["attempted_replicates"],
            "accepted_replicates": summary["accepted_replicates"],
            "acceptance_rate": summary["acceptance_rate"],
            "persistent_probability": summary["persistent_polymorphism"]["proportion"],
            "two_two_probability": summary["two_two_terminal_pattern"]["proportion"],
            "asymmetry_probability": summary["asymmetric_quartet_distribution"]["proportion"],
            "discordant_dominant_probability": summary["discordant_topology_dominant"]["proportion"],
        })
        aggregate.append(row)

    with (output_root / "parameter_grid_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(aggregate[0].keys()))
        writer.writeheader(); writer.writerows(aggregate)
    return aggregate, output_root
