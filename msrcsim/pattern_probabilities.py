from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from math import lgamma
from pathlib import Path
from time import perf_counter
from typing import Any, Mapping
import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from .analytic import compute_Hm
from .conditioning import terminal_pattern
from .experiments import _rearrangement_from_config, _tree_from_config
from .rearrangement import Rearrangement
from .species_tree import Node, SpeciesTree
from .wright_fisher import _selected_p, simulate_frequency_history

DEFAULT_MAX_STATES = 2001
PATTERNS = tuple("".join(map(str, z)) for z in product((0, 1), repeat=4))


@dataclass(frozen=True)
class PatternProbabilityResult:
    pattern_probabilities: dict[str, float]
    joint_persistent_pattern_probabilities: dict[str, float]
    conditional_persistent_pattern_probabilities: dict[str, float | None]
    summary: dict[str, Any]
    runtime_seconds: float
    max_state_count: int


def _assert_reasonable_state_count(Ne: int, max_states: int) -> int:
    total = 2 * int(Ne)
    states = total + 1
    if states > max_states:
        raise ValueError(
            f"Exact Wright-Fisher DP needs {states} count states for Ne={Ne}; "
            f"limit is {max_states}. Increase --max-states deliberately or use "
            "a smaller Ne for exact validation."
        )
    return total


def wf_transition_matrix(
    Ne: int,
    selection: float = 0.0,
    *,
    max_states: int = DEFAULT_MAX_STATES,
) -> np.ndarray:
    """Return the exact finite-population Wright-Fisher transition matrix.

    Rows and columns are chromosome-copy counts ``0..2Ne``. With non-zero
    selection this uses the same genic selected parental probability as the
    simulator.
    """
    total = _assert_reasonable_state_count(Ne, max_states)
    j = np.arange(total + 1, dtype=float)
    log_choose = (
        lgamma(total + 1)
        - np.array([lgamma(x + 1) for x in range(total + 1)])
        - np.array([lgamma(total - x + 1) for x in range(total + 1)])
    )
    matrix = np.zeros((total + 1, total + 1), dtype=float)
    matrix[0, 0] = 1.0
    matrix[total, total] = 1.0
    for i in range(1, total):
        p = _selected_p(i / total, float(selection))
        if p <= 0.0:
            matrix[i, 0] = 1.0
        elif p >= 1.0:
            matrix[i, total] = 1.0
        else:
            row = np.exp(log_choose + j * np.log(p) + (total - j) * np.log1p(-p))
            matrix[i, :] = row / row.sum()
    return matrix


def wf_distribution_after_t(
    Ne: int,
    initial_count: int,
    generations: int,
    selection: float = 0.0,
    *,
    max_states: int = DEFAULT_MAX_STATES,
) -> np.ndarray:
    """Propagate a point-mass count distribution for ``generations`` WF steps."""
    total = _assert_reasonable_state_count(Ne, max_states)
    initial_count = int(initial_count)
    if not 0 <= initial_count <= total:
        raise ValueError(f"initial_count must be between 0 and {total}")
    if int(generations) < 0:
        raise ValueError("generations must be non-negative")
    dist = np.zeros(total + 1, dtype=float)
    dist[initial_count] = 1.0
    if int(generations) == 0:
        return dist
    transition = wf_transition_matrix(Ne, selection, max_states=max_states)
    for _ in range(int(generations)):
        dist = dist @ transition
    return _renormalize(dist, "WF distribution")


def persistent_polymorphism_probability(dist: np.ndarray) -> float:
    """Probability that a count distribution is neither lost nor fixed."""
    if dist.size < 2:
        raise ValueError("count distribution must contain absorbing endpoints")
    return float(1.0 - dist[0] - dist[-1])


def equal_frequency_pattern_class_probabilities(p: float) -> dict[str, float]:
    """Closed-form pattern classes when all four terminal frequencies are p."""
    p = float(p)
    return {
        "P_2_2": 6.0 * p * p * (1.0 - p) * (1.0 - p),
        "P_3_1": 4.0 * p**3 * (1.0 - p) + 4.0 * p * (1.0 - p) ** 3,
        "P_4_0": p**4 + (1.0 - p) ** 4,
    }


def pattern_class_summary(pattern_probabilities: Mapping[str, float]) -> dict[str, float]:
    p22 = sum(prob for pattern, prob in pattern_probabilities.items() if pattern.count("1") == 2)
    p31 = sum(prob for pattern, prob in pattern_probabilities.items() if pattern.count("1") in {1, 3})
    p40 = float(pattern_probabilities["0000"] + pattern_probabilities["1111"])
    w1 = float(pattern_probabilities["0011"] + pattern_probabilities["1100"])
    w2 = float(pattern_probabilities["0101"] + pattern_probabilities["1010"])
    w3 = float(pattern_probabilities["0110"] + pattern_probabilities["1001"])
    return {
        "P_2_2": float(p22),
        "P_3_1": float(p31),
        "P_4_0": p40,
        "w1": w1,
        "w2": w2,
        "w3": w3,
    }


def conditional_pattern_probabilities(
    joint_pattern_probabilities: Mapping[str, float],
    persistent_probability: float,
) -> dict[str, float | None]:
    """Return P(pattern | persistent), or null values when persistence is impossible."""
    if persistent_probability <= 0.0:
        return {pattern: None for pattern in PATTERNS}
    return {
        pattern: float(joint_pattern_probabilities[pattern]) / float(persistent_probability)
        for pattern in PATTERNS
    }


def conditional_pattern_class_summary(
    conditional_probabilities: Mapping[str, float | None],
) -> dict[str, float | None]:
    if any(conditional_probabilities[pattern] is None for pattern in PATTERNS):
        return {
            "P_2_2_given_persistent": None,
            "P_3_1_given_persistent": None,
            "P_4_0_given_persistent": None,
            "w1_given_persistent": None,
            "w2_given_persistent": None,
            "w3_given_persistent": None,
            "discordant_2_2_given_persistent": None,
        }
    classes = pattern_class_summary({k: float(v) for k, v in conditional_probabilities.items()})
    return {
        "P_2_2_given_persistent": classes["P_2_2"],
        "P_3_1_given_persistent": classes["P_3_1"],
        "P_4_0_given_persistent": classes["P_4_0"],
        "w1_given_persistent": classes["w1"],
        "w2_given_persistent": classes["w2"],
        "w3_given_persistent": classes["w3"],
        "discordant_2_2_given_persistent": classes["w2"] + classes["w3"],
    }


def joint_persistent_pattern_class_summary(
    joint_pattern_probabilities: Mapping[str, float],
) -> dict[str, float]:
    classes = pattern_class_summary(joint_pattern_probabilities)
    return {
        "P_2_2_and_persistent": classes["P_2_2"],
        "P_3_1_and_persistent": classes["P_3_1"],
        "P_4_0_and_persistent": classes["P_4_0"],
        "w1_and_persistent": classes["w1"],
        "w2_and_persistent": classes["w2"],
        "w3_and_persistent": classes["w3"],
        "discordant_2_2_and_persistent": classes["w2"] + classes["w3"],
    }


def quartet_weights_to_zero_switching_q(
    weights: tuple[float, float, float],
    t: float,
    m01: float,
    m10: float,
    lambda0: float,
    lambda1: float,
) -> tuple[float, float, float]:
    """Map arrangement partition weights through the existing zero-switching kernel."""
    H, states = compute_Hm(t, m01, m10, lambda0, lambda1)
    lookup = {state: i for i, state in enumerate(states)}
    reps = ((0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 0))
    q = np.zeros(3, dtype=float)
    for weight, state in zip(weights, reps):
        q += float(weight) * H[lookup[state]]
    total = float(sum(weights))
    if total <= 0.0:
        return (float("nan"), float("nan"), float("nan"))
    q /= total
    return (float(q[0]), float(q[1]), float(q[2]))


def exact_pattern_probabilities(
    tree: SpeciesTree,
    rearrangement: Rearrangement,
    *,
    max_states: int = DEFAULT_MAX_STATES,
) -> PatternProbabilityResult:
    """Compute exact sampled-tip pattern probabilities for the four-taxon tree.

    Patterns are ordered by ``tree.taxa``. Quartet partition weights use that
    order: w1 = 0011 + 1100 (12|34), w2 = 0101 + 1010 (13|24), and
    w3 = 0110 + 1001 (14|23).

    ``persistent_at_all_required_speciation_events`` includes every internal
    species-tree node at or below the rearrangement origin age. For the standard
    quartet with origin on the root branch, these nodes are ROOT, A, and B.
    """
    start = perf_counter()
    origin_branch = tree.branches[rearrangement.origin_branch]
    origin_age = origin_branch.older_age - rearrangement.origin_time_from_branch_start
    if not (origin_branch.younger_age <= origin_age <= origin_branch.older_age):
        raise ValueError("Origin time falls outside origin branch")
    for branch in tree.branches.values():
        _assert_reasonable_state_count(branch.effective_population_size, max_states)

    children_by_name = {node.name: node.children for node in tree.nodes()}
    nodes_by_name = {node.name: node for node in tree.nodes()}
    transition_cache: dict[tuple[int, int, float], np.ndarray] = {}
    inherit_cache: dict[tuple[int, int, int], np.ndarray] = {}

    def transition(Ne: int, generations: int) -> np.ndarray:
        key = (int(Ne), int(generations), float(rearrangement.selection_coefficient))
        if key not in transition_cache:
            if generations == 0:
                transition_cache[key] = np.eye(2 * int(Ne) + 1)
            else:
                matrix = wf_transition_matrix(
                    Ne,
                    rearrangement.selection_coefficient,
                    max_states=max_states,
                )
                transition_cache[key] = np.linalg.matrix_power(matrix, generations)
        return transition_cache[key]

    def inherit(parent_total: int, child_total: int, parent_count: int) -> np.ndarray:
        key = (int(parent_total), int(child_total), int(parent_count))
        if key not in inherit_cache:
            p = parent_count / parent_total
            # Speciation inheritance is a fresh binomial draw from the parent
            # frequency, matching simulate_frequency_history exactly.
            j = np.arange(child_total + 1, dtype=float)
            if p <= 0.0:
                row = np.zeros(child_total + 1); row[0] = 1.0
            elif p >= 1.0:
                row = np.zeros(child_total + 1); row[child_total] = 1.0
            else:
                log_choose = (
                    lgamma(child_total + 1)
                    - np.array([lgamma(x + 1) for x in range(child_total + 1)])
                    - np.array([lgamma(child_total - x + 1) for x in range(child_total + 1)])
                )
                row = np.exp(log_choose + j * np.log(p) + (child_total - j) * np.log1p(-p))
                row = row / row.sum()
            inherit_cache[key] = row
        return inherit_cache[key]

    def branch_generations(branch_id: str, start_age: float) -> int:
        younger = tree.branches[branch_id].younger_age
        generations = round(start_age - younger)
        if abs((start_age - younger) - generations) > 1e-8:
            raise ValueError("Exact DP requires integer branch durations in generations")
        return int(generations)

    def combine(left: dict[tuple[str, bool], float], right: dict[tuple[str, bool], float]):
        out: dict[tuple[str, bool], float] = defaultdict(float)
        for (lp, lok), lprob in left.items():
            for (rp, rok), rprob in right.items():
                out[(lp + rp, lok and rok)] += lprob * rprob
        return out

    def subtree_from_branch_end(branch_id: str, end_count: int | None):
        node = nodes_by_name[tree.branches[branch_id].child_node]
        required_here = (not node.is_tip()) and (node.age <= origin_age + 1e-12)
        if end_count is None:
            current_ok = not required_here
        else:
            total = 2 * tree.branches[branch_id].effective_population_size
            current_ok = (not required_here) or (0 < end_count < total)

        if node.is_tip():
            p = 0.0 if end_count is None else end_count / (2 * tree.branches[branch_id].effective_population_size)
            return {
                ("0", current_ok): 1.0 - p,
                ("1", current_ok): p,
            }

        child_results = []
        for child in children_by_name[node.name]:
            child_branch = tree.branches[child.name]
            child_total = 2 * child_branch.effective_population_size
            child_start: dict[int | None, float] = defaultdict(float)
            if end_count is None:
                child_start[None] = 1.0
            else:
                parent_total = 2 * tree.branches[branch_id].effective_population_size
                row = inherit(parent_total, child_total, end_count)
                for count, prob in enumerate(row):
                    if prob:
                        child_start[count] += float(prob)
            child_results.append(propagate_branch(child.name, child_start))
        joined = combine(child_results[0], child_results[1])
        return {(pattern, current_ok and ok): prob for (pattern, ok), prob in joined.items()}

    def propagate_branch(branch_id: str, start_dist: Mapping[int | None, float]):
        branch = tree.branches[branch_id]
        total = 2 * branch.effective_population_size
        if branch_id == rearrangement.origin_branch:
            start_count = int(rearrangement.initial_copy_count)
            if not 0 <= start_count <= total:
                raise ValueError(
                    "initial_copy_count must be between 0 and "
                    f"{total} for origin branch {branch_id}"
                )
            dist = np.zeros(total + 1)
            dist[start_count] = 1.0
            gens = branch_generations(branch_id, origin_age)
        else:
            if set(start_dist) == {None}:
                return subtree_from_branch_end(branch_id, None)
            dist = np.zeros(total + 1)
            for count, prob in start_dist.items():
                if count is not None:
                    dist[int(count)] += float(prob)
            gens = branch_generations(branch_id, branch.older_age)
        end_dist = dist @ transition(branch.effective_population_size, gens)
        end_dist = _renormalize(end_dist, f"branch {branch_id} end distribution")
        out: dict[tuple[str, bool], float] = defaultdict(float)
        for count, prob in enumerate(end_dist):
            if prob:
                for key, subprob in subtree_from_branch_end(branch_id, count).items():
                    out[key] += float(prob) * subprob
        return out

    root_start = {None: 1.0}
    joint = propagate_branch(tree.root.name, root_start)
    pattern_probs = {pattern: 0.0 for pattern in PATTERNS}
    joint_persistent_probs = {pattern: 0.0 for pattern in PATTERNS}
    persistent = 0.0
    for (pattern, ok), prob in joint.items():
        pattern_probs[pattern] += prob
        if ok:
            joint_persistent_probs[pattern] += prob
            persistent += prob
    pattern_probs = _normalize_pattern_probabilities(pattern_probs)
    conditional_probs = conditional_pattern_probabilities(joint_persistent_probs, persistent)
    classes = pattern_class_summary(pattern_probs)
    joint_classes = joint_persistent_pattern_class_summary(joint_persistent_probs)
    conditional_classes = conditional_pattern_class_summary(conditional_probs)
    runtime = perf_counter() - start
    max_state_count = max(2 * b.effective_population_size + 1 for b in tree.branches.values())
    summary = {
        "version": "0.8.7",
        "statistic_definition": "persistent_at_all_required_speciation_events",
        "taxon_order": list(tree.taxa),
        "quartet_partition_weight_definitions": {
            "w1": "0011 + 1100 -> 12|34 relative to taxon_order",
            "w2": "0101 + 1010 -> 13|24 relative to taxon_order",
            "w3": "0110 + 1001 -> 14|23 relative to taxon_order",
        },
        "required_speciation_nodes": _required_speciation_nodes(tree, origin_age),
        "persistent_at_all_required_speciation_events": float(persistent),
        **classes,
        "discordant_2_2": classes["w2"] + classes["w3"],
        **joint_classes,
        **conditional_classes,
        "conditional_persistent_diagnostic": (
            None
            if persistent > 0.0
            else "P(persistent_at_all_required_speciation_events)=0; conditional values are null"
        ),
        "probability_sum": float(sum(pattern_probs.values())),
        "joint_persistent_probability_sum": float(sum(joint_persistent_probs.values())),
        "max_state_count": int(max_state_count),
        "runtime_seconds": float(runtime),
    }
    if abs(summary["probability_sum"] - 1.0) > 1e-8:
        raise RuntimeError(f"Pattern probabilities sum to {summary['probability_sum']}")
    if abs(summary["joint_persistent_probability_sum"] - persistent) > 1e-8:
        raise RuntimeError("Joint persistent pattern probabilities do not sum to P(persistent)")
    return PatternProbabilityResult(
        pattern_probs,
        joint_persistent_probs,
        conditional_probs,
        summary,
        runtime,
        max_state_count,
    )


def monte_carlo_pattern_probabilities(
    tree: SpeciesTree,
    rearrangement: Rearrangement,
    *,
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    counts = {pattern: 0 for pattern in PATTERNS}
    joint_persistent_counts = {pattern: 0 for pattern in PATTERNS}
    persistent = 0
    required_nodes = _required_speciation_nodes(
        tree,
        tree.branches[rearrangement.origin_branch].older_age
        - rearrangement.origin_time_from_branch_start,
    )
    required_branches = required_nodes
    for _ in range(int(replicates)):
        history = simulate_frequency_history(tree, rearrangement, rng)
        sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        pattern = terminal_pattern(sampled, tree.taxa)
        counts[pattern] += 1
        is_persistent = all(_branch_end_is_polymorphic(history, branch_id) for branch_id in required_branches)
        if is_persistent:
            joint_persistent_counts[pattern] += 1
            persistent += 1
    pattern_probs = {pattern: counts[pattern] / float(replicates) for pattern in PATTERNS}
    joint_persistent_probs = {
        pattern: joint_persistent_counts[pattern] / float(replicates)
        for pattern in PATTERNS
    }
    persistent_prob = persistent / float(replicates)
    conditional_probs = conditional_pattern_probabilities(joint_persistent_probs, persistent_prob)
    classes = pattern_class_summary(pattern_probs)
    joint_classes = joint_persistent_pattern_class_summary(joint_persistent_probs)
    conditional_classes = conditional_pattern_class_summary(conditional_probs)
    return {
        "pattern_probabilities": pattern_probs,
        "joint_persistent_pattern_probabilities": joint_persistent_probs,
        "conditional_persistent_pattern_probabilities": conditional_probs,
        "summary": {
            "taxon_order": list(tree.taxa),
            "persistent_at_all_required_speciation_events": persistent_prob,
            **classes,
            "discordant_2_2": classes["w2"] + classes["w3"],
            **joint_classes,
            **conditional_classes,
            "conditional_persistent_diagnostic": (
                None
                if persistent_prob > 0.0
                else "P(persistent_at_all_required_speciation_events)=0; conditional values are null"
            ),
        },
        "replicates": int(replicates),
        "seed": int(seed),
    }


def run_pattern_probability_analysis(
    config: Mapping[str, Any],
    output: str | Path,
    *,
    monte_carlo_replicates: int = 20000,
    seed: int | None = None,
    max_states: int = DEFAULT_MAX_STATES,
) -> dict[str, Any]:
    tree = _tree_from_config(config)
    rearrangement = _rearrangement_from_config(config)
    theory = exact_pattern_probabilities(tree, rearrangement, max_states=max_states)
    mc_seed = int(config.get("seed", 1) if seed is None else seed)
    simulation = monte_carlo_pattern_probabilities(
        tree,
        rearrangement,
        replicates=int(monte_carlo_replicates),
        seed=mc_seed,
    )

    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    _write_pattern_csv(
        out / "theoretical_pattern_probabilities.csv",
        theory.pattern_probabilities,
        theory.joint_persistent_pattern_probabilities,
        theory.conditional_persistent_pattern_probabilities,
    )
    with (out / "theoretical_summary.json").open("w") as handle:
        json.dump(theory.summary, handle, indent=2)
    comparison = _comparison_rows(theory, simulation)
    with (out / "theory_vs_simulation.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["statistic", "theory", "simulation", "abs_error"])
        writer.writeheader(); writer.writerows(comparison)
    _plot_theory_vs_simulation(out / "theory_vs_simulation_pattern_probabilities.pdf", comparison)
    return {
        "output": str(out),
        "theory": theory,
        "simulation": simulation,
        "comparison": comparison,
        "max_abs_error": max(
            float(row["abs_error"])
            for row in comparison
            if row["abs_error"] is not None
        ),
    }


def _comparison_rows(theory: PatternProbabilityResult, simulation: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = []
    sim_summary = simulation["summary"]
    mapping = {
        "persistent_polymorphism": "persistent_at_all_required_speciation_events",
        "P_2_2": "P_2_2",
        "P_3_1": "P_3_1",
        "P_4_0": "P_4_0",
        "w1": "w1",
        "w2": "w2",
        "w3": "w3",
        "P_2_2_and_persistent": "P_2_2_and_persistent",
        "P_3_1_and_persistent": "P_3_1_and_persistent",
        "P_4_0_and_persistent": "P_4_0_and_persistent",
        "discordant_2_2_and_persistent": "discordant_2_2_and_persistent",
        "P_2_2_given_persistent": "P_2_2_given_persistent",
        "P_3_1_given_persistent": "P_3_1_given_persistent",
        "P_4_0_given_persistent": "P_4_0_given_persistent",
        "w1_given_persistent": "w1_given_persistent",
        "w2_given_persistent": "w2_given_persistent",
        "w3_given_persistent": "w3_given_persistent",
        "discordant_2_2_given_persistent": "discordant_2_2_given_persistent",
    }
    for label, key in mapping.items():
        t = theory.summary[key]
        s = sim_summary[key]
        if t is None or s is None:
            rows.append({"statistic": label, "theory": t, "simulation": s, "abs_error": None})
        else:
            tf = float(t)
            sf = float(s)
            rows.append({"statistic": label, "theory": tf, "simulation": sf, "abs_error": abs(tf - sf)})
    return rows


def _write_pattern_csv(
    path: Path,
    pattern_probabilities: Mapping[str, float],
    joint_persistent_pattern_probabilities: Mapping[str, float],
    conditional_persistent_pattern_probabilities: Mapping[str, float | None],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "pattern",
                "P_pattern",
                "P_pattern_and_persistent",
                "P_pattern_given_persistent",
            ]
        )
        for pattern in PATTERNS:
            writer.writerow(
                [
                    pattern,
                    float(pattern_probabilities[pattern]),
                    float(joint_persistent_pattern_probabilities[pattern]),
                    conditional_persistent_pattern_probabilities[pattern],
                ]
            )


def _plot_theory_vs_simulation(path: Path, rows: list[Mapping[str, Any]]) -> None:
    finite_rows = [
        row for row in rows
        if row["theory"] is not None and row["simulation"] is not None
    ]
    if not finite_rows:
        return
    x = [float(row["theory"]) for row in finite_rows]
    y = [float(row["simulation"]) for row in finite_rows]
    labels = [str(row["statistic"]) for row in finite_rows]
    fig, ax = plt.subplots(figsize=(5.5, 5.0))
    ax.scatter(x, y, color="#2364aa")
    lo = min(x + y + [0.0])
    hi = max(x + y + [1.0])
    ax.plot([lo, hi], [lo, hi], color="#111111", linewidth=1.0)
    for xx, yy, label in zip(x, y, labels):
        ax.annotate(label, (xx, yy), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Theoretical probability")
    ax.set_ylabel("Simulated probability")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _required_speciation_nodes(tree: SpeciesTree, origin_age: float) -> list[str]:
    out = []
    for node in tree.nodes():
        if not node.is_tip() and node.age <= origin_age + 1e-12:
            out.append(node.name)
    return sorted(out, key=lambda name: tree.branches[name].older_age, reverse=True)


def _branch_end_is_polymorphic(history: Any, branch_id: str) -> bool:
    records = history.by_branch[branch_id]
    end = min(records, key=lambda row: row.absolute_age)
    return 0 < end.copy_count_A1 < end.population_chromosomes


def _renormalize(dist: np.ndarray, label: str) -> np.ndarray:
    total = float(dist.sum())
    if not np.isfinite(total) or total <= 0.0:
        raise RuntimeError(f"{label} has invalid probability mass {total}")
    if abs(total - 1.0) > 1e-10:
        dist = dist / total
    return dist


def _normalize_pattern_probabilities(pattern_probs: Mapping[str, float]) -> dict[str, float]:
    total = float(sum(pattern_probs.values()))
    if not np.isfinite(total) or total <= 0.0:
        raise RuntimeError(f"Pattern probabilities have invalid total {total}")
    return {pattern: float(prob) / total for pattern, prob in pattern_probs.items()}
