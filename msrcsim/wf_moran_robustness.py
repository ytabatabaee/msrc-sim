from __future__ import annotations

from dataclasses import dataclass
from math import exp, sqrt
from typing import Any, Mapping

import numpy as np

from .conditioning import branch_end_status, terminal_pattern
from .history_io import save_frozen_history
from .model_fitting import compare_models, msc_probabilities
from .population_process import simulate_population_history
from .rearrangement import Rearrangement
from .robustness import dominant_quartet_threshold, infer_with_strategy
from .species_tree import SpeciesTree
from .structured_coalescent import simulate_genealogy
from .wright_fisher import FrequencyHistory, FrequencyRecord


DISCORDANT_T2_PATTERNS = {"0101", "1010"}
DISCORDANT_T3_PATTERNS = {"0110", "1001"}
DISCORDANT_PATTERNS = DISCORDANT_T2_PATTERNS | DISCORDANT_T3_PATTERNS


@dataclass(frozen=True)
class QuartetSummary:
    n1: int
    n2: int
    n3: int
    q1: float
    q2: float
    q3: float
    favored_topology: int
    favored_margin: float
    unique_favored_dominant: bool
    beta: float
    se_favored_margin: float
    ci95_low: float
    ci95_high: float
    distance_to_nearest_msc_arm: float


class CachedFrequencyHistory:
    """Read-only frequency-history adapter for Monte Carlo replay experiments."""

    def __init__(self, history: FrequencyHistory):
        self.records = history.records
        self.by_branch = {
            branch_id: sorted(records, key=lambda r: r.absolute_age)
            for branch_id, records in history.by_branch.items()
        }
        self._ages = {
            branch_id: np.asarray([r.absolute_age for r in records], dtype=float)
            for branch_id, records in self.by_branch.items()
        }

    def frequency_at(self, branch_id: str, age: float) -> float:
        records = self.by_branch[branch_id]
        ages = self._ages[branch_id]
        idx = int(np.searchsorted(ages, float(age) + 1e-12, side="right") - 1)
        idx = min(max(idx, 0), len(records) - 1)
        return float(records[idx].frequency_A1)

    def next_frequency_boundary(self, branch_id: str, age: float) -> float:
        ages = self._ages[branch_id]
        idx = int(np.searchsorted(ages, float(age) + 1e-10, side="right"))
        if idx >= len(ages):
            return float("inf")
        return float(ages[idx])

    def terminal_frequency(self, taxon: str) -> float:
        return self.frequency_at(taxon, 0.0)


def four_taxon_tree(ne: int, tau: float, tip_time: float = 1.0, root_extension: float = 4.0) -> SpeciesTree:
    """Return the fixed quartet tree 12|34 with lengths in coalescent units."""
    scale = 2.0 * int(ne)
    tip = tip_time * scale
    internal = float(tau) * scale
    root_ext = root_extension * scale
    return SpeciesTree(f"((1:{tip:g},2:{tip:g})A:{internal:g},(3:{tip:g},4:{tip:g})B:{internal:g})ROOT;", ne, root_ext)


def root_rearrangement(tree: SpeciesTree, ne: int) -> Rearrangement:
    root = tree.branches["ROOT"]
    return Rearrangement("inv", "inversion", "ROOT", 1, int(ne), 0.0)


def favored_topology_for_pattern(pattern: str) -> int:
    pattern = str(pattern)
    if pattern in DISCORDANT_T2_PATTERNS:
        return 1
    if pattern in DISCORDANT_T3_PATTERNS:
        return 2
    if pattern in {"0011", "1100"}:
        return 0
    raise ValueError(f"Pattern {pattern!r} is not a 2:2 quartet pattern")


def is_two_two(pattern: str) -> bool:
    return str(pattern).count("1") == 2


def is_discordant_two_two(pattern: str) -> bool:
    return str(pattern) in DISCORDANT_PATTERNS


def is_persistent(history: FrequencyHistory) -> bool:
    return branch_end_status(history, "A") == "segregating" and branch_end_status(history, "B") == "segregating"


def _record(rearrangement_id: str, branch, age: float, frequency: float, index: int, status: str) -> FrequencyRecord:
    total = 2 * branch.effective_population_size
    k = int(round(float(frequency) * total))
    return FrequencyRecord(
        rearrangement_id,
        branch.branch_id,
        branch.parent_branch_id,
        index,
        float(age),
        k,
        total - k,
        total,
        k / total,
        1.0 - k / total,
        status,
        False,
        index == 0,
        index == 1,
        branch.selection_coefficient,
        branch.effective_population_size,
    )


def matched_discordant_history(tree: SpeciesTree, pattern: str, rearrangement_id: str = "inv") -> tuple[FrequencyHistory, dict[str, int]]:
    """Construct a deterministic frozen history for a chosen 2:2 sampled pattern."""
    if favored_topology_for_pattern(pattern) not in (1, 2):
        raise ValueError("matched history requires a discordant 2:2 pattern")
    sampled = {taxon: int(bit) for taxon, bit in zip(tree.taxa, str(pattern))}
    records: list[FrequencyRecord] = []
    by: dict[str, list[FrequencyRecord]] = {}
    for branch_id, branch in tree.branches.items():
        if branch_id == "ROOT":
            freq = 0.5
            status = "segregating"
        elif branch_id in {"A", "B"}:
            freq = 0.5
            status = "segregating"
        else:
            freq = float(sampled[branch_id])
            status = "fixed" if freq == 1.0 else "lost"
        rows = [
            _record(rearrangement_id, branch, branch.older_age, freq, 0, status),
            _record(rearrangement_id, branch, branch.younger_age, freq, 1, status),
        ]
        by[branch_id] = rows
        records.extend(rows)
    return FrequencyHistory(records, by), sampled


def simulate_quartet_summary(
    tree: SpeciesTree,
    history: FrequencyHistory,
    sampled: Mapping[str, int],
    num_loci: int,
    seed: int,
    recombination_rate: float,
    effective_fraction: float,
    favored_topology: int | None = None,
) -> QuartetSummary:
    rng = np.random.default_rng(int(seed))
    replay_history = CachedFrequencyHistory(history)
    results = [
        simulate_genealogy(i, tree, replay_history, sampled, recombination_rate, effective_fraction, rng, record_events=False)
        for i in range(int(num_loci))
    ]
    counts = np.bincount([r.topology_index for r in results], minlength=3)
    q = counts / int(num_loci)
    if favored_topology is None:
        favored_topology = favored_topology_for_pattern(terminal_pattern(sampled, tree.taxa))
    favored = float(q[favored_topology])
    species = float(q[0])
    other = max(float(q[i]) for i in range(3) if i != favored_topology)
    margin = favored - species
    # Conservative multinomial delta-method SE for q_favored - q_species.
    se = sqrt(max(0.0, favored + species - margin * margin) / int(num_loci))
    model = compare_models(counts)
    return QuartetSummary(
        n1=int(counts[0]),
        n2=int(counts[1]),
        n3=int(counts[2]),
        q1=float(q[0]),
        q2=float(q[1]),
        q3=float(q[2]),
        favored_topology=int(favored_topology),
        favored_margin=float(margin),
        unique_favored_dominant=bool(favored > other),
        beta=float(max(0.0, margin)),
        se_favored_margin=float(se),
        ci95_low=float(margin - 1.959963984540054 * se),
        ci95_high=float(margin + 1.959963984540054 * se),
        distance_to_nearest_msc_arm=float(model["distance_to_nearest_msc_arm"]),
    )


def forward_history_config(process: str, ne: int, tau: float, root_extension: float = 0.5) -> dict[str, Any]:
    tree = four_taxon_tree(ne, tau, root_extension=root_extension)
    return {
        "seed": 1,
        "population_process": {"model": process},
        "species_tree": {
            "newick": f"((1:{tree.branches['1'].older_age:g},2:{tree.branches['2'].older_age:g})A:{tree.branches['A'].older_age - tree.branches['A'].younger_age:g},(3:{tree.branches['3'].older_age:g},4:{tree.branches['4'].older_age:g})B:{tree.branches['B'].older_age - tree.branches['B'].younger_age:g})ROOT;",
            "default_effective_population_size": int(ne),
            "root_extension": float(tree.root_extension),
        },
        "rearrangement": {
            "id": "inv",
            "type": "inversion",
            "origin_branch": "ROOT",
            "origin_time_from_branch_start": 1,
            "initial_copy_count": int(ne),
            "selection": {"model": "genic", "coefficient": 0.0},
        },
    }


def simulate_forward_attempt(process: str, ne: int, tau: float, rng: np.random.Generator, root_extension: float = 0.5):
    tree = four_taxon_tree(ne, tau, root_extension=root_extension)
    history = simulate_population_history(tree, root_rearrangement(tree, ne), rng, {"population_process": {"model": process}})
    sampled = {t: int(rng.random() < history.terminal_frequency(t)) for t in tree.taxa}
    pattern = terminal_pattern(sampled, tree.taxa)
    return tree, history, sampled, pattern


def msc_delta(tau: float) -> float:
    return 1.0 - exp(-float(tau))


def threshold_from_beta(tau: float, beta: float) -> float:
    return dominant_quartet_threshold(float(tau), float(beta))


def mixture_probabilities(tau: float, msrc_q: tuple[float, float, float], epsilon: float) -> np.ndarray:
    return (1.0 - float(epsilon)) * msc_probabilities(0, float(tau)) + float(epsilon) * np.asarray(msrc_q, dtype=float)


def make_weighted_rows(q_msc: np.ndarray, q_msrc: np.ndarray, epsilon: float, windows: int, rng: np.random.Generator) -> list[dict[str, Any]]:
    rows = []
    for window_id in range(int(windows)):
        is_rearranged = bool(rng.random() < float(epsilon))
        q = q_msrc if is_rearranged else q_msc
        top = int(rng.choice(3, p=q))
        rows.append({
            "window_id": window_id,
            "block_id": window_id,
            "topology_index": top,
            "is_rearranged": is_rearranged,
            "rearrangement_id": "inv" if is_rearranged else "",
            "msrc_probability": 1.0 if is_rearranged else 0.0,
            "weight": 0.0 if is_rearranged else 1.0,
        })
    return rows


def correction_recovered(rows: list[dict[str, Any]], strategy: str, soft_ratio: float | None = None) -> bool:
    if strategy == "soft_weight":
        ratio = 0.5 if soft_ratio is None else float(soft_ratio)

        def soft(row: Mapping[str, Any]) -> float:
            return ratio if bool(row.get("is_rearranged", False)) else 1.0

        result = infer_with_strategy(rows, strategy, soft)
    else:
        result = infer_with_strategy(rows, strategy)
    return result.inferred_topology_index == 0


def save_matched_history(path, config: Mapping[str, Any], history: FrequencyHistory, sampled: Mapping[str, int], process_label: str, pattern: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    save_frozen_history(path, config, history, sampled, {
        "population_process": {"model": process_label},
        "terminal_pattern": pattern,
        "matched_history_control": True,
        "note": "The same deterministic frequency history is reused with different generator labels.",
    })
