from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any, Mapping

import numpy as np

from .species_tree import SpeciesTree


MODEL_RECIPROCAL_LOG_ODDS = "reciprocal_log_odds"
MODEL_MASS_CONSERVING = "mass_conserving"
ORIENTATION_RULE = "positive_child_override_else_newick_order"


@dataclass(frozen=True)
class SplitRule:
    eta: float
    positive_child: str | None
    positive_child_source: str


@dataclass(frozen=True)
class SpeciationSplitPolicy:
    model: str
    eta: float
    by_split: Mapping[str, Mapping[str, Any]]
    rules: Mapping[str, SplitRule]

    def rule_for(self, split_id: str) -> SplitRule:
        return self.rules[split_id]

    def provenance(self) -> dict[str, Any]:
        return {
            "enabled": True,
            "model": self.model,
            "eta": self.eta,
            "orientation_rule": ORIENTATION_RULE,
            "by_split": dict(self.by_split),
        }


@dataclass(frozen=True)
class SplitInitializationEvent:
    parent_branch_id: str | None
    child_branch_id: str
    child_side_or_sign: str
    positive_child: str | None
    parent_state: str
    parent_A1_count: int | None
    parent_total: int | None
    parent_frequency_A1: float | None
    eta: float
    log_odds_offset: float
    omega: float | None
    target_frequency_A1: float | None
    child_total: int
    realized_child_A1_count: int | None
    model: str


def _children_by_split(tree: SpeciesTree) -> dict[str, list[str]]:
    return {node.name: [child.name for child in node.children] for node in tree.nodes() if node.children}


def _finite_eta(value: Any, *, path: str) -> float:
    eta = float(value)
    if not math.isfinite(eta):
        raise ValueError(f"{path} must be finite")
    return eta


def resolve_speciation_split(config: Mapping[str, Any] | None, tree: SpeciesTree | None = None) -> SpeciationSplitPolicy | None:
    section = (config or {}).get("speciation_split")
    if section is None:
        return None
    if not isinstance(section, Mapping):
        raise ValueError("speciation_split must be a mapping")

    model = str(section.get("model", MODEL_RECIPROCAL_LOG_ODDS))
    if model not in {MODEL_RECIPROCAL_LOG_ODDS, MODEL_MASS_CONSERVING}:
        raise ValueError("speciation_split.model must be 'reciprocal_log_odds' or 'mass_conserving'")
    eta = _finite_eta(section.get("eta", 0.0), path="speciation_split.eta")
    by_split = section.get("by_split", {})
    if by_split is None:
        by_split = {}
    if not isinstance(by_split, Mapping):
        raise ValueError("speciation_split.by_split must be a mapping")

    rules: dict[str, SplitRule] = {}
    if tree is not None:
        children = _children_by_split(tree)
        for split_id, override in by_split.items():
            if split_id not in children:
                raise ValueError(f"speciation_split.by_split references unknown split {split_id!r}")
            if not isinstance(override, Mapping):
                raise ValueError(f"speciation_split.by_split.{split_id} must be a mapping")
            if "positive_child" in override and str(override["positive_child"]) not in children[split_id]:
                raise ValueError(
                    f"speciation_split.by_split.{split_id}.positive_child must be a direct child of {split_id!r}"
                )

        for split_id, child_ids in children.items():
            override = by_split.get(split_id, {}) or {}
            split_eta = _finite_eta(override.get("eta", eta), path=f"speciation_split.by_split.{split_id}.eta")
            if len(child_ids) != 2 and split_eta != 0.0:
                raise ValueError(
                    f"speciation_split {model} with nonzero eta requires a binary split; "
                    f"split {split_id!r} has {len(child_ids)} children"
                )
            if "positive_child" in override:
                positive_child = str(override["positive_child"])
                source = "override"
            elif len(child_ids) == 2:
                positive_child = child_ids[0]
                source = "newick_order_fallback"
            else:
                positive_child = None
                source = "not_applicable_nonbinary_eta_zero"
            rules[split_id] = SplitRule(split_eta, positive_child, source)

    return SpeciationSplitPolicy(model, eta, by_split, rules)


def daughter_target_frequency(parent_frequency: float, log_odds_offset: float) -> float:
    p = float(parent_frequency)
    offset = float(log_odds_offset)
    if not math.isfinite(offset):
        raise ValueError("log_odds_offset must be finite")
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    if offset == 0.0:
        return p
    logit = math.log(p) - math.log1p(-p) + offset
    if logit >= 0.0:
        z = math.exp(-logit)
        return 1.0 / (1.0 + z)
    z = math.exp(logit)
    return z / (1.0 + z)


def mass_conserving_daughter_frequencies(parent_frequency: float, eta: float) -> tuple[float, float]:
    p = float(parent_frequency)
    eta = float(eta)
    if not math.isfinite(eta):
        raise ValueError("eta must be finite")
    if p <= 0.0:
        return 0.0, 0.0
    if p >= 1.0:
        return 1.0, 1.0
    if eta == 0.0:
        return p, p
    delta = math.tanh(eta) * min(p, 1.0 - p)
    return p + delta, p - delta


def _omega(offset: float) -> float | None:
    try:
        return math.exp(offset)
    except OverflowError:
        return None


def _event(
    *,
    policy: SpeciationSplitPolicy,
    parent_branch_id: str | None,
    child_branch_id: str,
    child_side_or_sign: str,
    positive_child: str | None,
    parent_state: str,
    parent_count: int | None,
    parent_total: int | None,
    parent_frequency: float | None,
    eta: float,
    offset: float,
    target_frequency: float | None,
    child_total: int,
    realized_child_count: int | None,
) -> SplitInitializationEvent:
    return SplitInitializationEvent(
        parent_branch_id=parent_branch_id,
        child_branch_id=child_branch_id,
        child_side_or_sign=child_side_or_sign,
        positive_child=positive_child,
        parent_state=parent_state,
        parent_A1_count=parent_count,
        parent_total=parent_total,
        parent_frequency_A1=parent_frequency,
        eta=eta,
        log_odds_offset=offset,
        omega=_omega(offset) if policy.model == MODEL_RECIPROCAL_LOG_ODDS else None,
        target_frequency_A1=target_frequency,
        child_total=child_total,
        realized_child_A1_count=realized_child_count,
        model=policy.model,
    )


def sample_daughter_initialization(
    *,
    tree: SpeciesTree,
    parent_branch_id: str | None,
    child_branch_id: str,
    parent_count: int | None,
    parent_total: int | None,
    child_total: int,
    rng: np.random.Generator,
    policy: SpeciationSplitPolicy,
) -> tuple[int | None, SplitInitializationEvent]:
    if parent_branch_id is None:
        raise ValueError("Speciation split initialization requires a parent branch")
    rule = policy.rule_for(parent_branch_id)
    child_ids = _children_by_split(tree)[parent_branch_id]
    if child_branch_id not in child_ids:
        raise ValueError(f"{child_branch_id!r} is not a direct child of split {parent_branch_id!r}")
    if len(child_ids) != 2 and rule.eta != 0.0:
        raise ValueError(
            f"speciation_split {policy.model} with nonzero eta requires a binary split; "
            f"split {parent_branch_id!r} has {len(child_ids)} children"
        )

    if len(child_ids) == 2:
        offset = rule.eta if child_branch_id == rule.positive_child else -rule.eta
        side = "positive" if child_branch_id == rule.positive_child else "negative"
    else:
        offset = 0.0
        side = "neutral"

    if parent_count is None or parent_total is None:
        return None, _event(
            policy=policy,
            parent_branch_id=parent_branch_id,
            child_branch_id=child_branch_id,
            child_side_or_sign=side,
            positive_child=rule.positive_child,
            parent_state="not_present",
            parent_count=None,
            parent_total=None,
            parent_frequency=None,
            eta=rule.eta,
            offset=offset,
            target_frequency=None,
            child_total=child_total,
            realized_child_count=None,
        )

    parent_frequency = parent_count / parent_total
    if policy.model == MODEL_MASS_CONSERVING and len(child_ids) == 2:
        positive_target, negative_target = mass_conserving_daughter_frequencies(parent_frequency, rule.eta)
        target = positive_target if side == "positive" else negative_target
    else:
        target = daughter_target_frequency(parent_frequency, offset)
    realized = int(rng.binomial(child_total, target))
    if parent_count == 0:
        state = "lost"
    elif parent_count == parent_total:
        state = "fixed"
    else:
        state = "segregating"
    return realized, _event(
        policy=policy,
        parent_branch_id=parent_branch_id,
        child_branch_id=child_branch_id,
        child_side_or_sign=side,
        positive_child=rule.positive_child,
        parent_state=state,
        parent_count=int(parent_count),
        parent_total=int(parent_total),
        parent_frequency=parent_frequency,
        eta=rule.eta,
        offset=offset,
        target_frequency=target,
        child_total=child_total,
        realized_child_count=realized,
    )


def split_events_as_dicts(events: list[SplitInitializationEvent] | tuple[SplitInitializationEvent, ...]) -> list[dict[str, Any]]:
    return [dict(event) if isinstance(event, Mapping) else asdict(event) for event in events]
