from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from .wright_fisher import FrequencyHistory


@dataclass(frozen=True)
class ConditioningResult:
    accepted: bool
    reason: str


def terminal_pattern(sampled_arrangements: Mapping[str, int], taxa: tuple[str, ...]) -> str:
    return "".join(str(int(sampled_arrangements[t])) for t in taxa)


def branch_end_status(history: FrequencyHistory, branch_id: str) -> str:
    if branch_id not in history.by_branch:
        raise ValueError(f"Unknown branch in conditioning rule: {branch_id}")
    return history.by_branch[branch_id][-1].status


def evaluate_conditioning(
    spec: Mapping[str, Any] | None,
    history: FrequencyHistory,
    sampled_arrangements: Mapping[str, int],
    taxa: tuple[str, ...],
) -> ConditioningResult:
    """Evaluate an evolutionary replicate against a transparent conditioning rule.

    Supported modes:
      - none
      - persistent_polymorphism
      - terminal_pattern
      - persistent_and_pattern
    """
    spec = spec or {"mode": "none"}
    mode = str(spec.get("mode", "none"))

    if mode == "none":
        return ConditioningResult(True, "unconditional")

    pattern = terminal_pattern(sampled_arrangements, taxa)

    persistent_ok = True
    if mode in {"persistent_polymorphism", "persistent_and_pattern"}:
        branches = spec.get("require_segregating_at", [])
        if isinstance(branches, str):
            branches = [branches]
        if not branches:
            raise ValueError(
                "persistent_polymorphism conditioning requires require_segregating_at"
            )
        persistent_ok = all(branch_end_status(history, b) == "segregating" for b in branches)

    pattern_ok = True
    if mode in {"terminal_pattern", "persistent_and_pattern"}:
        accepted_patterns = [str(x) for x in spec.get("accepted_patterns", [])]
        if not accepted_patterns:
            raise ValueError("terminal_pattern conditioning requires accepted_patterns")
        pattern_ok = pattern in accepted_patterns

    accepted = persistent_ok and pattern_ok
    if accepted:
        return ConditioningResult(True, "accepted")
    reasons = []
    if not persistent_ok:
        reasons.append("not_persistent")
    if not pattern_ok:
        reasons.append(f"terminal_pattern={pattern}")
    return ConditioningResult(False, ";".join(reasons))
