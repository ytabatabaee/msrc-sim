from __future__ import annotations

from collections import Counter
from math import sqrt
from typing import Iterable, Mapping, Any


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return (float("nan"), float("nan"))
    p = successes / total
    denom = 1.0 + z * z / total
    center = (p + z * z / (2.0 * total)) / denom
    half = z * sqrt((p * (1.0 - p) + z * z / (4.0 * total)) / total) / denom
    return max(0.0, center - half), min(1.0, center + half)


def summarize_replicates(rows: Iterable[Mapping[str, Any]], asymmetry_threshold: float) -> dict[str, Any]:
    rows = list(rows)
    accepted = [r for r in rows if bool(r.get("accepted", True))]
    attempted = len(rows)
    n = len(accepted)

    patterns = Counter(str(r.get("terminal_pattern", "")) for r in accepted)
    persistent = sum(bool(r.get("persisted_to_target", False)) for r in accepted)
    two_two = sum(bool(r.get("is_2_2_pattern", False)) for r in accepted)
    asymmetric = sum(abs(float(r.get("q2_minus_q3", 0.0))) > asymmetry_threshold for r in accepted)
    discordant_dominant = sum(int(r.get("dominant_topology", -1)) in (1, 2) for r in accepted)

    def prop_with_ci(k: int) -> dict[str, float]:
        lo, hi = wilson_interval(k, n)
        return {"count": k, "proportion": (k / n if n else float("nan")), "ci95_low": lo, "ci95_high": hi}

    return {
        "attempted_replicates": attempted,
        "accepted_replicates": n,
        "acceptance_rate": n / attempted if attempted else float("nan"),
        "asymmetry_threshold": asymmetry_threshold,
        "persistent_polymorphism": prop_with_ci(persistent),
        "two_two_terminal_pattern": prop_with_ci(two_two),
        "asymmetric_quartet_distribution": prop_with_ci(asymmetric),
        "discordant_topology_dominant": prop_with_ci(discordant_dominant),
        "terminal_pattern_counts": dict(sorted(patterns.items())),
    }
