from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from typing import Any, Iterable, Mapping

import numpy as np

from .analytic import TOPOLOGY_NAMES
from .model_fitting import compare_models


@dataclass(frozen=True)
class RearrangementInterval:
    start: int
    end: int
    interval_id: str = "rearrangement"


def generate_locus_positions(length: int, count: int, placement: str, rng: np.random.Generator) -> np.ndarray:
    if length <= 0 or count <= 0:
        raise ValueError("genome length and locus count must be positive")
    if count > length:
        raise ValueError("locus count cannot exceed integer genome length")
    if placement == "evenly_spaced":
        pos = np.linspace(0, length - 1, count, dtype=int)
    elif placement == "uniform_random":
        pos = np.sort(rng.choice(length, size=count, replace=False).astype(np.int64))
    else:
        raise ValueError("placement must be 'evenly_spaced' or 'uniform_random'")
    if len(np.unique(pos)) != count:
        raise ValueError("locus positions must be sorted and unique")
    return pos


def classify_region(position: int, interval: RearrangementInterval) -> tuple[str, bool]:
    if position < interval.start:
        return "outside_left", False
    if position <= interval.end:
        return "rearrangement", True
    return "outside_right", False


def topology_counts(topologies: Iterable[int]) -> np.ndarray:
    return np.bincount(np.asarray(tuple(topologies), dtype=int), minlength=3)[:3]


def _q(counts: np.ndarray) -> np.ndarray:
    total = int(counts.sum())
    return counts / total if total else np.full(3, np.nan)


def quartet_summary(topologies: Iterable[int]) -> dict[str, Any]:
    counts = topology_counts(topologies)
    q = _q(counts)
    out: dict[str, Any] = {
        "n_loci": int(counts.sum()),
        "n1": int(counts[0]), "n2": int(counts[1]), "n3": int(counts[2]),
        "q1": float(q[0]), "q2": float(q[1]), "q3": float(q[2]),
        "dominant_topology": int(np.nanargmax(q)) if counts.sum() else -1,
    }
    if counts.sum():
        out.update(compare_models(counts))
    return out


def sliding_windows(locus_rows: list[Mapping[str, Any]], window_loci: int, step_loci: int) -> list[dict[str, Any]]:
    if window_loci <= 0 or step_loci <= 0:
        raise ValueError("window_loci and step_loci must be positive")
    rows = sorted(locus_rows, key=lambda r: int(r["position"]))
    out: list[dict[str, Any]] = []
    for window_id, start in enumerate(range(0, len(rows) - window_loci + 1, step_loci)):
        chunk = rows[start:start + window_loci]
        counts = topology_counts(int(r["topology_index"]) for r in chunk)
        q = _q(counts)
        model = compare_models(counts)
        first = int(chunk[0]["position"])
        last = int(chunk[-1]["position"])
        out.append({
            "window_id": window_id,
            "start_position": first,
            "end_position": last,
            "center_position": (first + last) / 2.0,
            "n_loci": len(chunk),
            "n1": int(counts[0]), "n2": int(counts[1]), "n3": int(counts[2]),
            "q1": float(q[0]), "q2": float(q[1]), "q3": float(q[2]),
            "dominant_topology": int(np.argmax(q)),
            "nearest_msc_topology": int(model["nearest_msc_topology"]),
            "distance_to_nearest_msc_arm": float(model["distance_to_nearest_msc_arm"]),
            "off_arm_difference": float(model["off_arm_difference"]),
            "off_arm_p_value": float(model["off_arm_p_value"]),
            "fraction_inside_rearrangement": float(np.mean([bool(r["inside_rearrangement"]) for r in chunk])),
        })
    return out


def breakpoint_jumps(locus_rows: list[Mapping[str, Any]], breakpoints: Iterable[int], bandwidth_loci: int) -> list[dict[str, Any]]:
    if bandwidth_loci <= 0:
        raise ValueError("breakpoint_bandwidth_loci must be positive")
    rows = sorted(locus_rows, key=lambda r: int(r["position"]))
    positions = np.asarray([int(r["position"]) for r in rows])
    out = []
    for bp in breakpoints:
        left_idx = np.where(positions < bp)[0][-bandwidth_loci:]
        right_idx = np.where(positions >= bp)[0][:bandwidth_loci]
        left = [rows[i] for i in left_idx]
        right = [rows[i] for i in right_idx]
        left_q = _q(topology_counts(int(r["topology_index"]) for r in left))
        right_q = _q(topology_counts(int(r["topology_index"]) for r in right))
        jump = float(np.linalg.norm(right_q - left_q)) if left and right else float("nan")
        out.append({
            "breakpoint_position": int(bp),
            "left_q1": float(left_q[0]), "left_q2": float(left_q[1]), "left_q3": float(left_q[2]),
            "right_q1": float(right_q[0]), "right_q2": float(right_q[1]), "right_q3": float(right_q[2]),
            "quartet_jump_l2": jump,
        })
    return out


def adjacent_window_jumps(windows: list[Mapping[str, Any]], breakpoints: Iterable[int], top_n: int = 5) -> list[dict[str, Any]]:
    bps = list(breakpoints)
    jumps = []
    for left, right in zip(windows, windows[1:]):
        ql = np.asarray([float(left["q1"]), float(left["q2"]), float(left["q3"])])
        qr = np.asarray([float(right["q1"]), float(right["q2"]), float(right["q3"])])
        pos = (float(left["center_position"]) + float(right["center_position"])) / 2.0
        nearest = min(abs(pos - bp) for bp in bps) if bps else float("nan")
        jumps.append({
            "left_window_id": int(left["window_id"]),
            "right_window_id": int(right["window_id"]),
            "between_center_position": pos,
            "quartet_jump_l2": float(np.linalg.norm(qr - ql)),
            "distance_to_nearest_breakpoint": float(nearest),
        })
    return sorted(jumps, key=lambda r: r["quartet_jump_l2"], reverse=True)[:top_n]


def spatial_summary(locus_rows: list[Mapping[str, Any]], windows: list[Mapping[str, Any]], interval: RearrangementInterval, bandwidth_loci: int) -> dict[str, Any]:
    inside = [r for r in locus_rows if bool(r["inside_rearrangement"])]
    outside = [r for r in locus_rows if not bool(r["inside_rearrangement"])]
    overall_s = quartet_summary(int(r["topology_index"]) for r in locus_rows)
    inside_s = quartet_summary(int(r["topology_index"]) for r in inside)
    outside_s = quartet_summary(int(r["topology_index"]) for r in outside)
    bps = [interval.start, interval.end]
    jumps = breakpoint_jumps(locus_rows, bps, bandwidth_loci)
    adjacent = adjacent_window_jumps(windows, bps)
    return {
        "linked_loci_model": False,
        "spatial_model_note": "Loci are conditionally independent given the realized rearrangement history and local region parameters; v0.7.0 does not simulate a multi-locus ARG.",
        "overall_q1": overall_s["q1"], "overall_q2": overall_s["q2"], "overall_q3": overall_s["q3"],
        "inside_n_loci": inside_s["n_loci"],
        "inside_q1": inside_s.get("q1"), "inside_q2": inside_s.get("q2"), "inside_q3": inside_s.get("q3"),
        "inside_dominant_topology": inside_s.get("dominant_topology"),
        "inside_distance_to_nearest_msc_arm": inside_s.get("distance_to_nearest_msc_arm"),
        "outside_n_loci": outside_s["n_loci"],
        "outside_q1": outside_s.get("q1"), "outside_q2": outside_s.get("q2"), "outside_q3": outside_s.get("q3"),
        "outside_dominant_topology": outside_s.get("dominant_topology"),
        "outside_distance_to_nearest_msc_arm": outside_s.get("distance_to_nearest_msc_arm"),
        "delta_q1_inside_minus_outside": inside_s.get("q1") - outside_s.get("q1"),
        "delta_q2_inside_minus_outside": inside_s.get("q2") - outside_s.get("q2"),
        "delta_q3_inside_minus_outside": inside_s.get("q3") - outside_s.get("q3"),
        "breakpoint_jumps": jumps,
        "strongest_adjacent_window_jumps": adjacent,
    }
