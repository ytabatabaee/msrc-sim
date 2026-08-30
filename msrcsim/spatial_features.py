from __future__ import annotations

from typing import Any
import math

import numpy as np

from .spatial_compare import SpatialModelRun


def _finite(values: list[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


def _topologies(run: SpatialModelRun) -> list[int]:
    if run.loci:
        return [int(row["topology_index"]) for row in run.loci]
    return [int(row["dominant_topology"]) for row in run.windows if int(row.get("dominant_topology", -1)) >= 0]


def _positions(run: SpatialModelRun) -> list[float]:
    if run.loci:
        return [float(row["position_bp"]) for row in run.loci]
    return [float(row["center_bp"]) for row in run.windows]


def infer_focal_topology(run: SpatialModelRun) -> int:
    q = np.asarray(run.marginal_q, dtype=float)
    return int(np.argmax(q))


def topology_same_probability(run: SpatialModelRun, lags_bp: list[float]) -> dict[str, float]:
    positions = np.asarray(_positions(run), dtype=float)
    tops = np.asarray(_topologies(run), dtype=int)
    order = np.argsort(positions)
    positions = positions[order]
    tops = tops[order]
    out: dict[str, float] = {}
    for lag in lags_bp:
        same = []
        for i, pos in enumerate(positions):
            j = int(np.searchsorted(positions, pos + float(lag), side="left"))
            if j < len(positions):
                same.append(tops[i] == tops[j])
        out[f"topology_same_prob_lag_{int(lag)}"] = float(np.mean(same)) if same else math.nan
    return out


def longest_contiguous_focal_interval(run: SpatialModelRun, focal_topology: int) -> float:
    positions = _positions(run)
    tops = _topologies(run)
    if not positions:
        return math.nan
    best = 0.0
    start_pos: float | None = None
    last_pos: float | None = None
    for pos, top in zip(positions, tops):
        if int(top) == int(focal_topology):
            if start_pos is None:
                start_pos = float(pos)
            last_pos = float(pos)
            best = max(best, last_pos - start_pos)
        else:
            start_pos = None
            last_pos = None
    return float(best)


def _inside_feature_interval(pos: float, run: SpatialModelRun) -> bool:
    for interval in run.feature_intervals:
        state = interval.get("state")
        label = str(interval.get("label", ""))
        if run.model_name == "Hybridization" and not (state == 1 or label == "introgressed"):
            continue
        if float(interval["start_bp"]) <= float(pos) <= float(interval["end_bp"]):
            return True
    return False


def focal_support_in_feature_fraction(run: SpatialModelRun, focal_topology: int) -> float:
    positions = _positions(run)
    tops = _topologies(run)
    focal_positions = [pos for pos, top in zip(positions, tops) if int(top) == int(focal_topology)]
    if not focal_positions:
        return math.nan
    return float(np.mean([_inside_feature_interval(pos, run) for pos in focal_positions]))


def _breakpoints(run: SpatialModelRun) -> list[float]:
    if run.model_name != "MSRC":
        return []
    bps = []
    for interval in run.feature_intervals:
        bps.extend([float(interval["start_bp"]), float(interval["end_bp"])])
    return bps


def change_point_distances_to_breakpoints(run: SpatialModelRun) -> tuple[float, float]:
    bps = _breakpoints(run)
    centers = [float(row["center_bp"]) for row in run.windows]
    dom = [int(row.get("dominant_topology", -1)) for row in run.windows]
    cps = [(a + b) / 2.0 for a, b, ta, tb in zip(centers, centers[1:], dom, dom[1:]) if ta >= 0 and tb >= 0 and ta != tb]
    if not cps or not bps:
        return (math.nan, math.nan)
    distances = [min(abs(cp - bp) for bp in bps) for cp in cps]
    return (float(np.mean(distances)), float(np.min(distances)))


def extract_spatial_features(
    run: SpatialModelRun,
    *,
    focal_topology: int | None = None,
    lags_bp: list[float] | None = None,
) -> dict[str, Any]:
    lags = [10000.0, 100000.0, 1000000.0] if lags_bp is None else [float(x) for x in lags_bp]
    focal = infer_focal_topology(run) if focal_topology is None else int(focal_topology)
    tops = _topologies(run)
    dom = [int(row.get("dominant_topology", -1)) for row in run.windows]
    distances = _finite([float(row.get("distance_to_nearest_msc_arm", math.nan)) for row in run.windows])
    intro_lengths = [
        float(interval["end_bp"]) - float(interval["start_bp"])
        for interval in run.feature_intervals
        if run.model_name == "Hybridization" and (interval.get("state") == 1 or str(interval.get("label", "")) == "introgressed")
    ]
    mean_cp, min_cp = change_point_distances_to_breakpoints(run)
    features: dict[str, Any] = {
        "marginal_q1": float(run.marginal_q[0]),
        "marginal_q2": float(run.marginal_q[1]),
        "marginal_q3": float(run.marginal_q[2]),
        "focal_topology": int(focal),
        "num_ancestry_tracts": len(run.feature_intervals) if run.model_name == "Hybridization" else 0,
        "mean_introgressed_tract_length": float(np.mean(intro_lengths)) if intro_lengths else math.nan,
        "median_introgressed_tract_length": float(np.median(intro_lengths)) if intro_lengths else math.nan,
        "num_dominant_topology_change_points": int(sum(a != b for a, b in zip(dom, dom[1:]) if a >= 0 and b >= 0)),
        "num_locus_topology_change_points": int(sum(a != b for a, b in zip(tops, tops[1:]))),
        "longest_contiguous_focal_topology_interval": longest_contiguous_focal_interval(run, focal),
        "focal_support_in_feature_fraction": focal_support_in_feature_fraction(run, focal),
        "mean_distance_to_nearest_msc_arm": float(np.mean(distances)) if distances.size else math.nan,
        "max_distance_to_nearest_msc_arm": float(np.max(distances)) if distances.size else math.nan,
        "mean_change_point_distance_to_msrc_breakpoint": mean_cp,
        "min_change_point_distance_to_msrc_breakpoint": min_cp,
    }
    features.update(topology_same_probability(run, lags))
    return features
