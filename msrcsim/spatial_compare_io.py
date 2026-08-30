from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable
import csv
import json
import math

import numpy as np

from .model_fitting import off_arm_statistics
from .spatial_compare import SpatialModelRun


MSRC_LOCI_FILES = ("msrc_loci.csv", "spatial_loci.csv")
MSRC_WINDOWS_FILES = ("msrc_windows.csv", "spatial_windows.csv")
MSRC_SUMMARY_FILES = ("msrc_summary.json", "spatial_summary.json", "summary.json")
HYB_LOCI_FILES = ("hybridization_loci.csv", "hyb_loci.csv")
HYB_WINDOWS_FILES = ("hybridization_windows.csv", "hyb_windows.csv")
HYB_SUMMARY_FILES = ("hybridization_summary.json", "hyb_summary.json", "summary.json")
HYB_TRACT_FILES = ("hybridization_tracts.csv", "hyb_tracts.csv")


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def _first_existing(run_dir: Path, names: Iterable[str], label: str, required: bool = True) -> Path | None:
    for name in names:
        path = run_dir / name
        if path.exists():
            return path
    if required:
        expected = ", ".join(f"`{name}`" for name in names)
        raise FileNotFoundError(f"{label} not found. Expected one of: {expected}")
    return None


def _as_float(value: Any, default: float = math.nan) -> float:
    if value is None or value == "":
        return default
    return float(value)


def _as_int(value: Any, default: int = 0) -> int:
    if value is None or value == "":
        return default
    return int(float(value))


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y", "inside", "rearrangement"}


def _pick(row: dict[str, Any], *keys: str, default: Any = None) -> Any:
    for key in keys:
        if key in row and row[key] != "":
            return row[key]
    return default


def _normalize_loci(rows: list[dict[str, Any]], model_type: str) -> list[dict[str, Any]]:
    loci: list[dict[str, Any]] = []
    for i, row in enumerate(rows):
        pos = _as_float(_pick(row, "position_bp", "position", "pos", default=math.nan))
        topo = _as_int(_pick(row, "topology_index", "topology_id", default=0))
        label = str(_pick(row, "topology_label", "topology", default=f"topology_{topo + 1}"))
        out = {
            "locus_id": _as_int(_pick(row, "locus_id", default=i), i),
            "position_bp": pos,
            "topology_index": topo,
            "topology_label": label,
            "model": str(_pick(row, "model", default=model_type)),
        }
        if model_type == "msrc":
            inside = _as_bool(_pick(row, "is_inside_rearranged_interval", "inside_rearrangement", "inside_rearranged_interval", default=False))
            out["is_inside_rearranged_interval"] = inside
            out["arrangement_state"] = _pick(row, "arrangement_state", "region", "local_model", default="rearranged" if inside else "background")
        else:
            ancestry_state = _as_int(_pick(row, "ancestry_state", default=0))
            out["ancestry_state"] = ancestry_state
            out["ancestry_label"] = str(_pick(row, "ancestry_label", default="introgressed" if ancestry_state == 1 else "major"))
            out["tract_id"] = _as_int(_pick(row, "tract_id", default=-1), -1)
        loci.append(out)
    return sorted(loci, key=lambda r: float(r["position_bp"]))


def _normalize_windows(rows: list[dict[str, Any]], model_type: str) -> list[dict[str, Any]]:
    windows: list[dict[str, Any]] = []
    for i, row in enumerate(rows):
        start = _as_float(_pick(row, "start_bp", "start_position", default=math.nan))
        end = _as_float(_pick(row, "end_bp", "end_position", default=math.nan))
        center = _as_float(_pick(row, "center_bp", "center_position", default=(start + end) / 2.0))
        out = {
            "window_id": _as_int(_pick(row, "window_id", default=i), i),
            "start_bp": start,
            "end_bp": end,
            "center_bp": center,
            "num_loci": _as_int(_pick(row, "num_loci", "n_loci", default=0)),
            "q1": _as_float(_pick(row, "q1", default=math.nan)),
            "q2": _as_float(_pick(row, "q2", default=math.nan)),
            "q3": _as_float(_pick(row, "q3", default=math.nan)),
            "dominant_topology": _as_int(_pick(row, "dominant_topology", default=-1), -1),
            "distance_to_nearest_msc_arm": _as_float(_pick(row, "distance_to_nearest_msc_arm", default=math.nan)),
            "off_arm_difference": _as_float(_pick(row, "off_arm_difference", default=math.nan)),
        }
        if model_type == "msrc":
            out["fraction_rearranged"] = _as_float(_pick(row, "fraction_rearranged", "fraction_inside_rearrangement", default=math.nan))
        else:
            out["introgressed_fraction"] = _as_float(_pick(row, "introgressed_fraction", default=math.nan))
        windows.append(out)
    return sorted(windows, key=lambda r: float(r["center_bp"]))


def _counts(topologies: Iterable[int]) -> np.ndarray:
    return np.bincount(np.asarray(list(topologies), dtype=int), minlength=3)[:3]


def compute_windows_from_loci(
    loci: list[dict[str, Any]],
    model_type: str = "msrc",
    window_size_loci: int = 50,
    step_loci: int = 10,
) -> list[dict[str, Any]]:
    if window_size_loci <= 0 or step_loci <= 0:
        raise ValueError("window_size_loci and step_loci must be positive")
    ordered = sorted(loci, key=lambda r: float(r["position_bp"]))
    if len(ordered) < window_size_loci:
        raise ValueError("not enough loci to compute one window")
    out: list[dict[str, Any]] = []
    for window_id, start_idx in enumerate(range(0, len(ordered) - window_size_loci + 1, step_loci)):
        chunk = ordered[start_idx:start_idx + window_size_loci]
        counts = _counts(_as_int(row["topology_index"]) for row in chunk)
        total = int(counts.sum())
        q = counts / total
        stats = off_arm_statistics(counts)
        first = float(chunk[0]["position_bp"])
        last = float(chunk[-1]["position_bp"])
        row: dict[str, Any] = {
            "window_id": window_id,
            "start_bp": first,
            "end_bp": last,
            "center_bp": (first + last) / 2.0,
            "num_loci": len(chunk),
            "q1": float(q[0]),
            "q2": float(q[1]),
            "q3": float(q[2]),
            "dominant_topology": int(np.argmax(q)),
            "distance_to_nearest_msc_arm": float(stats["distance_to_nearest_msc_arm"]),
            "off_arm_difference": float(stats["off_arm_difference"]),
        }
        if model_type == "msrc":
            row["fraction_rearranged"] = float(np.mean([bool(r.get("is_inside_rearranged_interval", False)) for r in chunk]))
        else:
            row["introgressed_fraction"] = float(np.mean([_as_int(r.get("ancestry_state", 0)) == 1 for r in chunk]))
        out.append(row)
    return out


def _fill_window_fraction_from_loci(windows: list[dict[str, Any]], loci: list[dict[str, Any]], model_type: str) -> None:
    for window in windows:
        if model_type == "msrc":
            key = "fraction_rearranged"
            if not math.isnan(float(window.get(key, math.nan))):
                continue
            chunk = [r for r in loci if float(window["start_bp"]) <= float(r["position_bp"]) <= float(window["end_bp"])]
            window[key] = float(np.mean([bool(r.get("is_inside_rearranged_interval", False)) for r in chunk])) if chunk else math.nan
        else:
            key = "introgressed_fraction"
            if not math.isnan(float(window.get(key, math.nan))):
                continue
            chunk = [r for r in loci if float(window["start_bp"]) <= float(r["position_bp"]) <= float(window["end_bp"])]
            window[key] = float(np.mean([_as_int(r.get("ancestry_state", 0)) == 1 for r in chunk])) if chunk else math.nan


def compute_marginal_q_from_loci_or_windows(
    run_or_loci: SpatialModelRun | list[dict[str, Any]],
    windows: list[dict[str, Any]] | None = None,
    summary: dict[str, Any] | None = None,
) -> tuple[float, float, float]:
    if isinstance(run_or_loci, SpatialModelRun):
        loci = run_or_loci.loci
        windows = run_or_loci.windows
        summary = run_or_loci.summary
    else:
        loci = run_or_loci
        windows = [] if windows is None else windows
        summary = {} if summary is None else summary
    for key in ("observed_marginal_q", "expected_marginal_q", "topology_frequencies"):
        value = summary.get(key)
        if isinstance(value, list) and len(value) >= 3:
            return (float(value[0]), float(value[1]), float(value[2]))
    if all(k in summary for k in ("overall_q1", "overall_q2", "overall_q3")):
        return (float(summary["overall_q1"]), float(summary["overall_q2"]), float(summary["overall_q3"]))
    if loci:
        counts = _counts(_as_int(row["topology_index"]) for row in loci)
        q = counts / counts.sum()
        return (float(q[0]), float(q[1]), float(q[2]))
    weighted = np.zeros(3, dtype=float)
    total = 0
    for row in windows:
        n = _as_int(row.get("num_loci", 1), 1)
        weighted += n * np.asarray([float(row["q1"]), float(row["q2"]), float(row["q3"])], dtype=float)
        total += n
    if total <= 0:
        raise ValueError("cannot compute marginal quartet vector from empty run")
    q = weighted / total
    return (float(q[0]), float(q[1]), float(q[2]))


def _chromosome_length(summary: dict[str, Any], loci: list[dict[str, Any]], windows: list[dict[str, Any]], intervals: list[dict[str, Any]]) -> float:
    for key in ("chromosome_length_bp", "genome_length", "chromosome_length"):
        if key in summary:
            return float(summary[key])
    maxima = [float(row["position_bp"]) for row in loci]
    maxima += [float(row["end_bp"]) for row in windows]
    maxima += [float(row["end_bp"]) for row in intervals]
    return max(maxima) if maxima else 1.0


def _msrc_intervals(summary: dict[str, Any]) -> list[dict[str, Any]]:
    raw = summary.get("rearranged_intervals") or summary.get("rearrangement_intervals")
    if raw is None and "rearrangement_interval" in summary:
        raw = [summary["rearrangement_interval"]]
    intervals = []
    for i, item in enumerate(raw or []):
        start = _as_float(_pick(item, "start_bp", "start", default=math.nan))
        end = _as_float(_pick(item, "end_bp", "end", default=math.nan))
        intervals.append({"interval_id": str(_pick(item, "interval_id", "id", default=f"rearrangement_{i + 1}")), "start_bp": start, "end_bp": end, "state": "rearranged"})
    return intervals


def _tract_intervals(path: Path | None) -> list[dict[str, Any]]:
    if path is None:
        return []
    intervals = []
    for row in _read_csv(path):
        intervals.append({
            "interval_id": str(_pick(row, "tract_id", default=len(intervals))),
            "start_bp": _as_float(_pick(row, "start_bp", "start", default=math.nan)),
            "end_bp": _as_float(_pick(row, "end_bp", "end", default=math.nan)),
            "state": _as_int(_pick(row, "ancestry_state", default=0)),
            "label": str(_pick(row, "ancestry_label", default="introgressed" if _as_int(_pick(row, "ancestry_state", default=0)) == 1 else "major")),
        })
    return intervals


def load_spatial_model_run(
    run_dir: str | Path,
    model_type: str,
    *,
    window_size_loci: int = 50,
    step_loci: int = 10,
) -> SpatialModelRun:
    run_path = Path(run_dir)
    normalized_model = model_type.lower()
    if normalized_model not in {"msrc", "hybridization", "hyb"}:
        raise ValueError("model_type must be 'msrc' or 'hybridization'")
    normalized_model = "hybridization" if normalized_model == "hyb" else normalized_model

    if normalized_model == "msrc":
        loci_path = _first_existing(run_path, MSRC_LOCI_FILES, "Spatial MSRC locus outputs", required=False)
        windows_path = _first_existing(run_path, MSRC_WINDOWS_FILES, "Spatial MSRC window outputs", required=False)
        summary_path = _first_existing(run_path, MSRC_SUMMARY_FILES, "Spatial MSRC summary", required=False)
    else:
        loci_path = _first_existing(run_path, HYB_LOCI_FILES, "Hybridization locus outputs", required=False)
        windows_path = _first_existing(run_path, HYB_WINDOWS_FILES, "Hybridization window outputs", required=False)
        summary_path = _first_existing(run_path, HYB_SUMMARY_FILES, "Hybridization summary", required=False)

    if loci_path is None and windows_path is None:
        if normalized_model == "msrc":
            expected = ", ".join(f"`{name}`" for name in (*MSRC_WINDOWS_FILES, *MSRC_LOCI_FILES))
            raise FileNotFoundError(f"Spatial MSRC outputs not found. Expected one of: {expected}")
        expected = ", ".join(f"`{name}`" for name in (*HYB_WINDOWS_FILES, *HYB_LOCI_FILES))
        raise FileNotFoundError(f"Hybridization outputs not found. Expected one of: {expected}")
    summary = _read_json(summary_path) if summary_path else {}
    if normalized_model == "msrc":
        feature_intervals = _msrc_intervals(summary)
    else:
        tract_path = _first_existing(run_path, HYB_TRACT_FILES, "Hybridization tract outputs", required=False)
        feature_intervals = _tract_intervals(tract_path)
    loci = _normalize_loci(_read_csv(loci_path), "msrc" if normalized_model == "msrc" else "hybridization") if loci_path else []
    if windows_path is not None:
        windows = _normalize_windows(_read_csv(windows_path), "msrc" if normalized_model == "msrc" else "hybridization")
    else:
        windows = compute_windows_from_loci(loci, "msrc" if normalized_model == "msrc" else "hybridization", window_size_loci, step_loci)
    _fill_window_fraction_from_loci(windows, loci, "msrc" if normalized_model == "msrc" else "hybridization")
    marginal_q = compute_marginal_q_from_loci_or_windows(loci, windows, summary)
    length = _chromosome_length(summary, loci, windows, feature_intervals)
    model_name = "MSRC" if normalized_model == "msrc" else "Hybridization"
    return SpatialModelRun(model_name, length, loci, windows, summary, marginal_q, feature_intervals, run_path)
