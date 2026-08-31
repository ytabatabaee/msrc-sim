from __future__ import annotations

from dataclasses import dataclass, asdict
from math import log, sqrt
from typing import Any

import numpy as np
from scipy.optimize import minimize

from .analytic import TOPOLOGY_NAMES
from .hybridization_match import mixture_q
from .model_fitting import off_arm_statistics
from .spatial_compare import SpatialModelRun
from .spatial_compare_io import compute_windows_from_loci


@dataclass(frozen=True)
class ExactGammaMatch:
    gamma: float
    t_major: float
    t_introgressed: float
    major_topology: int
    introgressed_topology: int
    shared_topology: int
    target_q: tuple[float, float, float]
    fitted_q: tuple[float, float, float]
    l1_error: float
    l2_error: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def feasible_gamma_bounds(
    target_q: tuple[float, float, float] | list[float] | np.ndarray,
    major_topology: int = 0,
    introgressed_topology: int = 1,
) -> tuple[float, float]:
    target = np.asarray(target_q, dtype=float)
    if target.shape != (3,) or np.any(target < 0.0) or not np.isclose(target.sum(), 1.0, atol=1e-6):
        raise ValueError("target_q must contain three nonnegative values summing to one")
    if major_topology == introgressed_topology:
        raise ValueError("major and introgressed topologies must differ")
    shared = ({0, 1, 2} - {int(major_topology), int(introgressed_topology)}).pop()
    d_major = float(target[int(major_topology)] - target[shared])
    d_intro = float(target[int(introgressed_topology)] - target[shared])
    return (max(0.0, d_intro), min(1.0, 1.0 - d_major))


def exact_match_for_gamma(
    target_q: tuple[float, float, float] | list[float] | np.ndarray,
    gamma: float,
    major_topology: int = 0,
    introgressed_topology: int = 1,
    *,
    max_branch_length: float = 50.0,
) -> ExactGammaMatch:
    target = np.asarray(target_q, dtype=float)
    low, high = feasible_gamma_bounds(target, major_topology, introgressed_topology)
    gamma = float(gamma)
    if gamma < low - 1e-12 or gamma > high + 1e-12:
        raise ValueError(f"gamma is outside the feasible interval [{low}, {high}]")
    shared = ({0, 1, 2} - {int(major_topology), int(introgressed_topology)}).pop()
    d_major = float(target[int(major_topology)] - target[shared])
    d_intro = float(target[int(introgressed_topology)] - target[shared])
    major_arg = max(0.0, 1.0 - d_major / max(1e-15, 1.0 - gamma))
    intro_arg = max(0.0, 1.0 - d_intro / max(1e-15, gamma))
    t_major = max_branch_length if major_arg <= 0.0 else min(max_branch_length, -log(major_arg))
    t_intro = max_branch_length if intro_arg <= 0.0 else min(max_branch_length, -log(intro_arg))
    fitted = mixture_q(int(major_topology), int(introgressed_topology), gamma, t_major, t_intro)
    diff = fitted - target
    return ExactGammaMatch(
        gamma=gamma,
        t_major=float(t_major),
        t_introgressed=float(t_intro),
        major_topology=int(major_topology),
        introgressed_topology=int(introgressed_topology),
        shared_topology=int(shared),
        target_q=(float(target[0]), float(target[1]), float(target[2])),
        fitted_q=(float(fitted[0]), float(fitted[1]), float(fitted[2])),
        l1_error=float(np.sum(np.abs(diff))),
        l2_error=float(sqrt(float(np.dot(diff, diff)))),
    )


def gamma_grid(
    target_q: tuple[float, float, float],
    num_gamma: int,
    major_topology: int = 0,
    introgressed_topology: int = 1,
) -> list[float]:
    if num_gamma <= 0:
        raise ValueError("num_gamma must be positive")
    low, high = feasible_gamma_bounds(target_q, major_topology, introgressed_topology)
    if low > high:
        raise ValueError("target q is not feasible for the selected topology pair")
    if num_gamma == 1:
        return [(low + high) / 2.0]
    width = high - low
    return [float(low + width * (i + 1) / (num_gamma + 1)) for i in range(num_gamma)]


def eta_to_hr(eta: float, rearranged_interval_length: float, gamma: float) -> float:
    eta = float(eta)
    length = float(rearranged_interval_length)
    gamma = float(gamma)
    if eta <= 0.0 or length <= 0.0:
        raise ValueError("eta and rearranged interval length must be positive")
    if gamma >= 1.0:
        raise ValueError("gamma must be less than one")
    return float(1.0 / (eta * length * max(1e-15, 1.0 - gamma)))


def rearranged_interval_length(run: SpatialModelRun) -> float:
    length = 0.0
    for interval in run.feature_intervals:
        length += max(0.0, float(interval["end_bp"]) - float(interval["start_bp"]))
    if length <= 0.0:
        raise ValueError("MSRC run does not contain a positive rearranged interval length")
    return float(length)


def _topology_counts(loci: list[dict[str, Any]]) -> np.ndarray:
    return np.bincount(np.asarray([int(row["topology_index"]) for row in loci], dtype=int), minlength=3)[:3]


def resample_msrc_chromosome(template: SpatialModelRun, rng: np.random.Generator, window_size_loci: int, step_loci: int) -> SpatialModelRun:
    inside = [row for row in template.loci if bool(row.get("is_inside_rearranged_interval", False))]
    outside = [row for row in template.loci if not bool(row.get("is_inside_rearranged_interval", False))]
    q_inside = _topology_counts(inside) / max(1, len(inside)) if inside else np.asarray(template.marginal_q, dtype=float)
    q_outside = _topology_counts(outside) / max(1, len(outside)) if outside else np.asarray(template.marginal_q, dtype=float)
    loci = []
    for row in template.loci:
        out = dict(row)
        q = q_inside if bool(row.get("is_inside_rearranged_interval", False)) else q_outside
        top = int(rng.choice(3, p=q / q.sum()))
        out["topology_index"] = top
        out["topology_label"] = TOPOLOGY_NAMES[top]
        loci.append(out)
    windows = compute_windows_from_loci(loci, "msrc", window_size_loci, step_loci)
    counts = _topology_counts(loci)
    observed = counts / counts.sum()
    summary = dict(template.summary)
    summary["observed_marginal_q"] = [float(x) for x in observed]
    return SpatialModelRun(
        "MSRC",
        template.chromosome_length_bp,
        loci,
        windows,
        summary,
        template.marginal_q,
        template.feature_intervals,
        template.run_dir,
    )


def _positions_and_topologies(run: SpatialModelRun) -> tuple[np.ndarray, np.ndarray]:
    positions = np.asarray([float(row["position_bp"]) for row in run.loci], dtype=float)
    topologies = np.asarray([int(row["topology_index"]) for row in run.loci], dtype=int)
    order = np.argsort(positions)
    return positions[order], topologies[order]


def _run_lengths(positions: np.ndarray, labels: np.ndarray) -> list[float]:
    if len(positions) == 0:
        return []
    out = []
    start = float(positions[0])
    last = float(positions[0])
    current = int(labels[0])
    for pos, label in zip(positions[1:], labels[1:]):
        if int(label) != current:
            out.append(max(0.0, last - start))
            start = float(pos)
            current = int(label)
        last = float(pos)
    out.append(max(0.0, last - start))
    return out


def _dominant_run_lengths(run: SpatialModelRun) -> list[float]:
    centers = np.asarray([float(row["center_bp"]) for row in run.windows], dtype=float)
    dom = np.asarray([int(row.get("dominant_topology", -1)) for row in run.windows], dtype=int)
    valid = dom >= 0
    return _run_lengths(centers[valid], dom[valid])


def _same_topology_probabilities(positions: np.ndarray, topologies: np.ndarray, lags_bp: list[float]) -> dict[str, float]:
    out = {}
    for lag in lags_bp:
        same = []
        for i, pos in enumerate(positions):
            j = int(np.searchsorted(positions, pos + float(lag), side="left"))
            if j < len(positions):
                same.append(topologies[i] == topologies[j])
        out[f"topology_same_prob_lag_{int(lag)}"] = float(np.mean(same)) if same else float("nan")
    return out


def _breakpoints(reference: SpatialModelRun) -> list[float]:
    bps = []
    for interval in reference.feature_intervals:
        bps.extend([float(interval["start_bp"]), float(interval["end_bp"])])
    return bps


def _inside_reference_interval(pos: float, reference: SpatialModelRun) -> bool:
    return any(float(interval["start_bp"]) <= pos <= float(interval["end_bp"]) for interval in reference.feature_intervals)


def _change_point_distances(run: SpatialModelRun, reference: SpatialModelRun) -> tuple[float, float]:
    bps = _breakpoints(reference)
    centers = [float(row["center_bp"]) for row in run.windows]
    dom = [int(row.get("dominant_topology", -1)) for row in run.windows]
    cps = [(left + right) / 2.0 for left, right, a, b in zip(centers, centers[1:], dom, dom[1:]) if a >= 0 and b >= 0 and a != b]
    if not cps or not bps:
        return (float("nan"), float("nan"))
    distances = [min(abs(cp - bp) for bp in bps) for cp in cps]
    return (float(np.mean(distances)), float(np.min(distances)))


def _longest_focal_run(positions: np.ndarray, topologies: np.ndarray, focal_topology: int) -> float:
    if len(positions) == 0:
        return float("nan")
    best = 0.0
    start: float | None = None
    last: float | None = None
    for pos, top in zip(positions, topologies):
        if int(top) == int(focal_topology):
            if start is None:
                start = float(pos)
            last = float(pos)
            best = max(best, last - start)
        else:
            start = None
            last = None
    return float(best) if best > 0.0 or np.any(topologies == int(focal_topology)) else float("nan")


def _observed_q(run: SpatialModelRun) -> tuple[float, float, float]:
    if run.loci:
        counts = _topology_counts(run.loci)
        q = counts / counts.sum()
        return (float(q[0]), float(q[1]), float(q[2]))
    return run.marginal_q


def identifiability_features(
    run: SpatialModelRun,
    reference_msrc: SpatialModelRun,
    *,
    focal_topology: int,
    lags_bp: list[float],
) -> dict[str, Any]:
    positions, topologies = _positions_and_topologies(run)
    observed_q = _observed_q(run)
    dominant_lengths = _dominant_run_lengths(run)
    distances = np.asarray([float(row.get("distance_to_nearest_msc_arm", float("nan"))) for row in run.windows], dtype=float)
    distances = distances[np.isfinite(distances)]
    focal_positions = positions[topologies == int(focal_topology)]
    mean_cp, min_cp = _change_point_distances(run, reference_msrc)
    intro_tracts = [
        interval for interval in run.feature_intervals
        if run.model_name == "Hybridization" and (interval.get("state") == 1 or str(interval.get("label", "")) == "introgressed")
    ]
    all_tracts = run.feature_intervals if run.model_name == "Hybridization" else []
    realized_intro = run.summary.get("realized_introgressed_fraction_by_length", run.summary.get("realized_introgressed_fraction_at_loci"))
    features = {
        "marginal_q1": float(run.marginal_q[0]),
        "marginal_q2": float(run.marginal_q[1]),
        "marginal_q3": float(run.marginal_q[2]),
        "observed_marginal_q1": observed_q[0],
        "observed_marginal_q2": observed_q[1],
        "observed_marginal_q3": observed_q[2],
        "num_dominant_topology_change_points": int(sum(int(a.get("dominant_topology", -1)) != int(b.get("dominant_topology", -1)) for a, b in zip(run.windows, run.windows[1:]))),
        "mean_dominant_topology_run_length": float(np.mean(dominant_lengths)) if dominant_lengths else float("nan"),
        "longest_contiguous_focal_topology_interval": _longest_focal_run(positions, topologies, focal_topology),
        "mean_moving_window_off_arm_distance": float(np.mean(distances)) if distances.size else float("nan"),
        "max_moving_window_off_arm_distance": float(np.max(distances)) if distances.size else float("nan"),
        "fraction_focal_support_inside_rearranged_interval": float(np.mean([_inside_reference_interval(float(pos), reference_msrc) for pos in focal_positions])) if len(focal_positions) else float("nan"),
        "mean_change_point_distance_to_msrc_breakpoint": mean_cp,
        "min_change_point_distance_to_msrc_breakpoint": min_cp,
        "num_ancestry_tracts": len(all_tracts),
        "realized_introgressed_fraction": float(realized_intro) if realized_intro is not None else float("nan"),
        "mean_introgressed_tract_length": float(np.mean([float(t["end_bp"]) - float(t["start_bp"]) for t in intro_tracts])) if intro_tracts else float("nan"),
    }
    features.update(_same_topology_probabilities(positions, topologies, lags_bp))
    return features


BAG_FEATURE_COLUMNS = ["marginal_q1", "marginal_q2", "marginal_q3"]


def spatial_classifier_columns(rows: list[dict[str, Any]]) -> list[str]:
    excluded = {
        "model",
        "gamma",
        "eta",
        "hr",
        "h",
        "r",
        "replicate",
        "target_q1",
        "target_q2",
        "target_q3",
        "t_major",
        "t_introgressed",
        "expected_introgressed_tract_length_bp",
        "marginal_q1",
        "marginal_q2",
        "marginal_q3",
        "observed_marginal_q1",
        "observed_marginal_q2",
        "observed_marginal_q3",
        "num_ancestry_tracts",
        "realized_introgressed_fraction",
        "mean_introgressed_tract_length",
    }
    if not rows:
        return []
    return [
        key for key, value in rows[0].items()
        if key not in excluded and isinstance(value, (int, float))
    ]


def _prepare_matrix(rows: list[dict[str, Any]], feature_cols: list[str]) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray([
        [float(row[col]) if row.get(col) not in (None, "") and np.isfinite(float(row[col])) else 0.0 for col in feature_cols]
        for row in rows
    ], dtype=float)
    y = np.asarray([1 if row["model"] == "hybridization" else 0 for row in rows], dtype=int)
    return x, y


def _fit_logistic(train_x: np.ndarray, train_y: np.ndarray) -> np.ndarray:
    x = np.column_stack([np.ones(train_x.shape[0]), train_x])

    def objective(beta: np.ndarray) -> float:
        z = x @ beta
        return float(np.sum(np.logaddexp(0.0, z) - train_y * z) + 0.5 * 1e-4 * np.dot(beta[1:], beta[1:]))

    result = minimize(objective, np.zeros(x.shape[1], dtype=float), method="BFGS")
    return np.asarray(result.x, dtype=float)


def _auc(y_true: np.ndarray, scores: np.ndarray) -> float:
    pos = scores[y_true == 1]
    neg = scores[y_true == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    wins = 0.0
    for score in pos:
        wins += float(np.sum(score > neg)) + 0.5 * float(np.sum(score == neg))
    return float(wins / (len(pos) * len(neg)))


def _sigmoid(values: np.ndarray) -> np.ndarray:
    out = np.empty_like(values, dtype=float)
    positive = values >= 0.0
    out[positive] = 1.0 / (1.0 + np.exp(-values[positive]))
    exp_values = np.exp(values[~positive])
    out[~positive] = exp_values / (1.0 + exp_values)
    return out


def _ci95(values: np.ndarray) -> tuple[float, float]:
    if len(values) == 0:
        return (float("nan"), float("nan"))
    if len(values) == 1:
        value = float(values[0])
        return (value, value)
    mean = float(np.mean(values))
    half_width = float(1.96 * np.std(values, ddof=1) / sqrt(len(values)))
    return (max(0.0, mean - half_width), min(1.0, mean + half_width))


def stratified_cv_metrics(
    rows: list[dict[str, Any]],
    feature_cols: list[str],
    *,
    folds: int = 5,
    seed: int = 1,
) -> dict[str, float | int]:
    if folds < 2:
        raise ValueError("folds must be at least 2")
    x, y = _prepare_matrix(rows, feature_cols)
    rng = np.random.default_rng(seed)
    fold_ids = np.zeros(len(y), dtype=int)
    for label in (0, 1):
        indices = np.where(y == label)[0]
        rng.shuffle(indices)
        for i, idx in enumerate(indices):
            fold_ids[idx] = i % min(folds, len(indices))
    n_folds = int(max(fold_ids) + 1)
    accuracies = []
    aucs = []
    for fold in range(n_folds):
        test = fold_ids == fold
        train = ~test
        mean = x[train].mean(axis=0)
        sd = x[train].std(axis=0)
        sd[sd < 1e-12] = 1.0
        train_x = (x[train] - mean) / sd
        test_x = (x[test] - mean) / sd
        beta = _fit_logistic(train_x, y[train])
        logits = np.column_stack([np.ones(test_x.shape[0]), test_x]) @ beta
        scores = _sigmoid(logits)
        pred = (scores >= 0.5).astype(int)
        accuracies.append(float(np.mean(pred == y[test])))
        aucs.append(_auc(y[test], scores))
    acc = np.asarray(accuracies, dtype=float)
    auc = np.asarray(aucs, dtype=float)
    auc = auc[np.isfinite(auc)]
    acc_low, acc_high = _ci95(acc)
    auc_low, auc_high = _ci95(auc)
    return {
        "accuracy": float(np.mean(acc)),
        "accuracy_ci95_low": acc_low,
        "accuracy_ci95_high": acc_high,
        "roc_auc": float(np.mean(auc)) if len(auc) else float("nan"),
        "roc_auc_ci95_low": auc_low,
        "roc_auc_ci95_high": auc_high,
        "folds": n_folds,
    }


def off_arm_distance_for_q(q: tuple[float, float, float]) -> float:
    counts = np.asarray(q, dtype=float) * 1000000
    return float(off_arm_statistics(counts.astype(int))["distance_to_nearest_msc_arm"])
