from __future__ import annotations

from pathlib import Path
from typing import Any
import math
import os
import tempfile

import numpy as np

from .spatial_compare import SpatialModelRun, marginal_l1_distance
from .simplex import barycentric_to_cartesian


QUARTET_COLORS = ["#1b9e77", "#d95f02", "#7570b3"]
QUARTET_LABELS = ["12|34", "13|24", "14|23"]
MSRC_FEATURE = "#3b6fb6"
HYB_MAJOR = "#8aa1b4"
HYB_INTRO = "#d95f02"


def configure_matplotlib_cache() -> None:
    if "MPLCONFIGDIR" not in os.environ:
        cache = Path(tempfile.gettempdir()) / "msrcsim-matplotlib"
        cache.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(cache)
    if "XDG_CACHE_HOME" not in os.environ:
        cache_home = Path(tempfile.gettempdir()) / "msrcsim-cache"
        cache_home.mkdir(parents=True, exist_ok=True)
        os.environ["XDG_CACHE_HOME"] = str(cache_home)


def _xmax(msrc_run: SpatialModelRun, hyb_run: SpatialModelRun) -> float:
    return max(float(msrc_run.chromosome_length_bp), float(hyb_run.chromosome_length_bp), 1.0)


def _centers(run: SpatialModelRun) -> list[float]:
    return [float(row["center_bp"]) for row in run.windows]


def _values(run: SpatialModelRun, key: str) -> list[float]:
    return [float(row.get(key, math.nan)) for row in run.windows]


def _setup_common_style(plt: Any) -> None:
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    })


def _plot_msrc_state(ax: Any, run: SpatialModelRun, xmax: float) -> None:
    ax.broken_barh([(0, xmax)], (0.38, 0.24), facecolors="#e6e6e6", edgecolors="0.55", linewidth=0.6)
    for interval in run.feature_intervals:
        start = float(interval["start_bp"])
        end = float(interval["end_bp"])
        ax.broken_barh([(start, end - start)], (0.32, 0.36), facecolors=MSRC_FEATURE, edgecolors="white", linewidth=0.6)
        ax.axvline(start, color="0.2", linewidth=0.65, alpha=0.8)
        ax.axvline(end, color="0.2", linewidth=0.65, alpha=0.8)
    ax.set_title("Structural state / rearranged interval")
    ax.set_yticks([])
    ax.set_ylim(0, 1)
    ax.set_xlim(0, xmax)
    ax.spines[["left", "bottom"]].set_visible(False)
    ax.tick_params(axis="x", length=0, labelbottom=False)


def _plot_hyb_state(ax: Any, run: SpatialModelRun, xmax: float) -> None:
    if run.feature_intervals:
        for interval in run.feature_intervals:
            start = float(interval["start_bp"])
            end = float(interval["end_bp"])
            label = str(interval.get("label", "major"))
            color = HYB_INTRO if label == "introgressed" or int(interval.get("state", 0)) == 1 else HYB_MAJOR
            ax.broken_barh([(start, end - start)], (0.32, 0.36), facecolors=color, edgecolors="white", linewidth=0.35)
            if start > 0:
                ax.axvline(start, color="0.2", linewidth=0.35, alpha=0.45)
    else:
        ax.broken_barh([(0, xmax)], (0.38, 0.24), facecolors=HYB_MAJOR, edgecolors="0.55", linewidth=0.6)
    ax.set_title("Ancestry state / introgressed tracts")
    ax.set_yticks([])
    ax.set_ylim(0, 1)
    ax.set_xlim(0, xmax)
    ax.spines[["left", "bottom"]].set_visible(False)
    ax.tick_params(axis="x", length=0, labelbottom=False)


def _plot_quartets(ax: Any, run: SpatialModelRun, *, show_ylabel: bool) -> None:
    x = _centers(run)
    labels = run.summary.get("topology_names", QUARTET_LABELS)
    if not isinstance(labels, list) or len(labels) < 3:
        labels = QUARTET_LABELS
    for key, label, color in zip(("q1", "q2", "q3"), labels, QUARTET_COLORS):
        ax.plot(x, _values(run, key), color=color, linewidth=1.55, label=str(label))
    ax.axhline(1 / 3, color="0.45", linestyle="--", linewidth=0.85)
    ax.set_ylim(0, 1.02)
    if show_ylabel:
        ax.set_ylabel("Moving-window quartet support")


def _plot_fraction(ax: Any, run: SpatialModelRun, key: str, color: str, *, show_ylabel: bool) -> None:
    ax.plot(_centers(run), _values(run, key), color=color, linewidth=1.55)
    ax.set_ylim(0, 1.02)
    if show_ylabel:
        ax.set_ylabel("Local feature fraction")


def _distance_summary(run: SpatialModelRun) -> tuple[float, float, float]:
    vals = np.asarray([float(row.get("distance_to_nearest_msc_arm", math.nan)) for row in run.windows], dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return (math.nan, math.nan, math.nan)
    threshold = 0.05
    return (float(np.mean(vals)), float(np.max(vals)), float(np.mean(vals > threshold)))


def _plot_bag_summary(ax: Any, run: SpatialModelRun, other: SpatialModelRun, *, show_ylabel: bool) -> None:
    x = np.arange(3)
    ax.bar(x, run.marginal_q, color=QUARTET_COLORS, width=0.65)
    ax.axhline(1 / 3, color="0.45", linestyle="--", linewidth=0.85)
    ax.set_xticks(x, ["q1", "q2", "q3"])
    ax.set_ylim(0, 1.02)
    if show_ylabel:
        ax.set_ylabel("Genome-wide support")
    l1 = marginal_l1_distance(run.marginal_q, other.marginal_q)
    mean_d, max_d, frac_d = _distance_summary(run)
    text = (
        f"L1 to other = {l1:.3f}\n"
        f"MSC-arm distance\n"
        f"mean {mean_d:.3f}  max {max_d:.3f}\n"
        f"frac > 0.05: {frac_d:.2f}"
    )
    ax.text(0.98, 0.94, text, transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "0.75", "boxstyle": "round,pad=0.25", "linewidth": 0.6})


def _simplex_xy(q1: float, q2: float, q3: float) -> tuple[float, float] | None:
    vals = np.asarray([q1, q2, q3], dtype=float)
    if not np.all(np.isfinite(vals)):
        return None
    total = float(vals.sum())
    if total <= 0:
        return None
    vals = vals / total
    return barycentric_to_cartesian(float(vals[0]), float(vals[1]), float(vals[2]))


def _plot_simplex(ax: Any, run: SpatialModelRun, *, title: str) -> None:
    verts = np.asarray([[0.0, 0.0], [1.0, 0.0], [0.5, math.sqrt(3.0) / 2.0], [0.0, 0.0]])
    centroid = np.asarray([0.5, math.sqrt(3.0) / 6.0])
    ax.plot(verts[:, 0], verts[:, 1], color="0.2", linewidth=0.9)
    for vertex in verts[:3]:
        ax.plot([centroid[0], vertex[0]], [centroid[1], vertex[1]], color="0.45", linestyle="--", linewidth=0.8)

    points = [_simplex_xy(float(row.get("q1", math.nan)), float(row.get("q2", math.nan)), float(row.get("q3", math.nan))) for row in run.windows]
    points = [p for p in points if p is not None]
    if points:
        arr = np.asarray(points)
        color = MSRC_FEATURE if str(run.model_name).lower() == "msrc" else HYB_INTRO
        ax.scatter(arr[:, 0], arr[:, 1], s=16, color=color, alpha=0.48, linewidths=0, label="windows")

    marginal = _simplex_xy(float(run.marginal_q[0]), float(run.marginal_q[1]), float(run.marginal_q[2]))
    if marginal is not None:
        ax.scatter([marginal[0]], [marginal[1]], marker="X", s=72, color="black", linewidths=0.6, label="genome-wide")

    for label, xy, offset in (
        ("q1", (0.0, 0.0), (-0.035, -0.035)),
        ("q2", (1.0, 0.0), (0.035, -0.035)),
        ("q3", (0.5, math.sqrt(3.0) / 2.0), (0.0, 0.035)),
    ):
        ax.text(xy[0] + offset[0], xy[1] + offset[1], label, ha="center", va="center", fontsize=8, color="0.25")
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.08, math.sqrt(3.0) / 2.0 + 0.08)
    ax.axis("off")
    ax.legend(frameon=False, loc="lower center", ncols=2, fontsize=7)


def make_spatial_compare_figure(msrc_run: SpatialModelRun, hyb_run: SpatialModelRun, *, title: str | None = None) -> Any:
    configure_matplotlib_cache()
    import matplotlib.pyplot as plt

    _setup_common_style(plt)
    xmax = _xmax(msrc_run, hyb_run)
    fig, axes = plt.subplots(
        5,
        2,
        figsize=(12.0, 10.7),
        constrained_layout=True,
        gridspec_kw={"height_ratios": [0.7, 2.25, 1.15, 1.45, 1.85], "hspace": 0.28, "wspace": 0.12},
    )
    axes[0, 0].set_ylabel("MSRC", rotation=0, labelpad=28, va="center", fontweight="bold")
    axes[0, 1].set_ylabel("Hybridization", rotation=0, labelpad=46, va="center", fontweight="bold")
    _plot_msrc_state(axes[0, 0], msrc_run, xmax)
    _plot_hyb_state(axes[0, 1], hyb_run, xmax)
    _plot_quartets(axes[1, 0], msrc_run, show_ylabel=True)
    _plot_quartets(axes[1, 1], hyb_run, show_ylabel=False)
    axes[1, 1].legend(frameon=False, loc="upper right", ncols=3)
    _plot_fraction(axes[2, 0], msrc_run, "fraction_rearranged", MSRC_FEATURE, show_ylabel=True)
    _plot_fraction(axes[2, 1], hyb_run, "introgressed_fraction", HYB_INTRO, show_ylabel=False)
    axes[2, 0].set_title("Window overlap with rearranged interval")
    axes[2, 1].set_title("Introgressed ancestry fraction")
    _plot_bag_summary(axes[3, 0], msrc_run, hyb_run, show_ylabel=True)
    _plot_bag_summary(axes[3, 1], hyb_run, msrc_run, show_ylabel=False)
    _plot_simplex(axes[4, 0], msrc_run, title="Window q-vectors against MSC arms")
    _plot_simplex(axes[4, 1], hyb_run, title="Window q-vectors against MSC arms")
    for ax in axes[:3, :].flat:
        ax.set_xlim(0, xmax)
    axes[2, 0].set_xlabel("Genomic position (bp)")
    axes[2, 1].set_xlabel("Genomic position (bp)")
    for i, ax in enumerate(axes.flat):
        ax.text(-0.04, 1.04, chr(ord("A") + i), transform=ax.transAxes, ha="right", va="bottom", fontweight="bold")
    q_text = "Same/matched genome-wide quartet frequencies, different spatial organization" if marginal_l1_distance(msrc_run.marginal_q, hyb_run.marginal_q) <= 0.05 else "Genome-wide quartet frequencies differ; inspect spatial profiles"
    params = _hybridization_parameter_text(hyb_run)
    suffix = f" ({params})" if params else ""
    fig.suptitle(title or f"{q_text}{suffix}", y=0.995)
    return fig


def _hybridization_parameter_text(run: SpatialModelRun) -> str:
    summary = run.summary
    pieces = []
    if "gamma" in summary:
        pieces.append(f"gamma={float(summary['gamma']):.3g}")
    h = summary.get("generations_since_pulse", summary.get("generations_since_hybridization"))
    if h is not None:
        pieces.append(f"h={float(h):.3g}")
    r = summary.get("recombination_rate_per_bp_per_generation", summary.get("recombination_rate"))
    if r is not None:
        pieces.append(f"r={float(r):.3g}")
    mean_len = summary.get("expected_mean_introgressed_tract_length_bp")
    if mean_len is None and h is not None and r is not None and "gamma" in summary and float(h) > 0 and float(r) > 0 and float(summary["gamma"]) < 1.0:
        mean_len = 1.0 / (float(h) * float(r) * max(1e-15, 1.0 - float(summary["gamma"])))
    if mean_len is not None:
        pieces.append(f"E[L_intro]={float(mean_len):.3g} bp")
    return ", ".join(pieces)


def make_spatial_overlay_figure(msrc_run: SpatialModelRun, hyb_run: SpatialModelRun, *, focal_topology: int | None = None) -> Any:
    configure_matplotlib_cache()
    import matplotlib.pyplot as plt

    _setup_common_style(plt)
    xmax = _xmax(msrc_run, hyb_run)
    if focal_topology is None:
        diffs = [abs(float(msrc_run.marginal_q[i]) - float(hyb_run.marginal_q[i])) for i in range(3)]
        focal_topology = int(np.argmax(diffs))
    focal_topology = min(max(int(focal_topology), 0), 2)
    qkey = f"q{focal_topology + 1}"
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 5.4), sharex=True)
    axes[0].plot(_centers(msrc_run), _values(msrc_run, "fraction_rearranged"), color=MSRC_FEATURE, linewidth=1.7, label="MSRC rearranged fraction")
    axes[0].plot(_centers(hyb_run), _values(hyb_run, "introgressed_fraction"), color=HYB_INTRO, linewidth=1.7, label="HYB introgressed fraction")
    axes[0].set_ylim(0, 1.02)
    axes[0].set_ylabel("Local feature fraction")
    axes[0].legend(frameon=False, loc="upper right")
    axes[1].plot(_centers(msrc_run), _values(msrc_run, qkey), color=MSRC_FEATURE, linewidth=1.7, label=f"MSRC {qkey}")
    axes[1].plot(_centers(hyb_run), _values(hyb_run, qkey), color=HYB_INTRO, linewidth=1.7, label=f"HYB {qkey}")
    axes[1].axhline(1 / 3, color="0.45", linestyle="--", linewidth=0.85)
    axes[1].set_ylim(0, 1.02)
    axes[1].set_xlim(0, xmax)
    axes[1].set_ylabel(f"Topology {focal_topology + 1} support")
    axes[1].set_xlabel("Genomic position (bp)")
    axes[1].legend(frameon=False, loc="upper right")
    fig.suptitle("Direct overlay of local spatial signal", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig
