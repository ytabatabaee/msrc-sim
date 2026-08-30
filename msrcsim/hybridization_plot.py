from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
import json
import os
import tempfile


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _f(row: dict[str, Any], key: str) -> float:
    return float(row[key])


def plot_hybridization_output(input_dir: str | Path, output_dir: str | Path | None = None, fmt: str = "png") -> list[Path]:
    inp = Path(input_dir)
    out = Path(output_dir) if output_dir else inp
    out.mkdir(parents=True, exist_ok=True)
    tracts = _read_csv(inp / "hybridization_tracts.csv")
    windows = _read_csv(inp / "hybridization_windows.csv")
    with (inp / "hybridization_summary.json").open() as handle:
        summary = json.load(handle)

    if "MPLCONFIGDIR" not in os.environ:
        cache = Path(tempfile.gettempdir()) / "msrcsim-matplotlib"
        cache.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(cache)
    if "XDG_CACHE_HOME" not in os.environ:
        cache_home = Path(tempfile.gettempdir()) / "msrcsim-cache"
        cache_home.mkdir(parents=True, exist_ok=True)
        os.environ["XDG_CACHE_HOME"] = str(cache_home)
    import matplotlib.pyplot as plt

    colors = {"major": "#4c78a8", "introgressed": "#f58518"}
    quartet_colors = ["#1b9e77", "#d95f02", "#7570b3"]
    length = float(summary["chromosome_length_bp"])
    x = [_f(r, "center_bp") for r in windows]

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(10.5, 6.8),
        sharex=True,
        gridspec_kw={"height_ratios": [0.9, 2.2, 1.2]},
    )

    ax = axes[0]
    for tract in tracts:
        start = _f(tract, "start_bp")
        end = _f(tract, "end_bp")
        label = tract["ancestry_label"]
        ax.broken_barh([(start, end - start)], (0.25, 0.5), facecolors=colors[label], edgecolors="none")
        if start > 0.0:
            ax.axvline(start, color="0.25", linewidth=0.45, alpha=0.45)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.5], ["ancestry"])
    ax.set_xlim(0, length)
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.tick_params(axis="y", length=0)

    ax = axes[1]
    for key, label, color in zip(("q1", "q2", "q3"), ("12|34", "13|24", "14|23"), quartet_colors):
        ax.plot(x, [_f(r, key) for r in windows], color=color, linewidth=1.8, label=label)
    ax.axhline(1 / 3, color="0.45", linestyle="--", linewidth=1.0)
    ax.set_ylabel("Moving-window quartet support")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, ncols=3, loc="upper right")

    ax = axes[2]
    ax.plot(x, [_f(r, "introgressed_fraction") for r in windows], color=colors["introgressed"], linewidth=1.8)
    ax.axhline(float(summary["gamma"]), color="0.45", linestyle="--", linewidth=1.0)
    ax.set_ylabel("Introgressed fraction")
    ax.set_xlabel("Genomic position (bp)")
    ax.set_ylim(0, 1.02)

    fig.suptitle("Pulse-hybridization ancestry mosaic and local quartet support", y=0.99)
    fig.tight_layout()
    profile = out / f"hybridization_spatial_profile.{fmt}"
    fig.savefig(profile, dpi=250)
    plt.close(fig)

    manifest = {"input": str(inp), "figures": [profile.name], "renderer": "matplotlib"}
    (out / "hybridization_figure_manifest.json").write_text(json.dumps(manifest, indent=2))
    return [profile]
