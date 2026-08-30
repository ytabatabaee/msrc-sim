from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
import json
import importlib.metadata


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _f(row: dict[str, Any], key: str) -> float:
    return float(row[key])


def _matplotlib_version_supported() -> bool:
    try:
        version = importlib.metadata.version("matplotlib")
    except importlib.metadata.PackageNotFoundError:
        return False
    parts = []
    for piece in version.split(".")[:2]:
        try:
            parts.append(int(piece))
        except ValueError:
            parts.append(0)
    return tuple(parts) >= (3, 7)


def _plot_with_matplotlib(inp: Path, out: Path, fmt: str, windows: list[dict[str, str]], loci: list[dict[str, str]], summary: dict[str, Any]) -> list[Path]:
    import matplotlib.pyplot as plt

    interval = summary["rearrangement_interval"]
    start = float(interval["start"])
    end = float(interval["end"])
    x = [_f(r, "center_position") for r in windows]
    figures: list[Path] = []

    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.axvspan(start, end, color="#d8b365", alpha=0.25, label="rearrangement interval")
    ax.axvline(start, color="#8c510a", linewidth=1.2)
    ax.axvline(end, color="#8c510a", linewidth=1.2)
    colors = ["#1b9e77", "#d95f02", "#7570b3"]
    for key, label, color in zip(("q1", "q2", "q3"), ("12|34", "13|24", "14|23"), colors):
        ax.plot(x, [_f(r, key) for r in windows], color=color, linewidth=1.8, label=label)
    if loci:
        y0 = -0.055
        for i, label in enumerate(("12|34", "13|24", "14|23")):
            xs = [float(r["position"]) for r in loci if r["topology"] == label]
            if xs:
                ax.scatter(xs, [y0 - i * 0.025] * len(xs), s=4, color=colors[i], alpha=0.35, linewidths=0)
    ax.axhline(1 / 3, color="0.45", linestyle="--", linewidth=1, label="1/3")
    ax.set_xlabel("Genomic position")
    ax.set_ylabel("Local quartet support")
    ax.set_ylim(-0.14 if loci else 0.0, 1.02)
    ax.legend(loc="upper right", ncols=2, frameon=False)
    fig.tight_layout()
    profile = out / f"spatial_quartet_profile.{fmt}"
    fig.savefig(profile, dpi=250)
    plt.close(fig)
    figures.append(profile)

    labels = ["overall", "inside", "outside"]
    values = [
        [summary["overall_q1"], summary["overall_q2"], summary["overall_q3"]],
        [summary["inside_q1"], summary["inside_q2"], summary["inside_q3"]],
        [summary["outside_q1"], summary["outside_q2"], summary["outside_q3"]],
    ]
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    width = 0.24
    offsets = [-width, 0, width]
    for i, topo in enumerate(("12|34", "13|24", "14|23")):
        ax.bar([j + offsets[i] for j in range(len(labels))], [v[i] for v in values], width=width, color=colors[i], label=topo)
    ax.axhline(1 / 3, color="0.45", linestyle="--", linewidth=1)
    ax.set_xticks(range(len(labels)), labels)
    ax.set_ylabel("Quartet support")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False)
    fig.tight_layout()
    contrast = out / f"spatial_bag_inside_outside.{fmt}"
    fig.savefig(contrast, dpi=250)
    plt.close(fig)
    figures.append(contrast)

    manifest = {"input": str(inp), "figures": [p.name for p in figures], "renderer": "matplotlib"}
    (out / "spatial_figure_manifest.json").write_text(json.dumps(manifest, indent=2))
    return figures


def _plot_with_pillow(inp: Path, out: Path, windows: list[dict[str, str]], summary: dict[str, Any]) -> list[Path]:
    from PIL import Image, ImageDraw, ImageFont

    width, height = 1200, 580
    left, right, top, bottom = 90, 40, 45, 95
    colors = [(27, 158, 119), (217, 95, 2), (117, 112, 179)]
    xmax = max(float(r["end_position"]) for r in windows)

    def sx(pos: float) -> int:
        return int(left + (pos / xmax) * (width - left - right))

    def sy(q: float) -> int:
        return int(top + (1.0 - q) * (height - top - bottom))

    interval = summary["rearrangement_interval"]
    figures: list[Path] = []
    font = ImageFont.load_default()
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.rectangle((left, top, width - right, height - bottom), outline=(30, 30, 30))
    draw.rectangle((sx(float(interval["start"])), top, sx(float(interval["end"])), height - bottom), fill=(239, 222, 185))
    for bp in (float(interval["start"]), float(interval["end"])):
        x = sx(bp)
        draw.line((x, top, x, height - bottom), fill=(140, 81, 10), width=2)
    draw.line((left, sy(1 / 3), width - right, sy(1 / 3)), fill=(120, 120, 120), width=1)
    xs = [sx(float(r["center_position"])) for r in windows]
    for key, color in zip(("q1", "q2", "q3"), colors):
        pts = list(zip(xs, [sy(float(r[key])) for r in windows]))
        if len(pts) > 1:
            draw.line(pts, fill=color, width=3)
    draw.text((left, 16), "Spatial quartet profile", fill=(0, 0, 0), font=font)
    draw.text((left, height - 40), "Genomic position", fill=(0, 0, 0), font=font)
    for i, (label, color) in enumerate(zip(("12|34", "13|24", "14|23"), colors)):
        x0 = width - 210
        y0 = 22 + 18 * i
        draw.line((x0, y0 + 5, x0 + 30, y0 + 5), fill=color, width=3)
        draw.text((x0 + 38, y0), label, fill=(0, 0, 0), font=font)
    profile = out / "spatial_quartet_profile.png"
    img.save(profile)
    figures.append(profile)

    img = Image.new("RGB", (780, 460), "white")
    draw = ImageDraw.Draw(img)
    labels = ["overall", "inside", "outside"]
    values = [
        [summary["overall_q1"], summary["overall_q2"], summary["overall_q3"]],
        [summary["inside_q1"], summary["inside_q2"], summary["inside_q3"]],
        [summary["outside_q1"], summary["outside_q2"], summary["outside_q3"]],
    ]
    base_y = 370
    scale = 300
    draw.line((70, base_y, 730, base_y), fill=(30, 30, 30))
    for group, vals in enumerate(values):
        gx = 150 + group * 220
        for i, val in enumerate(vals):
            h = int(float(val) * scale)
            x0 = gx + i * 42
            draw.rectangle((x0, base_y - h, x0 + 32, base_y), fill=colors[i])
        draw.text((gx, base_y + 18), labels[group], fill=(0, 0, 0), font=font)
    draw.text((70, 25), "Overall vs inside/outside quartet support", fill=(0, 0, 0), font=font)
    contrast = out / "spatial_bag_inside_outside.png"
    img.save(contrast)
    figures.append(contrast)

    manifest = {"input": str(inp), "figures": [p.name for p in figures], "renderer": "pillow_fallback"}
    (out / "spatial_figure_manifest.json").write_text(json.dumps(manifest, indent=2))
    return figures


def plot_spatial_output(input_dir: str | Path, output_dir: str | Path | None = None, fmt: str = "png") -> list[Path]:
    inp = Path(input_dir)
    out = Path(output_dir) if output_dir else inp
    out.mkdir(parents=True, exist_ok=True)
    windows = _read_csv(inp / "spatial_windows.csv")
    loci_path = inp / "spatial_loci.csv"
    loci = _read_csv(loci_path) if loci_path.exists() else []
    with (inp / "spatial_summary.json").open() as handle:
        summary = json.load(handle)
    if fmt == "png" and not _matplotlib_version_supported():
        return _plot_with_pillow(inp, out, windows, summary)
    try:
        return _plot_with_matplotlib(inp, out, fmt, windows, loci, summary)
    except Exception:
        if fmt != "png":
            raise
        return _plot_with_pillow(inp, out, windows, summary)
