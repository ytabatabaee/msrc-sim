from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
import csv
import math


A0_COLOR = "#2b6cb0"
A1_COLOR = "#d9483b"
TUBE_COLOR = "#b8c7d9"
ORIGIN_COLOR = "#f2b705"


@dataclass(frozen=True)
class InferredBranch:
    branch_id: str
    parent_branch_id: str | None
    children: tuple[str, ...]
    older_age: float
    younger_age: float
    effective_population_size: float | None


def _bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _optional_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def _natural_key(value: str) -> tuple[Any, ...]:
    parts: list[Any] = []
    buf = ""
    mode_digit: bool | None = None
    for ch in str(value):
        digit = ch.isdigit()
        if mode_digit is None or digit == mode_digit:
            buf += ch
            mode_digit = digit
            continue
        parts.append((0, int(buf)) if mode_digit else (1, buf))
        buf = ch
        mode_digit = digit
    if buf:
        parts.append((0, int(buf)) if mode_digit else (1, buf))
    return tuple(parts)


def read_frequency_history(path: str | Path) -> list[dict[str, Any]]:
    with Path(path).open(newline="") as handle:
        records = list(csv.DictReader(handle))
    for record in records:
        for key in (
            "forward_generation",
            "absolute_age",
            "copy_count_A1",
            "copy_count_A0",
            "population_chromosomes",
            "frequency_A1",
            "frequency_A0",
            "selection_coefficient",
            "effective_population_size",
        ):
            if key in record and record[key] != "":
                record[key] = float(record[key])
        for key in ("is_origin", "is_branch_start", "is_branch_end"):
            if key in record:
                record[key] = _bool(record[key])
        if record.get("parent_branch_id") == "":
            record["parent_branch_id"] = None
    return records


def read_sampled_arrangements(path: str | Path | None) -> dict[str, dict[str, str]]:
    if path is None or not Path(path).exists():
        return {}
    with Path(path).open(newline="") as handle:
        return {str(row["taxon"]): row for row in csv.DictReader(handle)}


def infer_population_tree(records: Iterable[dict[str, Any]]) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    parent_by_branch: dict[str, str | None] = {}
    children: dict[str, set[str]] = defaultdict(set)
    for record in records:
        branch_id = str(record["branch_id"])
        parent_id = record.get("parent_branch_id")
        parent_id = str(parent_id) if parent_id not in (None, "") else None
        grouped[branch_id].append(record)
        parent_by_branch.setdefault(branch_id, parent_id)
        if parent_id is not None:
            children[parent_id].add(branch_id)

    branches: dict[str, InferredBranch] = {}
    for branch_id, rows in grouped.items():
        ages = [float(r["absolute_age"]) for r in rows]
        nes = [_optional_float(r.get("effective_population_size")) for r in rows]
        ne_values = [x for x in nes if x is not None]
        branches[branch_id] = InferredBranch(
            branch_id=branch_id,
            parent_branch_id=parent_by_branch.get(branch_id),
            children=tuple(sorted(children.get(branch_id, ()), key=_natural_key)),
            older_age=max(ages),
            younger_age=min(ages),
            effective_population_size=ne_values[0] if ne_values else None,
        )
    roots = tuple(sorted((b for b, p in parent_by_branch.items() if p is None), key=_natural_key))
    tips = tuple(sorted((b for b in branches if not branches[b].children), key=_natural_key))
    return {"branches": branches, "roots": roots, "tips": tips, "children": {k: tuple(sorted(v, key=_natural_key)) for k, v in children.items()}}


def choose_display_records(
    records: Iterable[dict[str, Any]],
    max_rows_per_branch: int | None,
    preserve_events: bool = True,
) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[str(record["branch_id"])].append(record)

    selected: list[dict[str, Any]] = []
    for branch_id in sorted(grouped, key=_natural_key):
        rows = sorted(grouped[branch_id], key=lambda r: (float(r["absolute_age"]), float(r.get("forward_generation", 0))), reverse=True)
        if max_rows_per_branch is None or max_rows_per_branch <= 0 or len(rows) <= max_rows_per_branch:
            selected.extend(rows)
            continue
        keep: set[int] = set()
        if preserve_events:
            for i, row in enumerate(rows):
                if row.get("is_branch_start") or row.get("is_branch_end") or row.get("is_origin"):
                    keep.add(i)
                if i and row.get("status") != rows[i - 1].get("status"):
                    keep.add(i - 1)
                    keep.add(i)
        slots = max(max_rows_per_branch - len(keep), 0)
        if slots:
            if slots == 1:
                keep.add(len(rows) // 2)
            else:
                for j in range(slots):
                    keep.add(round(j * (len(rows) - 1) / (slots - 1)))
        if len(keep) > max_rows_per_branch:
            essential = {i for i in keep if rows[i].get("is_branch_start") or rows[i].get("is_branch_end") or rows[i].get("is_origin")}
            transitions = sorted(keep - essential)
            keep = set(sorted(essential))
            for i in transitions:
                if len(keep) >= max_rows_per_branch:
                    break
                keep.add(i)
        selected.extend(rows[i] for i in sorted(keep))
    return selected


def compute_tree_layout(records: Iterable[dict[str, Any]], tip_order: Iterable[str] | None = None) -> dict[str, Any]:
    tree = infer_population_tree(records)
    branches: dict[str, InferredBranch] = tree["branches"]
    tips = list(tip_order) if tip_order is not None else list(tree["tips"])
    tips = [str(t) for t in tips if str(t) in branches]
    for tip in tree["tips"]:
        if tip not in tips:
            tips.append(tip)
    x: dict[str, float] = {tip: float(i) * 2.2 for i, tip in enumerate(tips)}

    def assign_weighted(branch_id: str) -> tuple[float, int]:
        if branch_id in x:
            return x[branch_id], 1
        kids = branches[branch_id].children
        if not kids:
            x[branch_id] = float(len(x)) * 2.2
            return x[branch_id], 1
        child_values = [assign_weighted(child) for child in kids]
        total = sum(count for _value, count in child_values)
        x[branch_id] = sum(value * count for value, count in child_values) / total
        return x[branch_id], total

    for root in tree["roots"]:
        assign_weighted(root)
    return {"x": x, "tips": tuple(tips), "tree": tree}


def glyph_counts(frequency_a1: float, glyphs_per_row: int) -> tuple[int, int]:
    n_a1 = int(round(glyphs_per_row * float(frequency_a1)))
    n_a1 = max(0, min(glyphs_per_row, n_a1))
    return glyphs_per_row - n_a1, n_a1


def _branch_records(records: Iterable[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[str(record["branch_id"])].append(record)
    for rows in grouped.values():
        rows.sort(key=lambda r: float(r["absolute_age"]), reverse=True)
    return grouped


def _draw_arrow(ax: Any, x: float, y: float, direction: int, color: str, scale: float) -> None:
    from matplotlib.patches import FancyArrowPatch

    dx = 0.075 * scale * direction
    start = (x - dx, y)
    end = (x + dx, y)
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=6.5 * scale,
        linewidth=1.0,
        color=color,
        zorder=5,
    )
    ax.add_patch(patch)


def _draw_row(ax: Any, x: float, y: float, frequency_a1: float, glyphs_per_row: int, row_width: float) -> None:
    n_a0, n_a1 = glyph_counts(frequency_a1, glyphs_per_row)
    total = n_a0 + n_a1
    if total == 0:
        return
    step = row_width / max(total - 1, 1)
    x0 = x - row_width / 2
    scale = max(row_width, 0.8)
    for i in range(total):
        xi = x if total == 1 else x0 + i * step
        if i < n_a0:
            _draw_arrow(ax, xi, y, 1, A0_COLOR, scale)
        else:
            _draw_arrow(ax, xi, y, -1, A1_COLOR, scale)


def _tube_width(branch: InferredBranch, width_mode: str, max_ne: float | None) -> float:
    if width_mode == "ne" and branch.effective_population_size and max_ne:
        return 9.0 + 8.0 * math.sqrt(branch.effective_population_size / max_ne)
    return 15.0


def _terminal_frequency_label(record: dict[str, Any]) -> str:
    p = float(record.get("frequency_A1", 0.0))
    return f"p(A1)={p:.2f}"


def plot_frequency_history_tree(
    records: list[dict[str, Any]],
    *,
    sampled_arrangements: dict[str, dict[str, str]] | None = None,
    tip_order: Iterable[str] | None = None,
    glyphs_per_row: int = 12,
    max_rows_per_branch: int | None = 30,
    width_mode: str = "constant",
    show_internal_labels: bool = True,
    show_frequency_trace: bool = False,
    title: str | None = None,
):
    import matplotlib.lines as mlines
    import matplotlib.pyplot as plt

    if not records:
        raise ValueError("frequency history is empty")
    if width_mode not in {"constant", "ne"}:
        raise ValueError("width_mode must be 'constant' or 'ne'")

    selected = choose_display_records(records, max_rows_per_branch)
    layout = compute_tree_layout(records, tip_order)
    x = layout["x"]
    tree = layout["tree"]
    branches: dict[str, InferredBranch] = tree["branches"]
    grouped = _branch_records(records)
    displayed = _branch_records(selected)
    ages = [float(r["absolute_age"]) for r in records]
    age_min, age_max = min(ages), max(ages)
    max_ne = max((b.effective_population_size or 0 for b in branches.values()), default=0) or None
    span = max(age_max - age_min, 1.0)
    row_width = 0.88
    fig_w = max(7.5, 1.45 * max(len(layout["tips"]), 4) + 3.5)
    fig_h = max(6.0, min(12.0, 4.8 + 0.012 * len(selected)))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for branch_id, branch in branches.items():
        bx = x[branch_id]
        lw = _tube_width(branch, width_mode, max_ne)
        if branch.parent_branch_id and branch.parent_branch_id in x:
            px = x[branch.parent_branch_id]
            ax.plot([px, bx], [branch.older_age, branch.older_age], color=TUBE_COLOR, linewidth=lw, alpha=0.42, solid_capstyle="round", zorder=1)
        ax.plot([bx, bx], [branch.older_age, branch.younger_age], color=TUBE_COLOR, linewidth=lw, alpha=0.42, solid_capstyle="round", zorder=1)

    if show_frequency_trace:
        for branch_id, rows in grouped.items():
            bx = x[branch_id]
            trace_x = [bx - 0.52 + float(r.get("frequency_A1", 0.0)) * 1.04 for r in rows]
            trace_y = [float(r["absolute_age"]) for r in rows]
            ax.plot(trace_x, trace_y, color="#6b7280", linewidth=0.8, alpha=0.7, zorder=3)

    for branch_id, rows in displayed.items():
        for row in rows:
            _draw_row(ax, x[branch_id], float(row["absolute_age"]), float(row.get("frequency_A1", 0.0)), glyphs_per_row, row_width)

    origins = [r for r in records if r.get("is_origin")]
    if origins:
        origin = sorted(origins, key=lambda r: float(r["absolute_age"]), reverse=True)[0]
        ox = x[str(origin["branch_id"])]
        oy = float(origin["absolute_age"])
        side = 1 if ox <= (min(x.values()) + max(x.values())) / 2 else -1
        sx = ox + side * 0.78
        ax.scatter([sx], [oy], marker="*", s=180, color=ORIGIN_COLOR, edgecolor="#5f4300", linewidth=0.8, zorder=8)
        ax.plot([sx - side * 0.08, ox + side * 0.35], [oy, oy], color="#5f4300", linewidth=0.8, zorder=7)
        ax.annotate(
            f"Origin of inversion A1\nage = {oy:g} generations\ninitial p = {float(origin.get('frequency_A1', 0.0)):.2f}",
            xy=(sx, oy),
            xytext=(sx + side * 0.42, oy + 0.08 * span),
            ha="left" if side > 0 else "right",
            va="center",
            fontsize=9,
            color="#262626",
            arrowprops={"arrowstyle": "-", "color": "#5f4300", "linewidth": 0.8},
        )

    sampled_arrangements = sampled_arrangements or {}
    terminal_by_branch = {bid: min(rows, key=lambda r: float(r["absolute_age"])) for bid, rows in grouped.items()}
    for tip in layout["tips"]:
        row = terminal_by_branch.get(tip)
        if not row:
            continue
        arrangement = sampled_arrangements.get(tip, {}).get("arrangement")
        label = f"taxon {tip}"
        if arrangement in {"0", "A0", 0}:
            label += "   A0"
        elif arrangement in {"1", "A1", 1}:
            label += "   A1"
        else:
            label += f"\n{_terminal_frequency_label(row)}"
        ax.text(x[tip], age_min - 0.075 * span, label, ha="center", va="top", fontsize=9, color="#222222")

    if show_internal_labels:
        for branch_id, branch in branches.items():
            if branch_id in layout["tips"]:
                continue
            mid_y = (branch.older_age + branch.younger_age) / 2
            label = branch_id
            if branch.effective_population_size is not None:
                label += f"\nNe = {branch.effective_population_size:g}"
            ax.text(
                x[branch_id] + 0.48,
                mid_y,
                label,
                ha="left",
                va="center",
                fontsize=8,
                color="#30363d",
                zorder=6,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.2},
            )

    x_min, x_max = min(x.values()), max(x.values())
    arrow_x = x_min - 1.35
    ax.annotate("", xy=(arrow_x, age_max), xytext=(arrow_x, age_min), arrowprops={"arrowstyle": "->", "color": "#333333", "linewidth": 1.2})
    ax.text(arrow_x, age_max + 0.035 * span, "Past (ancestral)", ha="center", va="bottom", fontsize=9)
    ax.text(arrow_x, age_min - 0.035 * span, "Present", ha="center", va="top", fontsize=9)
    ax.text(arrow_x - 0.18, (age_min + age_max) / 2, "absolute age", rotation=90, ha="center", va="center", fontsize=9, color="#333333")

    a0 = mlines.Line2D([], [], color=A0_COLOR, marker=r"$\rightarrow$", linestyle="None", markersize=11, label="A0 original")
    a1 = mlines.Line2D([], [], color=A1_COLOR, marker=r"$\leftarrow$", linestyle="None", markersize=11, label="A1 rearranged")
    tube_note = "population tube (schematic width)" if width_mode == "constant" else "population tube (sqrt Ne-scaled)"
    tube = mlines.Line2D([], [], color=TUBE_COLOR, linewidth=9, alpha=0.5, label=tube_note)
    ax.legend(handles=[tube, a0, a1], loc="upper right", frameon=False, fontsize=9)
    fig.text(
        0.5,
        0.018,
        "Each displayed arrow represents a fixed fraction of the population; rows visualize the recorded arrangement frequency and do not correspond one-for-one to chromosome copies.",
        ha="center",
        va="bottom",
        fontsize=8,
        color="#4b5563",
    )

    if title:
        ax.set_title(title, fontsize=13, pad=14)
    ax.set_xlim(x_min - 1.85, x_max + 1.85)
    ax.set_ylim(age_min - 0.18 * span, age_max + 0.18 * span)
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelbottom=False)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig._msrc_history_metadata = {"num_records_displayed": len(selected), "tip_order": list(layout["tips"])}  # type: ignore[attr-defined]
    return fig, ax
