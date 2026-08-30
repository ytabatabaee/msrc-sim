import csv

import pytest

from msrcsim.frequency_history_plot import (
    choose_display_records,
    compute_tree_layout,
    glyph_counts,
    infer_population_tree,
    plot_frequency_history_tree,
    read_frequency_history,
)


FIELDNAMES = [
    "rearrangement_id",
    "branch_id",
    "parent_branch_id",
    "forward_generation",
    "absolute_age",
    "copy_count_A1",
    "copy_count_A0",
    "population_chromosomes",
    "frequency_A1",
    "frequency_A0",
    "status",
    "is_origin",
    "is_branch_start",
    "is_branch_end",
    "selection_coefficient",
    "effective_population_size",
]


def _record(branch, parent, generation, age, p, status="segregating", origin=False, start=False, end=False):
    n = 20
    a1 = round(n * p)
    return {
        "rearrangement_id": "inv_1",
        "branch_id": branch,
        "parent_branch_id": parent or "",
        "forward_generation": generation,
        "absolute_age": age,
        "copy_count_A1": a1,
        "copy_count_A0": n - a1,
        "population_chromosomes": n,
        "frequency_A1": p,
        "frequency_A0": 1 - p,
        "status": status,
        "is_origin": origin,
        "is_branch_start": start,
        "is_branch_end": end,
        "selection_coefficient": 0.0,
        "effective_population_size": 10,
    }


def _write_history(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
    except Exception as exc:
        pytest.skip(f"Matplotlib is not importable in this environment: {exc}")


def _balanced_rows():
    rows = []
    rows += [_record("ROOT", None, 0, 150, 0.2, "newly_originated", True, True), _record("ROOT", None, 1, 100, 0.5, end=True)]
    rows += [_record("A", "ROOT", 0, 100, 0.5, start=True), _record("A", "ROOT", 1, 50, 1.0, "fixed", end=True)]
    rows += [_record("B", "ROOT", 0, 100, 0.5, start=True), _record("B", "ROOT", 1, 50, 0.0, "lost", end=True)]
    for tip, parent, p in (("1", "A", 1.0), ("2", "A", 1.0), ("3", "B", 0.0), ("4", "B", 0.5)):
        rows += [_record(tip, parent, 0, 50, p, start=True), _record(tip, parent, 1, 0, p, end=True)]
    return rows


def test_parse_minimal_frequency_history(tmp_path):
    path = tmp_path / "frequency_history.csv"
    _write_history(path, [_record("ROOT", None, 0, 10, 0.25, origin=True, start=True)])
    rows = read_frequency_history(path)
    assert rows[0]["branch_id"] == "ROOT"
    assert rows[0]["parent_branch_id"] is None
    assert rows[0]["absolute_age"] == 10.0
    assert rows[0]["is_origin"] is True


def test_infer_parent_child_branch_structure():
    tree = infer_population_tree(_balanced_rows())
    assert tree["roots"] == ("ROOT",)
    assert tree["branches"]["ROOT"].children == ("A", "B")
    assert tree["branches"]["A"].children == ("1", "2")
    assert tree["tips"] == ("1", "2", "3", "4")


def test_balanced_tree_layout():
    layout = compute_tree_layout(_balanced_rows(), tip_order=["1", "2", "3", "4"])
    x = layout["x"]
    assert x["A"] == pytest.approx((x["1"] + x["2"]) / 2)
    assert x["B"] == pytest.approx((x["3"] + x["4"]) / 2)
    assert x["ROOT"] == pytest.approx((x["1"] + x["2"] + x["3"] + x["4"]) / 4)


def test_unbalanced_tree_layout():
    rows = []
    rows += [_record("ROOT", None, 0, 180, 0.2, origin=True, start=True), _record("ROOT", None, 1, 120, 0.2, end=True)]
    rows += [_record("1", "ROOT", 0, 120, 0.2, start=True), _record("1", "ROOT", 1, 0, 0.2, end=True)]
    rows += [_record("C", "ROOT", 0, 120, 0.2, start=True), _record("C", "ROOT", 1, 80, 0.2, end=True)]
    rows += [_record("2", "C", 0, 80, 0.2, start=True), _record("2", "C", 1, 0, 0.2, end=True)]
    rows += [_record("D", "C", 0, 80, 0.2, start=True), _record("D", "C", 1, 40, 0.2, end=True)]
    rows += [_record("3", "D", 0, 40, 0.2, start=True), _record("3", "D", 1, 0, 0.2, end=True)]
    rows += [_record("4", "D", 0, 40, 0.2, start=True), _record("4", "D", 1, 0, 0.2, end=True)]
    layout = compute_tree_layout(rows, tip_order=["1", "2", "3", "4"])
    x = layout["x"]
    assert x["D"] == pytest.approx((x["3"] + x["4"]) / 2)
    assert x["C"] == pytest.approx((x["2"] + x["3"] + x["4"]) / 3)
    assert x["ROOT"] == pytest.approx((x["1"] + x["2"] + x["3"] + x["4"]) / 4)


def test_deterministic_arrow_counts():
    assert glyph_counts(0.0, 12) == (12, 0)
    assert glyph_counts(0.5, 12) == (6, 6)
    assert glyph_counts(1.0, 12) == (0, 12)


def test_downsampling_preserves_events_and_transitions():
    rows = []
    statuses = ["not_present"] * 5 + ["newly_originated"] + ["segregating"] * 10 + ["fixed"] * 10 + ["lost"] * 10
    for i, status in enumerate(statuses):
        rows.append(
            _record(
                "ROOT",
                None,
                i,
                100 - i,
                0.0 if status in {"not_present", "lost"} else 1.0 if status == "fixed" else 0.5,
                status=status,
                origin=status == "newly_originated",
                start=i == 0,
                end=i == len(statuses) - 1,
            )
        )
    selected = choose_display_records(rows, max_rows_per_branch=12)
    assert any(r["is_origin"] for r in selected)
    assert any(r["is_branch_start"] for r in selected)
    assert any(r["is_branch_end"] for r in selected)
    assert {"not_present", "newly_originated", "segregating", "fixed", "lost"} <= {r["status"] for r in selected}


def test_plot_generation_headless(tmp_path):
    _require_matplotlib()
    fig, ax = plot_frequency_history_tree(_balanced_rows(), sampled_arrangements={"1": {"arrangement": "0"}}, title="History")
    out = tmp_path / "history.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0
    assert ax.get_title() == "History"
