from copy import deepcopy

from msrcsim.robustness import contribution_weights, infer_with_strategy


def _constructed_rows():
    rows = []
    for block in range(5):
        rows.append({"block_id": block, "topology_index": 0, "topology": "12|34", "is_rearranged": False, "msrc_probability": 0.0})
    for window in range(20):
        rows.append({"block_id": 5, "topology_index": 1, "topology": "13|24", "is_rearranged": True, "msrc_probability": 1.0})
    return rows


def test_block_collapse_gives_each_block_total_weight_one():
    rows = _constructed_rows()
    weights = contribution_weights(rows, "block_collapse")
    totals = {}
    for row, weight in zip(rows, weights):
        totals[row["block_id"]] = totals.get(row["block_id"], 0.0) + weight
    assert set(round(v, 12) for v in totals.values()) == {1.0}


def test_filtering_and_weighting_never_rewrite_topology():
    rows = _constructed_rows()
    original = deepcopy(rows)
    infer_with_strategy(rows, "oracle_filter")
    infer_with_strategy(rows, "soft_weight")
    infer_with_strategy(rows, "block_collapse")
    assert rows == original


def test_constructed_example_naive_flips_but_corrections_recover_t1():
    rows = _constructed_rows()
    assert infer_with_strategy(rows, "all_windows").inferred_topology_index == 1
    assert infer_with_strategy(rows, "oracle_filter").inferred_topology_index == 0
    assert infer_with_strategy(rows, "block_collapse").inferred_topology_index == 0
    assert infer_with_strategy(rows, "soft_weight").inferred_topology_index == 0
