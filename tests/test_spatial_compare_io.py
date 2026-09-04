from __future__ import annotations

import csv
import json

import pytest

from msrcsim.spatial_compare_io import compute_windows_from_loci, compute_windows_from_loci_bp, load_spatial_model_run


def _write_csv(path, rows, fieldnames=None):
    if fieldnames is None:
        fieldnames = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _msrc_dir(tmp_path, *, windows=True):
    out = tmp_path / "msrc"
    out.mkdir()
    loci = [
        {"locus_id": i, "position": i * 100, "inside_rearrangement": 2 <= i <= 5, "topology_index": 1 if 2 <= i <= 5 else 0, "topology": "13|24" if 2 <= i <= 5 else "12|34"}
        for i in range(10)
    ]
    _write_csv(out / "spatial_loci.csv", loci)
    if windows:
        _write_csv(out / "spatial_windows.csv", [
            {"window_id": 0, "start_position": 0, "end_position": 300, "center_position": 150, "n_loci": 4, "q1": 0.5, "q2": 0.5, "q3": 0.0, "dominant_topology": 0, "distance_to_nearest_msc_arm": 0.1, "off_arm_difference": 0.5},
            {"window_id": 1, "start_position": 400, "end_position": 700, "center_position": 550, "n_loci": 4, "q1": 0.5, "q2": 0.5, "q3": 0.0, "dominant_topology": 0, "distance_to_nearest_msc_arm": 0.1, "off_arm_difference": 0.5},
        ])
    (out / "spatial_summary.json").write_text(json.dumps({
        "genome_length": 1000,
        "overall_q1": 0.6,
        "overall_q2": 0.4,
        "overall_q3": 0.0,
        "rearrangement_interval": {"start": 200, "end": 500, "interval_id": "inv1"},
        "topology_names": ["12|34", "13|24", "14|23"],
    }))
    return out


def _hyb_dir(tmp_path, *, windows=True):
    out = tmp_path / "hyb"
    out.mkdir()
    loci = [
        {"locus_id": i, "position_bp": i * 100, "ancestry_state": 1 if i in {1, 2, 6, 7} else 0, "ancestry_label": "introgressed" if i in {1, 2, 6, 7} else "major", "tract_id": i // 2, "topology_index": 1 if i in {1, 2, 6, 7} else 0, "topology_label": "13|24" if i in {1, 2, 6, 7} else "12|34", "model": "pulse_hybridization"}
        for i in range(12)
    ]
    _write_csv(out / "hybridization_loci.csv", loci)
    _write_csv(out / "hybridization_tracts.csv", [
        {"tract_id": 0, "start_bp": 0, "end_bp": 100, "length_bp": 100, "ancestry_state": 0, "ancestry_label": "major"},
        {"tract_id": 1, "start_bp": 100, "end_bp": 300, "length_bp": 200, "ancestry_state": 1, "ancestry_label": "introgressed"},
        {"tract_id": 2, "start_bp": 300, "end_bp": 600, "length_bp": 300, "ancestry_state": 0, "ancestry_label": "major"},
        {"tract_id": 3, "start_bp": 600, "end_bp": 800, "length_bp": 200, "ancestry_state": 1, "ancestry_label": "introgressed"},
        {"tract_id": 4, "start_bp": 800, "end_bp": 1200, "length_bp": 400, "ancestry_state": 0, "ancestry_label": "major"},
    ])
    if windows:
        _write_csv(out / "hybridization_windows.csv", [
            {"window_id": 0, "start_bp": 0, "end_bp": 300, "center_bp": 150, "num_loci": 4, "introgressed_fraction": 0.5, "q1": 0.5, "q2": 0.5, "q3": 0.0, "dominant_topology": 0, "distance_to_nearest_msc_arm": 0.1, "off_arm_difference": 0.5},
            {"window_id": 1, "start_bp": 400, "end_bp": 900, "center_bp": 650, "num_loci": 6, "introgressed_fraction": 0.3333333333, "q1": 0.6666666667, "q2": 0.3333333333, "q3": 0.0, "dominant_topology": 0, "distance_to_nearest_msc_arm": 0.05, "off_arm_difference": 0.3333333333},
        ])
    (out / "hybridization_summary.json").write_text(json.dumps({
        "chromosome_length_bp": 1200,
        "observed_marginal_q": [2 / 3, 1 / 3, 0.0],
        "expected_marginal_q": [0.65, 0.35, 0.0],
    }))
    return out


def test_load_normalized_msrc_run(tmp_path):
    run = load_spatial_model_run(_msrc_dir(tmp_path), "msrc")
    assert run.model_name == "MSRC"
    assert run.chromosome_length_bp == 1000
    assert run.marginal_q == (0.6, 0.4, 0.0)
    assert run.feature_intervals[0]["start_bp"] == 200
    assert "fraction_rearranged" in run.windows[0]


def test_load_normalized_hybridization_run(tmp_path):
    run = load_spatial_model_run(_hyb_dir(tmp_path), "hybridization")
    assert run.model_name == "Hybridization"
    assert run.chromosome_length_bp == 1200
    assert run.marginal_q == (2 / 3, 1 / 3, 0.0)
    assert len(run.feature_intervals) == 5
    assert run.windows[0]["introgressed_fraction"] == 0.5


def test_fallback_window_computation_from_loci(tmp_path):
    run = load_spatial_model_run(_msrc_dir(tmp_path, windows=False), "msrc", window_size_loci=4, step_loci=2)
    assert len(run.windows) == 4
    assert run.windows[0]["fraction_rearranged"] == 0.5
    assert run.windows[0]["q1"] == 0.5
    assert run.windows[0]["q2"] == 0.5


def test_fixed_bp_window_computation_from_loci():
    loci = [
        {"position_bp": i * 100, "topology_index": 1 if i in {2, 3, 4} else 0, "is_inside_rearranged_interval": 2 <= i <= 4}
        for i in range(10)
    ]
    windows = compute_windows_from_loci_bp(loci, "msrc", window_size_bp=300, step_bp=300, chromosome_length_bp=900)
    assert len(windows) == 3
    assert [row["num_loci"] for row in windows] == [3, 3, 4]
    assert windows[0]["q1"] == pytest.approx(2 / 3)
    assert windows[1]["q2"] == pytest.approx(2 / 3)
    assert windows[1]["fraction_rearranged"] == pytest.approx(2 / 3)


def test_fixed_bp_window_grid_includes_last_exact_start():
    loci = [{"position_bp": i * 1_000_000, "topology_index": 0, "is_inside_rearranged_interval": False} for i in range(100)]
    windows = compute_windows_from_loci_bp(loci, "msrc", window_size_bp=10_000_000, step_bp=2_000_000, chromosome_length_bp=100_000_000)
    assert len(windows) == 46
    assert windows[-1]["start_bp"] == 90_000_000
    assert windows[-1]["end_bp"] == 100_000_000


def test_load_recomputes_fixed_bp_windows_even_when_saved_windows_exist(tmp_path):
    run = load_spatial_model_run(_msrc_dir(tmp_path), "msrc", window_size_bp=300, step_bp=300)
    assert len(run.windows) == 3
    assert run.windows[0]["start_bp"] == 0
    assert run.windows[0]["end_bp"] == 300
    assert run.windows[0]["num_loci"] == 3


def test_marginal_q_from_loci_when_summary_missing_q(tmp_path):
    out = _hyb_dir(tmp_path, windows=False)
    (out / "hybridization_summary.json").write_text(json.dumps({"chromosome_length_bp": 1200}))
    run = load_spatial_model_run(out, "hybridization", window_size_loci=4, step_loci=4)
    assert run.marginal_q == pytest.approx((2 / 3, 1 / 3, 0.0))


def test_compute_windows_validates_enough_loci():
    with pytest.raises(ValueError, match="not enough loci"):
        compute_windows_from_loci([{"position_bp": 0, "topology_index": 0}], "msrc", 2, 1)


def test_informative_error_when_required_files_missing(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError, match="Spatial MSRC outputs not found"):
        load_spatial_model_run(empty, "msrc")


def test_windows_only_run_loads_when_loci_absent(tmp_path):
    out = _msrc_dir(tmp_path)
    (out / "spatial_loci.csv").unlink()
    run = load_spatial_model_run(out, "msrc")
    assert run.loci == []
    assert len(run.windows) == 2
    assert run.marginal_q == (0.6, 0.4, 0.0)
