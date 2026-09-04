from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")

from msrcsim.spatial_compare_cli import main
from spatial_compare_helpers import hyb_dir, msrc_dir


def test_cli_writes_figure_and_summary_json(tmp_path):
    msrc = msrc_dir(tmp_path)
    hyb = hyb_dir(tmp_path)
    output = tmp_path / "comparison" / "spatial_compare.png"
    code = main([
        "--msrc-dir", str(msrc),
        "--hyb-dir", str(hyb),
        "--output", str(output),
        "--formats", "png,svg",
        "--no-overlay",
    ])
    assert code == 0
    assert output.exists()
    assert output.with_suffix(".svg").exists()
    summary = json.loads((output.parent / "spatial_compare_summary.json").read_text())
    assert summary["msrc_num_windows"] == 2
    assert summary["hyb_num_loci"] == 12
    assert summary["has_overlay"] is False
    assert summary["marginal_q_l1_difference"] >= 0


def test_cli_writes_overlay_and_is_deterministic_for_summary(tmp_path):
    msrc = msrc_dir(tmp_path)
    hyb = hyb_dir(tmp_path)
    output = tmp_path / "spatial_compare.pdf"
    args = [
        "--msrc-dir", str(msrc),
        "--hyb-dir", str(hyb),
        "--output", str(output),
        "--formats", "pdf",
    ]
    assert main(args) == 0
    summary1 = json.loads((tmp_path / "spatial_compare_summary.json").read_text())
    assert (tmp_path / "spatial_compare_overlay.pdf").exists()
    assert main(args) == 0
    summary2 = json.loads((tmp_path / "spatial_compare_summary.json").read_text())
    assert summary1 == summary2


def test_cli_returns_error_for_missing_hybridization_files(tmp_path, capsys):
    msrc = msrc_dir(tmp_path)
    empty = tmp_path / "empty_hyb"
    empty.mkdir()
    output = tmp_path / "spatial_compare.png"
    code = main(["--msrc-dir", str(msrc), "--hyb-dir", str(empty), "--output", str(output)])
    captured = capsys.readouterr()
    assert code == 2
    assert "Hybridization outputs not found" in captured.err


def test_cli_uses_matched_hybridization_json(tmp_path):
    msrc = msrc_dir(tmp_path)
    hyb = hyb_dir(tmp_path)
    matched = tmp_path / "matched.json"
    matched.write_text(json.dumps({"best_hyb_dir": str(hyb)}))
    output = tmp_path / "from_match" / "spatial_compare.png"
    code = main([
        "--msrc-dir", str(msrc),
        "--matched-hybridization-json", str(matched),
        "--output", str(output),
        "--no-overlay",
    ])
    assert code == 0
    summary = json.loads((output.parent / "spatial_compare_summary.json").read_text())
    assert summary["hyb_dir"] == str(hyb)


def test_cli_can_recompute_fixed_bp_windows(tmp_path):
    msrc = msrc_dir(tmp_path)
    hyb = hyb_dir(tmp_path)
    output = tmp_path / "fixed_bp" / "spatial_compare.png"
    code = main([
        "--msrc-dir", str(msrc),
        "--hyb-dir", str(hyb),
        "--output", str(output),
        "--window-size-bp", "300",
        "--step-bp", "300",
        "--no-overlay",
    ])
    assert code == 0
    summary = json.loads((output.parent / "spatial_compare_summary.json").read_text())
    assert summary["windowing"]["mode"] == "fixed_bp"
    assert summary["windowing"]["window_size_bp"] == 300.0
    assert summary["msrc_num_windows"] == 3
    assert summary["hyb_num_windows"] == 4
