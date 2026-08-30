from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from msrcsim.spatial_compare_plot import make_spatial_compare_figure, make_spatial_overlay_figure
from msrcsim.spatial_compare_io import load_spatial_model_run
from spatial_compare_helpers import hyb_dir, msrc_dir


def test_compare_figure_generation_headless(tmp_path):
    msrc = load_spatial_model_run(msrc_dir(tmp_path), "msrc")
    hyb = load_spatial_model_run(hyb_dir(tmp_path), "hybridization")
    fig = make_spatial_compare_figure(msrc, hyb)
    path = tmp_path / "compare.png"
    fig.savefig(path)
    assert path.exists()
    assert path.stat().st_size > 0


def test_overlay_figure_generation_headless_with_different_window_counts(tmp_path):
    msrc = load_spatial_model_run(msrc_dir(tmp_path), "msrc")
    hyb = load_spatial_model_run(hyb_dir(tmp_path), "hybridization")
    hyb.windows.append({**hyb.windows[-1], "window_id": 99, "center_bp": 950, "start_bp": 900, "end_bp": 1000})
    fig = make_spatial_overlay_figure(msrc, hyb, focal_topology=1)
    path = tmp_path / "overlay.svg"
    fig.savefig(path)
    assert path.exists()
    assert path.read_text().startswith("<?xml") or "<svg" in path.read_text()[:500]
