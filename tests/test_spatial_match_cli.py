from __future__ import annotations

import json

from msrcsim.spatial_match_cli import main
from spatial_compare_helpers import hyb_dir, msrc_dir


def _set_hyb_q(run_dir, q):
    summary = json.loads((run_dir / "hybridization_summary.json").read_text())
    summary["observed_marginal_q"] = list(q)
    (run_dir / "hybridization_summary.json").write_text(json.dumps(summary))


def test_find_matched_hybridization_picks_nearest_run(tmp_path):
    msrc = msrc_dir(tmp_path)
    grid = tmp_path / "grid"
    grid.mkdir()
    close = grid / "close"
    far = grid / "far"
    close.mkdir()
    far.mkdir()
    close_src = hyb_dir(tmp_path / "close_src")
    far_src = hyb_dir(tmp_path / "far_src")
    for src, dst in ((close_src, close), (far_src, far)):
        for path in src.iterdir():
            (dst / path.name).write_text(path.read_text())
    _set_hyb_q(close, (0.61, 0.39, 0.0))
    _set_hyb_q(far, (0.2, 0.2, 0.6))

    output = tmp_path / "matched.json"
    code = main([
        "--msrc-dir", str(msrc),
        "--hyb-grid-dir", str(grid),
        "--metric", "l1",
        "--output", str(output),
    ])

    assert code == 0
    data = json.loads(output.read_text())
    assert data["best_hyb_dir"] == str(close)
    assert data["best_distance"] < 0.03
    assert data["is_close_match"] is True
    assert data["ranked_candidates"][0]["hyb_dir"] == str(close)


def test_find_matched_hybridization_errors_when_no_candidates(tmp_path):
    output = tmp_path / "matched.json"
    code = main([
        "--msrc-dir", str(msrc_dir(tmp_path)),
        "--hyb-grid-dir", str(tmp_path / "missing"),
        "--output", str(output),
    ])
    assert code == 2
    assert not output.exists()
