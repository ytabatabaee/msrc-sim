from __future__ import annotations

import json

import matplotlib
import yaml

matplotlib.use("Agg")

from msrcsim.hybridization_match import fit_hybridization_to_q
from msrcsim.spatial_distinguishability_cli import main
from spatial_compare_helpers import msrc_dir


def _matched_yaml(path):
    match = fit_hybridization_to_q((0.6, 0.4, 0.0), "12|34", "13|24")
    path.write_text(yaml.safe_dump({
        "mode": "matched_pulse_hybridization",
        "target_q": list(match.target_q),
        "fitted_q": list(match.fitted_q),
        "l1_error": match.l1_error,
        "hybridization": {
            "gamma": match.gamma,
            "t_major": match.t_major,
            "t_introgressed": match.t_introgressed,
            "major_topology": match.major_topology,
            "introgressed_topology": match.introgressed_topology,
        },
        "simulation_template": {
            "chromosome": {"length_bp": 1000, "num_loci": 60, "locus_positions": {"mode": "evenly_spaced"}},
            "windows": {"loci_per_window": 10, "step_loci": 5},
        },
    }, sort_keys=False))


def test_spatial_distinguishability_cli_small_reproducible_run(tmp_path):
    msrc = msrc_dir(tmp_path)
    matched = tmp_path / "matched.yaml"
    _matched_yaml(matched)
    out = tmp_path / "experiment"
    args = [
        "--msrc-dir", str(msrc),
        "--matched-hybridization-yaml", str(matched),
        "--output-dir", str(out),
        "--h-values", "10,20",
        "--r-multipliers", "1",
        "--baseline-r", "1e-6",
        "--replicates", "2",
        "--num-loci", "60",
        "--window-size-loci", "10",
        "--step-loci", "5",
        "--lags-bp", "100,200",
        "--seed", "9",
    ]
    assert main(args) == 0
    rows1 = (out / "spatial_distinguishability.csv").read_text()
    metrics1 = (out / "spatial_distinguishability_metrics.csv").read_text()
    assert (out / "spatial_distinguishability.png").exists()
    summary = json.loads((out / "spatial_distinguishability_summary.json").read_text())
    assert summary["num_rows"] == 8
    assert "topology_same_prob_lag_100" in summary["feature_columns"]
    assert main(args) == 0
    assert rows1 == (out / "spatial_distinguishability.csv").read_text()
    assert metrics1 == (out / "spatial_distinguishability_metrics.csv").read_text()
