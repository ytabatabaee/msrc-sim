from __future__ import annotations

import csv
import json

import matplotlib
import numpy as np

matplotlib.use("Agg")

from msrcsim.spatial_identifiability import (
    eta_to_hr,
    exact_match_for_gamma,
    feasible_gamma_bounds,
    gamma_grid,
    spatial_classifier_columns,
)
from msrcsim.spatial_identifiability_cli import main
from spatial_compare_helpers import msrc_dir


def _rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def test_exact_quartet_matching_over_feasible_gamma_grid():
    target = (0.45, 0.40, 0.15)
    gammas = gamma_grid(target, 5, 0, 1)
    for gamma in gammas:
        match = exact_match_for_gamma(target, gamma, 0, 1)
        assert match.l1_error < 1e-8
        assert np.max(np.abs(np.asarray(match.fitted_q) - np.asarray(target))) < 1e-8


def test_feasible_gamma_boundaries_for_t1_t2_mixture():
    low, high = feasible_gamma_bounds((0.45, 0.40, 0.15), 0, 1)
    assert low == 0.25
    assert high == 0.70


def test_gamma_grid_uses_interior_points_not_endpoints():
    low, high = feasible_gamma_bounds((0.45, 0.40, 0.15), 0, 1)
    gammas = gamma_grid((0.45, 0.40, 0.15), 2, 0, 1)
    assert gammas == [0.4, 0.55]
    assert all(low < gamma < high for gamma in gammas)


def test_eta_to_hr_conversion():
    hr = eta_to_hr(0.5, 1000.0, 0.2)
    assert hr == 0.0025


def test_spatial_classifier_columns_exclude_latent_and_model_parameters():
    rows = [{
        "model": "hybridization",
        "gamma": 0.4,
        "eta": 1.0,
        "hr": 1e-8,
        "replicate": 0,
        "marginal_q1": 0.4,
        "observed_marginal_q1": 0.42,
        "num_ancestry_tracts": 10,
        "realized_introgressed_fraction": 0.5,
        "mean_introgressed_tract_length": 100.0,
        "t_major": 0.5,
        "t_introgressed": 1.0,
        "expected_introgressed_tract_length_bp": 1000.0,
        "num_dominant_topology_change_points": 3,
        "mean_dominant_topology_run_length": 25.0,
    }]
    cols = spatial_classifier_columns(rows)
    assert "num_dominant_topology_change_points" in cols
    assert "mean_dominant_topology_run_length" in cols
    assert "gamma" not in cols
    assert "eta" not in cols
    assert "t_major" not in cols
    assert "t_introgressed" not in cols
    assert "expected_introgressed_tract_length_bp" not in cols
    assert "num_ancestry_tracts" not in cols
    assert "realized_introgressed_fraction" not in cols


def test_spatial_identifiability_cli_outputs_are_reproducible(tmp_path):
    out = tmp_path / "ident"
    args = [
        "--msrc-dir", str(msrc_dir(tmp_path)),
        "--output-dir", str(out),
        "--major-topology", "12|34",
        "--introgressed-topology", "13|24",
        "--num-gamma", "1",
        "--eta-values", "0.5,1",
        "--replicates", "2",
        "--allow-small-sample",
        "--folds", "2",
        "--window-size-loci", "4",
        "--step-loci", "2",
        "--lags-bp", "100",
        "--seed", "22",
    ]
    assert main(args) == 0
    reps1 = (out / "spatial_identifiability_replicates.csv").read_text()
    summary1 = (out / "spatial_identifiability_summary.csv").read_text()
    assert len(_rows(out / "spatial_identifiability_replicates.csv")) == 8
    summary_rows = _rows(out / "spatial_identifiability_summary.csv")
    assert len(summary_rows) == 2
    assert {float(row["bag_roc_auc"]) for row in summary_rows} == {0.5}
    for name in ("identifiability_heatmap.pdf", "auc_vs_eta.pdf", "bag_vs_spatial_auc.pdf", "example_profiles.pdf"):
        assert (out / name).exists()
        assert (out / name).stat().st_size > 0
    metadata = json.loads((out / "spatial_identifiability_metadata.json").read_text())
    assert "num_ancestry_tracts" not in metadata["spatial_feature_columns"]
    assert "gamma" not in metadata["spatial_feature_columns"]
    assert main(args) == 0
    assert reps1 == (out / "spatial_identifiability_replicates.csv").read_text()
    assert summary1 == (out / "spatial_identifiability_summary.csv").read_text()


def test_spatial_identifiability_cli_rejects_small_analysis_without_override(tmp_path):
    out = tmp_path / "ident_small"
    code = main([
        "--msrc-dir", str(msrc_dir(tmp_path)),
        "--output-dir", str(out),
        "--num-gamma", "1",
        "--eta-values", "1",
        "--replicates", "2",
        "--folds", "2",
        "--window-size-loci", "4",
        "--step-loci", "2",
    ])
    assert code == 2
    assert not (out / "spatial_identifiability_summary.csv").exists()
