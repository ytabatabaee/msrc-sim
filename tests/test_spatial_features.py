from __future__ import annotations

from msrcsim.spatial_compare_io import load_spatial_model_run
from msrcsim.spatial_features import extract_spatial_features
from spatial_compare_helpers import hyb_dir, msrc_dir


def test_spatial_features_for_msrc_include_breakpoint_distances(tmp_path):
    run = load_spatial_model_run(msrc_dir(tmp_path), "msrc")
    features = extract_spatial_features(run, focal_topology=1, lags_bp=[100])
    assert features["num_ancestry_tracts"] == 0
    assert features["focal_support_in_feature_fraction"] > 0.0
    assert "topology_same_prob_lag_100" in features
    assert features["num_dominant_topology_change_points"] >= 0


def test_spatial_features_for_hybridization_include_tract_lengths(tmp_path):
    run = load_spatial_model_run(hyb_dir(tmp_path), "hybridization")
    features = extract_spatial_features(run, focal_topology=1, lags_bp=[100])
    assert features["num_ancestry_tracts"] == 5
    assert features["mean_introgressed_tract_length"] == 200.0
    assert features["focal_support_in_feature_fraction"] == 1.0
