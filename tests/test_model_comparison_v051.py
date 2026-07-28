import csv
from pathlib import Path

from msrcsim.model_comparison_cli import _normalize_terminal_pattern
from msrcsim.model_fitting import compare_models, fit_network_pair


def test_terminal_pattern_normalization():
    assert _normalize_terminal_pattern("101") == "0101"
    assert _normalize_terminal_pattern("101.0") == "0101"
    assert _normalize_terminal_pattern("0101") == "0101"


def test_comparison_separates_geometry_and_evidence():
    result = compare_models((169, 758, 73))
    assert result["network_representable"]
    assert result["off_arm_supported"]
    assert "model_geometry_classification" in result
    assert "model_evidence_classification" in result


def test_network_boundary_flags_at_optimizer_ceiling():
    fit = fit_network_pair((169, 758, 73), 0, 1, max_branch=20.0)
    assert isinstance(fit.boundary_warning, bool)
    assert isinstance(fit.well_interior, bool)
