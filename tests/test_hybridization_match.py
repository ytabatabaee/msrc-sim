from __future__ import annotations

import json

import numpy as np
import yaml

from msrcsim.hybridization_match import fit_hybridization_to_q, mixture_q
from msrcsim.hybridization_match_cli import main
from spatial_compare_helpers import msrc_dir


def test_exact_hybridization_match_for_representable_target():
    target = mixture_q(0, 1, 0.35, 0.7, 1.4)
    match = fit_hybridization_to_q(target, "12|34", "13|24")
    assert match.l1_error < 1e-5
    assert np.max(np.abs(np.asarray(match.fitted_q) - target)) < 1e-5


def test_exact_hybridization_match_major_parent_limit():
    target = mixture_q(0, 1, 0.0, 0.5, 1.0)
    match = fit_hybridization_to_q(target, 0, 1)
    assert match.l1_error < 1e-5
    assert match.gamma < 1e-3 or np.linalg.norm(np.asarray(match.fitted_q) - target) < 1e-5


def test_match_hybridization_cli_writes_yaml(tmp_path):
    output = tmp_path / "matched.yaml"
    code = main([
        "--msrc-dir", str(msrc_dir(tmp_path)),
        "--major-topology", "12|34",
        "--introgressed-topology", "13|24",
        "--output", str(output),
    ])
    assert code == 0
    data = yaml.safe_load(output.read_text())
    assert data["hybridization"]["gamma"] >= 0.0
    assert data["hybridization"]["t_major"] >= 0.0
    assert data["target_q"] == [0.6, 0.4, 0.0]
    assert "fitted_q" in data
