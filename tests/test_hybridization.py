import csv
import json
import subprocess
import sys

import numpy as np
import yaml

from msrcsim.ancestry_tracts import introgressed_fraction_by_length, simulate_ancestry_tracts
from msrcsim.hybridization import generate_hybridization_locus_positions, simulate_hybridization
from msrcsim.model_fitting import msc_probabilities


def _config(tmp_path, count=4000, length=20_000_000, make_plots=False):
    return {
        "mode": "pulse_hybridization",
        "seed": 123,
        "chromosome": {
            "length_bp": length,
            "num_loci": count,
            "locus_positions": {"mode": "evenly_spaced"},
        },
        "hybridization": {
            "gamma": 0.25,
            "generations_since_pulse": 200,
            "major": {"topology": "12|34", "internal_branch_length": 0.5},
            "introgressed": {"topology": "13|24", "internal_branch_length": 1.5},
            "donor": "3",
            "recipient": "2",
        },
        "recombination": {"rate_per_bp_per_generation": 1e-8},
        "windows": {"loci_per_window": 50, "step_loci": 25},
        "spatial_statistics": {"lags_bp": [10_000, 100_000, 1_000_000]},
        "output": {
            "directory": str(tmp_path / "hyb"),
            "record_tracts": True,
            "record_loci": True,
            "record_windows": True,
            "make_plots": make_plots,
        },
    }


def _read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def test_marginal_ancestry_fraction_and_tract_lengths():
    rng = np.random.default_rng(10)
    gamma = 0.3
    h = 200
    r = 1e-8
    tracts = simulate_ancestry_tracts(2_000_000_000, gamma, h, r, rng)
    assert abs(introgressed_fraction_by_length(tracts) - gamma) < 0.04

    intro = [t.length_bp for t in tracts if t.ancestry_state == 1 and t.end_bp < 2_000_000_000]
    major = [t.length_bp for t in tracts if t.ancestry_state == 0 and t.end_bp < 2_000_000_000]
    assert abs(np.mean(intro) - 1 / (h * r * (1 - gamma))) / (1 / (h * r * (1 - gamma))) < 0.20
    assert abs(np.mean(major) - 1 / (h * r * gamma)) / (1 / (h * r * gamma)) < 0.20


def test_no_recombination_and_gamma_boundaries():
    rng = np.random.default_rng(11)
    assert len(simulate_ancestry_tracts(1000, 0.2, 0, 1e-8, rng)) == 1
    assert len(simulate_ancestry_tracts(1000, 0.2, 10, 0, rng)) == 1
    zero = simulate_ancestry_tracts(1000, 0, 10, 1e-8, rng)
    one = simulate_ancestry_tracts(1000, 1, 10, 1e-8, rng)
    assert len(zero) == 1 and zero[0].ancestry_state == 0
    assert len(one) == 1 and one[0].ancestry_state == 1


def test_marginal_quartet_mixture_and_parent_limits(tmp_path):
    cfg = _config(tmp_path, count=12_000, length=100_000_000)
    cfg["output"]["record_windows"] = False
    out = simulate_hybridization(cfg)
    summary = json.loads((out / "hybridization_summary.json").read_text())
    expected = np.asarray(summary["expected_marginal_q"])
    observed = np.asarray(summary["observed_marginal_q"])
    assert np.max(np.abs(observed - expected)) < 0.04

    cfg["hybridization"]["gamma"] = 0.0
    cfg["output"]["directory"] = str(tmp_path / "major")
    major = json.loads((simulate_hybridization(cfg) / "hybridization_summary.json").read_text())
    assert np.max(np.abs(np.asarray(major["observed_marginal_q"]) - msc_probabilities(0, 0.5))) < 0.04

    cfg["hybridization"]["gamma"] = 1.0
    cfg["output"]["directory"] = str(tmp_path / "intro")
    intro = json.loads((simulate_hybridization(cfg) / "hybridization_summary.json").read_text())
    assert np.max(np.abs(np.asarray(intro["observed_marginal_q"]) - msc_probabilities(1, 1.5))) < 0.04


def test_spatial_correlation_decays_with_lag(tmp_path):
    cfg = _config(tmp_path, count=6000, length=100_000_000)
    cfg["output"]["record_windows"] = False
    cfg["spatial_statistics"]["lags_bp"] = [10_000, 10_000_000]
    summary = json.loads((simulate_hybridization(cfg) / "hybridization_summary.json").read_text())
    lag_stats = summary["lag_statistics"]
    assert lag_stats[0]["prob_same_ancestry"] > lag_stats[1]["prob_same_ancestry"]


def test_determinism_fixed_config(tmp_path):
    cfg = _config(tmp_path, count=1000)
    cfg["output"]["record_windows"] = False
    out1 = simulate_hybridization(cfg)
    tracts1 = (out1 / "hybridization_tracts.csv").read_text()
    loci1 = (out1 / "hybridization_loci.csv").read_text()
    cfg["output"]["directory"] = str(tmp_path / "hyb2")
    out2 = simulate_hybridization(cfg)
    assert tracts1 == (out2 / "hybridization_tracts.csv").read_text()
    assert loci1 == (out2 / "hybridization_loci.csv").read_text()


def test_file_locus_positions_are_sorted_and_ids_preserved(tmp_path):
    path = tmp_path / "loci.tsv"
    path.write_text("locus_id\tposition_bp\n7\t90\n3\t10\n")
    chromosome = {
        "length_bp": 100,
        "num_loci": 2,
        "locus_positions": {"mode": "file", "path": str(path), "column": "position_bp"},
    }
    positions, locus_ids = generate_hybridization_locus_positions(chromosome, np.random.default_rng(1))
    assert positions.tolist() == [10.0, 90.0]
    assert locus_ids == [3, 7]


def test_outputs_cli_and_plot_smoke(tmp_path):
    cfg = _config(tmp_path, count=50, length=1_000_000, make_plots=True)
    cfg["windows"] = {"loci_per_window": 50, "step_loci": 50}
    cfg_path = tmp_path / "pulse.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))
    result = subprocess.run(
        [sys.executable, "-m", "msrcsim.hybridization_cli", "--config", str(cfg_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    out = tmp_path / "hyb"
    assert "Wrote" in result.stdout
    for name in (
        "hybridization_tracts.csv",
        "hybridization_loci.csv",
        "hybridization_windows.csv",
        "hybridization_spatial_autocorrelation.csv",
        "hybridization_summary.json",
        "hybridization_spatial_profile.png",
    ):
        assert (out / name).exists()
    loci = _read_csv(out / "hybridization_loci.csv")
    windows = _read_csv(out / "hybridization_windows.csv")
    assert loci[0]["model"] == "pulse_hybridization"
    assert windows[0]["model"] == "pulse_hybridization"
