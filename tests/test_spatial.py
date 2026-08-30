import csv
import json
import subprocess
import sys

import numpy as np
import pytest

from msrcsim.history_io import save_frozen_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.spatial import simulate_spatial
from msrcsim.spatial_plotting import plot_spatial_output
from msrcsim.spatial_statistics import (
    RearrangementInterval,
    breakpoint_jumps,
    classify_region,
    generate_locus_positions,
    sliding_windows,
)
from msrcsim.species_tree import SpeciesTree
from msrcsim.structured_coalescent import simulate_msc_genealogy
from msrcsim.wright_fisher import simulate_frequency_history


def _config(tmp_path, count=60, placement="evenly_spaced"):
    return {
        "mode": "spatial",
        "seed": 7,
        "species_tree": {
            "newick": "((1:40,2:40)A:20,(3:40,4:40)B:20)ROOT;",
            "root_extension": 100,
            "default_effective_population_size": 50,
        },
        "rearrangement": {
            "id": "inv_1",
            "type": "inversion",
            "origin_branch": "ROOT",
            "origin_time_from_branch_start": 40,
            "initial_copy_count": 50,
            "selection": {"model": "genic", "coefficient": 0.0},
        },
        "recombination": {"baseline_rate": 0.02, "effective_cross_arrangement_fraction": 0.1},
        "genome": {
            "length": 1000,
            "loci": {"count": count, "placement": placement},
            "rearrangement_interval": {"start": 250, "end": 650},
            "inside_model": {"type": "msrc", "effective_cross_arrangement_fraction": 0.05},
            "outside_model": {"type": "msc"},
        },
        "spatial_summary": {"window_loci": 10, "step_loci": 5, "breakpoint_bandwidth_loci": 10},
        "output": {
            "directory": str(tmp_path / "spatial"),
            "record_gene_trees": True,
            "record_frequency_history": True,
            "make_plots": False,
        },
    }


def _read_csv(path, delimiter=","):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def test_region_classification():
    interval = RearrangementInterval(25, 65)
    assert classify_region(24, interval) == ("outside_left", False)
    assert classify_region(25, interval) == ("rearrangement", True)
    assert classify_region(65, interval) == ("rearrangement", True)
    assert classify_region(66, interval) == ("outside_right", False)


def test_locus_placement_reproducible_and_unique():
    a = generate_locus_positions(1000, 25, "evenly_spaced", np.random.default_rng(1))
    b = generate_locus_positions(1000, 25, "evenly_spaced", np.random.default_rng(2))
    assert np.array_equal(a, b)
    c = generate_locus_positions(1000, 25, "uniform_random", np.random.default_rng(3))
    d = generate_locus_positions(1000, 25, "uniform_random", np.random.default_rng(3))
    assert np.array_equal(c, d)
    assert np.all(np.diff(c) > 0)


def test_sliding_window_counts_and_breakpoint_fraction():
    rows = []
    interval = RearrangementInterval(50, 149)
    for i, pos in enumerate(range(0, 200, 10)):
        region, inside = classify_region(pos, interval)
        rows.append({"position": pos, "topology_index": i % 3, "inside_rearrangement": inside, "region": region})
    windows = sliding_windows(rows, window_loci=10, step_loci=5)
    for row in windows:
        assert row["n1"] + row["n2"] + row["n3"] == row["n_loci"]
        assert abs(row["q1"] + row["q2"] + row["q3"] - 1.0) < 1e-12
    crossing = [w for w in windows if w["start_position"] < 50 < w["end_position"]][0]
    assert 0.0 < crossing["fraction_inside_rearrangement"] < 1.0


def test_outside_msc_equal_discordants():
    tree = SpeciesTree("((1:80,2:80)A:20,(3:80,4:80)B:20)ROOT;", 80, 200)
    rng = np.random.default_rng(11)
    tops = [simulate_msc_genealogy(i, tree, rng, False).topology_index for i in range(2500)]
    q = np.bincount(tops, minlength=3) / len(tops)
    assert q[0] > q[1]
    assert abs(q[1] - q[2]) < 0.04


def test_spatial_outputs_and_gene_tree_positions(tmp_path):
    out = simulate_spatial(_config(tmp_path, count=40))
    loci = _read_csv(out / "spatial_loci.csv")
    trees = _read_csv(out / "spatial_gene_trees.tsv", delimiter="\t")
    windows = _read_csv(out / "spatial_windows.csv")
    summary = json.loads((out / "spatial_summary.json").read_text())
    assert len(loci) == 40
    assert len(trees) == 40
    assert trees[0]["position"] == loci[0]["position"]
    assert windows
    assert summary["linked_loci_model"] is False
    assert {r["local_model"] for r in loci} == {"msc", "msrc"}


def test_frozen_history_replay_uses_same_terminal_pattern(tmp_path):
    cfg = _config(tmp_path, count=20)
    tree = SpeciesTree(cfg["species_tree"]["newick"], 50, 100)
    rr = Rearrangement("inv_1", "inversion", "ROOT", 40, 50, 0.0)
    hist = simulate_frequency_history(tree, rr, np.random.default_rng(12))
    frozen = tmp_path / "frozen.yaml"
    save_frozen_history(frozen, cfg, hist, {"1": 1, "2": 0, "3": 1, "4": 0}, {"label": "fixed"})
    cfg["history"] = {"frozen_history": str(frozen)}
    out1 = simulate_spatial(cfg)
    cfg2 = _config(tmp_path, count=20)
    cfg2["output"]["directory"] = str(tmp_path / "spatial2")
    cfg2["genome"]["inside_model"]["effective_cross_arrangement_fraction"] = 0.5
    cfg2["history"] = {"frozen_history": str(frozen)}
    out2 = simulate_spatial(cfg2)
    s1 = json.loads((out1 / "spatial_summary.json").read_text())
    s2 = json.loads((out2 / "spatial_summary.json").read_text())
    assert s1["terminal_pattern"] == s2["terminal_pattern"] == "1010"
    assert s1["history_metadata"]["label"] == s2["history_metadata"]["label"] == "fixed"


def test_breakpoint_jump_synthetic_sequence():
    rows = []
    for i, pos in enumerate(range(100)):
        rows.append({"position": pos, "topology_index": 0 if pos < 50 else 1, "inside_rearrangement": pos >= 50})
    jump = breakpoint_jumps(rows, [50], bandwidth_loci=10)[0]
    assert jump["left_q1"] == 1.0
    assert jump["right_q2"] == 1.0
    assert abs(jump["quartet_jump_l2"] - np.sqrt(2.0)) < 1e-12


def test_cli_smoke_spatial(tmp_path):
    cfg = _config(tmp_path, count=20)
    cfg_path = tmp_path / "spatial.yaml"
    import yaml
    cfg_path.write_text(yaml.safe_dump(cfg))
    result = subprocess.run(
        [sys.executable, "-m", "msrcsim.spatial_cli", "--config", str(cfg_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "Wrote" in result.stdout
    assert (tmp_path / "spatial" / "spatial_loci.csv").exists()


def test_plot_generation_smoke(tmp_path):
    cfg = _config(tmp_path, count=40)
    out = simulate_spatial(cfg)
    figures = plot_spatial_output(out, fmt="png")
    assert figures
    assert all(p.exists() and p.stat().st_size > 0 for p in figures)
