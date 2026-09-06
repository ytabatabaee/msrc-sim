import csv
import json

import numpy as np

from msrcsim.genomic import GenomicInterval
from msrcsim.linked_spatial import generate_genealogy_breakpoints, simulate_linked_spatial
from msrcsim.spatial import simulate_spatial


def _config(tmp_path, directory="linked", seed=23, kappa=0.2):
    return {
        "mode": "spatial",
        "seed": seed,
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
            "loci": {"count": 40, "placement": "evenly_spaced"},
            "rearrangement_interval": {"id": "inv_1", "start": 250, "end": 650},
            "inside_model": {"type": "msrc", "effective_cross_arrangement_fraction": 0.05},
            "outside_model": {"type": "msc"},
        },
        "linked_spatial": {
            "enabled": True,
            "kappa": kappa,
            "breakpoint_rate_per_bp": 0.02,
            "windows": {"count": 40},
        },
        "output": {
            "directory": str(tmp_path / directory),
            "record_resolved_config": True,
            "make_plots": False,
        },
    }


def _legacy_config(tmp_path, directory="spatial", seed=7):
    cfg = _config(tmp_path, directory=directory, seed=seed)
    cfg.pop("linked_spatial")
    cfg["spatial_summary"] = {"window_loci": 10, "step_loci": 5, "breakpoint_bandwidth_loci": 10}
    cfg["output"].update({"record_gene_trees": True, "record_frequency_history": True})
    return cfg


def _read_rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def _mean_block_lengths(kappa, replicates=400):
    intervals = [GenomicInterval("chr1", 250.0, 750.0, "inv")]
    inside = []
    outside = []
    for seed in range(replicates):
        bps = generate_genealogy_breakpoints(1000.0, intervals, 0.03, kappa, np.random.default_rng(seed))
        for start, end in zip(bps, bps[1:]):
            mid = (start + end) / 2.0
            if 250.0 <= mid < 750.0:
                inside.append(end - start)
            else:
                outside.append(end - start)
    return float(np.mean(inside)), float(np.mean(outside))


def test_legacy_spatial_unchanged_when_linked_mode_off(tmp_path):
    cfg1 = _legacy_config(tmp_path, "a")
    cfg2 = _legacy_config(tmp_path, "b")
    cfg2["linked_spatial"] = {"enabled": False}
    out1 = simulate_spatial(cfg1)
    out2 = simulate_spatial(cfg2)
    assert (out1 / "spatial_loci.csv").read_text() == (out2 / "spatial_loci.csv").read_text()
    assert not (out2 / "spatial_genealogies.csv").exists()
    assert json.loads((out2 / "spatial_summary.json").read_text())["linked_loci_model"] is False


def test_positions_and_block_ids_are_valid_and_contiguous(tmp_path):
    out = simulate_linked_spatial(_config(tmp_path))
    rows = _read_rows(out / "spatial_genealogies.csv")
    diagnostics = _read_rows(out / "spatial_linkage_diagnostics.csv")
    ids = [int(row["block_id"]) for row in rows]
    assert ids == sorted(ids)
    assert sorted(set(ids)) == list(range(max(ids) + 1))
    assert {row["region"] for row in diagnostics} == {"inside", "outside"}
    previous_midpoint = -1.0
    for row in rows:
        start = float(row["start"])
        end = float(row["end"])
        midpoint = float(row["midpoint"])
        assert 0.0 <= start < end <= 1000.0
        assert start <= midpoint <= end
        assert midpoint > previous_midpoint
        previous_midpoint = midpoint


def test_kappa_one_has_no_rearrangement_specific_suppression():
    interval = [GenomicInterval("chr1", 250.0, 750.0, "inv")]
    inside_counts = []
    outside_counts = []
    for seed in range(100):
        bps = generate_genealogy_breakpoints(1000.0, interval, 0.03, 1.0, np.random.default_rng(seed))
        inside_counts.append(sum(250.0 < bp < 750.0 for bp in bps))
        outside_counts.append(sum((0.0 < bp < 250.0) or (750.0 < bp < 1000.0) for bp in bps))
    assert abs(np.mean(inside_counts) - np.mean(outside_counts)) < 3.0


def test_decreasing_kappa_increases_expected_rearranged_block_persistence():
    interval = [GenomicInterval("chr1", 250.0, 750.0, "inv")]
    def mean_inside_breakpoints(kappa):
        vals = []
        for seed in range(150):
            bps = generate_genealogy_breakpoints(1000.0, interval, 0.04, kappa, np.random.default_rng(seed))
            vals.append(sum(250.0 < bp < 750.0 for bp in bps))
        return float(np.mean(vals))
    assert mean_inside_breakpoints(0.1) < mean_inside_breakpoints(0.8)


def test_kappa_controls_expected_block_lengths_statistically():
    inside_1, outside_1 = _mean_block_lengths(1.0)
    inside_05, outside_05 = _mean_block_lengths(0.5)
    inside_02, outside_02 = _mean_block_lengths(0.2)
    assert abs(inside_1 - outside_1) / outside_1 < 0.15
    assert inside_05 > outside_05
    assert inside_02 > inside_05
    assert outside_02 > 0.0


def test_kappa_zero_has_no_internal_rearrangement_breakpoints():
    intervals = [GenomicInterval("chr1", 250.0, 750.0, "inv")]
    bps = generate_genealogy_breakpoints(1000.0, intervals, 0.05, 0.0, np.random.default_rng(2))
    assert not any(250.0 < bp < 750.0 for bp in bps)
    assert 250.0 in bps and 750.0 in bps


def test_reproducible_by_seed(tmp_path):
    out1 = simulate_linked_spatial(_config(tmp_path, "one", seed=91))
    out2 = simulate_linked_spatial(_config(tmp_path, "two", seed=91))
    assert (out1 / "spatial_genealogies.csv").read_text() == (out2 / "spatial_genealogies.csv").read_text()
