from __future__ import annotations

import csv
import json
import math

import numpy as np
import pytest
import yaml

from msrcsim.arbitrary import run_arbitrary_simulation
from msrcsim.experiments import simulate_parameter_grid
from msrcsim.history_io import load_frozen_history, save_frozen_history
from msrcsim.population_process import simulate_population_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.speciation_split import (
    daughter_target_frequency,
    mass_conserving_daughter_frequencies,
    resolve_speciation_split,
    sample_daughter_initialization,
)
from msrcsim.species_tree import SpeciesTree


def _binary_tree(ne: int = 20) -> SpeciesTree:
    return SpeciesTree("(L:1,R:1)ROOT;", ne, 0)


def _root_rearrangement(k: int = 10) -> Rearrangement:
    return Rearrangement("x", "inversion", "ROOT", 0, k, 0.0)


def _rows(history):
    return [
        (
            r.rearrangement_id,
            r.branch_id,
            r.parent_branch_id,
            r.forward_generation,
            r.absolute_age,
            r.copy_count_A1,
            r.copy_count_A0,
            r.population_chromosomes,
            r.frequency_A1,
            r.status,
            r.is_origin,
            r.is_branch_start,
            r.is_branch_end,
        )
        for r in history.records
    ]


def test_reciprocal_log_odds_transform_examples_and_identity():
    eta = math.log(3.0)
    assert daughter_target_frequency(0.25, eta) == pytest.approx(0.5)
    assert daughter_target_frequency(0.25, -eta) == pytest.approx(0.1)
    for p in (1e-12, 0.01, 0.25, 0.5, 0.99, 1.0 - 1e-12):
        assert daughter_target_frequency(p, 0.0) == pytest.approx(p)


def test_reciprocal_log_odds_state_preserving_and_stable():
    for eta in (-100.0, -50.0, 0.0, 50.0, 100.0):
        assert daughter_target_frequency(0.0, eta) == 0.0
        assert daughter_target_frequency(1.0, eta) == 1.0
        value = daughter_target_frequency(0.25, eta)
        assert math.isfinite(value)
        assert 0.0 <= value <= 1.0
    assert daughter_target_frequency(0.25, 2.0) > 0.25
    assert daughter_target_frequency(0.25, -2.0) < 0.25
    assert daughter_target_frequency(0.25, -2.0) == pytest.approx(
        1.0 - daughter_target_frequency(0.75, 2.0)
    )
    for bad in (math.inf, -math.inf, math.nan):
        with pytest.raises(ValueError, match="finite"):
            daughter_target_frequency(0.25, bad)


def test_mass_conserving_transform_properties():
    for p in (0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0):
        for eta in (-4.0, -2.0, 0.0, 2.0, 4.0):
            positive, negative = mass_conserving_daughter_frequencies(p, eta)
            assert 0.0 <= positive <= 1.0
            assert 0.0 <= negative <= 1.0
            assert (positive + negative) / 2.0 == pytest.approx(p)

    assert mass_conserving_daughter_frequencies(0.0, 4.0) == (0.0, 0.0)
    assert mass_conserving_daughter_frequencies(1.0, -4.0) == (1.0, 1.0)
    assert mass_conserving_daughter_frequencies(0.25, 0.0) == (0.25, 0.25)

    positive_eta = mass_conserving_daughter_frequencies(0.25, 2.0)
    negative_eta = mass_conserving_daughter_frequencies(0.25, -2.0)
    assert positive_eta[0] > 0.25
    assert positive_eta[1] < 0.25
    assert negative_eta[0] < 0.25
    assert negative_eta[1] > 0.25
    assert positive_eta == pytest.approx((negative_eta[1], negative_eta[0]))


def test_resolve_speciation_split_validation_and_orientation():
    tree = SpeciesTree("(B:1,A:1)ROOT;", 10, 0)
    policy = resolve_speciation_split({"speciation_split": {"model": "reciprocal_log_odds", "eta": 1.5}}, tree)
    assert policy.rule_for("ROOT").positive_child == "B"
    assert policy.rule_for("ROOT").positive_child_source == "newick_order_fallback"

    cfg = {
        "speciation_split": {
            "model": "reciprocal_log_odds",
            "eta": 0.0,
            "by_split": {"ROOT": {"eta": 2.0, "positive_child": "A"}},
        }
    }
    policy = resolve_speciation_split(cfg, tree)
    assert policy.rule_for("ROOT").eta == 2.0
    assert policy.rule_for("ROOT").positive_child == "A"

    mass_policy = resolve_speciation_split({"speciation_split": {"model": "mass_conserving", "eta": 1.5}}, tree)
    assert mass_policy.model == "mass_conserving"

    with pytest.raises(ValueError, match="unknown split"):
        resolve_speciation_split({"speciation_split": {"by_split": {"NOPE": {"eta": 1.0}}}}, tree)
    with pytest.raises(ValueError, match="direct child"):
        resolve_speciation_split({"speciation_split": {"by_split": {"ROOT": {"positive_child": "NOPE"}}}}, tree)
    with pytest.raises(ValueError, match="model"):
        resolve_speciation_split({"speciation_split": {"model": "other"}}, tree)
    with pytest.raises(ValueError, match="by_split must be a mapping"):
        resolve_speciation_split({"speciation_split": {"by_split": []}}, tree)
    with pytest.raises(ValueError, match="finite"):
        resolve_speciation_split({"speciation_split": {"eta": math.nan}}, tree)


def test_multifurcation_eta_zero_allowed_nonzero_rejected():
    tree = SpeciesTree("(A:1,B:1,C:1)ROOT;", 10, 0)
    policy = resolve_speciation_split({"speciation_split": {"model": "reciprocal_log_odds", "eta": 0.0}}, tree)
    assert policy.rule_for("ROOT").eta == 0.0
    with pytest.raises(ValueError, match="requires a binary split"):
        resolve_speciation_split({"speciation_split": {"model": "reciprocal_log_odds", "eta": 1.0}}, tree)


def test_no_recurrent_creation_and_boundary_states():
    tree = _binary_tree()
    policy = resolve_speciation_split({"speciation_split": {"eta": 100.0, "positive_child": "L"}}, tree)

    class NoDrawRng:
        def binomial(self, *_args, **_kwargs):
            raise AssertionError("binomial should not be called")

    start, event = sample_daughter_initialization(
        tree=tree,
        parent_branch_id="ROOT",
        child_branch_id="L",
        parent_count=None,
        parent_total=None,
        child_total=40,
        rng=NoDrawRng(),
        policy=policy,
    )
    assert start is None
    assert event.parent_state == "not_present"

    rng = np.random.default_rng(1)
    start, event = sample_daughter_initialization(
        tree=tree,
        parent_branch_id="ROOT",
        child_branch_id="L",
        parent_count=0,
        parent_total=40,
        child_total=40,
        rng=rng,
        policy=policy,
    )
    assert start == 0
    assert event.target_frequency_A1 == 0.0
    start, event = sample_daughter_initialization(
        tree=tree,
        parent_branch_id="ROOT",
        child_branch_id="R",
        parent_count=40,
        parent_total=40,
        child_total=40,
        rng=rng,
        policy=policy,
    )
    assert start == 40
    assert event.target_frequency_A1 == 1.0


def test_stochastic_daughter_sampling_means_match_targets():
    tree = _binary_tree()
    eta = math.log(3.0)
    policy = resolve_speciation_split({"speciation_split": {"eta": eta, "by_split": {"ROOT": {"positive_child": "L"}}}}, tree)
    rng = np.random.default_rng(7)
    pos = []
    neg = []
    for _ in range(2000):
        k, _event = sample_daughter_initialization(
            tree=tree, parent_branch_id="ROOT", child_branch_id="L",
            parent_count=10, parent_total=40, child_total=40, rng=rng, policy=policy
        )
        pos.append(k / 40)
        k, _event = sample_daughter_initialization(
            tree=tree, parent_branch_id="ROOT", child_branch_id="R",
            parent_count=10, parent_total=40, child_total=40, rng=rng, policy=policy
        )
        neg.append(k / 40)
    assert float(np.mean(pos)) == pytest.approx(0.5, abs=0.02)
    assert float(np.mean(neg)) == pytest.approx(0.1, abs=0.02)


def test_wright_fisher_and_moran_share_split_targets():
    cfg = {"speciation_split": {"model": "reciprocal_log_odds", "eta": math.log(3.0), "by_split": {"ROOT": {"positive_child": "L"}}}}
    for model in ("wright_fisher", "moran"):
        tree = _binary_tree()
        config = dict(cfg)
        config["population_process"] = {"model": model}
        history = simulate_population_history(tree, _root_rearrangement(), np.random.default_rng(11), config)
        by_child = {event.child_branch_id: event for event in history.split_events}
        assert by_child["L"].parent_frequency_A1 == pytest.approx(0.25)
        assert by_child["L"].target_frequency_A1 == pytest.approx(0.5)
        assert by_child["R"].target_frequency_A1 == pytest.approx(0.1)
        assert by_child["L"].eta == pytest.approx(math.log(3.0))


def test_absent_and_explicit_eta_zero_histories_match():
    tree = _binary_tree()
    rr = _root_rearrangement()
    absent = simulate_population_history(tree, rr, np.random.default_rng(12), {})
    explicit = simulate_population_history(
        tree,
        rr,
        np.random.default_rng(12),
        {"speciation_split": {"model": "reciprocal_log_odds", "eta": 0.0}},
    )
    assert _rows(absent) == _rows(explicit)
    assert absent.split_events == []
    assert len(explicit.split_events) == 2


def test_eta_nonzero_reproducibility():
    tree = _binary_tree()
    cfg = {"speciation_split": {"model": "reciprocal_log_odds", "eta": 1.25}}
    h1 = simulate_population_history(tree, _root_rearrangement(), np.random.default_rng(13), cfg)
    h2 = simulate_population_history(tree, _root_rearrangement(), np.random.default_rng(13), cfg)
    assert _rows(h1) == _rows(h2)
    assert h1.split_events == h2.split_events


def test_arbitrary_tree_records_split_events_and_valid_gene_trees(tmp_path):
    cfg = {
        "mode": "arbitrary_tree",
        "seed": 14,
        "num_loci": 4,
        "species_tree": {
            "newick": "((1:1,2:1)A:1,(3:1,4:1)B:1)ROOT;",
            "root_extension": 0,
            "default_effective_population_size": 20,
        },
        "rearrangement": {"id": "inv", "type": "inversion", "origin_branch": "ROOT", "origin_time_from_branch_start": 0, "initial_copy_count": 10},
        "recombination": {"baseline_rate": 0.01, "effective_cross_arrangement_fraction": 0.1},
        "genome": {"length_bp": 1000, "breakpoint_rate_per_bp": 0.0, "kappa": 0.2, "rearrangement_interval": {"start": 100, "end": 900}},
        "quartets": [{"taxa": ["1", "2", "3", "4"]}],
        "speciation_split": {"model": "reciprocal_log_odds", "eta": math.log(3.0), "by_split": {"ROOT": {"positive_child": "A"}}},
        "output": {"directory": str(tmp_path / "arbitrary")},
    }
    out = run_arbitrary_simulation(cfg)
    summary = json.loads((out / "simulation_summary.json").read_text())
    assert summary["gene_tree_validity"]["all_expected_taxa_once"]
    assert summary["speciation_split"]["eta"] == pytest.approx(math.log(3.0))
    rows = list(csv.DictReader((out / "speciation_split_events.tsv").open(), delimiter="\t"))
    root_rows = [row for row in rows if row["parent_branch_id"] == "ROOT"]
    assert {row["child_branch_id"] for row in root_rows} == {"A", "B"}
    assert float(next(row for row in root_rows if row["child_branch_id"] == "A")["target_frequency_A1"]) == pytest.approx(0.5)


def test_parameter_grid_can_vary_split_eta(tmp_path):
    cfg = {
        "seed": 15,
        "species_tree": {"newick": "((1:1,2:1)A:1,(3:1,4:1)B:1)ROOT;", "root_extension": 0, "default_effective_population_size": 20},
        "rearrangement": {"id": "inv", "type": "inversion", "origin_branch": "ROOT", "origin_time_from_branch_start": 0, "initial_copy_count": 10},
        "recombination": {"baseline_rate": 0.01, "effective_cross_arrangement_fraction": 0.1},
        "conditioning": {"mode": "none", "max_attempts": 10},
        "experiment": {"replicates": 1, "loci_per_replicate": 2, "asymmetry_threshold": 0.1},
        "speciation_split": {"model": "reciprocal_log_odds", "eta": 0.0, "by_split": {"ROOT": {"eta": 0.0, "positive_child": "A"}}},
        "parameter_grid": {"speciation_split.by_split.ROOT.eta": [-4.0, 0.0, 4.0]},
        "output": {"directory": str(tmp_path / "grid")},
    }
    rows, out = simulate_parameter_grid(cfg)
    assert [row["speciation_split.by_split.ROOT.eta"] for row in rows] == [-4.0, 0.0, 4.0]
    assert (out / "cell_0002" / "prevalence_summary.json").exists()


def test_frozen_history_roundtrip_preserves_split_events_and_old_files_load(tmp_path):
    tree = _binary_tree()
    cfg = {"species_tree": {"newick": "(L:1,R:1)ROOT;", "default_effective_population_size": 20, "root_extension": 0}, "speciation_split": {"eta": 1.0}}
    history = simulate_population_history(tree, _root_rearrangement(), np.random.default_rng(16), cfg)
    path = save_frozen_history(tmp_path / "new.yaml", cfg, history, {"L": 1, "R": 0})
    _config, loaded, _sampled, metadata = load_frozen_history(path)
    assert len(loaded.split_events) == 2
    assert metadata["speciation_split"]["enabled"] is True

    old_path = tmp_path / "old.yaml"
    payload = yaml.safe_load(path.read_text())
    payload.pop("split_events")
    payload["metadata"].pop("speciation_split")
    old_path.write_text(yaml.safe_dump(payload, sort_keys=False))
    _config, old_loaded, _sampled, metadata = load_frozen_history(old_path)
    assert old_loaded.split_events == []
    assert "speciation_split" not in metadata
