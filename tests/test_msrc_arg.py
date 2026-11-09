import numpy as np
import pytest

from msrcsim.arg.backgrounds.msrc_two_state import MSRCTwoBackground
from msrcsim.arg.engine import ARGConfig, simulate_arg
from msrcsim.arg.lineages import ARGLineage
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap
from msrcsim.arg.segments import Segment
from msrcsim.species_tree import SpeciesTree


class ConstantHistory:
    structured_origin_age = 1_000.0

    def frequency_at(self, branch_id, age):
        return 0.25

    def next_frequency_boundary(self, branch_id, age):
        return float("inf")

    def terminal_frequency(self, taxon):
        return 0.25


def _model(kappa=0.1):
    return MSRCTwoBackground(ConstantHistory(), 50.0, (25.0, 75.0), kappa)


def test_background_dependent_coalescence_uses_finite_copy_convention():
    model = _model()
    assert model.coalescence_rate("ROOT", 0, 0, 100, 10) == pytest.approx(1 / 150)
    assert model.coalescence_rate("ROOT", 1, 1, 100, 10) == pytest.approx(1 / 50)
    assert model.coalescence_rate("ROOT", 0, 1, 100, 10) == 0


def test_crossover_background_inheritance_is_marker_oriented():
    model = _model()
    assert model.recombination_states("ROOT", 0, (), (), 10, position=40, parent_background=1) == (1, 0)
    assert model.recombination_states("ROOT", 0, (), (), 10, position=60, parent_background=1) == (0, 1)


def test_marker_transfer_measure_includes_stored_marker_without_empty_parent():
    model = _model(kappa=0.0)
    rec_map = PiecewiseConstantRecombinationMap.constant(100, 1.0)
    # Neutral active ancestry is only left of the marker, but the stored
    # segment retains the marker and therefore permits marker-only transfer.
    lineage = ARGLineage(1, "ROOT", (Segment(0, 100, 1),), 0, active=((0, 10),))
    hazard = model.recombination_hazard("ROOT", 0, lineage.active, lineage.segments, 10, rec_map)
    assert hazard == pytest.approx(25 * 1.0 + 25 * 0.75)
    sampled = next(model.sample_recombination("ROOT", 0, lineage.active, lineage.segments, 10, rec_map, np.random.default_rng(seed))
                    for seed in range(100)
                    if model.sample_recombination("ROOT", 0, lineage.active, lineage.segments, 10, rec_map, np.random.default_rng(seed))["kind"] == "background_transfer")
    assert sampled["kind"] == "background_transfer"
    assert not sampled["neutral_region"]


def test_msrc_arg_produces_rooted_valid_tree_sequence_and_reproducible_events():
    t = SpeciesTree("((1:10,2:10)A:5,(3:10,4:10)B:5)ROOT;", 20, 100)
    rec_map = PiecewiseConstantRecombinationMap.constant(1000, 2e-5)
    model = _model()
    kwargs = dict(tree=t, config=ARGConfig(1000, rec_map, model, True, None, {"1": 1, "2": 1, "3": 0, "4": 0}))
    one = simulate_arg(rng=np.random.default_rng(44), **kwargs)
    two = simulate_arg(rng=np.random.default_rng(44), **kwargs)
    assert one.edges == two.edges and one.events == two.events
    tskit = pytest.importorskip("tskit")
    from msrcsim.arg.output.tree_sequence import build_tree_sequence
    ts = build_tree_sequence(one)
    assert ts.num_samples == 4
    assert all(len(tree.roots) == 1 for tree in ts.trees())
    assert all(tree.mrca(0, 1) != tskit.NULL for tree in ts.trees())
