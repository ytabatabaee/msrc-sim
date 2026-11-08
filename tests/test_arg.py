import numpy as np
import pytest

from msrcsim.arg.engine import ARGConfig, simulate_arg
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap
from msrcsim.species_tree import SpeciesTree


def _tree():
    return SpeciesTree("((1:10,2:10)A:5,(3:10,4:10)B:5)ROOT;", 20, 100)


def test_segments_are_sorted_and_no_recombination_has_one_shared_process():
    result = simulate_arg(_tree(), ARGConfig(1000, PiecewiseConstantRecombinationMap.constant(1000, 0.0)), np.random.default_rng(4))
    assert sum(e.event_type == "recombination" for e in result.events) == 0
    assert len(result.nodes) == 7
    assert all(0 <= left < right <= 1000 for left, right, *_ in result.edges)


def test_recombination_creates_correlated_segment_ancestry():
    result = simulate_arg(_tree(), ARGConfig(1000, PiecewiseConstantRecombinationMap.constant(1000, 2e-4)), np.random.default_rng(7))
    recombinations = [e for e in result.events if e.event_type == "recombination"]
    assert recombinations
    assert all(e.rate >= 0 for e in result.events if e.rate is not None)
    # The process tracks the full chromosome, rather than independent blocks.
    assert max(right for _, right, *_ in result.edges) == 1000


def test_tskit_output_is_valid_when_installed(tmp_path):
    tskit = pytest.importorskip("tskit")
    from msrcsim.arg.output.tree_sequence import build_tree_sequence
    result = simulate_arg(_tree(), ARGConfig(1000, PiecewiseConstantRecombinationMap.constant(1000, 1e-5)), np.random.default_rng(8))
    ts = build_tree_sequence(result)
    assert ts.num_samples == 4
    assert ts.num_trees >= 1
