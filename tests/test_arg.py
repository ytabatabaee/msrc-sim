import numpy as np
import pytest
import csv
import re

from msrcsim.arg.engine import ARGConfig, simulate_arg
from msrcsim.arg.coalescence import union_lineage
from msrcsim.arg.lineages import ARGLineage
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap
from msrcsim.species_tree import SpeciesTree
from msrcsim.arg.segments import Segment


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


def test_gap_between_segments_is_a_valid_crossover_cut():
    rec_map = PiecewiseConstantRecombinationMap.constant(10, 1.0)
    segments = ((0.0, 1.0), (3.0, 4.0))
    assert rec_map.breakpoint_hazard(segments) == pytest.approx(4.0)
    positions = [rec_map.sample_position(segments, np.random.default_rng(i)) for i in range(100)]
    assert any(1.0 < x < 3.0 for x in positions)


def test_disjoint_coalescence_does_not_change_local_roots():
    a = ARGLineage(1, "ROOT", (Segment(0, 1, 10),), 0.0)
    b = ARGLineage(2, "ROOT", (Segment(2, 3, 11),), 0.0)
    parent = union_lineage(a, b, None, 5.0, 3)
    assert [(s.left, s.right, s.node) for s in parent.segments] == [(0, 1, 10), (2, 3, 11)]


def test_exported_gene_trees_have_original_taxa(tmp_path):
    pytest.importorskip("tskit")
    from msrcsim.arg.output.genealogy_export import export_genealogies
    from msrcsim.arg.output.tree_sequence import build_tree_sequence
    result = simulate_arg(_tree(), ARGConfig(1000, PiecewiseConstantRecombinationMap.constant(1000, 1e-5)), np.random.default_rng(9))
    ts = build_tree_sequence(result)
    newick = tmp_path / "trees.nwk"
    csv_path = tmp_path / "trees.csv"
    export_genealogies(ts, newick, csv_path)
    lines = newick.read_text().splitlines()
    assert lines and all(not line.startswith("[") and line.endswith(";") for line in lines)
    for line in lines:
        assert set(re.findall(r"(?<![0-9.])([1-4])(?=[:),])", line)) == {"1", "2", "3", "4"}
    with csv_path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert rows and list(rows[0])[:6] == ["chrom", "left", "right", "midpoint", "tree_id", "newick"]
