import numpy as np
import pytest
from msrcsim.species_tree import SpeciesTree
from msrcsim.smc import PiecewiseConstantMap, SMCConfig, simulate_msc_smc

def tree():
    return SpeciesTree("((A:100,B:100)AB:50,(C:100,D:100)CD:50)ROOT;", 100, 1000)

def run(rate, seed=4, recmap=None):
    m=recmap or PiecewiseConstantMap.constant(10000,rate)
    return simulate_msc_smc(tree(),SMCConfig(10000,m,seed))

def test_zero_recombination_one_genealogy_all_positions():
    r=run(0.0)
    assert len(r.trees)==1 and r.trees[0].start==0 and r.trees[0].end==10000
    assert not r.events

def test_seed_reproduces_all_outputs():
    a,b=run(2e-6,77),run(2e-6,77)
    assert a.trees==b.trees and a.events==b.events and a.changes==b.changes

def test_equal_interval_map_is_identical_to_constant_map():
    a=run(1e-7,21)
    b=run(1e-7,21,PiecewiseConstantMap(10000,(3000.,7000.),(1e-7,1e-7,1e-7)))
    assert len(a.trees)==len(b.trees) and len(a.events)==len(b.events)
    assert [t.topology for t in a.trees]==[t.topology for t in b.trees]
    assert [e.event_type for e in a.events]==[e.event_type for e in b.events]
    assert [t.start for t in a.trees]==pytest.approx([t.start for t in b.trees])

def test_map_crosses_boundaries_and_suppression_is_less_dense():
    m=PiecewiseConstantMap(10000,(2000.,8000.),(1e-5,0.,1e-5))
    assert m.advance_hazard(0,0.03)==pytest.approx(9000)
    assert len(run(2e-5,11).events)>len(run(1e-6,11).events)

def test_zero_rate_suppressed_interval_has_no_recombination_breakpoints():
    m=PiecewiseConstantMap(10000,(2000.,8000.),(2e-5,0.,2e-5))
    r=run(2e-5,8,m)
    assert all(not (2000 < e.position < 8000) for e in r.events)

def test_partial_suppression_reduces_event_density_inside_region():
    m=PiecewiseConstantMap(10000,(2000.,8000.),(2e-5,1e-6,2e-5))
    r=run(2e-5,11,m)
    inside=sum(2000<=e.position<8000 for e in r.events)
    outside=len(r.events)-inside
    assert inside/6000 < outside/4000

def test_higher_recombination_increases_tree_turnover():
    low,high=run(2e-7,120),run(2e-5,120)
    low_changes=sum(e.event_type!="no_genealogy_change" for e in low.events)
    high_changes=sum(e.event_type!="no_genealogy_change" for e in high.events)
    assert high_changes>low_changes

def test_topology_autocorrelation_decreases_with_recombination_rate():
    low,high=run(2e-7,305),run(2e-5,305)
    def same_fraction(result, lag=500):
        xs=np.linspace(0,9500-lag,200)
        def at(x):
            return next(t.topology for t in result.trees if t.start<=x<t.end)
        return np.mean([at(x)==at(x+lag) for x in xs])
    assert same_fraction(high)<same_fraction(low)

def test_local_tree_labels_and_nonnegative_branches():
    from msrcsim.species_tree import _Parser
    r=run(1e-5,9)
    for local in r.trees:
        root=_Parser(local.newick).parse()
        tips=[n.name for n in _walk(root) if n.is_tip()]
        assert sorted(tips)==["A","B","C","D"]
        assert all(n.length>=-1e-9 for n in _walk(root))
        assert local.end>local.start

def test_recombination_outcomes_preserve_no_change_intervals():
    r=run(1e-5,9)
    kinds={e.event_type for e in r.events}
    assert {"no_genealogy_change","branch_time_change","topology_change"}.issubset(kinds)
    assert len(r.trees)==len(r.events)+1
    for left,right,event in zip(r.trees,r.trees[1:],r.events):
        assert left.end==pytest.approx(right.start)
        if event.event_type=="no_genealogy_change":
            assert left.newick==right.newick

def test_smc_topology_marginal_is_checked_against_independent_msc():
    # Deliberately statistical rather than assuming a Markov transition
    # preserves the MSC marginal. One test seed keeps this fast and repeatable.
    from msrcsim.structured_coalescent import simulate_msc_genealogy
    tr=tree(); rng=np.random.default_rng(9142)
    independent_draws=[simulate_msc_genealogy(i,tr,rng,False) for i in range(600)]
    independent=[g.topology_index for g in independent_draws]
    linked=simulate_msc_smc(tr,SMCConfig(30000,PiecewiseConstantMap.constant(30000,2e-6),9143),np.random.default_rng(9143))
    # Compare quartet split frequencies on a tolerant Monte Carlo scale.
    def freq(values, weights=None):
        weights=np.ones(len(values)) if weights is None else np.asarray(weights)
        return np.array([weights[np.asarray(values)==k].sum()/weights.sum() for k in range(3)])
    from msrcsim.species_tree import _Parser
    pairs=[frozenset(("A","B")),frozenset(("A","C")),frozenset(("A","D"))]
    def linked_code(nw):
        root=_Parser(nw).parse(); clades=[]
        def visit(n):
            if n.is_tip(): return frozenset((n.name,))
            d=frozenset().union(*(visit(c) for c in n.children))
            if len(d)==2: clades.append(d)
            return d
        visit(root)
        pair=clades[0]
        return pairs.index(pair if pair in pairs else frozenset(set(("A","B","C","D"))-set(pair)))
    widths=[t.end-t.start for t in linked.trees]
    a,b=freq(independent),freq([linked_code(t.newick) for t in linked.trees],widths)
    assert np.max(np.abs(a-b)) < 0.25, (a,b,len(linked.trees))
    independent_root=np.mean([max(g.coalescence_times) for g in independent_draws])
    from msrcsim.species_tree import _Parser
    linked_root=np.average([SpeciesTree(t.newick,100,1000,ultrametric_tolerance=1e-5).root.age for t in linked.trees],weights=widths)
    assert abs(linked_root-independent_root) < 0.25*independent_root, (linked_root,independent_root)

def _walk(n):
    yield n
    for c in n.children: yield from _walk(c)
