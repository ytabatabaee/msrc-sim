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

def test_map_inverse_has_exact_boundary_convention():
    m=PiecewiseConstantMap(10,(2.,4.,6.),(1.,0.,2.,0.))
    assert m.advance_hazard(0,2.0)==pytest.approx(2.0)
    assert m.advance_hazard(0,2.5)==pytest.approx(4.25)
    assert m.advance_hazard(0,6.0)==pytest.approx(6.0)
    assert m.advance_hazard(2,0.0)==pytest.approx(2.0)

def test_large_generation_times_preserve_ultrametricity():
    tr=SpeciesTree("((T1:40000,T2:40000)I:40000,T3:80000)ROOT;",10000,1e7)
    r=simulate_msc_smc(tr,SMCConfig(1000,PiecewiseConstantMap.constant(1000,1.25e-7),520000),np.random.default_rng(520000))
    from msrcsim.species_tree import _Parser
    for local in r.trees:
        root=_Parser(local.newick).parse()
        def height(n):
            if n.is_tip(): return 0.0
            hs=[height(c)+c.length for c in n.children]
            assert max(hs)-min(hs)<1e-8
            return hs[0]
        height(root)

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

class _ForcedRng:
    def __init__(self, edge_draw=0.0, cut_fraction=0.5, wait=0.1, target_index=0):
        self.draws=iter((edge_draw,cut_fraction)); self.wait=wait; self.target_index=target_index
    def random(self): return next(self.draws)
    def exponential(self, scale): return self.wait
    def integers(self, high):
        assert self.target_index < high
        return self.target_index

def test_invisible_recoalescence_uses_original_edge_identity():
    from msrcsim.smc.engine import _transition, _newick
    species=SpeciesTree("(A:0,B:0)ROOT;",100,1000)
    guide=SpeciesTree("(A:2,B:2)G;",100,1000).root
    original=_newick(guide)
    result,tr,tcoal,detached,target,unchanged,target_kind,pop=_transition(guide,species,_ForcedRng(target_index=0))
    assert unchanged
    assert target_kind=="original_edge"
    assert detached==frozenset({"A"}) and target==detached
    assert tcoal>tr and _newick(result)==original
    assert pop=="ROOT"

def test_recoalescence_on_ancestral_stem_changes_time_not_topology():
    from msrcsim.smc.engine import _transition, _newick, _topology
    species=SpeciesTree("(A:0,B:0)ROOT;",100,1000)
    guide=SpeciesTree("(A:2,B:2)G;",100,1000).root
    original_topology=_topology(guide); original_newick=_newick(guide)
    result,tr,tcoal,detached,target,unchanged,target_kind,pop=_transition(guide,species,_ForcedRng(target_index=1))
    assert not unchanged and target==frozenset({"B"})
    assert target_kind=="ancestral_stem"
    assert pop=="ROOT" and tcoal>tr
    assert _topology(result)==original_topology
    assert _newick(result)!=original_newick
    assert all(c.age<=result.age for c in result.children)

def test_original_target_expires_at_its_parent_coalescence():
    from msrcsim.smc.engine import _transition
    species=SpeciesTree("(A:0,B:0)ROOT;",100,1000)
    guide=SpeciesTree("(A:10,B:10)G;",100,1000).root
    class BoundaryRng:
        def __init__(self): self.waits=iter((6.0,.1)); self.draws=iter((0.0,.5))
        def random(self): return next(self.draws)
        def exponential(self,scale): return next(self.waits)
        def integers(self,high): return 0
    result,tr,tcoal,detached,target,unchanged,target_kind,pop=_transition(guide,species,BoundaryRng())
    assert not unchanged
    assert target_kind=="ancestral_stem"
    assert tcoal==pytest.approx(10.1)
    assert pop=="ROOT"

def test_fixed_two_lineage_recoalescence_matches_analytic_kernel():
    import copy
    from msrcsim.smc.engine import _transition
    species=SpeciesTree("(A:0,B:0)ROOT;",100,1000)
    guide=SpeciesTree("(A:10,B:10)G;",100,1000).root
    n=12000; seed=98321; base=np.random.default_rng(seed); original=0; before_boundary=0; times=[]
    class FixedCutRng:
        def random(self): return next(draws)
        def exponential(self,scale): return base.exponential(scale)
        def integers(self,high): return base.integers(high)
    for _ in range(n):
        draws=iter((0.0,.5))
        _,tr,tcoal,_,_,unchanged,target_kind,_=_transition(copy.deepcopy(guide),species,FixedCutRng())
        assert tr==pytest.approx(5.0)
        original+=target_kind=="original_edge"
        before_boundary+=tcoal<10.0
        times.append(tcoal)
    c=1/(2*species.default_ne); span=5.0
    p_original=.5*(1-np.exp(-2*c*span))
    p_before=1-np.exp(-2*c*span)
    assert abs(original/n-p_original)<.008
    assert abs(before_boundary/n-p_before)<.012
    # Survival at age 10 is followed by the single-lineage ancestral-stem
    # hazard c; compare an unconditional CDF point after that boundary.
    point=10.5
    expected=1-np.exp(-2*c*span)*np.exp(-c*(point-10))
    assert abs(np.mean(np.asarray(times)<=point)-expected)<.012

def test_recoalescence_population_matches_both_lineages():
    from msrcsim.smc.engine import _pop_at
    tr=tree(); r=run(1e-5,31)
    for e in r.events:
        detached=frozenset(e.detached_descendants.split(","))
        target=frozenset(e.target_descendants.split(","))
        assert _pop_at(tr,detached,e.recoal_time)==e.recoal_population
        assert _pop_at(tr,target,e.recoal_time)==e.recoal_population

def _walk(n):
    yield n
    for c in n.children: yield from _walk(c)
