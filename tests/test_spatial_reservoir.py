import numpy as np
import pytest

from msrcsim.spatial_reservoir import ReservoirBackground, simulate_reservoir_arg
from msrcsim.arg.engine import ARGConfig, simulate_arg
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.arg.segments import Segment
from msrcsim.arg.output.tree_sequence import build_tree_sequence
from msrcsim.species_tree import SpeciesTree


class History:
    structured_origin_age = 100.0
    def __init__(self, p=.5): self.p = p
    def frequency_at(self, population, time): return self.p
    def next_frequency_boundary(self, population, time): return float('inf')
    def terminal_frequency(self, taxon): return self.p


@pytest.mark.parametrize('p', [.2, .5, .8])
@pytest.mark.parametrize('h', [0, 1])
def test_single_focal_generator_includes_exchange_after_marker_separation(p, h):
    model = ReservoirBackground(History(p), 0, (0, 100), .1)
    rec_map = Map.constant(100, .01)
    # No stored marker: this case lost exchange in the previous ARG.
    active = ((80., 90.),); segments = (Segment(80, 90, 0),)
    same = p if h else 1-p; cross = 1-same
    assert model.recombination_hazard('ROOT', h, active, segments, 0, rec_map) == pytest.approx(.9*(same+.1*cross))
    # Of the event measure, cross-parent events in [0,85) switch the focal 85.
    assert sum((b-a)*r*k for a,b,r,k in model._pieces(rec_map,0,85))*cross == pytest.approx(.085*cross)


def test_origin_and_impossible_state_boundary():
    model = ReservoirBackground(History(0), 0, (0, 100), 1)
    with pytest.raises(ValueError, match='absent'): model.coalescence_rate('ROOT',1,1,20,99.99)
    assert model.coalescence_rate('ROOT',1,0,20,100) == .025
    assert model.coalescence_rate('ROOT',0,0,20,99.99) == .025
    history = History(); history.structured_origin_age = None
    with pytest.raises(ValueError, match='origin'): ReservoirBackground(history,0,(0,100))


def test_zero_recombination_single_tree_and_reproducible():
    tree=SpeciesTree('((A:2,B:2)AB:2,(C:2,D:2)CD:2)ROOT;',20,200)
    kwargs=dict(tree=tree,rec_map=Map.constant(100,0),history=History(),marker=0,
                suppressed_region=(0,100),sampled_backgrounds=dict(A=0,B=1,C=0,D=1),seed=41)
    a=simulate_reservoir_arg(**kwargs); b=simulate_reservoir_arg(**kwargs)
    assert a.events==b.events and a.edges==b.edges
    assert build_tree_sequence(a).simplify().num_trees==1
    assert not any(e.event_type in ('background_transfer','recombination') for e in a.events)


def test_fixed_present_background_has_neutral_two_tip_root_distribution():
    tree=SpeciesTree('(A:0,B:0)ROOT;',20,200); rec_map=Map.constant(100,.0001)
    model=ReservoirBackground(History(1),0,(0,100),1,dict(A=1,B=1))
    n=1200; older=0
    for seed in range(n):
        a=simulate_arg(tree,ARGConfig(100,rec_map,model),np.random.default_rng(5000+seed))
        t=build_tree_sequence(a).at(75)
        older += t.time(t.root)>40
    target=np.exp(-1)
    assert abs(older/n-target)<4*np.sqrt(target*(1-target)/n)


@pytest.mark.parametrize('p', [.2, .5, .8])
def test_sampled_joint_marks_produce_target_focal_switch_probability(p):
    model=ReservoirBackground(History(p),0,(0,100),.1)
    rec_map=Map.constant(100,.01); active=((80.,90.),)
    segs=(Segment(80,90,0),);rng=np.random.default_rng(4102)
    n=12000; switches=0
    for _ in range(n):
        e=model.sample_recombination('ROOT',0,active,segs,0,rec_map,rng)
        # A focal point at 85 inherits the distal parent only for y<85.
        switches += e['position']<85 and e['parent_background']==1
    target=(.85*.1*p)/(.9*((1-p)+.1*p))
    assert abs(switches/n-target)<4.5*np.sqrt(target*(1-target)/n)


def test_complete_cross_suppression_still_allows_within_background_split():
    tree=SpeciesTree('(A:0,B:0)ROOT;',20,200)
    result=simulate_reservoir_arg(tree,Map.constant(100,.001),History(),0,(0,100),
                                  kappa=0,seed=17,sampled_backgrounds=dict(A=0,B=0))
    assert any(e.event_type=='recombination' for e in result.events)
    assert not any(e.background_changed for e in result.events)
    assert all(t.num_roots==1 for t in build_tree_sequence(result).trees())


def test_frequency_jump_stranding_lineage_raises_without_resampling():
    class JumpHistory(History):
        def frequency_at(self, population, time): return .5 if time<1 else 1.
        def next_frequency_boundary(self, population, time):
            return 1. if time<1 else float('inf')
    tree=SpeciesTree('(A:2,B:2)ROOT;',20,200)
    with pytest.raises(ValueError,match='absent'):
        simulate_reservoir_arg(tree,Map.constant(100,0),JumpHistory(),0,(0,100),
                               sampled_backgrounds=dict(A=0,B=0))


def test_nonfinite_map_is_rejected_at_prototype_boundary():
    tree=SpeciesTree('(A:2,B:2)ROOT;',20,200)
    with pytest.raises(ValueError,match='map'):
        simulate_reservoir_arg(tree,Map.constant(100,float('nan')),History(),0,(0,100))
