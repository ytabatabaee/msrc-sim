import copy
import numpy as np
import pytest
from msrcsim.species_tree import SpeciesTree, Node
from msrcsim.smc.structured import (Kernel, ConstantHistory, StructuredConfig,
    simulate_structured_smc, background_signature, genealogy_signature)
from msrcsim.smc.engine import SMCConfig, simulate_msc_smc
from msrcsim.smc.recombination_map import PiecewiseConstantMap as Map


def setup(p=.5,k=.1,rate=.0002,seed=33):
    tree=SpeciesTree('((A:2,B:2)AB:2,(C:2,D:2)CD:2)ROOT;',20,200)
    history=ConstantHistory(p,100)
    cfg=StructuredConfig(Map.constant(100,rate),0,(20,80),k,seed)
    return tree,history,cfg,dict(A=0,B=1,C=0,D=1)


def test_neutral_delegation_is_seed_identical():
    t,h,c,s=setup();c.neutral=True
    a=simulate_structured_smc(t,h,c,s)
    b=simulate_msc_smc(t,SMCConfig(100,c.recombination_map,c.seed))
    assert a.trees==b.trees and a.events==b.events


def test_zero_map_shares_genealogy_and_histories():
    t,h,c,s=setup(rate=0)
    r=simulate_structured_smc(t,h,c,s)
    assert len(r.trees)==1 and not r.events
    assert r.metadata['rho_reference']==0


def test_reproducibility_and_tree_population_validity():
    t,h,c,s=setup(rate=.001)
    a=simulate_structured_smc(t,h,c,s);b=simulate_structured_smc(t,h,c,s)
    assert a.trees==b.trees and a.events==b.events
    assert a.trees[0]['start']==0 and a.trees[-1]['end']==100
    for old,new in zip(a.trees,a.trees[1:]):assert old['end']==new['start']
    for row in a.trees:
        gene=SpeciesTree(row['newick'],20,200,ultrametric_tolerance=1e-7)
        assert gene.taxa==t.taxa
        for n in gene.nodes():
            assert n.age>=0
            if n.parent:assert n.parent.age>=n.age
    assert all(np.isfinite(e['log_acceptance_ratio']) for e in a.events)


def test_within_background_recombination_survives_complete_cross_suppression():
    t,h,c,s=setup(k=0,rate=.001);c.suppressed_region=(0,100)
    r=simulate_structured_smc(t,h,c,s)
    assert r.events
    assert all(e['proposed_path_switches']==0 for e in r.events)
    assert not any(e['background_changed'] for e in r.events)
    assert any(e['genealogy_changed'] for e in r.events)


def test_two_tip_augmented_density_matches_manual_survival():
    tree=SpeciesTree('(A:0,B:0)ROOT;',20,200)
    cfg=StructuredConfig(Map.constant(100,0),0,(0,100))
    model=Kernel(tree,ConstantHistory(.5,100),cfg)
    a=Node('A',age=0);b=Node('B',age=0);r=Node('R',children=[a,b],age=10)
    a.parent=r;b.parent=r
    a.bh=b.bh=[(0,0)];r.bh=[(10,0)]
    # Two same-background lineages have pair hazard 1/20 until coalescence.
    assert model.log_density(r)==pytest.approx(-10/20+np.log(1/20))


def test_original_edge_expires_at_parent_and_background_history_boundaries():
    t,h,c,s=setup();model=Kernel(t,h,c)
    rng=np.random.default_rng(9);root=model.initial(s,rng)
    from msrcsim.smc.engine import _edge_rows,_detach
    child,parent,lo,hi,desc=_edge_rows(root)[0]
    sibling=next(n for n in parent.children if n is not child)
    from msrcsim.smc.structured import _join,_state
    sibling.bh=_join(sibling.bh,parent.bh);guide=_detach(root,child,parent)
    before,_=model.eligible(guide,(child.bh,hi),desc,hi-1e-7,_state(child.bh,hi-1e-7,model.origin))
    after,_=model.eligible(guide,(child.bh,hi),desc,hi,_state(child.bh,hi,model.origin))
    assert any(e[2]=='original_edge' for e in before)
    assert not any(e[2]=='original_edge' for e in after)


@pytest.mark.parametrize('p',[.2,.5,.8])
def test_extended_reverse_density_balance_with_background_paths(p):
    t,h,c,s=setup(p=p,rate=.001);c.audit_density=True;kernel=Kernel(t,h,c);rng=np.random.default_rng(505)
    root=kernel.initial(s,rng);seen=set();changed_bg=False
    for _ in range(60):
        root,event=kernel.transition(root,50,rng)
        assert np.isfinite(kernel.log_density(root))
        # For this inherited-state proposal the structured density and
        # conditional path terms analytically cancel. This checks the actual
        # forward/reverse densities, not just empirical marginal agreement.
        assert abs(event['log_acceptance_ratio'])<1e-7
        seen.add(event['target_edge_kind'])
        changed_bg |= event['background_history_changed']
    assert 'original_edge' in seen and changed_bg


def test_analytic_cancellation_fast_path_preserves_seeded_states():
    t,h,c,s=setup(rate=.001)
    fast=simulate_structured_smc(t,h,c,s)
    c.audit_density=True
    audited=simulate_structured_smc(t,h,c,s)
    assert fast.trees==audited.trees
    assert [(e['position'],e['accepted'],e['genealogy_changed'],e['background_changed']) for e in fast.events]==[(e['position'],e['accepted'],e['genealogy_changed'],e['background_changed']) for e in audited.events]


def test_origin_and_absent_background_policy():
    t,h,c,s=setup();kernel=Kernel(t,ConstantHistory(0,100),c)
    assert kernel.pair('ROOT',None,100)==pytest.approx(1/40)
    assert kernel.switch('ROOT',None,100)==0
    with pytest.raises(ValueError,match='absent'):kernel.pair('ROOT',1,99)
    with pytest.raises(ValueError,match='origin'):Kernel(t,ConstantHistory(.5,None),c)


def test_piecewise_history_boundaries_and_rate_integration():
    class History(ConstantHistory):
        def frequency_at(self,pop,t):return .2 if t<10 else .8
        def next_frequency_boundary(self,pop,t):return 10 if t<10 else float('inf')
    t,h,c,s=setup();kernel=Kernel(t,History(.2,100),c)
    assert kernel.next_time(9,{100})==10
    assert kernel.switch('ROOT',0,9)==pytest.approx(kernel.rho*.2)
    assert kernel.switch('ROOT',0,10)==pytest.approx(kernel.rho*.8)
    assert kernel.effective_rate(50)==pytest.approx(.0002*(20+30*.1))
    root=kernel.initial(s,np.random.default_rng(16))
    for _ in range(10):root,e=kernel.transition(root,50,np.random.default_rng(_))
    assert np.isfinite(kernel.log_density(root))


def test_zero_map_plateau_has_no_proposal_events():
    t,h,c,s=setup(rate=.001)
    c.recombination_map=Map(100,(20.,80.),(.001,0.,.001))
    r=simulate_structured_smc(t,h,c,s)
    assert r.events
    assert not any(20<=e['position']<80 for e in r.events)


def test_origin_nextafter_uses_legacy_numeric_boundary_convention():
    t,h,c,s=setup();kernel=Kernel(t,h,c)
    assert kernel.switch('ROOT',None,np.nextafter(100.,0.))==0
    assert kernel.pair('ROOT',None,np.nextafter(100.,0.))==pytest.approx(1/40)


@pytest.mark.parametrize('model',['moran','wright_fisher'])
def test_realized_history_is_shared_without_mutation(model):
    from msrcsim.population_process import simulate_population_history
    from msrcsim.rearrangement import Rearrangement
    t=SpeciesTree('((A:2,B:2)AB:2,(C:2,D:2)CD:2)ROOT;',20,10)
    h=simulate_population_history(t,Rearrangement('b','suppressed_recombination','ROOT',6,20,0),
                                  np.random.default_rng(93),{'population_process':{'model':model}})
    before=repr(h.records)
    c=StructuredConfig(Map.constant(100,.0002),0,(20,80),.1,19500)
    r=simulate_structured_smc(t,h,c,dict(A=0,B=1,C=0,D=1))
    assert r.trees and before==repr(h.records)


def test_single_lineage_switch_count_has_analytic_poisson_mean():
    t,h,c,s=setup(k=1,rate=.001);c.suppressed_region=(0,100)
    kernel=Kernel(t,h,c);rng=np.random.default_rng(89)
    counts=np.array([len(kernel.free_path(frozenset(t.taxa),0,0,rng))-1 for _ in range(2000)])
    target=100*kernel.rho*.5
    assert abs(counts.mean()-target)<4*np.sqrt(target/len(counts))
