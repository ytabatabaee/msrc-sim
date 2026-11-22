"""Version consolidation regression tests captured before moving adapters."""
import copy
import csv
import dataclasses
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from msrcsim.arg import simulate_msrc_arg
from msrcsim.arg.backgrounds import create_background_model, ReservoirBackground
from msrcsim.arg.backgrounds.msrc_two_state import MSRCTwoBackground
from msrcsim.arg.engine import ARGConfig, simulate_arg, run_arg_config
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.arg.segments import Segment
from msrcsim.spatial_reservoir import simulate_reservoir_arg
from msrcsim.spatial_reservoir.model import ReservoirBackground as CompatibilityBackground
from msrcsim.species_tree import SpeciesTree

FIXTURES=Path(__file__).parent/'fixtures'
BASELINE=json.loads((FIXTURES/'arg_consolidation_baseline.json').read_text())


class History:
    def __init__(self,p=.5,origin=100):self.p=p;self.structured_origin_age=origin
    def frequency_at(self,pop,t):return self.p
    def terminal_frequency(self,taxon):return self.p
    def next_frequency_boundary(self,pop,t):return float('inf')


def tree(n=4):
    return SpeciesTree('(A:0,B:0)ROOT;' if n==2 else '((A:2,B:2)AB:2,(C:2,D:2)CD:2)ROOT;',20,200)


def digest(result):
    value=dict(nodes=result.nodes,edges=result.edges,events=[dataclasses.asdict(e) for e in result.events],
               peak_lineages=result.peak_lineages,peak_segments=result.peak_segments)
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


@pytest.mark.parametrize('case',BASELINE['cases'],ids=lambda c:f"{c['version']}-{c['taxa']}-{c['p']}-{c['kappa']}-{c['sampling']}")
def test_exact_preconsolidation_regression(case):
    t=tree(case['taxa']);states={taxon:i%2 for i,taxon in enumerate(t.taxa)}
    origin=None if case['sampling']=='legacy_threshold_missing_origin' else 100
    hist=History(case['p'],origin)
    rec=Map(100,(25.,75.),(.0002,0.,.0005)) if case['sampling']=='explicit' else Map.constant(100,.0002)
    if case['version']=='neutral':r=simulate_arg(t,ARGConfig(100,rec),np.random.default_rng(case['seed']))
    elif case['sampling']=='legacy_threshold_missing_origin':
        model=MSRCTwoBackground(hist,50,(20,80),case['kappa'])
        r=simulate_arg(t,ARGConfig(100,rec,model),np.random.default_rng(case['seed']))
        # New API permits exact reproduction by making the old implicit
        # threshold states explicit; it then consumes no terminal RNG draws.
        threshold={taxon:int(case['p']>=.5) for taxon in t.taxa}
        unified=simulate_msrc_arg(t,rec,hist,50,(20,80),case['kappa'],'legacy',case['seed'],threshold)
        assert digest(unified)==case['digest']
    else:
        r=simulate_msrc_arg(t,rec,hist,50,(20,80),case['kappa'],case['version'],case['seed'],
                            None if case['sampling']=='bernoulli' else states)
    assert digest(r)==case['digest']
    assert len(r.events)==case['events']
    assert r.metadata['background_model_version']==case['version']


@pytest.mark.parametrize('mode',['legacy','reservoir_v1'])
def test_unified_matches_direct_model_with_explicit_states(mode):
    t=tree();hist=History(.37);rec=Map.constant(100,.001);states=dict(A=0,B=1,C=0,D=1)
    model=create_background_model(mode,hist,0,(20,80),.1,states)
    raw=simulate_arg(t,ARGConfig(100,rec,model,sampled_backgrounds=states),np.random.default_rng(83))
    unified=simulate_msrc_arg(t,rec,hist,0,(20,80),.1,mode,83,states)
    assert digest(raw)==digest(unified)
    assert unified.metadata['seed']==83


@pytest.mark.parametrize('states',[None,dict(A=0,B=1,C=0,D=1)])
def test_compatibility_wrapper_and_default_unified_are_exact(states):
    assert CompatibilityBackground is ReservoirBackground
    t=tree();rec=Map.constant(100,.001);hist=History(.37)
    old=simulate_reservoir_arg(t,rec,hist,0,(20,80),.1,83,states)
    new=simulate_msrc_arg(t,rec,hist,0,(20,80),kappa=.1,seed=83,sampled_backgrounds=states)
    assert digest(old)==digest(new) and old.metadata==new.metadata


def test_caller_rng_is_continued_without_reseeding():
    t=tree();rec=Map.constant(100,.0002);hist=History(.5);states=dict(A=0,B=1,C=0,D=1)
    a=np.random.default_rng(83);b=np.random.default_rng(83);a.random(7);b.random(7)
    r=simulate_msrc_arg(t,rec,hist,0,(20,80),sampled_backgrounds=states,rng=a)
    model=ReservoirBackground(hist,0,(20,80),1.,states)
    direct=simulate_arg(t,ARGConfig(100,rec,model,sampled_backgrounds=states),b)
    assert digest(r)==digest(direct)
    assert a.random()==b.random()
    assert r.metadata['seed'] is None and r.metadata['rng_source']=='caller_generator'


@pytest.mark.parametrize('p',[.2,.5,.8])
@pytest.mark.parametrize('h',[0,1])
def test_matching_occupied_rates_and_origin(p,h):
    for mode in ('legacy','reservoir_v1'):
        model=create_background_model(mode,History(p),0,(0,100))
        ph=p if h else 1-p
        assert model.coalescence_rate('ROOT',h,h,20,50)==pytest.approx(1/(40*ph))
        assert model.coalescence_rate('ROOT',h,1-h,20,50)==0
        assert model.coalescence_rate('ROOT',h,1-h,20,100)==pytest.approx(1/40)


def test_marker_absence_policy_is_versioned_not_reintroduced():
    active=((80.,90.),);segments=(Segment(80,90,0),);rec=Map.constant(100,.01)
    legacy=create_background_model('legacy',History(.5),0,(0,100),.1)
    res=create_background_model('reservoir_v1',History(.5),0,(0,100),.1)
    assert legacy._region(active,segments)[1]==((80.,90.),)
    assert res._region(active,segments)[1]==((0.,90.),)
    assert legacy.recombination_hazard('ROOT',0,active,segments,0,rec)==pytest.approx(.1*(.5+.05))
    assert res.recombination_hazard('ROOT',0,active,segments,0,rec)==pytest.approx(.9*(.5+.05))


@pytest.mark.parametrize('p,h',[(0,1),(1,0)])
def test_absent_states_legacy_floor_and_reservoir_error(p,h):
    legacy=create_background_model('legacy',History(p),0,(0,100))
    res=create_background_model('reservoir_v1',History(p),0,(0,100))
    assert legacy.coalescence_rate('ROOT',h,h,20,1)==1
    with pytest.raises(ValueError,match='absent'):res.coalescence_rate('ROOT',h,h,20,1)
    assert res.coalescence_rate('ROOT',h,1-h,20,100)==1/40


def test_positive_rare_background_difference_remains_explicit():
    legacy=create_background_model('legacy',History(.001),0,(0,100))
    res=create_background_model('reservoir_v1',History(.001),0,(0,100))
    assert legacy.coalescence_rate('ROOT',1,1,20,1)==1
    assert res.coalescence_rate('ROOT',1,1,20,1)==25


def test_reservoir_explicit_origin_required_legacy_missing_origin_unchanged():
    h=History(.5,None)
    legacy=create_background_model('legacy',h,0,(0,100))
    assert legacy.coalescence_rate('ROOT',0,1,20,1e6)==0
    with pytest.raises(ValueError,match='origin'):create_background_model('reservoir_v1',h,0,(0,100))


def test_legacy_cli_artifacts_match_preconsolidation_hashes(tmp_path):
    cfg=copy.deepcopy(BASELINE['cli_configuration']);cfg['output']['directory']=str(tmp_path)
    cfg['arg']['frequency_history_path']=str(FIXTURES/'arg_consolidation_frequency_history.csv')
    run_arg_config(cfg)
    for name,expected in BASELINE['cli_hashes'].items():
        assert hashlib.sha256((tmp_path/name).read_bytes()).hexdigest()==expected
    assert json.loads((tmp_path/'metadata.json').read_text())['background_model_version']=='legacy'


@pytest.mark.parametrize('mode,version',[('msrc',None),('msrc','legacy'),('msrc_two_state','reservoir_v1')])
def test_cli_version_selection_and_metadata(tmp_path,mode,version):
    cfg=copy.deepcopy(BASELINE['cli_configuration']);cfg['output']['directory']=str(tmp_path)
    cfg['arg']['background_model']=mode
    cfg['arg']['frequency_history_path']=str(FIXTURES/'arg_consolidation_frequency_history.csv')
    # A frozen history needs no redundant forward rearrangement configuration.
    cfg.pop('rearrangement')
    if version:cfg['arg']['background_model_version']=version
    run_arg_config(cfg)
    expected=version or 'reservoir_v1'
    assert json.loads((tmp_path/'metadata.json').read_text())['background_model_version']==expected
    assert json.loads((tmp_path/'arg_summary.json').read_text())['background_model_version']==expected


def test_unknown_versions_are_explicit_errors():
    with pytest.raises(ValueError,match='version'):
        simulate_msrc_arg(tree(),Map.constant(100,0),History(),0,(0,100),background_model='typo')
    cfg=copy.deepcopy(BASELINE['cli_configuration']);cfg['arg']['background_model']='typo'
    with pytest.raises(ValueError,match='background_model'):run_arg_config(cfg)


@pytest.mark.parametrize('make_seed',[
    lambda:np.int64(51),lambda:[1,2],lambda:np.array([1,2]),
    lambda:np.random.SeedSequence(51),lambda:np.random.default_rng(51),lambda:np.random.PCG64(51)])
def test_seed_metadata_preserves_numpy_seed_inputs(make_seed):
    t=tree();hist=History();rec=Map.constant(100,.0002);states=dict(A=0,B=1,C=0,D=1)
    r=simulate_msrc_arg(t,rec,hist,0,(0,100),seed=make_seed(),sampled_backgrounds=states)
    model=ReservoirBackground(hist,0,(0,100),1.,states)
    direct=simulate_arg(t,ARGConfig(100,rec,model,sampled_backgrounds=states),np.random.default_rng(make_seed()))
    assert digest(r)==digest(direct)
    json.dumps(r.metadata)


@pytest.mark.parametrize('case',BASELINE['legacy_edge_inputs'])
def test_accepted_legacy_map_and_marker_inputs_remain_reproducible(case):
    t=tree();rec=Map(100,tuple(case['breakpoints']),tuple(case['rates']));states=dict(A=0,B=1,C=0,D=1)
    r=simulate_msrc_arg(t,rec,History(),case['marker'],(20,80),.1,'legacy',case['seed'],states)
    assert digest(r)==case['digest']
    # The same bounds/duplicate-map guards used by the original reservoir
    # wrapper remain active for that version, rather than leaking into legacy.
    with pytest.raises(ValueError):
        simulate_msrc_arg(t,rec,History(),case['marker'],(20,80),.1,'reservoir_v1',case['seed'],states)
