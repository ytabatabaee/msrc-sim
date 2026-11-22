"""Tests of the thin theory adapter, not another implementation of the theorem."""
import copy
import json
import os
from pathlib import Path
import numpy as np
import pytest
from validation.moran_spatial.theory import load_theory,make_tree,predict,rho_effective,freeze_history
from validation.moran_spatial.run import prescribed,moran_history,assess,tree_topology
from msrcsim.moran import moran_rates
from msrcsim.species_tree import SpeciesTree

CONFIG=Path(__file__).resolve().parents[1]/'validation/moran_spatial/cases.json'
CASES=json.loads(CONFIG.read_text())['cases']

@pytest.fixture(scope='module')
def theory():
    source=os.environ.get('MSRC_THEORY_SOURCE')
    if not source:pytest.skip('set MSRC_THEORY_SOURCE to the external exact-theory checkout')
    return load_theory(source)[0]


def test_piecewise_marker_path_integral():
    c=copy.deepcopy(CASES[0]);c['kappa']=.1
    assert rho_effective(c,50)==pytest.approx(20*.0002+30*.0003*.1)
    assert rho_effective(c,95)==pytest.approx(20*.0002+60*.0003*.1+15*.0001)
    c['marker']=50
    assert rho_effective(c,5)==pytest.approx(15*.0002+30*.0003*.1)
    assert rho_effective(c,50)==0
    c['M']=9
    with pytest.raises(ValueError,match='even'):make_tree(c)
    c['M']=40;c['origin_age']=1
    with pytest.raises(ValueError,match='origins'):make_tree(c)


@pytest.mark.parametrize('index',[0,1,2,3,6,7])
def test_joint_normalization_and_weighted_conditionals(theory,index):
    c=CASES[index];h=prescribed(c,make_tree(c));r=predict(theory,c,.01,h)
    assert np.all(r.joint>=-1e-12)
    assert r.joint.sum()==pytest.approx(1,abs=1e-10)
    assert np.allclose(r.joint.sum(axis=0),r.quartet_probabilities,atol=1e-10)
    good=r.pattern_probabilities>0
    assert np.allclose((r.conditional_quartet_probabilities[good]*r.pattern_probabilities[good,None]).sum(axis=0),r.quartet_probabilities,atol=1e-10)


@pytest.mark.parametrize('index',[8,9,10])
def test_unconditional_full_moran_recursion(theory,index):
    c=CASES[index];r=predict(theory,c,.017)
    assert r.joint.sum()==pytest.approx(1,abs=1e-10)
    assert np.allclose(r.joint.sum(axis=1),r.pattern_probabilities,atol=1e-10)
    if c['shape']=='pectinate':
        assert r.early_topology_joint[:,1:].sum()>0 # A3 absorption before taxon D enters.


@pytest.mark.parametrize('selection',[-.2,0,.2])
def test_moran_rates_identical(theory,selection):
    for M in (8,20,40):
        for k in range(M+1):
            a=moran_rates(k,M,selection);b=theory['moran_quartet_operator'].moran_rates(k,M,selection)
            assert (a.q_plus,a.q_minus)==(b.q_plus,b.q_minus)


def test_freeze_history_preserves_intervals():
    c=CASES[6];t=make_tree(c);h=moran_history(c,t,np.random.default_rng(c['history_seed']));f=freeze_history(t,h)
    for bid,b in t.branches.items():
        for u in np.linspace(b.younger_age+1e-8,b.older_age-1e-8,31):
            assert f.frequency_at(bid,u)==h.frequency_at(bid,u)
        assert sum(d for p,d in f.segments[bid])==pytest.approx(b.older_age-b.younger_age)


def test_failure_cannot_produce_survivor_only_equivalence():
    rows=assess([1000,1000,1000],np.ones(3)/3,.05,failed=1)
    assert all(r['status']=='FAIL_incompatible' and r['observed'] is None for r in rows)
    assert all(r['status']=='PASS' for r in assess([1000,1000,1000],np.ones(3)/3,.05,0))
    assert any(r['status']=='FAIL_uncertainty' for r in assess([10,10,10],np.ones(3)/3,.05,0))
    assert all(r['status']=='NOT_DEFINED_zero_probability_pattern' for r in assess([0,0,0],np.full(3,np.nan),.05,0))


@pytest.mark.parametrize('newick,expected',[('((A:1,B:1):2,(C:2,D:2):1);',0),('((A:1,C:1):2,(B:2,D:2):1);',1),('(((B:1,C:1):1,A:2):1,D:3);',2)])
def test_unrooted_classification_including_pectinate(newick,expected):
    assert tree_topology(SpeciesTree(newick,20,1))==expected


def test_joint_output_does_not_treat_fixed_patterns_as_prior_samples(tmp_path):
    from validation.moran_spatial.aggregate import aggregate
    from validation.moran_spatial.run import write_csv
    import csv
    dest=tmp_path/'case';dest.mkdir()
    (dest/'summary.json').write_text(json.dumps(dict(case='case',n=4,failures={'ARG':0})))
    (dest/'case.json').write_text(json.dumps(dict(sampling='1010')))
    write_csv(dest/'draws.csv',[dict(method='ARG',position=5,pattern='1010',topology=q) for q in (0,0,1,2)])
    write_csv(dest/'comparisons.csv',[dict(method='ARG',pattern='all',status='PASS')])
    (dest/'theory.json').write_text(json.dumps([dict(position=5,joint=(np.ones((16,3))/48).tolist())]))
    aggregate(tmp_path)
    rows=list(csv.DictReader((tmp_path/'joint_pattern_comparisons.csv').open()))
    r=next(r for r in rows if r['pattern']=='1010' and r['topology']=='1')
    assert r['joint_empirical']=='' and float(r['conditional_empirical'])==.5
    (dest/'case.json').write_text(json.dumps(dict(sampling='unconditional')))
    aggregate(tmp_path)
    rows=list(csv.DictReader((tmp_path/'joint_pattern_comparisons.csv').open()))
    r=next(r for r in rows if r['pattern']=='1010' and r['topology']=='1')
    assert float(r['joint_empirical'])==.5
