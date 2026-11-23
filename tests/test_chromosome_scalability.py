import csv,json
import numpy as np
import pytest
from validation.chromosome_scalability.analysis import windows,tracts,span_ecdf,binomial_interval
from validation.chromosome_scalability.worker import species,validate
from msrcsim.arg.engine import simulate_arg,ARGConfig
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.arg.output.tree_sequence import build_tree_sequence
from msrcsim.arg.output.genealogy_export import export_genealogies
from msrcsim.species_tree import SpeciesTree


def rows():
    return [dict(left=a,right=b,topology=q,root_age=root,pair_mean=pair,quartet_internal_generations=internal,genealogy_id=gid) for a,b,q,root,pair,internal,gid in [(0,20,0,10,5,2,'a'),(20,60,1,20,10,4,'b'),(60,100,1,30,20,6,'c')]]


def test_windows_use_genomic_span_not_number_of_correlated_trees():
    a,b,v=windows(rows(),100,50)
    assert np.allclose(v[:,:3],[[.4,.6,0],[0,1,0]])
    assert np.allclose(v[:,3:],[ [16,8,3.2],[28,18,5.6]])
    lo,hi=binomial_interval(0,1)
    assert lo==0 and hi>.9 # one identical observation cannot imply precision


def test_tract_censoring_and_span_weighted_distributions():
    tr=tracts(rows(),'genealogy_id')
    assert [(r['length'],r['censored']) for r in tr]==[(20,True),(40,False),(40,True)]
    assert len(tracts(rows(),'topology'))==2
    x,y=span_ecdf(rows(),'root_age',True,30,70)
    assert np.allclose(x,[20,30]) and np.allclose(y,[.75,1])


def test_full_tree_and_export_validation_checks_labels_times_and_species(tmp_path):
    import tskit
    t=species(4,20,2,2);r=simulate_arg(t,ARGConfig(1000,Map.constant(1000,1e-5)),np.random.default_rng(17));ts=build_tree_sequence(r).simplify()
    assert ts.time_units=='generations' and validate(ts,t)['valid']
    export_genealogies(ts,tmp_path/'trees.nwk',tmp_path/'trees.csv','chr1')
    rr=list(csv.DictReader((tmp_path/'trees.csv').open()));assert len(rr)==ts.num_trees
    assert float(rr[0]['left'])==0 and float(rr[-1]['right'])==1000
    for row in rr:assert SpeciesTree(row['newick'],20,1).taxa==t.taxa
    tables=tskit.TableCollection(1000)
    for taxon in t.taxa:tables.nodes.add_row(flags=1,time=0,metadata=json.dumps(dict(taxon=taxon)).encode())
    for age in (.5,1,1.5):tables.nodes.add_row(time=age)
    for parent,child in ((4,0),(4,2),(5,1),(5,3),(6,4),(6,5)):tables.edges.add_row(0,1000,parent,child)
    tables.sort()
    with pytest.raises(AssertionError):validate(tables.tree_sequence(),t)
