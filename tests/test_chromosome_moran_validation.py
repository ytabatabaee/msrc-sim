import numpy as np
import pytest
from validation.chromosome_moran_validation.windows import summarize
from msrcsim.moran import moran_rates,simulate_frequency_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.species_tree import SpeciesTree

def test_window_intersections_and_resolution_conservation():
 for width in (250000,50000,25000,10000):
  x,mean=summarize([0,3123456],[3123456,10000000],np.array([[1.,0.],[0.,1.]]),10000000,width)
  assert len(x)==10000000//width
  np.testing.assert_allclose(mean.sum(axis=0),1)
  np.testing.assert_allclose(mean.mean(axis=1),[.3123456,.6876544])
 with pytest.raises(ValueError):summarize([0,5],[4,10],[1,2],10,2)

def test_native_selection_clock_and_absorbing_endpoints():
 M=20;k=7;s=.02;p=k/M;ps=p*(1+s)/(1+s*p);r=moran_rates(k,M,s)
 assert r.q_plus==M/2*ps*(1-p) and r.q_minus==M/2*(1-ps)*p
 assert moran_rates(0,M,s).total==moran_rates(M,M,s).total==0

def test_shared_history_seed_replay_and_origin():
 t=SpeciesTree('((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;',10,10)
 def draw():return simulate_frequency_history(t,Rearrangement('test','suppressed_recombination','ROOT',0,10,0),np.random.default_rng(np.random.SeedSequence([20261123,1,0])))
 a,b=draw(),draw();assert a.records==b.records;assert a.structured_origin_age==12
 assert a.frequency_at('ROOT',12-1e-7)==.5
 assert a.frequency_at('ROOT',12)==0 # older-side absence, not Moran loss
 for pop in t.branches:
  for x in np.linspace(t.branches[pop].younger_age,t.branches[pop].older_age,7):assert a.frequency_at(pop,x)==b.frequency_at(pop,x)
