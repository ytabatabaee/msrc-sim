import numpy as np
from msrcsim.species_tree import SpeciesTree
from msrcsim.rearrangement import Rearrangement
from msrcsim.wright_fisher import simulate_frequency_history
from msrcsim.structured_coalescent import simulate_genealogy

def test_full_event_logging_and_three_coalescences():
    t=SpeciesTree('(((1:10,2:10)A:5,3:15)B:5,4:20)ROOT;',50,100)
    rng=np.random.default_rng(3); h=simulate_frequency_history(t,Rearrangement('x','inversion','ROOT',5,20,0),rng)
    z={x:0 for x in t.taxa}; g=simulate_genealogy(0,t,h,z,0.01,0.1,rng,True)
    assert len(g.coalescences)==3
    assert list(g.coalescence_times)==sorted(g.coalescence_times)
    assert any(e.event_type=='species_boundary' for e in g.events)
