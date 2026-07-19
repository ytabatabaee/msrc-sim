import numpy as np
from msrcsim.species_tree import SpeciesTree
from msrcsim.rearrangement import Rearrangement
from msrcsim.wright_fisher import simulate_frequency_history

def test_status_and_columns():
    t=SpeciesTree('((1:10,2:10)A:5,(3:10,4:10)B:5)ROOT;',100,20)
    h=simulate_frequency_history(t,Rearrangement('x','inversion','ROOT',5,10,0),np.random.default_rng(1))
    assert {'not_present','newly_originated','segregating','fixed','lost'} & {r.status for r in h.records}
    assert all(0<=r.frequency_A1<=1 for r in h.records)
