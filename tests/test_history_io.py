import numpy as np
from msrcsim.species_tree import SpeciesTree
from msrcsim.rearrangement import Rearrangement
from msrcsim.wright_fisher import simulate_frequency_history
from msrcsim.history_io import save_frozen_history,load_frozen_history

def test_history_roundtrip(tmp_path):
    tree=SpeciesTree('((1:10,2:10)A:5,(3:10,4:10)B:5)ROOT;',20,20)
    rr=Rearrangement('x','inversion','ROOT',2,10,0.0)
    h=simulate_frequency_history(tree,rr,np.random.default_rng(2))
    path=save_frozen_history(tmp_path/'h.yaml',{'species_tree':{'newick':'((1:10,2:10)A:5,(3:10,4:10)B:5)ROOT;','default_effective_population_size':20,'root_extension':20}},h,{'1':0,'2':1,'3':0,'4':1})
    c,h2,s,m=load_frozen_history(path)
    assert len(h.records)==len(h2.records)
    assert s['2']==1
