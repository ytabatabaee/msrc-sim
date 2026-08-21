import numpy as np
from msrcsim.species_tree import SpeciesTree
from msrcsim.rearrangement import Rearrangement
from msrcsim.wright_fisher import simulate_frequency_history
from msrcsim.history_summary import summarize_frequency_history

def test_branch_summaries_have_integrals():
    tree=SpeciesTree('((1:10,2:10)A:5,(3:10,4:10)B:5)ROOT;',20,20)
    h=simulate_frequency_history(tree,Rearrangement('x','inversion','ROOT',2,10,0.0),np.random.default_rng(3))
    rows=summarize_frequency_history(h)
    assert rows and all('integrated_A1_frequency' in x for x in rows)
