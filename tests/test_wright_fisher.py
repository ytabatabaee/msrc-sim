import numpy as np
from msrcsim.species_tree import BalancedQuartetTree
from msrcsim.rearrangement import Rearrangement
from msrcsim.wright_fisher import simulate_frequency_history

def test_absorbing_boundaries():
    tree=BalancedQuartetTree.create(2,4,8,20)
    rng=np.random.default_rng(1)
    h=simulate_frequency_history(tree,Rearrangement(initial_copies=0),rng)
    assert all(np.all(b.counts_forward==0) for b in h.branches.values())

def test_neutral_one_step_mean():
    rng=np.random.default_rng(2); n=100; k=60
    draws=rng.binomial(2*n,k/(2*n),size=200000)
    assert abs(draws.mean()-k)<0.1
