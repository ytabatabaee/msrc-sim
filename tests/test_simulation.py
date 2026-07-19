import numpy as np
from msrcsim.structured_coalescent import simulate_conditional_quartets

def test_simulation_matches_exact():
    r=simulate_conditional_quartets("1010",1,0.05,0.05,1,1,50000,12345)
    assert np.all(r.absolute_error<0.015)

def test_reproducibility():
    kw=dict(configuration="1010",t=1,m01=0.05,m10=0.05,lambda0=1,lambda1=1,num_loci=1000,seed=7)
    a=simulate_conditional_quartets(**kw); b=simulate_conditional_quartets(**kw)
    assert np.array_equal(a.counts,b.counts)
