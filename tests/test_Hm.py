import numpy as np
from msrcsim.analytic import compute_Hm

def test_Hm_shape_and_probabilities():
    H,states=compute_Hm(1.0,0.05,0.05,1.0,1.0)
    assert H.shape==(16,3); assert len(states)==16
    assert np.all(H>=0) and np.all(H<=1); assert np.allclose(H.sum(axis=1),1)

def test_complements_match_symmetric_rates():
    H,states=compute_Hm(1.2,0.1,0.1,0.8,0.8)
    for z in states:
        c=tuple(1-x for x in z)
        assert np.allclose(H[states.index(z)],H[states.index(c)])
