import numpy as np
from msrcsim.analytic import compute_Hm,asymmetry_1010

def test_formula():
    H,s=compute_Hm(1,.05,.05,1,1); row=H[s.index((1,0,1,0))]
    assert np.isclose(row[1]-row[2],asymmetry_1010(1,.05,1))
