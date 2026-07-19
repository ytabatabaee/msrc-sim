import numpy as np, pytest
from msrcsim.analytic import asymmetry_1010, compute_Hm
@pytest.mark.parametrize("t,m,lam",[(0.1,0,1),(0.5,0.05,1),(1,0.1,1),(2,1,0.5),(5,5,2)])
def test_formula(t,m,lam):
    H,s=compute_Hm(t,m,m,lam,lam); row=H[s.index((1,0,1,0))]
    assert np.isclose(row[1]-row[2],asymmetry_1010(t,m,lam),atol=1e-10)
