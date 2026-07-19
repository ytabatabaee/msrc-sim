import numpy as np
from .analytic import asymmetry_1010, compute_Hm

def validate_Hm(H, atol=1e-10):
    if H.shape!=(16,3): raise ValueError(f"Expected (16,3), got {H.shape}")
    if np.any(H<-atol) or np.any(H>1+atol): raise ValueError("Invalid probabilities")
    if not np.allclose(H.sum(axis=1),1,atol=atol): raise ValueError("Rows do not sum to one")

def validate_analytic_asymmetry(t,m,lam,atol=1e-10):
    H,states=compute_Hm(t,m,m,lam,lam); row=H[states.index((1,0,1,0))]
    if not np.isclose(row[1]-row[2],asymmetry_1010(t,m,lam),atol=atol): raise ValueError("Analytic mismatch")
