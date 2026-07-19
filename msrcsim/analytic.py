from __future__ import annotations
from itertools import product
import numpy as np
from scipy.linalg import expm

TOPOLOGY_NAMES=("12|34","13|24","14|23")
TOPOLOGY_PAIRS={0:((0,1),(2,3)),1:((0,2),(1,3)),2:((0,3),(1,2))}

def all_configurations(): return list(product((0,1), repeat=4))

def build_generator_and_absorption(m01,m10,lambda0,lambda1):
    states=all_configurations(); idx={z:i for i,z in enumerate(states)}
    Q=np.zeros((16,16)); R=np.zeros((16,3))
    for row,z in enumerate(states):
        sw=0.0
        for i in range(4):
            rate=m01 if z[i]==0 else m10
            if rate:
                zz=list(z); zz[i]=1-zz[i]; Q[row,idx[tuple(zz)]]+=rate; sw+=rate
        coal=0.0
        for k,pairs in TOPOLOGY_PAIRS.items():
            rr=0.0
            for i,j in pairs:
                if z[i]==z[j]: rr += lambda0 if z[i]==0 else lambda1
            R[row,k]=rr; coal+=rr
        Q[row,row]=-(sw+coal)
    return Q,R,states

def compute_Hm(t,m01,m10,lambda0,lambda1):
    Q,R,states=build_generator_and_absorption(m01,m10,lambda0,lambda1)
    B=np.zeros((19,19)); B[:16,:16]=Q; B[:16,16:]=R
    E=expm(B*t); H=E[:16,16:]+E[:16,:16]@np.ones((16,1))@np.full((1,3),1/3)
    return H,states

def asymmetry_1010(t,m,lam):
    return lam/(lam+2*m)*(1-np.exp(-2*(lam+2*m)*t))
