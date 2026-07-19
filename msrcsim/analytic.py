from __future__ import annotations
from itertools import product
from typing import Dict, List, Tuple
import numpy as np
from scipy.linalg import expm

Config = Tuple[int, int, int, int]
TOPOLOGY_NAMES = ("12|34", "13|24", "14|23")
TOPOLOGY_PAIRS: Dict[int, Tuple[Tuple[int, int], Tuple[int, int]]] = {
    0: ((0, 1), (2, 3)),
    1: ((0, 2), (1, 3)),
    2: ((0, 3), (1, 2)),
}

def all_configurations() -> List[Config]:
    return list(product((0, 1), repeat=4))

def parse_configuration(value: str | Config) -> Config:
    if isinstance(value, tuple): z=value
    else:
        if len(value)!=4 or any(c not in "01" for c in value):
            raise ValueError("configuration must be a four-character binary string")
        z=tuple(int(c) for c in value)
    if len(z)!=4 or any(x not in (0,1) for x in z):
        raise ValueError(f"invalid configuration {z}")
    return z  # type: ignore[return-value]

def _flip(z: Config, i: int) -> Config:
    x=list(z); x[i]=1-x[i]; return tuple(x)  # type: ignore[return-value]

def _pair_rate(z: Config, i: int, j: int, l0: float, l1: float) -> float:
    if z[i]!=z[j]: return 0.0
    return l0 if z[i]==0 else l1

def build_generator_and_absorption(m01: float,m10: float,lambda0: float,lambda1: float):
    for n,v in {"m01":m01,"m10":m10,"lambda0":lambda0,"lambda1":lambda1}.items():
        if v<0: raise ValueError(f"{n} must be nonnegative")
    states=all_configurations(); idx={z:i for i,z in enumerate(states)}
    Q=np.zeros((16,16)); R=np.zeros((16,3))
    for row,z in enumerate(states):
        st=0.0
        for i in range(4):
            rate=m01 if z[i]==0 else m10
            if rate>0: Q[row,idx[_flip(z,i)]]+=rate; st+=rate
        ct=0.0
        for top,pairs in TOPOLOGY_PAIRS.items():
            rate=sum(_pair_rate(z,i,j,lambda0,lambda1) for i,j in pairs)
            R[row,top]=rate; ct+=rate
        Q[row,row]=-(st+ct)
    return Q,R,states

def compute_Hm(t: float,m01: float,m10: float,lambda0: float,lambda1: float):
    if t<0: raise ValueError("t must be nonnegative")
    Q,R,states=build_generator_and_absorption(m01,m10,lambda0,lambda1)
    block=np.zeros((19,19)); block[:16,:16]=Q; block[:16,16:]=R
    E=expm(block*t); survival=E[:16,:16]; absorbed=E[:16,16:]
    H=absorbed + survival@np.ones((16,1))@np.full((1,3),1/3)
    H[np.abs(H)<1e-14]=0; H=np.clip(H,0,1)
    if not np.allclose(H.sum(axis=1),1,atol=1e-10):
        raise RuntimeError("rows of H_m do not sum to one")
    return H,states

def conditional_distribution(configuration: str|Config, **kwargs) -> np.ndarray:
    z=parse_configuration(configuration); H,states=compute_Hm(**kwargs)
    return H[states.index(z)].copy()

def asymmetry_1010(t: float,m: float,lam: float) -> float:
    if t<0 or m<0 or lam<=0: raise ValueError("require t>=0,m>=0,lambda>0")
    return lam/(lam+2*m)*(1-np.exp(-2*(lam+2*m)*t))
