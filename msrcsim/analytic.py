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
    if isinstance(value, tuple):
        z = value
    else:
        if len(value) != 4 or any(c not in "01" for c in value):
            raise ValueError(f"Configuration must be a four-character binary string, got {value!r}")
        z = tuple(int(c) for c in value)
    if len(z) != 4 or any(x not in (0, 1) for x in z):
        raise ValueError(f"Invalid configuration: {z!r}")
    return z  # type: ignore[return-value]

def _flip(z: Config, i: int) -> Config:
    out = list(z); out[i] = 1 - out[i]
    return tuple(out)  # type: ignore[return-value]

def _pair_rate(z: Config, i: int, j: int, lambda0: float, lambda1: float) -> float:
    if z[i] != z[j]: return 0.0
    return lambda0 if z[i] == 0 else lambda1

def build_generator_and_absorption(m01: float, m10: float, lambda0: float, lambda1: float):
    for name, value in {"m01":m01,"m10":m10,"lambda0":lambda0,"lambda1":lambda1}.items():
        if value < 0: raise ValueError(f"{name} must be nonnegative, got {value}")
    states = all_configurations(); index = {z:i for i,z in enumerate(states)}
    Q = np.zeros((16,16)); R = np.zeros((16,3))
    for row,z in enumerate(states):
        switch_total = 0.0
        for i in range(4):
            rate = m01 if z[i] == 0 else m10
            if rate > 0:
                Q[row,index[_flip(z,i)]] += rate; switch_total += rate
        coal_total = 0.0
        for topo,pairs in TOPOLOGY_PAIRS.items():
            rate = sum(_pair_rate(z,i,j,lambda0,lambda1) for i,j in pairs)
            R[row,topo] = rate; coal_total += rate
        Q[row,row] = -(switch_total + coal_total)
    return Q,R,states

def compute_Hm(t: float, m01: float, m10: float, lambda0: float, lambda1: float):
    if t < 0: raise ValueError("t must be nonnegative")
    Q,R,states = build_generator_and_absorption(m01,m10,lambda0,lambda1)
    block = np.zeros((19,19)); block[:16,:16] = Q; block[:16,16:] = R
    E = expm(block*t); survival = E[:16,:16]; absorbed = E[:16,16:]
    H = absorbed + survival @ np.ones((16,1)) @ np.full((1,3),1/3)
    H[np.abs(H)<1e-14]=0.0; H=np.clip(H,0,1)
    if not np.allclose(H.sum(axis=1),1.0,atol=1e-10): raise RuntimeError("Rows of H_m do not sum to one")
    return H,states

def conditional_distribution(configuration, t, m01, m10, lambda0, lambda1):
    z = parse_configuration(configuration)
    H,states = compute_Hm(t,m01,m10,lambda0,lambda1)
    return H[states.index(z)].copy()

def asymmetry_1010(t: float, m: float, lam: float) -> float:
    if t < 0 or m < 0 or lam <= 0: raise ValueError("Require t>=0, m>=0, lam>0")
    return lam/(lam+2*m)*(1-np.exp(-2*(lam+2*m)*t))
