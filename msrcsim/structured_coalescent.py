from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .analytic import Config, TOPOLOGY_PAIRS, conditional_distribution, parse_configuration

@dataclass(frozen=True)
class ConditionalSimulationResult:
    configuration: Config
    num_loci: int
    counts: np.ndarray
    empirical_probabilities: np.ndarray
    exact_probabilities: np.ndarray
    @property
    def absolute_error(self): return np.abs(self.empirical_probabilities-self.exact_probabilities)

def _events(z, m01, m10, lambda0, lambda1):
    events=[]
    for i in range(4):
        rate = m01 if z[i]==0 else m10
        if rate>0: events.append(("switch",(i,),rate))
    for topo,pairs in TOPOLOGY_PAIRS.items():
        for i,j in pairs:
            if z[i]==z[j]:
                rate = lambda0 if z[i]==0 else lambda1
                if rate>0: events.append(("coalescence",(i,j,topo),rate))
    return events

def simulate_one_quartet(configuration, t, m01, m10, lambda0, lambda1, rng):
    if t<0: raise ValueError("t must be nonnegative")
    z=parse_configuration(configuration); now=0.0
    while now<t:
        events=_events(z,m01,m10,lambda0,lambda1); total=sum(e[2] for e in events)
        if total<=0: break
        wait=rng.exponential(1/total)
        if now+wait>=t: break
        now += wait; draw=rng.random()*total; cumulative=0.0
        for kind,payload,rate in events:
            cumulative += rate
            if draw<=cumulative:
                if kind=="switch":
                    i=payload[0]; zl=list(z); zl[i]=1-zl[i]; z=tuple(zl)
                else:
                    return int(payload[2])
                break
    return int(rng.integers(0,3))

def simulate_conditional_quartets(configuration, t, m01, m10, lambda0, lambda1, num_loci, seed=None):
    if num_loci<=0: raise ValueError("num_loci must be positive")
    z=parse_configuration(configuration); rng=np.random.default_rng(seed); counts=np.zeros(3,dtype=np.int64)
    for _ in range(num_loci):
        counts[simulate_one_quartet(z,t,m01,m10,lambda0,lambda1,rng)] += 1
    empirical=counts/num_loci
    exact=conditional_distribution(z,t,m01,m10,lambda0,lambda1)
    return ConditionalSimulationResult(z,num_loci,counts,empirical,exact)
