from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np
from .analytic import Config,TOPOLOGY_PAIRS,conditional_distribution,parse_configuration
from .lineages import Lineage
from .species_tree import BalancedQuartetTree
from .wright_fisher import FrequencyHistory

@dataclass(frozen=True)
class ConditionalSimulationResult:
    configuration: Config
    num_loci: int
    counts: np.ndarray
    empirical_probabilities: np.ndarray
    exact_probabilities: np.ndarray
    @property
    def absolute_error(self): return np.abs(self.empirical_probabilities-self.exact_probabilities)

@dataclass(frozen=True)
class GenealogyResult:
    newick: str
    topology_index: int
    coalescence_times: Tuple[float,...]


def _eligible_conditional(z,m01,m10,l0,l1):
    ev=[]
    for i in range(4):
        r=m01 if z[i]==0 else m10
        if r>0: ev.append(("switch",(i,),r))
    for top,pairs in TOPOLOGY_PAIRS.items():
        for i,j in pairs:
            if z[i]==z[j]:
                r=l0 if z[i]==0 else l1
                if r>0: ev.append(("coal",(i,j,top),r))
    return ev

def simulate_one_quartet(configuration,t,m01,m10,lambda0,lambda1,rng):
    z=parse_configuration(configuration); now=0.0
    while now<t:
        ev=_eligible_conditional(z,m01,m10,lambda0,lambda1); total=sum(x[2] for x in ev)
        if total<=0: break
        wait=rng.exponential(1/total)
        if now+wait>=t: break
        now+=wait; x=rng.random()*total; c=0
        for kind,payload,rate in ev:
            c+=rate
            if x<=c:
                if kind=="switch":
                    q=list(z); q[payload[0]]=1-q[payload[0]]; z=tuple(q)
                else: return int(payload[2])
                break
    return int(rng.integers(0,3))

def simulate_conditional_quartets(configuration,t,m01,m10,lambda0,lambda1,num_loci,seed=None):
    if num_loci<=0: raise ValueError("num_loci must be positive")
    z=parse_configuration(configuration); rng=np.random.default_rng(seed); counts=np.zeros(3,dtype=int)
    for _ in range(num_loci): counts[simulate_one_quartet(z,t,m01,m10,lambda0,lambda1,rng)]+=1
    exact=conditional_distribution(z,t=t,m01=m01,m10=m10,lambda0=lambda0,lambda1=lambda1)
    return ConditionalSimulationResult(z,num_loci,counts,counts/num_loci,exact)


def _topology_from_first_cherry(desc: frozenset[str]) -> int:
    s=frozenset(desc)
    if s in (frozenset({"1","2"}),frozenset({"3","4"})): return 0
    if s in (frozenset({"1","3"}),frozenset({"2","4"})): return 1
    if s in (frozenset({"1","4"}),frozenset({"2","3"})): return 2
    raise ValueError(f"cannot identify quartet from first cherry {sorted(s)}")


def _p(history: FrequencyHistory,pop: str,age: float) -> float:
    if pop=="ancestral": return 0.0
    return history.frequency(pop,age)


def _rates(lineages: List[Lineage],tree: BalancedQuartetTree,history: FrequencyHistory,
           age: float,rho: float):
    events=[]
    for i,L in enumerate(lineages):
        p=_p(history,L.population_id,age)
        rate=rho*p if L.arrangement==0 else rho*(1-p)
        if rate>0: events.append(("switch",(i,),rate))
    for i in range(len(lineages)):
        for j in range(i+1,len(lineages)):
            a,b=lineages[i],lineages[j]
            if a.population_id!=b.population_id or a.arrangement!=b.arrangement: continue
            N=tree.effective_size(a.population_id); p=_p(history,a.population_id,age)
            freq=(1-p) if a.arrangement==0 else p
            if freq<=0: continue
            rate=1/(2*N*freq)
            events.append(("coal",(i,j),rate))
    return events


def _next_boundary(lineages: List[Lineage],tree: BalancedQuartetTree,age: float) -> float:
    vals=[]
    for L in lineages:
        if L.population_id=="ancestral": continue
        vals.append(float(tree.boundary_age(L.population_id)))
        # generation boundaries make the WF trajectory piecewise constant
        next_gen=np.floor(age)+1.0
        if next_gen < tree.boundary_age(L.population_id)-1e-12: vals.append(next_gen)
    return min(vals) if vals else float("inf")


def _cross_boundaries(lineages: List[Lineage],tree: BalancedQuartetTree,age: float):
    for L in lineages:
        if L.population_id=="ancestral": continue
        b=tree.branches[L.population_id]
        if age >= b.older_age-1e-10:
            parent=b.parent_id
            L.population_id="ancestral" if parent is None else parent
            if L.population_id=="ancestral": L.arrangement=0
            else:
                p=0.0 if L.population_id=="ancestral" else None
    # If a segregating class has disappeared at the new time, collapse labels.
    for L in lineages:
        if L.population_id=="ancestral": L.arrangement=0


def simulate_genealogy(tree: BalancedQuartetTree,history: FrequencyHistory,
                       terminal_arrangements: Dict[str,int],rho: float,
                       rng: np.random.Generator) -> GenealogyResult:
    lineages=[]
    for i,taxon in enumerate(tree.taxa):
        lineages.append(Lineage(i,frozenset({taxon}),tree.terminal_by_taxon[taxon],
                                int(terminal_arrangements[taxon]),0.0,taxon))
    age=0.0; next_id=4; coal_times=[]; first_topology=None
    while len(lineages)>1:
        boundary=_next_boundary(lineages,tree,age)
        events=_rates(lineages,tree,history,age,rho)
        total=sum(e[2] for e in events)
        wait=rng.exponential(1/total) if total>0 else float("inf")
        if age+wait >= boundary:
            age=boundary; _cross_boundaries(lineages,tree,age); continue
        age+=wait; x=rng.random()*total; c=0.0
        chosen=None
        for e in events:
            c+=e[2]
            if x<=c: chosen=e; break
        if chosen is None: raise RuntimeError("failed to choose event")
        kind,payload,_=chosen
        if kind=="switch":
            L=lineages[payload[0]]; L.arrangement=1-L.arrangement
        else:
            i,j=payload; a=lineages[i]; b=lineages[j]
            if first_topology is None and len(a.descendants)==1 and len(b.descendants)==1:
                first_topology=_topology_from_first_cherry(a.descendants|b.descendants)
            bl_a=age-a.node_age; bl_b=age-b.node_age
            nw=f"({a.newick}:{bl_a:.8f},{b.newick}:{bl_b:.8f})"
            parent=Lineage(next_id,a.descendants|b.descendants,a.population_id,a.arrangement,age,nw)
            next_id+=1; coal_times.append(age)
            for k in sorted((i,j),reverse=True): lineages.pop(k)
            lineages.append(parent)
    root=lineages[0]
    if first_topology is None: raise RuntimeError("no identifiable first cherry")
    return GenealogyResult(root.newick+";",first_topology,tuple(coal_times))
