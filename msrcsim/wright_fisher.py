from __future__ import annotations
from dataclasses import dataclass, field
from bisect import bisect_right
import numpy as np
from .species_tree import SpeciesTree
from .rearrangement import Rearrangement
from .speciation_split import SpeciationSplitPolicy, sample_daughter_initialization

@dataclass(frozen=True)
class FrequencyRecord:
    rearrangement_id:str; branch_id:str; parent_branch_id:str|None; forward_generation:int; absolute_age:float
    copy_count_A1:float; copy_count_A0:float; population_chromosomes:int; frequency_A1:float; frequency_A0:float
    status:str; is_origin:bool; is_branch_start:bool; is_branch_end:bool; selection_coefficient:float; effective_population_size:int

@dataclass
class FrequencyHistory:
    records:list[FrequencyRecord]
    by_branch:dict[str,list[FrequencyRecord]]
    split_events:list = field(default_factory=list)
    _sorted_by_branch:dict[str,list[FrequencyRecord]] = field(default_factory=dict, init=False, repr=False)
    _ages_by_branch:dict[str,list[float]] = field(default_factory=dict, init=False, repr=False)

    def _ascending(self, branch_id):
        if branch_id not in self._sorted_by_branch:
            self._sorted_by_branch[branch_id] = sorted(self.by_branch[branch_id], key=lambda r: r.absolute_age)
            self._ages_by_branch[branch_id] = [r.absolute_age for r in self._sorted_by_branch[branch_id]]
        return self._sorted_by_branch[branch_id]

    def frequency_at(self, branch_id, age):
        """Return the piecewise-constant frequency at backward age ``age``.

        A Wright-Fisher record at integer age g describes the population during
        the backward-time interval [g, g+1). Beyond the oldest recorded root
        generation, the oldest state is held constant.
        """
        recs = self._ascending(branch_id)
        ages = self._ages_by_branch[branch_id]
        return recs[max(0, bisect_right(ages, age + 1e-12) - 1)].frequency_A1

    def next_frequency_boundary(self, branch_id, age):
        """Smallest recorded age strictly greater than ``age``."""
        self._ascending(branch_id)
        ages = self._ages_by_branch[branch_id]
        index = bisect_right(ages, age + 1e-10)
        return ages[index] if index < len(ages) else float('inf')

    def terminal_frequency(self,taxon): return self.frequency_at(taxon,0.0)

def _status(k,total,ever):
    if not ever: return 'not_present'
    if k==0: return 'lost'
    if k==total: return 'fixed'
    return 'segregating'

def _selected_p(p,s):
    if s==0: return p
    return p*(1+s)/(1+s*p)

def simulate_frequency_history(tree:SpeciesTree,rearrangement:Rearrangement,rng:np.random.Generator, split_policy:SpeciationSplitPolicy|None=None):
    if rearrangement.origin_branch not in tree.branches: raise ValueError("Unknown origin branch")
    records=[]; by={k:[] for k in tree.branches}
    # branch state at its older end; None means inversion absent
    starts={k:None for k in tree.branches}
    ob=tree.branches[rearrangement.origin_branch]
    origin_age=ob.older_age-rearrangement.origin_time_from_branch_start
    if not (ob.younger_age<=origin_age<=ob.older_age): raise ValueError("Origin time falls outside origin branch")
    # process branches from older to younger
    order=sorted(tree.branches.values(), key=lambda b:b.older_age, reverse=True)
    end_counts={}
    split_events=[]
    for b in order:
        total=2*b.effective_population_size
        if b.branch_id==rearrangement.origin_branch:
            start_count=None
        else:
            parent=end_counts.get(b.parent_branch_id,None)
            if split_policy is None:
                if parent is None: start_count=None
                else: start_count=int(rng.binomial(total,parent[0]/parent[1]))
            else:
                parent_count = None if parent is None else parent[0]
                parent_total = None if parent is None else parent[1]
                start_count,event=sample_daughter_initialization(tree=tree,parent_branch_id=b.parent_branch_id,child_branch_id=b.branch_id,parent_count=parent_count,parent_total=parent_total,child_total=total,rng=rng,policy=split_policy)
                split_events.append(event)
        ages=list(np.arange(b.older_age,b.younger_age-1e-9,-1.0))
        if ages[-1]!=b.younger_age: ages.append(b.younger_age)
        k=start_count; ever=k is not None and k>0
        for gi,age in enumerate(ages):
            s = b.selection_coefficient if b.selection_coefficient != 0.0 else rearrangement.selection_coefficient
            is_origin=False
            if b.branch_id==rearrangement.origin_branch and age<=origin_age and k is None:
                k=min(rearrangement.initial_copy_count,total); ever=True; is_origin=True
            if k is None: kk=0; status='not_present'
            else: kk=k; status='newly_originated' if is_origin else _status(kk,total,ever)
            rec=FrequencyRecord(rearrangement.rearrangement_id,b.branch_id,b.parent_branch_id,gi,float(age),kk,total-kk,total,kk/total,1-kk/total,status,is_origin,gi==0,gi==len(ages)-1,s,b.effective_population_size)
            records.append(rec); by[b.branch_id].append(rec)
            if gi<len(ages)-1 and k is not None and 0<k<total:
                p=_selected_p(k/total,s); k=int(rng.binomial(total,p))
        end_counts[b.branch_id]=None if k is None else (k,total)
    return FrequencyHistory(records,by,split_events)
