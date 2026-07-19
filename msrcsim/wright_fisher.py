from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .species_tree import SpeciesTree
from .rearrangement import Rearrangement

@dataclass(frozen=True)
class FrequencyRecord:
    rearrangement_id:str; branch_id:str; parent_branch_id:str|None; forward_generation:int; absolute_age:float
    copy_count_A1:int; copy_count_A0:int; population_chromosomes:int; frequency_A1:float; frequency_A0:float
    status:str; is_origin:bool; is_branch_start:bool; is_branch_end:bool; selection_coefficient:float; effective_population_size:int

@dataclass
class FrequencyHistory:
    records:list[FrequencyRecord]
    by_branch:dict[str,list[FrequencyRecord]]
    def frequency_at(self,branch_id,age):
        recs=self.by_branch[branch_id]
        return min(recs,key=lambda r:abs(r.absolute_age-age)).frequency_A1
    def terminal_frequency(self,taxon): return self.frequency_at(taxon,0.0)

def _status(k,total,ever):
    if not ever: return 'not_present'
    if k==0: return 'lost'
    if k==total: return 'fixed'
    return 'segregating'

def _selected_p(p,s):
    if s==0: return p
    return p*(1+s)/(1+s*p)

def simulate_frequency_history(tree:SpeciesTree,rearrangement:Rearrangement,rng:np.random.Generator):
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
    for b in order:
        total=2*b.effective_population_size
        if b.branch_id==rearrangement.origin_branch:
            start_count=None
        else:
            parent=end_counts.get(b.parent_branch_id,None)
            if parent is None: start_count=None
            else: start_count=int(rng.binomial(total,parent[0]/parent[1]))
        ages=list(np.arange(b.older_age,b.younger_age-1e-9,-1.0))
        if ages[-1]!=b.younger_age: ages.append(b.younger_age)
        k=start_count; ever=k is not None and k>0
        for gi,age in enumerate(ages):
            is_origin=False
            if b.branch_id==rearrangement.origin_branch and age<=origin_age and k is None:
                k=min(rearrangement.initial_copy_count,total); ever=True; is_origin=True
            if k is None: kk=0; status='not_present'
            else: kk=k; status='newly_originated' if is_origin else _status(kk,total,ever)
            rec=FrequencyRecord(rearrangement.rearrangement_id,b.branch_id,b.parent_branch_id,gi,float(age),kk,total-kk,total,kk/total,1-kk/total,status,is_origin,gi==0,gi==len(ages)-1,rearrangement.selection_coefficient,b.effective_population_size)
            records.append(rec); by[b.branch_id].append(rec)
            if gi<len(ages)-1 and k is not None and 0<k<total:
                p=_selected_p(k/total,rearrangement.selection_coefficient); k=int(rng.binomial(total,p))
        end_counts[b.branch_id]=None if k is None else (k,total)
    return FrequencyHistory(records,by)
