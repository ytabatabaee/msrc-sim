from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List
import numpy as np
from .species_tree import BalancedQuartetTree
from .rearrangement import Rearrangement

@dataclass(frozen=True)
class BranchFrequencyHistory:
    branch_id: str
    younger_age: int
    older_age: int
    counts_forward: np.ndarray
    effective_size: int

    @property
    def frequencies_forward(self) -> np.ndarray:
        return self.counts_forward/(2*self.effective_size)

    @property
    def end_frequency(self) -> float:
        return float(self.frequencies_forward[-1])

    def frequency_at_age(self, age: float) -> float:
        if not (self.younger_age <= age <= self.older_age):
            raise ValueError(f"age {age} outside branch {self.branch_id}")
        # Forward index 0 is at older boundary; age decreases forward.
        elapsed=self.older_age-age
        idx=int(np.floor(elapsed))
        idx=max(0,min(idx,len(self.counts_forward)-1))
        return float(self.counts_forward[idx]/(2*self.effective_size))

@dataclass(frozen=True)
class FrequencyHistory:
    branches: Dict[str, BranchFrequencyHistory]
    terminal_frequencies: Dict[str,float]
    origin_age: int

    def frequency(self, branch_id: str, age: float) -> float:
        return self.branches[branch_id].frequency_at_age(age)


def _selected_frequency(p: float, s: float) -> float:
    if s==0: return p
    denom=1+s*p
    if denom<=0: raise ValueError("selection coefficient produces invalid fitness")
    return p*(1+s)/denom


def _simulate_branch(branch_id: str,younger: int,older: int,Ne: int,initial_count: int,
                     s: float,rng: np.random.Generator) -> BranchFrequencyHistory:
    duration=older-younger
    counts=np.empty(duration+1,dtype=np.int64); counts[0]=initial_count
    total=2*Ne
    for g in range(duration):
        k=int(counts[g])
        if k==0 or k==total:
            counts[g+1]=k
        else:
            p=_selected_frequency(k/total,s)
            counts[g+1]=rng.binomial(total,p)
    return BranchFrequencyHistory(branch_id,younger,older,counts,Ne)


def simulate_frequency_history(tree: BalancedQuartetTree,rearrangement: Rearrangement,
                               rng: np.random.Generator) -> FrequencyHistory:
    if rearrangement.origin_branch!="root":
        raise NotImplementedError("milestone 2 currently supports origin_branch='root'")
    root=tree.branches["root"]
    if not (0 <= rearrangement.initial_copies <= 2*root.effective_size):
        raise ValueError("invalid initial copy count")
    hist: Dict[str,BranchFrequencyHistory]={}
    hist["root"]=_simulate_branch("root",root.younger_age,root.older_age,root.effective_size,
                                  rearrangement.initial_copies,rearrangement.selection_coefficient,rng)
    # At each split, descendants are independently founded by binomial draws from parent frequency.
    for parent_id,child_ids in [("root",("left","right")),("left",("t1","t2")),("right",("t3","t4"))]:
        p=hist[parent_id].end_frequency
        for cid in child_ids:
            b=tree.branches[cid]
            initial=rng.binomial(2*b.effective_size,p)
            hist[cid]=_simulate_branch(cid,b.younger_age,b.older_age,b.effective_size,initial,
                                       rearrangement.selection_coefficient,rng)
    terminal={taxon:hist[bid].end_frequency for taxon,bid in tree.terminal_by_taxon.items()}
    return FrequencyHistory(hist,terminal,tree.origin_age)


def sample_terminal_arrangements(history: FrequencyHistory,tree: BalancedQuartetTree,
                                 rng: np.random.Generator) -> Dict[str,int]:
    return {taxon:int(rng.random()<history.terminal_frequencies[taxon]) for taxon in tree.taxa}
