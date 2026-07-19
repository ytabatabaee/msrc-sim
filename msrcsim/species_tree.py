from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Tuple

@dataclass(frozen=True)
class PopulationBranch:
    branch_id: str
    parent_id: str | None
    children: Tuple[str, ...]
    younger_age: int
    older_age: int
    effective_size: int
    taxon: str | None = None

    @property
    def duration(self) -> int:
        return self.older_age-self.younger_age

    def contains_age(self, age: float) -> bool:
        return self.younger_age <= age < self.older_age

@dataclass(frozen=True)
class BalancedQuartetTree:
    branches: Dict[str, PopulationBranch]
    terminal_by_taxon: Dict[str, str]
    ancestral_population_size: int

    @classmethod
    def create(cls, cherry_age: int, root_age: int, origin_age: int,
               effective_size: int, ancestral_population_size: int | None=None):
        if not (0 < cherry_age < root_age < origin_age):
            raise ValueError("require 0 < cherry_age < root_age < origin_age")
        if effective_size<=0: raise ValueError("effective_size must be positive")
        na=ancestral_population_size or effective_size
        b={
          "root": PopulationBranch("root",None,("left","right"),root_age,origin_age,effective_size),
          "left": PopulationBranch("left","root",("t1","t2"),cherry_age,root_age,effective_size),
          "right": PopulationBranch("right","root",("t3","t4"),cherry_age,root_age,effective_size),
          "t1": PopulationBranch("t1","left",(),0,cherry_age,effective_size,"1"),
          "t2": PopulationBranch("t2","left",(),0,cherry_age,effective_size,"2"),
          "t3": PopulationBranch("t3","right",(),0,cherry_age,effective_size,"3"),
          "t4": PopulationBranch("t4","right",(),0,cherry_age,effective_size,"4"),
        }
        return cls(b,{str(i):f"t{i}" for i in range(1,5)},na)

    @property
    def origin_age(self)->int: return self.branches["root"].older_age
    @property
    def root_age(self)->int: return self.branches["root"].younger_age
    @property
    def taxa(self): return ("1","2","3","4")

    def parent_of(self, branch_id: str) -> str | None:
        return self.branches[branch_id].parent_id

    def boundary_age(self, branch_id: str) -> int:
        return self.branches[branch_id].older_age

    def effective_size(self, branch_id: str | None) -> int:
        if branch_id is None or branch_id=="ancestral": return self.ancestral_population_size
        return self.branches[branch_id].effective_size
