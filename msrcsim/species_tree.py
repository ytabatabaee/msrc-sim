from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

@dataclass
class Node:
    name: str
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)
    parent: Optional["Node"] = None
    age: float = 0.0
    def is_tip(self): return not self.children

@dataclass(frozen=True)
class PopulationBranch:
    branch_id: str
    parent_branch_id: str|None
    younger_age: float
    older_age: float
    effective_population_size: int
    child_node: str

class _Parser:
    def __init__(self,s): self.s=s.strip().rstrip(';'); self.i=0; self.auto=0
    def parse(self):
        n=self.node();
        if self.i!=len(self.s): raise ValueError(f"Unexpected Newick text at {self.s[self.i:]}")
        return n
    def node(self):
        if self.s[self.i]=='(':
            self.i+=1; kids=[self.node()]
            while self.s[self.i]==',': self.i+=1; kids.append(self.node())
            if self.s[self.i]!=')': raise ValueError("Malformed Newick")
            self.i+=1; name=self.label()
            if not name: self.auto+=1; name=f"N{self.auto}"
            length=self.branch_length(); n=Node(name,length,kids)
            for c in kids: c.parent=n
            return n
        name=self.label()
        if not name: raise ValueError("Tip names are required")
        return Node(name,self.branch_length())
    def label(self):
        j=self.i
        while self.i<len(self.s) and self.s[self.i] not in ':,()': self.i+=1
        return self.s[j:self.i].strip()
    def branch_length(self):
        if self.i<len(self.s) and self.s[self.i]==':':
            self.i+=1; j=self.i
            while self.i<len(self.s) and self.s[self.i] not in ',()': self.i+=1
            return float(self.s[j:self.i])
        return 0.0

class SpeciesTree:
    def __init__(self,newick,default_ne,root_extension,branch_parameters=None):
        self.root=_Parser(newick).parse(); self.default_ne=int(default_ne); self.root_extension=float(root_extension)
        self.branch_parameters=branch_parameters or {}; self._assign_ages(); self._name_nodes(); self.branches=self._make_branches()
        self.taxa=tuple(sorted(n.name for n in self.nodes() if n.is_tip()))
        if len(self.taxa)!=4: raise ValueError("This release supports exactly four sampled taxa")
    def nodes(self):
        out=[]
        def rec(n): out.append(n); [rec(c) for c in n.children]
        rec(self.root); return out
    def _assign_ages(self):
        def height(n):
            if n.is_tip(): n.age=0; return 0.0
            vals=[height(c)+c.length for c in n.children]
            if max(vals)-min(vals)>1e-8: raise ValueError("Species tree must be ultrametric")
            n.age=vals[0]; return n.age
        height(self.root)
    def _name_nodes(self):
        seen=set()
        for n in self.nodes():
            if n.name in seen: raise ValueError(f"Duplicate node name {n.name}")
            seen.add(n.name)
    def _ne(self,bid): return int(self.branch_parameters.get(bid,{}).get('effective_population_size',self.default_ne))
    def _make_branches(self):
        d={}
        for n in self.nodes():
            if n is self.root:
                d[n.name]=PopulationBranch(n.name,None,n.age,n.age+self.root_extension,self._ne(n.name),n.name)
            else:
                parent_id=n.parent.name
                d[n.name]=PopulationBranch(n.name,parent_id,n.age,n.parent.age,self._ne(n.name),n.name)
        return d
    def branch_for_tip(self,taxon): return taxon
    def parent_branch(self,bid): return self.branches[bid].parent_branch_id
    def next_boundary(self,bid): return self.branches[bid].older_age
    def descendants(self,node_name):
        n=next(x for x in self.nodes() if x.name==node_name); out=[]
        def rec(x):
            if x.is_tip(): out.append(x.name)
            else:
                for c in x.children: rec(c)
        rec(n); return tuple(out)
