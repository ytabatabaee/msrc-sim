"""A dated-tree SMC' transition constrained by a rooted species tree.

Times are generations before present; pairwise coalescence rate is 1/(2Ne).
The transition is the standard SMC' prune/recoalesce kernel, with target edges
eligible only while they share a species-tree population with the detached edge.
"""
from dataclasses import dataclass, field
import copy, json
from pathlib import Path
import numpy as np
from ..species_tree import SpeciesTree, Node
from .recombination_map import PiecewiseConstantMap


@dataclass
class SMCConfig:
    chromosome_length: float
    recombination_map: PiecewiseConstantMap
    seed: int = 1
    sampling: dict = field(default_factory=lambda: {"samples_per_species": 1, "time_generations": 0.0})

@dataclass
class Event:
    position: float
    recombination_time: float
    event_type: str
    detached_descendants: str
    recoal_time: float
    target_descendants: str
    target_edge_kind: str
    recoal_population: str

@dataclass
class LocalTree:
    start: float
    end: float
    newick: str
    topology: str

@dataclass
class SMCResult:
    trees: list[LocalTree] = field(default_factory=list)
    events: list[Event] = field(default_factory=list)
    changes: list[dict] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)


def _clone(n, parent=None):
    x=Node(n.name,n.length,[],parent,n.age)
    x.children=[_clone(c,x) for c in n.children]
    return x

def _desc(n):
    return frozenset([n.name]) if n.is_tip() else frozenset().union(*(_desc(c) for c in n.children))

def _newick(n, parent_age=None):
    if n.is_tip(): body=n.name
    else: body="("+",".join(_newick(c,n.age) for c in n.children)+")"+n.name
    if parent_age is None: return body+";"
    length=parent_age-n.age
    if length < -1e-8:
        raise ValueError(f"negative gene-tree branch length on edge ending at {n.name}: {length}")
    # Preserve dated-tree ultrametricity at large generation scales. Ten-digit
    # formatting accumulated visibly across branches for large Ne simulations.
    return body+":"+format(max(0.0,length),".17g")

def _topology(n):
    if n.is_tip(): return n.name
    return "("+",".join(sorted(_topology(c) for c in n.children))+ ")"

def _nodes(root):
    out=[]
    def rec(n): out.append(n); [rec(c) for c in n.children]
    rec(root); return out

def _species_lca(tree, taxa):
    candidates=[n for n in tree.nodes() if set(taxa).issubset(set(tree.descendants(n.name)))]
    return min(candidates,key=lambda n:n.age).name

def _pop_at(tree, taxa, t):
    node=_species_lca(tree,taxa)
    while tree.branches[node].older_age <= t + 1e-10 and tree.branches[node].parent_branch_id is not None:
        node=tree.branches[node].parent_branch_id
    return node

def _edge_rows(root):
    rows=[]
    def rec(n):
        for c in n.children:
            rows.append((c,n,c.age,n.age,_desc(c)))
            rec(c)
    rec(root)
    return rows

def _detach(root, child, parent):
    parent.children.remove(child); child.parent=None
    if parent.parent is None:
        if len(parent.children)==1:
            new=parent.children[0]; new.parent=None; return new
        return parent
    if len(parent.children)==1:
        gp=parent.parent; sib=parent.children[0]; idx=gp.children.index(parent); gp.children[idx]=sib; sib.parent=gp; parent.parent=None
    return root

def _insert(root, target_child, target_parent, detached, t):
    # New coalescence node lies on target edge at t and receives detached lineage.
    names={n.name for n in _nodes(root)} | {n.name for n in _nodes(detached)}
    i=1
    while f"_smc{i}" in names: i+=1
    join=Node(f"_smc{i}",0.0,[target_child,detached],target_parent,t)
    target_child.parent=join; detached.parent=join
    if target_parent is None:
        root=join
    else:
        idx=target_parent.children.index(target_child); target_parent.children[idx]=join; join.parent=target_parent
    return root

def _transition(root, tree, rng):
    edges=_edge_rows(root)
    lengths=np.array([hi-lo for _,_,lo,hi,_ in edges],float)
    total=float(lengths.sum())
    if total <= 0: raise RuntimeError("genealogy has no eligible branches")
    q=float(rng.random()*total); ei=int(np.searchsorted(np.cumsum(lengths),q,side="right")); ei=min(ei,len(edges)-1)
    child,parent,lo,hi,desc=edges[ei]; tr=lo+float(rng.random())*(hi-lo)
    old=copy.deepcopy(root)
    root=_detach(root,child,parent)
    detached_edge=(child,parent,desc,hi)
    # Track detached lineage and extant guide edges through demographic epochs.
    t=tr; recoal=None; target_desc=frozenset()
    species_root_age=tree.root.age
    while True:
        guide_edges=_edge_rows(root)
        boundaries={v for _,_,lo2,hi2,_ in guide_edges for v in (lo2,hi2) if v>t+1e-10}
        # The original-edge target expires when it reaches its old parent
        # coalescence. If this endpoint is omitted, the current target set is
        # incorrectly held constant beyond that genealogy event, producing
        # excess invisible events and a downward root-age bias.
        if detached_edge[3]>t+1e-10:
            boundaries.add(detached_edge[3])
        boundaries.update(tree.branches[b].older_age for b in tree.branches if tree.branches[b].older_age>t+1e-10)
        boundaries=sorted(boundaries)
        end=boundaries[0] if boundaries else float("inf")
        pop=_pop_at(tree,desc,t)
        eligible=[]
        for ch,pa,l,h,d in guide_edges:
            if l-1e-10 <= t < h-1e-10 and _pop_at(tree,d,t)==pop:
                eligible.append((ch,pa,d,"guide"))
        # SMC' retains the detached branch as a possible re-coalescence target.
        # This is the invisible-event/no-change outcome absent from SMC.
        dch,dpa,ddesc,dend=detached_edge
        if t < dend-1e-10 and _pop_at(tree,ddesc,t)==pop:
            eligible.append((dch,dpa,ddesc,"original_edge"))
        # The single lineage ancestral to the guide-tree root persists in the
        # ancestral species population. It is a recoalescence target, but is
        # not part of finite L_G used to place recombinations.
        if t >= max(species_root_age,root.age)-1e-10 and pop==tree.root.name:
            eligible.append((root,None,_desc(root),"ancestral_stem"))
        rate=len(eligible)/(2*tree.branches[pop].effective_population_size)
        if rate>0:
            wait=float(rng.exponential(1/rate))
            if t+wait<end:
                t+=wait; ch,pa,target_desc,target_kind=eligible[int(rng.integers(len(eligible)))]
                if target_kind=="original_edge":
                    return old,tr,t,desc,target_desc,True,target_kind,pop
                # The cut is on the ancestral edge above this subtree; its
                # root age remains the original child coalescence time.
                root=_insert(root,ch,pa,child,t); recoal=t; break
        t=end
    # Rebuild branch lengths from absolute node ages.
    def lengths_rec(n):
        for c in n.children: c.length=n.age-c.age; lengths_rec(c)
    lengths_rec(root)
    return root,tr,recoal,desc,target_desc,False,target_kind,pop

def simulate_msc_smc(species_tree: SpeciesTree, config: SMCConfig, rng=None):
    """Sample the first MSC tree, then advance under spatially inhomogeneous SMC'."""
    if abs(config.chromosome_length-config.recombination_map.length)>1e-9: raise ValueError("map length must equal chromosome length")
    if any(branch.effective_population_size <= 0 for branch in species_tree.branches.values()):
        raise ValueError("all effective population sizes must be positive")
    if config.sampling.get("samples_per_species",1)!=1 or config.sampling.get("time_generations",0.0)!=0.0:
        raise NotImplementedError("the initial SMC module supports one contemporaneous sample per species")
    rng=np.random.default_rng(config.seed) if rng is None else rng
    from ..structured_coalescent import simulate_msc_genealogy
    g=simulate_msc_genealogy(0,species_tree,rng,False)
    # Reuse the project's ultrametric Newick parser; the sampled MSC tree has
    # one lineage per extant species and hence is ultrametric in generations.
    root=SpeciesTree(g.newick, species_tree.default_ne, species_tree.root_extension,
                     ultrametric_tolerance=1e-5).root
    x=0.0; start=0.0; result=SMCResult()
    while x<config.chromosome_length:
        hazard=float(rng.exponential())
        # The map gives r(x) per bp per generation. Across a genealogy with
        # total eligible branch length L_G the spatial hazard is r(x)*L_G.
        gene_length=sum(hi-lo for _,_,lo,hi,_ in _edge_rows(root))
        if gene_length<=0: raise RuntimeError("genealogy has no eligible branch length")
        pos=config.recombination_map.advance_hazard(x,hazard/gene_length)
        if pos>=config.chromosome_length:
            result.trees.append(LocalTree(start,config.chromosome_length,_newick(root),_topology(root))); break
        result.trees.append(LocalTree(start,pos,_newick(root),_topology(root)))
        oldtop=_topology(root)
        root,tr,tcoal,desc,target,unchanged,target_kind,recoal_population=_transition(root,species_tree,rng)
        newtop=_topology(root)
        kind="no_genealogy_change" if unchanged else ("topology_change" if oldtop!=newtop else "branch_time_change")
        ev=Event(pos,tr,kind,",".join(sorted(desc)),tcoal,",".join(sorted(target)),target_kind,recoal_population)
        result.events.append(ev)
        result.changes.append({"position":pos,"event_type":kind,"before_topology":oldtop,"after_topology":newtop})
        x=pos; start=pos
    result.metadata={"model":"multispecies_SMC_prime","seed":config.seed,"time_unit":"generations","pairwise_coalescence_rate":"1/(2Ne)","sampling":config.sampling,"species_tree":species_tree.root.name,"species_tree_newick":_newick(species_tree.root),"population_sizes_by_species_branch":{bid:branch.effective_population_size for bid,branch in species_tree.branches.items()},"root_extension_generations":species_tree.root_extension,"chromosome_length":config.chromosome_length,"recombination_map":config.recombination_map.to_json(),"transition":"SMC' detachment and conditional recoalescence, target edges constrained to shared species-tree population","limitations":["SMC' approximation, not exact ARG","no migration, serial sampling, or background state","one sample per extant species","genealogy edge embedding is inferred from descendant species and divergence times"]}
    return result

def write_result(result, output_dir):
    out=Path(output_dir); out.mkdir(parents=True,exist_ok=True)
    import csv
    with (out/"local_genealogies.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["start","end","newick","topology"]);w.writeheader();w.writerows([t.__dict__ for t in result.trees])
    (out/"local_gene_trees.nwk").write_text("".join(f"[&start={t.start:.10g},end={t.end:.10g}] {t.newick}\n" for t in result.trees))
    with (out/"recombination_events.csv").open("w",newline="") as f:
        rows=[e.__dict__ for e in result.events]; w=csv.DictWriter(f,fieldnames=list(Event.__dataclass_fields__));w.writeheader();w.writerows(rows)
    with (out/"genealogy_changes.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["position","event_type","before_topology","after_topology"]);w.writeheader();w.writerows(result.changes)
    (out/"metadata.json").write_text(json.dumps(result.metadata,indent=2)+"\n")
