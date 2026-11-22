"""Thin adapter to msrc-study's existing operators; no quartet algebra copied."""
from pathlib import Path
import hashlib
import importlib
import sys
import types
import numpy as np
from msrcsim.species_tree import SpeciesTree

SOURCE_FILES = ('moran_quartet_operator.py', 'full_species_tree.py',
                'unbalanced_species_tree.py', 'prescribed_frequency_quartet.py',
                'matrix_free_quartet.py')


def load_theory(root):
    path = Path(root).resolve()
    if (path/'simulations/quartet/moran_theory').is_dir():
        path = path/'simulations/quartet/moran_theory'
    if any(not (path/f).is_file() for f in SOURCE_FILES):
        raise FileNotFoundError('theory source must contain the five original operator modules')
    name = '_msrc_theory_' + hashlib.sha256(str(path).encode()).hexdigest()[:12]
    if name not in sys.modules:
        package = types.ModuleType(name); package.__path__ = [str(path)]
        sys.modules[name] = package
    modules = {f[:-3]: importlib.import_module(name+'.'+f[:-3]) for f in SOURCE_FILES}
    provenance = {f: hashlib.sha256((path/f).read_bytes()).hexdigest() for f in SOURCE_FILES}
    return modules, provenance


def make_tree(case):
    sizes=[case['M'],*case.get('branch_M',{}).values()]
    if any(M<=0 or int(M)!=M or M%2 for M in sizes):
        raise ValueError('simulator requires positive even chromosome counts M=2Ne')
    tree=SpeciesTree(case['newick'],case['M']//2,case['root_duration'],
                     {k:{'effective_population_size':v//2} for k,v in case.get('branch_M',{}).items()})
    if tree.taxa!=('A','B','C','D') or tree.root.name!='ROOT':
        raise ValueError('map theory taxa 1,2,3,4 to A,B,C,D and name the root ROOT')
    expected={'AB':{'A','B'}}
    expected.update({'CD':{'C','D'}} if case['shape']=='balanced' else {'ABC':{'A','B','C'}})
    for bid,desc in expected.items():
        if bid not in tree.branches or set(tree.descendants(bid))!=desc:
            raise ValueError('internal branch labels/clades do not match the declared theory shape')
    if 'origin_age' in case and not np.isclose(case['origin_age'],tree.branches['ROOT'].older_age):
        raise ValueError('this interface compares root-older-boundary origins; set root_duration accordingly')
    return tree


def rho_effective(case, x):
    """Integral of physical base map times cross retention, in generations^-1."""
    m = case['map']; lo,hi = sorted((x,case['marker']))
    cuts = sorted({lo,hi,*[v for v in m['breakpoints']+case['suppressed_region'] if lo<v<hi]})
    def rate(y): return m['rates'][int(np.searchsorted(m['breakpoints'],y,side='right'))]
    a,b = case['suppressed_region']
    return sum((v-u)*rate((u+v)/2)*(case['kappa'] if a <= (u+v)/2 < b else 1)
               for u,v in zip(cuts,cuts[1:]))


class PrescribedHistory:
    """One externally specified young-to-old trajectory on every population."""
    def __init__(self, tree, segments):
        self.tree = tree; self.segments = segments
        self.structured_origin_age = tree.branches[tree.root.name].older_age
        for bid,b in tree.branches.items():
            if not np.isclose(sum(d for p,d in segments[bid]),b.older_age-b.younger_age):
                raise ValueError('trajectory duration differs from population branch: '+bid)
    def frequency_at(self, bid, t):
        end = self.tree.branches[bid].younger_age
        for p,d in self.segments[bid]:
            end += d
            if t < end-1e-12: return p
        return self.segments[bid][-1][0]
    def terminal_frequency(self, taxon): return self.segments[taxon][0][0]
    def next_frequency_boundary(self, bid, t):
        end = self.tree.branches[bid].younger_age
        for p,d in self.segments[bid]:
            end += d
            if end > t+1e-10: return end
        return float('inf')


def freeze_history(tree, history):
    """Condition on every recorded branch interval; no path rejection/resampling."""
    segments = {}
    for bid,b in tree.branches.items():
        cuts = sorted({b.younger_age,b.older_age,*[r.absolute_age for r in history.by_branch[bid]
                         if b.younger_age+1e-10 < r.absolute_age < b.older_age-1e-10]})
        segments[bid] = [(float(history.frequency_at(bid,(a+c)/2)),c-a) for a,c in zip(cuts,cuts[1:])]
    return PrescribedHistory(tree,segments)


def predict(modules, case, rho, history=None):
    """Return original theory result for full balanced/pectinate recursion."""
    t = make_tree(case); b=t.branches
    size=lambda bid:2*b[bid].effective_population_size
    duration=lambda bid:b[bid].older_age-b[bid].younger_age
    common=dict(M_root=size('ROOT'),M_terminal=tuple(size(x) for x in t.taxa),rho=rho)
    if history is None:
        common.update(root_length=duration('ROOT'),terminal_lengths=tuple(duration(x) for x in t.taxa),
                      pi_old=np.array(case['pi_old']),selection=case.get('selection',0.))
        if case['shape']=='balanced':
            return modules['full_species_tree'].full_species_tree_quartet_probabilities(
                **common,M_left=size('AB'),M_right=size('CD'),sister_lengths=(duration('AB'),duration('CD')))
        return modules['unbalanced_species_tree'].unbalanced_species_tree_conditional_quartet_probabilities(
            **common,M_intermediate=size('ABC'),M_sister=size('AB'),intermediate_length=duration('ABC'),sister_length=duration('AB'))
    p=modules['prescribed_frequency_quartet']; trajectory=lambda bid:p.PrescribedTrajectory(tuple(history.segments[bid]))
    common.update(terminal_trajectories=tuple(trajectory(x) for x in t.taxa),root_trajectory=trajectory('ROOT'),with_joint=True)
    if case['shape']=='balanced':
        return p.balanced_prescribed_quartet_probabilities(**common,M_left=size('AB'),M_right=size('CD'),sister_trajectories=(trajectory('AB'),trajectory('CD')))
    return p.unbalanced_prescribed_quartet_probabilities(**common,M_intermediate=size('ABC'),M_sister=size('AB'),sister_trajectory=trajectory('AB'),intermediate_trajectory=trajectory('ABC'))
