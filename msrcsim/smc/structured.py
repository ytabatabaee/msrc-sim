"""Experimental marginal-corrected structured SMC' on a declared window.

Retains a dated tree and branchwise background paths, never ARG material.
The target marginal is frozen at reference_position; spatial approximation
errors away from this coordinate must be measured. See the model document.
"""
from dataclasses import dataclass, field
import copy
import csv
import json
import math
from pathlib import Path

import numpy as np
from ..species_tree import Node
from ..structured_coalescent import simulate_genealogy_reference
from .engine import (_desc, _edge_rows, _nodes, _newick, _topology, _pop_at,
                     _detach, _insert, SMCConfig, simulate_msc_smc)
from .recombination_map import PiecewiseConstantMap

EPS = 1e-10


@dataclass
class StructuredConfig:
    recombination_map: PiecewiseConstantMap
    marker: float
    suppressed_region: tuple[float, float]
    cross_retention: float = 1.0
    seed: int = 1
    start: float = 0.0
    end: float | None = None
    reference_position: float | None = None
    neutral: bool = False
    max_events: int = 100000
    audit_density: bool = False


@dataclass
class StructuredResult:
    trees: list = field(default_factory=list)
    events: list = field(default_factory=list)
    metadata: dict = field(default_factory=dict)


@dataclass
class ConstantHistory:
    p: float
    structured_origin_age: float

    def frequency_at(self, population, time): return self.p
    def terminal_frequency(self, taxon): return self.p
    def next_frequency_boundary(self, population, time): return float('inf')


def _state(hist, t, origin):
    if t >= origin-EPS: return None
    return next(h for age, h in reversed(hist) if age <= t + 1e-12)


def _segment(hist, lo, hi, origin):
    return [(lo, _state(hist, lo, origin))] + [(t, h) for t, h in hist if lo < t < min(hi, origin)]


def _join(a, b):
    out = list(a)
    for t, h in b:
        if not out or out[-1][1] != h:
            out.append((t, h))
    return out


def _rows(root):
    return [(n, p, lo, hi, desc) for n, p, lo, hi, desc in _edge_rows(root)] + [
        (root, None, root.age, float('inf'), _desc(root))]


def genealogy_signature(root):
    return tuple(sorted((tuple(sorted(_desc(n))), n.age) for n in _nodes(root)))


def background_signature(root):
    return tuple(sorted((tuple(sorted(_desc(n))), tuple(n.bh)) for n in _nodes(root)))


def background_jump_signature(root):
    """Actual stochastic background marks, excluding edge-birth annotations."""
    return tuple(sorted((tuple(sorted(_desc(n))),t,h) for n in _nodes(root)
                        for t,h in n.bh[1:] if h in (0,1)))


class Kernel:
    def __init__(self, tree, history, config):
        self.tree, self.history, self.config = tree, history, config
        self._population_lca = {}
        self.origin = float(getattr(history, 'structured_origin_age', float('nan')) or float('nan'))
        if not np.isfinite(self.origin) or self.origin <= 0:
            raise ValueError('an explicit positive finite origin age is required')
        if any(b.effective_population_size <= 0 for b in tree.branches.values()):
            raise ValueError('effective sizes must be positive')
        L = config.recombination_map.length
        end = L if config.end is None else config.end
        if not 0 <= config.start < end <= L: raise ValueError('invalid scan window')
        if not 0 <= config.marker <= L: raise ValueError('marker must lie on map')
        a, b = config.suppressed_region
        if not 0 <= a < b <= L or not 0 <= config.cross_retention <= 1:
            raise ValueError('invalid region of suppressed recombination or cross retention')
        self.reference = (config.start + end) / 2 if config.reference_position is None else config.reference_position
        if not config.start <= self.reference <= end: raise ValueError('reference must lie in scan window')
        self.rho = self.effective_rate(self.reference)
        self.end = end

    def pop(self, desc, t):
        # Species LCA depends on the clade, not on genealogy/background time.
        # Only the ascent through divergence boundaries depends on t.
        if desc not in self._population_lca:
            self._population_lca[desc] = _pop_at(self.tree,desc,0.)
        pop = self._population_lca[desc]
        while (self.tree.branches[pop].older_age <= t+EPS
               and self.tree.branches[pop].parent_branch_id is not None):
            pop = self.tree.branches[pop].parent_branch_id
        return pop

    def kappa(self, x):
        a, b = self.config.suppressed_region
        return self.config.cross_retention if a <= x < b else 1.0

    def effective_rate(self, u):
        lo, hi = sorted((u, self.config.marker))
        cuts = sorted({lo, hi, *[x for x in self.config.recombination_map.breakpoints if lo < x < hi],
                       *[x for x in self.config.suppressed_region if lo < x < hi]})
        return sum((b-a) * self.config.recombination_map.rate_at((a+b)/2) * self.kappa((a+b)/2)
                   for a,b in zip(cuts,cuts[1:]))

    def p(self, pop, t):
        p = float(self.history.frequency_at(pop, t))
        if not np.isfinite(p) or not 0 <= p <= 1: raise ValueError('invalid frequency history')
        return p

    def pair(self, pop, h, t):
        ne = self.tree.branches[pop].effective_population_size
        if t >= self.origin-EPS: return 1/(2*ne)
        p = self.p(pop,t); ph = p if h == 1 else 1-p
        if h not in (0,1) or ph <= 0:
            raise ValueError(f'occupied absent background: strict_v1 needs a boundary bridge (population={pop}, time={t!r}, state={h}, p={p})')
        if ph < 1/(2*ne) - 1e-12:
            raise ValueError('strict_v1 initialization requires occupied frequency >= one copy')
        return 1/(2*ne*ph)

    def switch(self, pop, h, t):
        if t >= self.origin-EPS: return 0.0
        self.pair(pop,h,t)
        p = self.p(pop,t)
        return self.rho * (p if h == 0 else 1-p)

    def next_time(self, t, static, end=float('inf')):
        values = [end] + [x for x in static if x > t + EPS]
        if t < self.origin - EPS:
            values += [self.origin]
            for bid, b in self.tree.branches.items():
                if b.younger_age <= t + EPS and (t < b.older_age - EPS or b.parent_branch_id is None):
                    value = self.history.next_frequency_boundary(bid,t)
                    if value > t + EPS: values.append(value)
        return min(values)

    def static_times(self, rows):
        return {0., self.origin, *[b.older_age for b in self.tree.branches.values()],
                *[v for n,p,lo,hi,d in rows for v in (lo,hi)],
                *[t for n,p,lo,hi,d in rows for t,h in n.bh]}

    def log_density(self, root):
        """MSRC density wrt labelled-tree node ages and background jump times."""
        rows = _rows(root); static = self.static_times(rows)
        t = 0.; stop = max(root.age,self.origin); logp = 0.
        while t < stop - EPS:
            end = self.next_time(t,static,stop); mid = (t+end)/2
            active = [(n,d) for n,p,lo,hi,d in rows if lo <= mid < hi]
            rate = 0.; labels = []
            for n,d in active:
                pop = self.pop(d,mid); h = _state(n.bh,mid,self.origin)
                rate += self.switch(pop,h,mid); labels.append((pop,h))
            for i,(pop,h) in enumerate(labels):
                for p2,h2 in labels[:i]:
                    if pop==p2 and h==h2: rate += self.pair(pop,h,mid)
            logp -= (end-t)*rate; t=end
        for n,p,lo,hi,d in rows:
            for t,h in n.bh[1:]:
                if lo < t < min(hi,self.origin):
                    old = _state(n.bh,np.nextafter(t,-np.inf)-1e-11,self.origin)
                    q = self.switch(self.pop(d,t),old,t)
                    if q <= 0: return -float('inf')
                    logp += math.log(q)
        for n in _nodes(root):
            if n.is_tip(): continue
            h = _state(n.children[0].bh,n.age-1e-11,self.origin)
            h2 = _state(n.children[1].bh,n.age-1e-11,self.origin)
            pop = self.pop(_desc(n.children[0]),n.age)
            if pop != self.pop(_desc(n.children[1]),n.age) or h != h2:
                return -float('inf')
            if n.age < self.origin and _state(n.bh,n.age,self.origin) != h: return -float('inf')
            logp += math.log(self.pair(pop,h,n.age))
        return logp

    def initial(self, states, rng):
        if set(states) != set(self.tree.taxa) or any(h not in (0,1) for h in states.values()):
            raise ValueError('one binary terminal background is required per taxon')
        for taxon,h in states.items(): self.pair(taxon,h,0.)
        g = simulate_genealogy_reference(0,self.tree,self.history,states,self.rho,1.,rng,True)
        nodes = {}
        for i,taxon in enumerate(self.tree.taxa):
            n = Node(taxon,age=0.); n.bh=[(0.,states[taxon])]; nodes[i]=n
        root = None
        for e in g.events:
            if e.event_type == 'switch':
                nodes[int(e.lineage_ids_before)].bh.append((e.backward_time,int(e.arrangements_after)))
            elif e.event_type == 'coalescence':
                a,b = [nodes[int(x)] for x in e.lineage_ids_before.split(';')]
                h = int(e.arrangements_after) if e.backward_time < self.origin else None
                root = Node(f'_s{e.lineage_id_after}',children=[a,b],age=e.backward_time)
                root.bh=[(root.age,h)]; a.parent=root; b.parent=root
                nodes[int(e.lineage_id_after)]=root
        # The augmented root stem matters when a later cut regrafts above it.
        root.bh = self.free_path(_desc(root),root.age,root.bh[0][1],rng)
        root._log_pi = self.log_density(root)
        if not np.isfinite(root._log_pi):
            raise ValueError('initial genealogy/history outside strict_v1 target support')
        return root

    def free_path(self, desc, t, h, rng):
        hist=[(t,h)]
        static={self.origin,*[b.older_age for b in self.tree.branches.values()]}
        while t < self.origin - EPS:
            pop=self.pop(desc,t); q=self.switch(pop,h,t)
            end=self.next_time(t,static,self.origin)
            w=rng.exponential(1/q) if q else float('inf')
            if t+w < end:
                t+=w; h=1-h; hist.append((t,h))
            else: t=end
        return hist

    def cut_pieces(self, root, x):
        rows=_edge_rows(root); static=self.static_times(rows); pieces=[]
        for n,p,lo,hi,d in rows:
            t=lo
            while t < hi-EPS:
                end=self.next_time(t,static,hi); mid=(t+end)/2
                if mid >= self.origin-EPS: factor=1.
                else:
                    h=_state(n.bh,mid,self.origin); pop=self.pop(d,mid)
                    self.pair(pop,h,mid)
                    freq=self.p(pop,mid); same=freq if h else 1-freq
                    factor=same+self.kappa(x)*(1-same)
                if factor>0: pieces.append((n,p,t,end,d,factor))
                t=end
        return pieces

    def eligible(self, guide, original, desc, t, h):
        pop=self.pop(desc,t); found=[]
        for n,p,lo,hi,d in _rows(guide):
            if lo <= t < hi and self.pop(d,t)==pop and _state(n.bh,t,self.origin)==h:
                found.append((n,p,'ancestral_stem' if p is None else 'guide'))
        hist, upper=original
        if t < upper-EPS and _state(hist,t,self.origin)==h:
            found.append((None,None,'original_edge'))
        return found, pop

    def path_log_density(self, guide, original, desc, path, tr, end):
        rows=_rows(guide);static=self.static_times(rows)|{original[1], *[t for t,h in original[0]], *[t for t,h in path],end}
        t=tr; logq=0.
        while t < end-EPS:
            hi=self.next_time(t,static,end);mid=(t+hi)/2
            h=_state(path,mid,self.origin); eligible,pop=self.eligible(guide,original,desc,mid,h)
            rate=self.switch(pop,h,mid)
            if eligible: rate+=len(eligible)*self.pair(pop,h,mid)
            logq-=(hi-t)*rate;t=hi
        for t,h in path[1:]:
            if tr<t<min(end,self.origin):
                old=_state(path,t-1e-11,self.origin)
                q=self.switch(self.pop(desc,t),old,t)
                if q<=0:return -float('inf')
                logq+=math.log(q)
        h=_state(path,end-1e-11,self.origin);pop=self.pop(desc,end)
        return logq+math.log(self.pair(pop,h,end))

    def sample_path(self, guide, original, desc, tr, h, rng):
        path=[(tr,h)];t=tr
        static=self.static_times(_rows(guide))|{original[1], *[t for t,h in original[0]]}
        while True:
            h=_state(path,t,self.origin);eligible,pop=self.eligible(guide,original,desc,t,h)
            q=self.switch(pop,h,t); c=self.pair(pop,h,t) if eligible else 0.
            total=q+len(eligible)*c; end=self.next_time(t,static)
            if total==0 and not np.isfinite(end): raise RuntimeError('no re-coalescence opportunity')
            w=rng.exponential(1/total) if total else float('inf')
            if t+w>=end: t=end;continue
            t+=w
            if rng.random()*total<q:
                path.append((t,1-h));continue
            target=eligible[int(rng.integers(len(eligible)))]
            return path,t,target,pop

    def transition(self, root, x, rng):
        old=copy.deepcopy(root); work=copy.deepcopy(root)
        pieces=self.cut_pieces(work,x); weights=np.array([(hi-lo)*f for n,p,lo,hi,d,f in pieces])
        i=int(rng.choice(len(pieces),p=weights/weights.sum()))
        child,parent,lo,hi,desc,f=pieces[i]; tr=float(rng.uniform(lo,hi))
        old_upper=parent.age; old_hist=list(child.bh); parent_hist=list(parent.bh)
        sibling=next(n for n in parent.children if n is not child)
        sibling.bh=_join(sibling.bh,parent_hist)
        guide=_detach(work,child,parent)
        original=(old_hist,old_upper)
        path,tcoal,(target,target_parent,kind),pop=self.sample_path(guide,original,desc,tr,_state(old_hist,tr,self.origin),rng)
        forward=self.path_log_density(guide,original,desc,path,tr,tcoal) if self.config.audit_density else None
        if kind=='original_edge':
            proposal=copy.deepcopy(old)
            newchild=next(n for n in _nodes(proposal) if _desc(n)==desc)
            newchild.bh=_join(_segment(old_hist,newchild.age,tr,self.origin),
                              _join(path,_segment(old_hist,tcoal,old_upper,self.origin)))
            reverse_path=_segment(old_hist,tr,tcoal,self.origin)
            reverse_original=(newchild.bh,old_upper)
        else:
            target_hist=list(target.bh)
            child.bh=_join(_segment(old_hist,child.age,tr,self.origin),path)
            target.bh=_segment(target_hist,target.age,tcoal,self.origin)
            proposal=_insert(guide,target,target_parent,child,tcoal)
            join=child.parent
            join.bh=_segment(target_hist,tcoal,float('inf'),self.origin)
            reverse_path=_segment(old_hist,tr,old_upper,self.origin)
            reverse_original=(child.bh,tcoal)
        # Recreate the common guide after constructing the proposal. The
        # accepted state never retains a second genealogy or ancestral material.
        reverse=None
        if self.config.audit_density:
            revwork=copy.deepcopy(proposal)
            rc=next(n for n in _nodes(revwork) if _desc(n)==desc); rp=rc.parent
            rs=next(n for n in rp.children if n is not rc)
            rs.bh=_join(rs.bh,rp.bh); reverse_guide=_detach(revwork,rc,rp)
            reverse_end=tcoal if kind=='original_edge' else old_upper
            reverse=self.path_log_density(reverse_guide,reverse_original,desc,reverse_path,tr,reverse_end)
            lp_old=old._log_pi if getattr(old,'_log_pi',None) is not None else self.log_density(old)
            lp_new=self.log_density(proposal);proposal._log_pi=lp_new
            log_ratio=lp_new-lp_old+reverse-forward
        else:
            # Analytic cancellation for this inherited-state closure only.
            # Preserve RNG semantics: the acceptance uniform is still drawn.
            log_ratio=0.;proposal._log_pi=None
        # Spatial proposal intensity is r(x)*factor(h_at_cut,t,x) dt.
        # The retained subtree and its state below tr are identical in the
        # reverse move, so these unnormalised cut measures cancel. Do NOT
        # insert total-tree-length normalisers in this continuous-position MH.
        accepted=bool(np.log(max(rng.random(),np.finfo(float).tiny))<min(0.,log_ratio))
        result=proposal if accepted else old
        changed=genealogy_signature(result)!=genealogy_signature(old)
        bg_changed=background_signature(result)!=background_signature(old)
        event=dict(position=x,recombination_time=tr,recoal_time=tcoal,target_edge_kind=kind,
                   population=pop,accepted=accepted,genealogy_changed=changed,
                   topology_changed=_topology(result)!=_topology(old),background_history_changed=bg_changed,
                   background_changed=background_jump_signature(result)!=background_jump_signature(old),
                   old_background=_state(old_hist,tr,self.origin),new_background=_state(path,tcoal-1e-11,self.origin),
                   proposed_path_switches=max(0,len(path)-1),log_acceptance_ratio=float(log_ratio),
                   detached_descendants=','.join(sorted(desc)),log_forward_path=forward,log_reverse_path=reverse)
        return result,event


def simulate_structured_smc(tree, history, config, sampled_backgrounds=None):
    """Generate sequential local genealogies with a frozen reference marginal.

    Neutral mode delegates directly to the validated neutral engine. In the
    structured mode, each candidate is one SPR/path proposal, never an ARG.
    """
    if config.neutral:
        if config.start!=0 or config.end not in (None,config.recombination_map.length):
            raise ValueError('neutral delegation currently requires a full-map scan')
        return simulate_msc_smc(tree,SMCConfig(config.recombination_map.length,config.recombination_map,config.seed))
    kernel=Kernel(tree,history,config);rng=np.random.default_rng(config.seed)
    if sampled_backgrounds is None:
        sampled_backgrounds={taxon:int(rng.random()<kernel.p(taxon,0)) for taxon in tree.taxa}
    root=kernel.initial(sampled_backgrounds,rng);result=StructuredResult();x=config.start
    boundaries=sorted({config.recombination_map.length,*config.recombination_map.breakpoints,*config.suppressed_region})
    while x<kernel.end-EPS:
        start=x;budget=float(rng.exponential());event_position=kernel.end
        while x<kernel.end-EPS:
            end=min([b for b in boundaries if b>x+EPS]+[kernel.end]);end=min(end,kernel.end)
            pieces=kernel.cut_pieces(root,x)
            length=sum((hi-lo)*f for n,p,lo,hi,d,f in pieces)
            rate=config.recombination_map.rate_at(x)*length
            if rate>0 and budget < (end-x)*rate:
                event_position=x+budget/rate;break
            budget-=(end-x)*rate;x=end
        result.trees.append(dict(start=start,end=event_position,newick=_newick(root),topology=_topology(root),
                                 histories=[dict(descendants=','.join(sorted(_desc(n))),birth=n.age,path=n.bh) for n in _nodes(root)]))
        if event_position>=kernel.end-EPS:break
        root,event=kernel.transition(root,event_position,rng);result.events.append(event);x=event_position
        if len(result.events)>=config.max_events: raise RuntimeError('max_events reached; output would be incomplete')
    result.metadata=dict(model='structured_SMC_prime_MH_window_v1',seed=config.seed,origin_age=kernel.origin,
                         reference_position=kernel.reference,rho_reference=kernel.rho,scan_start=config.start,scan_end=kernel.end,
                         marker=config.marker,suppressed_region=config.suppressed_region,cross_retention=config.cross_retention,
                         sampled_backgrounds=sampled_backgrounds,boundary_policy='strict_v1',
                         density_audit=config.audit_density,
                         target='frozen_reference_MSRC_path_density',
                         event_semantics='candidate SMC SPR/path proposals; not a physical ARG event census',
                         approximation='frozen reference marginal across window; MH-corrected conditional path closure')
    check_positions=[config.start,kernel.end]
    if config.start<=config.marker<=kernel.end:check_positions.append(config.marker)
    rates=[kernel.effective_rate(u) for u in check_positions]
    result.metadata['rho_window_range']=[min(rates),max(rates)]
    result.metadata['target_rate_constant_on_window']=max(rates)==min(rates)
    result.metadata['rho_range_is_genealogy_error_bound']=False
    return result


def write_structured_result(result, output):
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    if not isinstance(result,StructuredResult):
        from .engine import write_result
        return write_result(result,out)
    def write(name,rows,fields):
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
    write('local_genealogies.csv',[{k:t[k] for k in ('start','end','newick','topology')} for t in result.trees],['start','end','newick','topology'])
    (out/'local_gene_trees.nwk').write_text(''.join(f"[&start={t['start']},end={t['end']}] {t['newick']}\n" for t in result.trees))
    fields=list(result.events[0]) if result.events else ['position','recombination_time','accepted','genealogy_changed','topology_changed','background_history_changed']
    write('recombination_events.csv',result.events,fields)
    write('genealogy_changes.csv',[e for e in result.events if e['genealogy_changed']],fields)
    write('background_transitions.csv',[e for e in result.events if e['background_changed']],fields)
    (out/'background_histories.jsonl').write_text(''.join(json.dumps(dict(start=t['start'],end=t['end'],histories=t['histories']))+'\n' for t in result.trees))
    (out/'metadata.json').write_text(json.dumps(result.metadata,indent=2)+'\n')
    tracks=[]
    from ..species_tree import SpeciesTree
    for row in result.trees:
        root=SpeciesTree(row['newick'],1,1,ultrametric_tolerance=1e-7).root
        taxa=sorted(_desc(root));quartet=-1
        if len(taxa)==4:
            pair=next(_desc(n) for n in _nodes(root) if len(_desc(n))==2)
            quartet=next(i for i,x in enumerate(({taxa[0],taxa[1]},{taxa[0],taxa[2]},{taxa[0],taxa[3]})) if pair==x or set(taxa)-pair==x)
        tracks.append(dict(start=row['start'],end=row['end'],quartet_index=quartet,root_age=root.age,topology=row['topology']))
    write('quartet_tracks.csv',tracks,['start','end','quartet_index','root_age','topology'])
    write('sample_backgrounds.csv',[dict(taxon=t,background=h) for t,h in result.metadata['sampled_backgrounds'].items()],['taxon','background'])


def main(argv=None):
    import argparse,yaml
    from ..experiments import _tree_from_config
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args(argv);cfg=yaml.safe_load(Path(args.config).read_text());tree=_tree_from_config(cfg)
    c=cfg['structured_smc'];m=c['recombination_map']
    rec=PiecewiseConstantMap(float(c['chromosome_length']),tuple(m.get('breakpoints',[])),tuple(m['rates']))
    config=StructuredConfig(rec,float(c['marker']),tuple(c['suppressed_region']),float(c.get('cross_retention',1)),
                            int(cfg.get('seed',1)),float(c.get('start',0)),c.get('end'),c.get('reference_position'),bool(c.get('neutral',False)),
                            audit_density=bool(c.get('audit_density',False)))
    hc=c['frequency_history'];kind=hc.get('model','constant')
    if kind=='constant':history=ConstantHistory(float(hc['p']),float(hc['origin_age']))
    elif kind=='frozen':
        from ..arg.history import load_frequency_history
        history=load_frequency_history(hc['path'],float(hc['origin_age']))
    elif kind in ('moran','wright_fisher'):
        from ..population_process import simulate_population_history
        from ..rearrangement import Rearrangement
        rr=cfg['rearrangement']
        rearrangement=Rearrangement(rr.get('id','background_1'),rr.get('type','suppressed_recombination'),rr['origin_branch'],
                                    int(rr['origin_time_from_branch_start']),int(rr.get('initial_copy_count',1)),float(rr.get('selection_coefficient',0)))
        config_copy=copy.deepcopy(cfg);config_copy['population_process']={'model':kind}
        history=simulate_population_history(tree,rearrangement,np.random.default_rng(config.seed+1),config_copy)
    else:raise ValueError('history model must be constant, frozen, moran, or wright_fisher')
    states=c.get('sampled_backgrounds')
    if states is not None: states={str(t):h for t,h in states.items()}
    result=simulate_structured_smc(tree,history,config,states)
    result.metadata['configuration']=cfg
    write_structured_result(result,args.output)
    records=getattr(history,'records',None)
    if records is None:
        from ..wright_fisher import FrequencyRecord
        records=[]
        for bid,b in tree.branches.items():
            ages=sorted({b.younger_age,b.older_age,*([history.structured_origin_age] if b.younger_age<=history.structured_origin_age<=b.older_age else [])},reverse=True)
            for i,age in enumerate(ages):
                p=history.frequency_at(bid,age);M=2*b.effective_population_size
                records.append(FrequencyRecord('frozen_background',bid,b.parent_branch_id,i,age,M*p,M*(1-p),M,p,1-p,'externally_fixed',age==history.structured_origin_age,i==0,i==len(ages)-1,0.,b.effective_population_size))
    if records:
        from dataclasses import asdict
        with (Path(args.output)/'frequency_history.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(asdict(records[0])),lineterminator='\n');w.writeheader();w.writerows(asdict(r) for r in records)
    return 0


if __name__=='__main__':raise SystemExit(main())
