"""Marginal and spatial audit of frozen-reference structured SMC'.

Run from the repository root. No windows are treated as independent samples.
"""
import argparse,csv,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from scipy.stats import norm
from scipy.linalg import expm
from scipy.integrate import quad
from msrcsim.smc.structured import (StructuredConfig,ConstantHistory,simulate_structured_smc)
from msrcsim.smc.engine import _nodes,_desc
from msrcsim.smc.recombination_map import PiecewiseConstantMap as Map
from msrcsim.species_tree import SpeciesTree
from msrcsim.structured_coalescent import simulate_genealogy_reference
from msrcsim.spatial_reservoir import simulate_reservoir_arg
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as AMap
from msrcsim.arg.output.tree_sequence import build_tree_sequence


def gene_metrics(root):
    nodes=_nodes(root);taxa=sorted(_desc(root));pairs=[]
    for i,a in enumerate(taxa):
        for b in taxa[:i]:pairs.append((min(n.age for n in nodes if {a,b}<=_desc(n)),a,b))
    top=-1
    if len(taxa)==4:
        _,a,b=min(pairs);pair={a,b}
        top=next(i for i,x in enumerate(({taxa[0],taxa[1]},{taxa[0],taxa[2]},{taxa[0],taxa[3]})) if pair==x or set(taxa)-pair==x)
    internal=sum(n.parent.age-n.age for n in nodes if n.parent and not n.is_tip())
    return [root.age,np.mean([p[0] for p in pairs]),internal,top]


def legacy_metrics(g,nt):
    taxa=('A','B') if nt==2 else ('A','B','C','D');pairs=[]
    for i,a in enumerate(taxa):
        for b in taxa[:i]:pairs.append(min(c.time for c in g.coalescences if {a,b}<=set(c.descendant_set_1.split(','))|set(c.descendant_set_2.split(','))))
    # Build dated tree from full-precision coalescence records, never rounded Newick.
    from msrcsim.species_tree import Node
    live={i:Node(taxon,age=0) for i,taxon in enumerate(taxa)}
    for c in g.coalescences:
        a,b=live[c.child_lineage_1],live[c.child_lineage_2]
        parent=Node(str(c.parent_lineage),children=[a,b],age=c.time);a.parent=parent;b.parent=parent
        live[c.parent_lineage]=parent
    return gene_metrics(parent)


def ts_metrics(tree,nt):
    pairs=[(tree.time(tree.mrca(i,j)),i,j) for i in range(nt) for j in range(i)]
    top=-1
    if nt==4:
        _,a,b=min(pairs);pair={a,b};top=next(k for k,x in enumerate(({0,1},{0,2},{0,3})) if pair==x or set(range(4))-pair==x)
    internal=sum(tree.branch_length(n) for n in tree.nodes() if tree.parent(n)!=-1 and not tree.is_sample(n))
    return [tree.time(tree.root),np.mean([p[0] for p in pairs]),internal,top]


def focal_row(result,x):
    return next(t for t in result.trees if t['start']<=x<t['end'])


def path_summary(row,origin):
    paths=row['histories'];total0=total1=nswitch=0.
    clades={frozenset(p['descendants'].split(',')):p for p in paths}
    for desc,p in clades.items():
        parent=min([q for d,q in clades.items() if desc<d],key=lambda q:len(q['descendants'].split(',')),default=None)
        if parent is None:continue
        end=min(parent['birth'],origin)
        events=[(t,h) for t,h in p['path'] if t<end]
        for j,(t,h) in enumerate(events):
            stop=events[j+1][0] if j+1<len(events) else end
            if h==0:total0+=max(0,stop-t)
            elif h==1:total1+=max(0,stop-t)
        nswitch+=max(0,len(events)-1)
    return [nswitch,total0,total1]


def metric_vectors(values,nt):
    a=np.asarray(values);out={name:a[:,i] for i,name in enumerate(('root_age','pair_mean','internal_length'))}
    if nt==4:
        out.update({f'q{i+1}':(a[:,3]==i).astype(float) for i in range(3)})
        out['D']=(a[:,3]==1).astype(float)-(a[:,3]==0).astype(float)
    return out


def compare(v,ref,case,method):
    out=[]
    for key,a in v.items():
        b=ref[key];diff=float(a.mean()-b.mean());se=float(np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b)))
        tol=.05 if key.startswith('q') or key=='D' else .1*float(b.mean())
        out.append(dict(case=case,method=method,metric=key,observed=float(a.mean()),reference=float(b.mean()),difference=diff,se=se,ci_low=diff-1.96*se,ci_high=diff+1.96*se,tolerance=tol,pvalue=float(2*norm.sf(abs(diff/se))) if se else float(diff==0)))
    return out


def run(out,nrep,nt,p,k,seed):
    out.mkdir(parents=True,exist_ok=True)
    t=SpeciesTree('(A:0,B:0)ROOT;' if nt==2 else '((A:2,B:2)AB:2,(C:2,D:2)CD:2)ROOT;',20,200)
    states={taxon:i%2 for i,taxon in enumerate(t.taxa)};hist=ConstantHistory(p,100)
    cfg=StructuredConfig(Map.constant(100,.0002),0,(0,100),k,seed,reference_position=50,audit_density=True)
    rmap=AMap.constant(100,.0002);rho=.0002*50*k
    values={m:[] for m in ('SMC','legacy','reservoir')};spatial=[];draws=[];paths=[];events=[]
    start=time.perf_counter()
    for i in range(nrep):
        seeds=np.random.SeedSequence([seed,nt,int(p*100),int(k*100),i]).generate_state(3)
        cfg.seed=int(seeds[0]);r=simulate_structured_smc(t,hist,cfg,states)
        row=focal_row(r,50);root=SpeciesTree(row['newick'],20,200,ultrametric_tolerance=1e-7).root
        values['SMC'].append(gene_metrics(root));paths.append(path_summary(row,100))
        g=simulate_genealogy_reference(0,t,hist,states,rho,1,np.random.default_rng(seeds[1]),False)
        values['legacy'].append(legacy_metrics(g,nt))
        arg=simulate_reservoir_arg(t,rmap,hist,0,(0,100),k,int(seeds[2]),states)
        ts=build_tree_sequence(arg).simplify();values['reservoir'].append(ts_metrics(ts.at(50),nt))
        for x in (5,25,50,75,95):
            row=focal_row(r,x);root=SpeciesTree(row['newick'],20,200,ultrametric_tolerance=1e-7).root
            sm=gene_metrics(root);am=ts_metrics(ts.at(x),nt)
            spatial.append(dict(replicate=i,position=x,smc_root=sm[0],smc_pair=sm[1],smc_top=sm[3],arg_root=am[0],arg_pair=am[1],arg_top=am[3],
                                smc_newick=row['newick'],arg_signature=';'.join(f'{ts.at(x).time(ts.at(x).mrca(a,b)):.15g}' for a in range(nt) for b in range(a))))
        for m in values:
            a=values[m][-1];draws.append(dict(replicate=i,method=m,root_age=a[0],pair_mean=a[1],internal_length=a[2],topology=a[3]))
        for e in r.events:events.append(dict(replicate=i,**e))
        if (i+1)%100==0:print(nt,p,k,i+1,flush=True)
    case=f'{nt}tax_p{p}_k{k}'
    tests=[];ref=metric_vectors(values['legacy'],nt)
    for method in ('SMC','reservoir'):tests.extend(compare(metric_vectors(values[method],nt),ref,case,method))
    if nt==2:
        c0=1/(40*(1-p));c1=1/(40*p)
        Q=np.array([[-(2*rho*p+c0),2*rho*p,0],[rho*(1-p),-rho,rho*p],[0,2*rho*(1-p),-(2*rho*(1-p)+c1)]])
        survival=lambda x:expm(Q*x)[1].sum()
        target=quad(survival,0,100)[0]+survival(100)*40
        a=np.array(values['SMC'])[:,0];diff=float(a.mean()-target);se=float(a.std(ddof=1)/np.sqrt(nrep))
        tests.append(dict(case=case,method='SMC_analytic',metric='root_age',observed=float(a.mean()),reference=float(target),difference=diff,se=se,ci_low=diff-1.96*se,ci_high=diff+1.96*se,tolerance=.1*target,pvalue=float(2*norm.sf(abs(diff/se)))))
    def write(name,rows):
        if not rows:return
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    write('marginal_tests.csv',tests);write('draws.csv',draws);write('spatial.csv',spatial);write('events.csv',events)
    write('path_summaries.csv',[dict(replicate=i,switches=a[0],dwell0=a[1],dwell1=a[2]) for i,a in enumerate(paths)])
    summary=dict(case=case,replicates=nrep,seed=seed,runtime_all_methods=time.perf_counter()-start,smc_events=len(events),
                 accepted=sum(e['accepted'] for e in events),topology_changes=sum(e['topology_changed'] for e in events),
                 genealogy_changes=sum(e['genealogy_changed'] for e in events),background_changes=sum(e['background_history_changed'] for e in events),
                 max_abs_log_acceptance=max([abs(e['log_acceptance_ratio']) for e in events],default=0),
                 approximation='frozen reference 50; spatial marginals away from reference differ from marker-distance target')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);a.add_argument('--replicates',type=int,default=1000)
    a.add_argument('--taxa',type=int,choices=(2,4),default=4);a.add_argument('--p',type=float,default=.5);a.add_argument('--kappa',type=float,default=.1);a.add_argument('--seed',type=int,default=20261121)
    args=a.parse_args();print(json.dumps(run(args.output,args.replicates,args.taxa,args.p,args.kappa,args.seed)),flush=True)
