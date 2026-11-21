#!/usr/bin/env python3
"""Seeded simple-case audit. Run from repo root; see the implementation plan."""
import argparse, csv, json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
import platform, scipy, tskit
from scipy.linalg import expm
from scipy.stats import norm
from msrcsim.species_tree import SpeciesTree
from msrcsim.spatial_reservoir import simulate_reservoir_arg
from msrcsim.arg.engine import ARGConfig, simulate_arg
from msrcsim.arg.backgrounds.msrc_two_state import MSRCTwoBackground
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.arg.output.tree_sequence import build_tree_sequence
from msrcsim.structured_coalescent import simulate_genealogy_reference


class Frozen:
    structured_origin_age = 100.
    def __init__(self, p): self.p=p
    def frequency_at(self, population, time): return self.p
    def next_frequency_boundary(self, population, time): return float('inf')
    def terminal_frequency(self, taxon): return self.p


def metrics(tree, n):
    pairs=[(tree.time(tree.mrca(i,j)), i,j) for i in range(n) for j in range(i)]
    top=-1
    if n==4:
        _,i,j=min(pairs)
        pair={i,j}; top=next(k for k,x in enumerate(({0,1},{0,2},{0,3})) if pair==x or set(range(4))-pair==x)
    return [tree.time(tree.root), np.mean([x[0] for x in pairs]), top]


def legacy_metrics(g,n):
    taxa=('A','B') if n==2 else ('A','B','C','D')
    vals=[]
    for i,a in enumerate(taxa):
        for b in taxa[:i]:
            vals.append(min(c.time for c in g.coalescences if {a,b} <= set(c.descendant_set_1.split(',')) | set(c.descendant_set_2.split(','))))
    return [max(g.coalescence_times),np.mean(vals),g.topology_index]


def run(out, nrep, case=None):
    out.mkdir(parents=True,exist_ok=True); rows=[]; spatial=[]; timing=[]; selected=None
    for nt in (2,4):
        newick='(A:0,B:0)ROOT;' if nt==2 else '((A:2,B:2)AB:2,(C:2,D:2)CD:2)ROOT;'
        tree=SpeciesTree(newick,20,200)
        states={taxon: i%2 for i,taxon in enumerate(tree.taxa)}
        for p in (.2,.5,.8):
            for k in (0., .1, 1.):
                if case is not None and (nt,p,k)!=case:
                    continue
                hist=Frozen(p); rate=.0002; length=100.; focal=75.; rho=rate*k*focal
                rec_map=Map.constant(length,rate); values={key:[] for key in ('reservoir','legacy','current_ARG')}
                start=time.perf_counter()
                for i in range(nrep):
                    # SeedSequence streams separate all methods and cases.
                    seq=np.random.SeedSequence([20261121,nt,int(p*100),int(k*100),i])
                    s=seq.generate_state(3)
                    result=simulate_reservoir_arg(tree,rec_map,hist,0,(0,length),k,int(s[0]),states)
                    ts=build_tree_sequence(result).simplify()
                    values['reservoir'].append(metrics(ts.at(focal),nt))
                    left,right=metrics(ts.at(25),nt),metrics(ts.at(75),nt)
                    spatial.append(dict(taxa=nt,p=p,kappa=k,replicate=i,root_left=left[0],root_right=right[0],same_topology=int(left[2]==right[2])))
                    g=simulate_genealogy_reference(0,tree,hist,states,rho,1,np.random.default_rng(s[1]),False)
                    values['legacy'].append(legacy_metrics(g,nt))
                    old=MSRCTwoBackground(hist,0,(0,length),k,states)
                    arg=simulate_arg(tree,ARGConfig(length,rec_map,old),np.random.default_rng(s[2]))
                    values['current_ARG'].append(metrics(build_tree_sequence(arg).simplify().at(focal),nt))
                timing.append(dict(taxa=nt,p=p,kappa=k,replicates=nrep,all_three_methods_seconds=time.perf_counter()-start))
                arrays={key:np.asarray(val) for key,val in values.items()}
                if nt==4 and p==.5 and k==1:
                    selected=arrays
                for method in ('reservoir','current_ARG'):
                    for name,col in [('root_age',0),('pair_mean',1)]+([(f'q{j+1}',j+2) for j in range(3)] if nt==4 else []):
                        a=arrays[method][:,col] if col<2 else (arrays[method][:,2]==col-2).astype(float)
                        b=arrays['legacy'][:,col] if col<2 else (arrays['legacy'][:,2]==col-2).astype(float)
                        diff=a.mean()-b.mean(); se=np.sqrt(a.var(ddof=1)/nrep+b.var(ddof=1)/nrep)
                        rows.append(dict(taxa=nt,p=p,kappa=k,method=method,metric=name,observed=a.mean(),reference=b.mean(),difference=diff,se=se,ci_low=diff-1.96*se,ci_high=diff+1.96*se,pvalue=float(2*norm.sf(abs(diff/se))) if se else float(diff==0)))
                if nt==2:
                    c0=1/(40*(1-p));c1=1/(40*p)
                    Q=np.array([[-(2*rho*p+c0),2*rho*p,0],[rho*(1-p),-rho,rho*p],[0,2*rho*(1-p),-(2*rho*(1-p)+c1)]])
                    surv=expm(Q*hist.structured_origin_age)[1].sum()
                    # Q is invertible except in degenerate no-switch regimes.
                    from scipy.integrate import quad
                    target=quad(lambda t:expm(Q*t)[1].sum(),0,hist.structured_origin_age)[0]+surv*40
                    a=arrays['reservoir'][:,0]; se=a.std(ddof=1)/np.sqrt(nrep);diff=a.mean()-target
                    rows.append(dict(taxa=nt,p=p,kappa=k,method='reservoir_analytic',metric='root_age',observed=a.mean(),reference=target,difference=diff,se=se,ci_low=diff-1.96*se,ci_high=diff+1.96*se,pvalue=float(2*norm.sf(abs(diff/se)))))
                print(nt,p,k,'completed',flush=True)
    summarize(out,nrep,rows,spatial,timing,selected)


def summarize(out,nrep,rows,spatial,timing,selected=None):
    # Holm correction across the complete primary marginal family.
    order=sorted(range(len(rows)),key=lambda i:rows[i]['pvalue']); prev=0
    for rank,i in enumerate(order):
        prev=max(prev,min(1.,(len(rows)-rank)*rows[i]['pvalue'])); rows[i]['holm_p']=prev
    def write(name,data):
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)
    write('marginals.csv',rows);write('spatial_replicates.csv',spatial);write('runtime.csv',timing)
    rng=np.random.default_rng(271828); space=[]
    for p in (.2,.5,.8):
        for k in (0.,.1,1.):
            a=[r for r in spatial if r['taxa']==4 and r['p']==p and r['kappa']==k]
            if not a: continue
            x=np.array([r['root_left'] for r in a]);y=np.array([r['root_right'] for r in a]);n=len(a);bs=[]
            for _ in range(500):
                ix=rng.integers(n,size=n);bs.append(np.corrcoef(x[ix],y[ix])[0,1])
            top=np.mean([r['same_topology'] for r in a]);z=1.96
            den=1+z*z/n;center=(top+z*z/(2*n))/den;half=z*np.sqrt(top*(1-top)/n+z*z/(4*n*n))/den
            space.append(dict(p=p,kappa=k,n_chromosomes=n,separation_bp=50,root_age_correlation=np.corrcoef(x,y)[0,1],correlation_bootstrap_low=np.quantile(bs,.025),correlation_bootstrap_high=np.quantile(bs,.975),topology_agreement=top,topology_Wilson_low=center-half,topology_Wilson_high=center+half))
    if space: write('spatial_summary.csv',space)
    summary={'replicates_per_case':nrep,'seed':20261121,'cases':len(timing),
             'versions':dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,tskit=tskit.__version__),
             'criteria':{'holm_family_alpha':.05,'quartet_practical_tolerance':.05,'root_mean_relative_tolerance':.10},
             'holm_failures':{m:sum(r['holm_p']<.05 for r in rows if r['method']==m) for m in ('reservoir','current_ARG','reservoir_analytic')},
             'scope':'fixed interior p, explicit origin, two and four taxa; no whole-genome scalability claim'}
    summary['practical_equivalence']={}
    for method in ('reservoir','current_ARG','reservoir_analytic'):
        subset=[r for r in rows if r['method']==method]
        def within(r):
            tol=.05 if r['metric'].startswith('q') else .1*r['reference']
            return r['ci_low']>=-tol and r['ci_high']<=tol
        summary['practical_equivalence'][method]={'comparisons':len(subset),'95pct_difference_CI_within_tolerance':int(sum(within(r) for r in subset))}
    (out/'validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for m,c in [('reservoir','tab:blue'),('current_ARG','tab:orange')]:
        rr=[r for r in rows if r['taxa']==4 and r['metric']=='q1' and r['method']==m]
        axes[0].errorbar(np.arange(len(rr))+(0 if m=='reservoir' else .15),[r['difference'] for r in rr],yerr=[1.96*r['se'] for r in rr],fmt='o',label=m,color=c)
        rr=[r for r in rows if r['metric']=='root_age' and r['method']==m]
        axes[1].plot([r['reference'] for r in rr],[r['observed'] for r in rr],'o',label=m,color=c)
    axes[0].axhline(0,color='k',lw=.5);axes[0].set(title='Four-taxon q1 difference vs legacy',xlabel='p × κ case',ylabel='Difference with 95% CI');axes[0].legend()
    axes[1].plot([0,180],[0,180],'k--');axes[1].set(title='Root-age mean',xlabel='Legacy',ylabel='Spatial process')
    for p in (.2,.5,.8):
        cor=[]; ks=[]
        for k in (0.,.1,1.):
            a=[r for r in spatial if r['taxa']==4 and r['p']==p and r['kappa']==k]
            if not a: continue
            ks.append(k)
            cor.append(np.corrcoef([r['root_left'] for r in a],[r['root_right'] for r in a])[0,1])
        if ks: axes[2].plot(ks,cor,'o-',label=f'p={p}')
    axes[2].set(title='Across-chromosome root correlation',xlabel='κ cross',ylabel='Corr(T25, T75)');axes[2].legend()
    fig.tight_layout();fig.savefig(out/'diagnostics.png',dpi=160);plt.close(fig)
    if selected is not None:
        selected_rows=[dict(method=name,replicate=i,root_age=val[0],pair_mean=val[1],topology=int(val[2])) for name,vals in selected.items() for i,val in enumerate(vals)]
        write('selected_case_draws.csv',selected_rows)
        fig,ax=plt.subplots(figsize=(6,4))
        for name,vals in selected.items():
            t=np.sort(vals[:,0]);ax.step(t,np.arange(1,len(t)+1)/len(t),where='post',label=name)
        ax.set(title='Four taxa, p=0.5, κ=1',xlabel='Root age (generations)',ylabel='Empirical CDF');ax.legend();fig.tight_layout()
        fig.savefig(out/'root_time_distributions.png',dpi=160);plt.close(fig)
    print(json.dumps(summary),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--replicates',type=int,default=1000)
    p.add_argument('--case',help='Optional n,p,k subset, e.g. 4,0.2,1')
    p.add_argument('--summarize-only',action='store_true',help='Regenerate summaries/figures from saved CSVs')
    args=p.parse_args();case=None
    if args.case:
        n,f,k=args.case.split(',');case=(int(n),float(f),float(k))
    if args.replicates<2: p.error('at least two independent replicates are required')
    if args.summarize_only:
        def read(name):
            with (args.output/name).open() as f:
                return [{k:(v if k in ('method','metric') else float(v)) for k,v in row.items()} for row in csv.DictReader(f)]
        timing=read('runtime.csv')
        summarize(args.output,int(timing[0]['replicates']),read('marginals.csv'),read('spatial_replicates.csv'),timing)
    else:
        run(args.output,args.replicates,case)
