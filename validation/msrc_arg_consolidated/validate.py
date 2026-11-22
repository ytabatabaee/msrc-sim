"""Unified-API marginal, D, switching-exposure and autocorrelation audit.

Existing validators are rerun independently; this supplies missing focal
background diagnostics and routes both versions through the primary API.
"""
import argparse,csv,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from scipy.stats import norm
from msrcsim.arg import simulate_msrc_arg
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.arg.output.tree_sequence import build_tree_sequence
from msrcsim.species_tree import SpeciesTree
from msrcsim.structured_coalescent import simulate_genealogy_reference
from validation.spatial_reservoir.validate import Frozen,metrics,legacy_metrics


def focal_background_exposure(result,states,focal,origin):
    # Track only carriers of focal material through the recorded ARG. A split
    # selects one parent; a background transfer preserves the carrier ID.
    live={i:h for i,h in enumerate(states.values())};last=0.;count=dwell0=dwell1=0.
    for e in result.events:
        stop=min(e.time,origin);dt=max(0.,stop-last)
        dwell0+=dt*sum(h==0 for h in live.values());dwell1+=dt*sum(h==1 for h in live.values());last=max(last,stop)
        ids=e.lineage_ids
        if e.event_type=='background_transfer' and ids[0] in live:
            h=int(e.background_after);count+=int(e.time<origin and h!=live[ids[0]]);live[ids[0]]=h
        elif e.event_type=='recombination' and ids[0] in live:
            old=live.pop(ids[0]);side=0 if focal<e.position else 1
            h=int(e.background_after.split(',')[side]);live[ids[side+1]]=h
            count+=int(e.time<origin and h!=old)
        elif e.event_type=='coalescence':
            h=live.pop(ids[0],None);g=live.pop(ids[1],None)
            if h is not None or g is not None:live[ids[2]]=h if h is not None else g
        if len(live)==1:break
    return count,dwell0,dwell1


def vectors(draws):
    a=np.array(draws)
    return dict(root_age=a[:,0],pair_mean=a[:,1],q1=(a[:,2]==0).astype(float),
                q2=(a[:,2]==1).astype(float),q3=(a[:,2]==2).astype(float),
                D=(a[:,2]==1).astype(float)-(a[:,2]==0).astype(float))


def write(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def main(out,n):
    if n<1000:raise ValueError('primary audit requires at least 1000 independent chromosomes')
    out.mkdir(parents=True,exist_ok=True);tests=[];rates=[];spatial=[];raw=[]
    tree=SpeciesTree('((A:2,B:2)AB:2,(C:2,D:2)CD:2)ROOT;',20,200);states=dict(A=0,B=1,C=0,D=1)
    rec=Map.constant(100,.0002);focal=75;points=(10,25,50,75,90);rng=np.random.default_rng(72822)
    for p in (.2,.5,.8):
        for k in (0.,.01,.1,.5,1.):
            hist=Frozen(p);rho=.0002*focal*k;ref=[];values={v:[] for v in ('legacy','reservoir_v1')};exposures={v:[] for v in values};pairs={v:[] for v in values}
            for i in range(n):
                seeds=np.random.SeedSequence([20261122,int(p*100),int(k*100),i]).generate_state(3)
                g=simulate_genealogy_reference(0,tree,hist,states,rho,1,np.random.default_rng(seeds[0]),False)
                ref.append(legacy_metrics(g,4))
                for j,version in enumerate(values):
                    result=simulate_msrc_arg(tree,rec,hist,0,(0,100),k,version,int(seeds[j+1]),states)
                    assert result.metadata['background_model_version']==version
                    ts=build_tree_sequence(result).simplify();val=metrics(ts.at(focal),4);values[version].append(val)
                    exposure=focal_background_exposure(result,states,focal,100);exposures[version].append(exposure)
                    pairs[version].append([metrics(ts.at(x),4)[1] for x in points])
                    raw.append(dict(p=p,kappa=k,version=version,replicate=i,root_age=val[0],pair_mean=val[1],topology=val[2],switches=exposure[0],dwell0=exposure[1],dwell1=exposure[2]))
            ref=vectors(ref)
            for version in values:
                for key,a in vectors(values[version]).items():
                    b=ref[key];diff=a.mean()-b.mean();se=np.sqrt(a.var(ddof=1)/n+b.var(ddof=1)/n)
                    tol=.05 if key.startswith('q') or key=='D' else .1*b.mean()
                    tests.append(dict(p=p,kappa=k,version=version,metric=key,observed=a.mean(),reference=b.mean(),difference=diff,se=se,ci_low=diff-1.96*se,ci_high=diff+1.96*se,tolerance=tol,pvalue=float(2*norm.sf(abs(diff/se))) if se else float(diff==0)))
                ex=np.array(exposures[version]);target=rho*(p*ex[:,1]+(1-p)*ex[:,2]);residual=ex[:,0]-target;se=residual.std(ddof=1)/np.sqrt(n)
                rates.append(dict(p=p,kappa=k,version=version,mean_switches=ex[:,0].mean(),mean_integrated_target=target.mean(),difference=residual.mean(),se=se,ci_low=residual.mean()-1.96*se,ci_high=residual.mean()+1.96*se,pvalue=float(2*norm.sf(abs(residual.mean()/se))) if se else 1.))
                xy=np.array(pairs[version])
                for j in (1,2,3,4):
                    r=np.corrcoef(xy[:,0],xy[:,j])[0,1];bs=[]
                    for _ in range(250):
                        ix=rng.integers(n,size=n);bs.append(np.corrcoef(xy[ix,0],xy[ix,j])[0,1])
                    lo,hi=np.quantile(bs,[.025,.975])
                    spatial.append(dict(p=p,kappa=k,version=version,start=10,end=points[j],distance=points[j]-10,correlation=r,ci_low=lo,ci_high=hi,n_chromosomes=n))
            print(p,k,'complete',flush=True)
    for rows in (tests,rates):
        order=sorted(range(len(rows)),key=lambda i:rows[i]['pvalue']);prev=0.
        for rank,i in enumerate(order):
            prev=max(prev,min(1.,(len(rows)-rank)*rows[i]['pvalue']));rows[i]['holm_p']=prev
    write(out/'marginal_tests.csv',tests);write(out/'switching_rate_tests.csv',rates);write(out/'autocorrelation.csv',spatial);write(out/'draws.csv',raw)
    summarize(out,n)


def summarize(out,n):
    def read(name):
        with (out/name).open() as f:
            return [{k:(v if k in ('version','metric') else float(v)) for k,v in row.items()} for row in csv.DictReader(f)]
    tests=read('marginal_tests.csv');rates=read('switching_rate_tests.csv')
    summary=dict(seed=20261122,independent_chromosomes_per_version_per_setting=n,settings=15,
                 marginal_comparisons={v:dict(n=sum(r['version']==v for r in tests),holm_significant=sum(r['version']==v and r['holm_p']<.05 for r in tests),practical_CI_pass=sum(r['version']==v and r['ci_low']>=-r['tolerance'] and r['ci_high']<=r['tolerance'] for r in tests)) for v in ('legacy','reservoir_v1')},
                 rate_comparisons={v:dict(n=sum(r['version']==v for r in rates),holm_significant=sum(r['version']==v and r['holm_p']<.05 for r in rates)) for v in ('legacy','reservoir_v1')},
                 scope='constant interior p; explicit origin; no claim about arbitrary fixation/loss paths')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)



if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);a.add_argument('--replicates',type=int,default=1000);a.add_argument('--summarize-only',action='store_true');cfg=a.parse_args();summarize(cfg.output,cfg.replicates) if cfg.summarize_only else main(cfg.output,cfg.replicates)
