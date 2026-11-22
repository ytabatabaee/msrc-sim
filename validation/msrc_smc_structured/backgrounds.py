"""Independent focal switching/dwell comparison using carrier event histories."""
import argparse,csv,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from msrcsim.smc.structured import ConstantHistory
from msrcsim.structured_coalescent import simulate_genealogy_reference
from msrcsim.spatial_reservoir import simulate_reservoir_arg
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.species_tree import SpeciesTree
from validate import compare


def exposures(events,states,focal=50,origin=100,legacy=False):
    live={i:h for i,h in enumerate(states.values())};last=0.;n=zero=one=0.
    for e in events:
        t=e.backward_time if legacy else e.time
        stop=min(t,origin);dt=max(0,stop-last)
        zero+=dt*sum(h==0 for h in live.values());one+=dt*sum(h==1 for h in live.values());last=max(last,stop)
        kind=e.event_type
        if kind=='species_boundary':continue
        if legacy:
            ids=[int(x) for x in e.lineage_ids_before.split(';')]
            if kind=='switch':
                h=live[ids[0]];live[ids[0]]=1-h;n+=int(t<origin)
            elif kind=='coalescence':
                h=live.pop(ids[0]);live.pop(ids[1]);live[int(e.lineage_id_after)]=h
        else:
            ids=e.lineage_ids
            if kind=='background_transfer' and ids[0] in live:
                h=int(e.background_after);n+=int(t<origin and live[ids[0]]!=h);live[ids[0]]=h
            elif kind=='recombination' and ids[0] in live:
                old=live.pop(ids[0]);side=0 if focal<e.position else 1
                h=int(e.background_after.split(',')[side]);live[ids[side+1]]=h
                n+=int(t<origin and old!=h)
            elif kind=='coalescence':
                h=live.pop(ids[0],None);h2=live.pop(ids[1],None)
                if h is not None or h2 is not None:live[ids[2]]=h if h is not None else h2
        if len(live)==1:break
    return np.array([n,zero,one])


def main(src,out,nrep):
    out.mkdir(parents=True,exist_ok=True);tests=[];raw=[]
    t=SpeciesTree('((A:2,B:2)AB:2,(C:2,D:2)CD:2)ROOT;',20,200);states=dict(A=0,B=1,C=0,D=1)
    for p in (.2,.5,.8):
        for k in (0.,.01,.1,.5,1.):
            case=f'4tax_p{p}_k{k}';file=src/case/'path_summaries.csv'
            with file.open() as f:sdata=list(csv.DictReader(f))[:nrep]
            smc=np.array([[float(r[x]) for x in ('switches','dwell0','dwell1')] for r in sdata])
            values={'SMC':smc,'legacy':[],'reservoir':[]};hist=ConstantHistory(p,100)
            for i in range(nrep):
                seeds=np.random.SeedSequence([88151,int(p*100),int(k*100),i]).generate_state(2)
                g=simulate_genealogy_reference(0,t,hist,states,.0002*50*k,1,np.random.default_rng(seeds[0]),True)
                values['legacy'].append(exposures(g.events,states,legacy=True))
                arg=simulate_reservoir_arg(t,Map.constant(100,.0002),hist,0,(0,100),k,int(seeds[1]),states)
                values['reservoir'].append(exposures(arg.events,states))
            arrays={m:np.asarray(v) for m,v in values.items()}
            ref={key:arrays['legacy'][:,j] for j,key in enumerate(('switches','dwell0','dwell1'))}
            for method in ('SMC','reservoir'):
                vals={key:arrays[method][:,j] for j,key in enumerate(('switches','dwell0','dwell1'))}
                tests.extend(compare(vals,ref,case,method))
            for method,a in arrays.items():
                for i,v in enumerate(a):raw.append(dict(case=case,method=method,replicate=i,switches=v[0],dwell0=v[1],dwell1=v[2]))
            print(case,'done',flush=True)
    order=sorted(range(len(tests)),key=lambda i:tests[i]['pvalue']);prev=0
    for rank,i in enumerate(order):
        prev=max(prev,min(1,(len(tests)-rank)*tests[i]['pvalue']));tests[i]['holm_p']=prev
    for name,rows in [('background_tests.csv',tests),('background_draws.csv',raw)]:
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--input',type=Path,required=True);a.add_argument('--output',type=Path,required=True);a.add_argument('--replicates',type=int,default=1000)
    cfg=a.parse_args();main(cfg.input,cfg.output,cfg.replicates)
