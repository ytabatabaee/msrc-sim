"""Suppression-boundary spatial audit with a nonzero constant marker gap."""
import argparse,csv,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from msrcsim.smc.structured import ConstantHistory,StructuredConfig,simulate_structured_smc,Kernel
from msrcsim.smc.recombination_map import PiecewiseConstantMap as Map
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as AMap
from msrcsim.spatial_reservoir import simulate_reservoir_arg
from msrcsim.arg.output.tree_sequence import build_tree_sequence
from msrcsim.species_tree import SpeciesTree
from validate import gene_metrics,ts_metrics,focal_row


def main(out,nrep):
    out.mkdir(parents=True,exist_ok=True);rows=[];jumps=[]
    tree=SpeciesTree('((A:2,B:2)AB:2,(C:2,D:2)CD:2)ROOT;',20,200);states=dict(A=0,B=1,C=0,D=1);history=ConstantHistory(.5,100)
    points=(275,325,500,675,725);rec=Map(1000,(300.,700.),(1e-5,.0002,1e-5));arec=AMap(1000,(300.,700.),(1e-5,.0002,1e-5))
    for k in (0.,.01,.1,.5,1.):
        cfg=StructuredConfig(rec,0,(300,700),k,start=250,end=750,reference_position=500)
        for i in range(nrep):
            seeds=np.random.SeedSequence([71821,int(k*100),i]).generate_state(2);cfg.seed=int(seeds[0])
            result=simulate_structured_smc(tree,history,cfg,states)
            arg=build_tree_sequence(simulate_reservoir_arg(tree,arec,history,0,(300,700),k,int(seeds[1]),states)).simplify()
            for x in points:
                root=SpeciesTree(focal_row(result,x)['newick'],20,200,ultrametric_tolerance=1e-7).root
                a=gene_metrics(root);b=ts_metrics(arg.at(x),4)
                rows.append(dict(kappa=k,replicate=i,position=x,smc_root=a[0],smc_pair=a[1],smc_top=a[3],arg_root=b[0],arg_pair=b[1],arg_top=b[3],rho_actual=Kernel(tree,history,cfg).effective_rate(x),rho_reference=result.metadata['rho_reference']))
            jumps.append(dict(kappa=k,replicate=i,proposals=len(result.events),genealogy_changes=sum(e['genealogy_changed'] for e in result.events),background_changes=sum(e['background_changed'] for e in result.events)))
        print(k,nrep,'complete',flush=True)
    summary=[];rng=np.random.default_rng(5321)
    for k in (0.,.01,.1,.5,1.):
        for x,y in ((275,325),(325,500),(325,675),(275,725)):
            left=[r for r in rows if r['kappa']==k and r['position']==x];right=[r for r in rows if r['kappa']==k and r['position']==y]
            a=np.array([(l['smc_pair'],r['smc_pair']) for l,r in zip(left,right)]);b=np.array([(l['arg_pair'],r['arg_pair']) for l,r in zip(left,right)])
            av=np.corrcoef(a.T)[0,1];bv=np.corrcoef(b.T)[0,1];bs=[]
            for _ in range(250):
                ix=rng.integers(nrep,size=nrep);iy=rng.integers(nrep,size=nrep);bs.append(np.corrcoef(a[ix].T)[0,1]-np.corrcoef(b[iy].T)[0,1])
            summary.append(dict(kappa=k,start=x,end=y,distance=y-x,smc_correlation=av,arg_correlation=bv,difference=av-bv,ci_low=np.quantile(bs,.025),ci_high=np.quantile(bs,.975),n_chromosomes=nrep))
    for name,data in [('window_spatial_summary.csv',summary),('window_draws.csv',rows),('window_events.csv',jumps)]:
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);a.add_argument('--replicates',type=int,default=500)
    cfg=a.parse_args();main(cfg.output,cfg.replicates)
