"""Compare the optimized engine to pinned unoptimized Git source, not a second maintained engine."""
import argparse,dataclasses,hashlib,json,pickle,subprocess,sys,time,types
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_scalability.worker import species,Map
from msrcsim.arg.engine import simulate_arg,ARGConfig
from msrcsim.arg.backgrounds.reservoir import ReservoirBackground


def digest(r):
    x=dict(nodes=r.nodes,edges=r.edges,events=[dataclasses.asdict(e) for e in r.events],peak_lineages=r.peak_lineages,peak_segments=r.peak_segments)
    return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--baseline-revision',default='5138370aa455eedc54b694c13f5aba91c9de5fa6');p.add_argument('--history',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    source=subprocess.check_output(['git','show',a.baseline_revision+':msrcsim/arg/engine.py'],text=True)
    module=types.ModuleType('msrcsim.arg._scalability_baseline');module.__package__='msrcsim.arg';sys.modules[module.__name__]=module
    exec(compile(source,'/baseline/arg/engine.py','exec'),module.__dict__)
    rows=[]
    for n,length,kind in ((4,500000,'neutral_uniform'),(8,100000,'neutral_uniform'),(4,500000,'neutral_low'),(4,100000,'MSRC_ARG')):
        t=species(n,1000,100,50);rec=Map(length,(.3*length,.7*length),(1e-8,1e-10,1e-8)) if kind=='neutral_low' else Map.constant(length,1e-8)
        values=[];seconds=[]
        for fn,configclass in ((module.simulate_arg,module.ARGConfig),(simulate_arg,ARGConfig)):
            states=None;bg='neutral'
            if kind=='MSRC_ARG':
                with a.history.open('rb') as f:h,states=pickle.load(f)
                bg=ReservoirBackground(h,0,(.3*length,.7*length),.01,states)
            cfg=configclass(length,rec,bg,sampled_backgrounds=states);start=time.perf_counter()
            r=fn(t,cfg,np.random.default_rng(np.random.SeedSequence([101,n,1])));seconds.append(time.perf_counter()-start);values.append(digest(r))
        row=dict(taxa=n,length=length,model=kind,baseline_seconds=seconds[0],optimized_seconds=seconds[1],baseline_digest=values[0],optimized_digest=values[1],identical=values[0]==values[1]);rows.append(row);print(row,flush=True)
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(dict(baseline=a.baseline_revision,comparisons=rows),indent=2)+'\n')
        assert values[0]==values[1], 'scientific result changed; do not weaken this test'
