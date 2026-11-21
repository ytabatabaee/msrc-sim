"""Small diagnostic benchmark; no genome-scale speed claim."""
import csv,sys,time,tracemalloc,cProfile,pstats
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from msrcsim.species_tree import SpeciesTree
from msrcsim.spatial_reservoir import simulate_reservoir_arg
from msrcsim.arg.engine import ARGConfig,simulate_arg
from msrcsim.arg.backgrounds.msrc_two_state import MSRCTwoBackground
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.smc.engine import SMCConfig,simulate_msc_smc
from msrcsim.smc.recombination_map import PiecewiseConstantMap as SMap
from validate import Frozen


def species(n):
    nodes=[(f'T{i}',0) for i in range(n)]
    while len(nodes)>1:
        age=nodes[0][1]+2
        nodes=[(f'({nodes[i][0]}:{age-nodes[i][1]},{nodes[i+1][0]}:{age-nodes[i+1][1]})',age) for i in range(0,len(nodes),2)]
    return SpeciesTree(nodes[0][0]+'ROOT;',20,200)


def main(out):
    out.mkdir(exist_ok=True,parents=True); rows=[]
    for n in (4,8,16,32):
        tree=species(n);hist=Frozen(.5);states={t:i%2 for i,t in enumerate(tree.taxa)}
        for length in (100,1000):
            rate=1e-5;rec=Map.constant(length,rate)
            for rep in range(2):
                seed=31415+rep+n
                for name in ('reservoir','current_ARG','neutral_SMC_prime'):
                    tracemalloc.start();start=time.perf_counter()
                    if name=='reservoir':
                        result=simulate_reservoir_arg(tree,rec,hist,0,(0,length),.1,seed,states)
                    elif name=='current_ARG':
                        model=MSRCTwoBackground(hist,0,(0,length),.1,states)
                        result=simulate_arg(tree,ARGConfig(length,rec,model),np.random.default_rng(seed))
                    else:
                        result=simulate_msc_smc(tree,SMCConfig(length,SMap.constant(length,rate),seed))
                    duration=time.perf_counter()-start;_,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
                    events=len(result.events)
                    rows.append(dict(taxa=n,length=length,rate=rate,method=name,seed=seed,runtime_seconds=duration,peak_python_bytes=peak,logged_events=events))
            print(n,length,'done',flush=True)
    with (out/'benchmark.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    profiler=cProfile.Profile()
    profiler.runcall(simulate_reservoir_arg,species(32),Map.constant(1000,1e-5),Frozen(.5),0,(0,1000),.1,31447)
    with (out/'profile.txt').open('w') as f:
        pstats.Stats(profiler,stream=f).strip_dirs().sort_stats('cumulative').print_stats(18)


if __name__=='__main__': main(Path(sys.argv[1]))
