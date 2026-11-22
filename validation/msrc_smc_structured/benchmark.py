"""Diagnostic costs for four distinct simulator kernels; no speed assumption."""
import argparse,csv,cProfile,pstats,sys,time,tracemalloc
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from msrcsim.species_tree import SpeciesTree
from msrcsim.smc.structured import StructuredConfig,ConstantHistory,simulate_structured_smc
from msrcsim.smc.engine import SMCConfig,simulate_msc_smc
from msrcsim.smc.recombination_map import PiecewiseConstantMap as Map
from msrcsim.spatial_reservoir import simulate_reservoir_arg
from msrcsim.arg.engine import ARGConfig,simulate_arg
from msrcsim.arg.backgrounds.msrc_two_state import MSRCTwoBackground
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as AMap


def species(n):
    nodes=[(f'T{i}',0.) for i in range(n)]
    while len(nodes)>1:
        age=nodes[0][1]+2
        nodes=[(f'({nodes[i][0]}:{age-nodes[i][1]},{nodes[i+1][0]}:{age-nodes[i+1][1]})',age) for i in range(0,len(nodes),2)]
    return SpeciesTree(nodes[0][0]+'ROOT;',20,200)


class VariableHistory(ConstantHistory):
    def frequency_at(self,pop,t):return (.2,.5,.8,.5,.2)[min(int(t//20),4)]
    def next_frequency_boundary(self,pop,t):
        return 20*(int(t//20)+1) if t<80 else float('inf')


def main(out):
    out.mkdir(parents=True,exist_ok=True);rows=[]
    for n in (4,8,16,32):
        t=species(n);states={taxon:i%2 for i,taxon in enumerate(t.taxa)}
        for L in (100,1000):
            for name in ('structured_SMC','neutral_SMC','original_ARG','reservoir_ARG'):
                rec=Map.constant(L,1e-5);ar=AMap.constant(L,1e-5);h=ConstantHistory(.5,100)
                rng=np.random.default_rng(991+n);tracemalloc.start();start=time.perf_counter()
                if name=='structured_SMC':r=simulate_structured_smc(t,h,StructuredConfig(rec,0,(0,L),.1,991+n),states)
                elif name=='neutral_SMC':r=simulate_msc_smc(t,SMCConfig(L,rec,991+n))
                elif name=='original_ARG':r=simulate_arg(t,ARGConfig(L,ar,MSRCTwoBackground(h,0,(0,L),.1,states)),rng)
                else:r=simulate_reservoir_arg(t,ar,h,0,(0,L),.1,991+n,states)
                duration=time.perf_counter()-start;_,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
                if name=='structured_SMC':
                    marks=len(r.events);changed=sum(e['genealogy_changed'] for e in r.events);distinct=1+changed
                elif name=='neutral_SMC':
                    marks=len(r.events);changed=sum(e.event_type!='no_genealogy_change' for e in r.events);distinct=1+changed
                else:
                    marks=sum(e.event_type in ('recombination','background_transfer') for e in r.events)
                    from msrcsim.arg.output.tree_sequence import build_tree_sequence
                    distinct=build_tree_sequence(r).simplify().num_trees;changed=distinct-1
                rows.append(dict(taxa=n,length=L,method=name,history_changes=0,seconds=duration,peak_python_bytes=peak,
                                 marks=marks,seconds_per_mark=duration/marks if marks else None,adjacent_distinct_genealogies=distinct,genealogy_changes=changed))
                print(n,L,name,round(duration,3),flush=True)
                with (out/'benchmark.csv').open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    t=species(8);states={taxon:i%2 for i,taxon in enumerate(t.taxa)}
    for h in (ConstantHistory(.5,100),VariableHistory(.5,100)):
        start=time.perf_counter();r=simulate_structured_smc(t,h,StructuredConfig(Map.constant(1000,1e-5),0,(0,1000),.1,53),states)
        rows.append(dict(taxa=8,length=1000,method='structured_history_complexity',history_changes=4 if isinstance(h,VariableHistory) else 0,
                         seconds=time.perf_counter()-start,peak_python_bytes=None,marks=len(r.events),seconds_per_mark=None,adjacent_distinct_genealogies=1+sum(e['genealogy_changed'] for e in r.events),genealogy_changes=sum(e['genealogy_changed'] for e in r.events)))
    with (out/'benchmark.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    prof=cProfile.Profile();prof.runcall(simulate_structured_smc,species(16),ConstantHistory(.5,100),StructuredConfig(Map.constant(1000,1e-5),0,(0,1000),.1,53))
    with (out/'profile.txt').open('w') as f:pstats.Stats(prof,stream=f).strip_dirs().sort_stats('cumtime').print_stats(20)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);cfg=a.parse_args();main(cfg.output)
