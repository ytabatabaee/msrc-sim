"""Shared realized Moran/WF histories: admitted-path temporal comparison."""
import argparse,csv,sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from msrcsim.smc.structured import StructuredConfig,simulate_structured_smc
from msrcsim.smc.recombination_map import PiecewiseConstantMap as Map
from msrcsim.species_tree import SpeciesTree
from msrcsim.population_process import simulate_population_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.structured_coalescent import simulate_genealogy_reference
from validate import gene_metrics,legacy_metrics,metric_vectors,compare,focal_row


def main(out,nrep):
    out.mkdir(parents=True,exist_ok=True);tests=[];info=[]
    t=SpeciesTree('((A:2,B:2)AB:2,(C:2,D:2)CD:2)ROOT;',20,10)
    rr=Rearrangement('background_1','suppressed_recombination','ROOT',6,20,0.)
    states=dict(A=0,B=1,C=0,D=1)
    for name in ('moran','wright_fisher'):
        history=simulate_population_history(t,rr,np.random.default_rng(93),{'population_process':{'model':name}})
        records=repr(history.records);config=StructuredConfig(Map.constant(100,.0002),0,(20,80),.1,reference_position=50)
        rho=.0002*(20+30*.1);values=[];ref=[]
        for i in range(nrep):
            config.seed=19500+i
            result=simulate_structured_smc(t,history,config,states)
            root=SpeciesTree(focal_row(result,50)['newick'],20,200,ultrametric_tolerance=1e-7).root
            values.append(gene_metrics(root))
            g=simulate_genealogy_reference(0,t,history,states,rho,1,np.random.default_rng(39500+i),False)
            ref.append(legacy_metrics(g,4))
        assert repr(history.records)==records
        tests.extend(compare(metric_vectors(values,4),metric_vectors(ref,4),name,'SMC'))
        info.append(dict(history=name,origin=history.structured_origin_age,records=len(history.records),replicates=nrep,history_seed=93))
        # Preserve the exact shared trajectory for independent replay.
        from msrcsim.arg.engine import _write_frequency_history
        _write_frequency_history(history,out/f'{name}_frequency_history.csv')
        print(name,nrep,'complete',flush=True)
    with (out/'temporal_tests.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(tests[0]),lineterminator='\n');w.writeheader();w.writerows(tests)
    (out/'histories.json').write_text(json.dumps(info,indent=2)+'\n')


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);a.add_argument('--replicates',type=int,default=500)
    cfg=a.parse_args();main(cfg.output,cfg.replicates)
