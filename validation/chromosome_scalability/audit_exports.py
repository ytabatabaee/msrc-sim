"""Post-run export audit; this cost is separate from benchmark timings."""
import argparse,csv,json,sys,time
from pathlib import Path
import numpy as np
import tskit
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from msrcsim.species_tree import SpeciesTree
from validation.chromosome_scalability.worker import species,validate

def audit(roots):
    results=[]
    for root in roots:
        for p in sorted(root.glob('*/summary.json')):
            s=json.loads(p.read_text());j=s.get('job',{})
            if s['status']!='complete' or j.get('kind')=='history':continue
            d=p.parent;ts=tskit.load(d/'chromosome.trees')
            valid=validate(ts,species(j['taxa'],j['ne'],j['divergence_step'],j['root_extension']))
            with (d/'local_trees.csv').open() as f:rows=list(csv.DictReader(f))
            lines=(d/'local_trees.nwk').read_text().splitlines()
            assert len(rows)==len(lines)==ts.num_trees
            labels={json.loads(ts.node(u).metadata)['taxon'] for u in ts.samples()}
            prev=0.
            for row,line,gt in zip(rows,lines,ts.trees()):
                left,right=float(row['left']),float(row['right'])
                assert left==prev and right>left and row['chrom']==j.get('chrom','chr1')
                assert (left,right)==tuple(gt.interval) and line==row['newick']
                parsed=SpeciesTree(line,j['ne'],1.,ultrametric_tolerance=1e-7)
                assert set(parsed.taxa)==labels
                for node in parsed.nodes():
                    if node.parent is not None:assert node.parent.age>=node.age-1e-7
                prev=right
            assert prev==j['length']
            # Some early raw files predate the additive time_units metadata fix.
            assert ts.time_units in ('unknown','generations')
            results.append(dict(case=d.name,root=str(root),intervals=ts.num_trees,taxa=len(labels),valid=True,raw_time_units=ts.time_units))
    return results

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--benchmarks',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();start=time.perf_counter();rows=audit(a.benchmarks);a.output.write_text(json.dumps(dict(cases=rows,wall_seconds=time.perf_counter()-start,consumer_integration='not tested; no Positional ASTRAL input contract available'),indent=2)+'\n');print(f'{len(rows)} complete exports passed')
