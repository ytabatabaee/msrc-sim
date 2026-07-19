from __future__ import annotations
import argparse, numpy as np
from .config import load_config
from .species_tree import SpeciesTree
from .rearrangement import Rearrangement
from .wright_fisher import simulate_frequency_history
from .structured_coalescent import simulate_genealogy
from .io import write_outputs

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',required=True); a=ap.parse_args(); c=load_config(a.config); rng=np.random.default_rng(c['seed'])
    st=c['species_tree']; tree=SpeciesTree(st['newick'],st['default_effective_population_size'],st['root_extension'],st.get('branch_parameters'))
    rr=c['rearrangement']; sel=rr.get('selection',{}).get('coefficient',rr.get('selection_coefficient',0.0))
    rearr=Rearrangement(rr.get('id','inv_1'),rr['type'],rr['origin_branch'],int(rr['origin_time_from_branch_start']),int(rr.get('initial_copy_count',1)),float(sel))
    hist=simulate_frequency_history(tree,rearr,rng)
    sampled={t:int(rng.random()<hist.terminal_frequency(t)) for t in tree.taxa}
    rec=c['recombination']; base=rec.get('baseline_rate',rec.get('rate')); frac=rec.get('effective_cross_arrangement_fraction',rec.get('suppression_factor'))
    logopt=c['output']['event_log_loci']; first_n=logopt.get('first_n',c['num_loci']) if isinstance(logopt,dict) else c['num_loci']
    results=[simulate_genealogy(i,tree,hist,sampled,float(base),float(frac),rng,c['output']['record_backward_events'] and i<first_n) for i in range(c['num_loci'])]
    out=write_outputs(c,tree,rearr,hist,sampled,results); print(f"Wrote {out}")
if __name__=='__main__': main()
