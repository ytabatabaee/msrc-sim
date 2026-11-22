"""Replay the exact prior seeds to document frequency dynamics and absorption."""
import argparse,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.moran_spatial.run import make_tree,moran_history

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--replicates',type=int,default=4000);a=p.parse_args();rows=[]
    for file in a.config:
        cfg=json.loads(file.read_text())
        for c in cfg['cases']:
            if c['history']!='unconditional_moran':continue
            t=make_tree(c);loss=fix=dynamic=0
            for i in range(a.replicates):
                seed=np.random.SeedSequence([cfg['seed'],c['case_number'],i]).generate_state(5)[0]
                h=moran_history(c,t,np.random.default_rng(seed))
                records=[r for r in h.records if r.absolute_age<h.structured_origin_age-1e-10]
                loss+=int(any(r.copy_count_A1==0 for r in records))
                fix+=int(any(r.copy_count_A1==r.population_chromosomes for r in records))
                dynamic+=int(any(len({r.copy_count_A1 for r in records if r.branch_id==bid})>1 for bid in t.branches))
            rows.append(dict(case=c['id'],n=a.replicates,histories_with_loss=loss,histories_with_fixation=fix,histories_with_within_branch_count_changes=dynamic,scope='genuine structured intervals only; pre-origin annotation excluded; exact prior seeds replayed'))
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
