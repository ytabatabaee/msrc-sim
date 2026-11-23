"""Progressive bounded cohorts; rejection occurs before any genealogy, without failure replacement."""
import argparse,json,subprocess,sys,time,gc
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_biological_validation.theory_grid import case
from validation.moran_spatial.theory import make_tree
from msrcsim.moran import simulate_frequency_history
from msrcsim.rearrangement import Rearrangement
root=Path(__file__).parent
p=argparse.ArgumentParser();p.add_argument('--replicates',type=int,default=2);p.add_argument('--seconds',type=int,default=600);p.add_argument('--max-extra-priors',type=int,default=6400);p.add_argument('--history-seconds',type=int,default=3600);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
if a.replicates<1:raise SystemExit('replicates must be positive')
a.output.mkdir(parents=True,exist_ok=True)
# These exact draws were already audited. For larger cohorts continue the same
# prior seed stream, before any genealogy, and retain every rejection.
selected=[30,64][:a.replicates];rows=[];t=make_tree(case());start=time.perf_counter()
for i in range(65,65+a.max_extra_priors):
 if len(selected)>=a.replicates or time.perf_counter()-start>a.history_seconds:break
 rng=np.random.default_rng(np.random.SeedSequence([20261123,i,0]));h=simulate_frequency_history(t,Rearrangement('MSRC','suppressed_recombination','ROOT',0,1000,0),rng);states={x:int(rng.random()<h.terminal_frequency(x)) for x in t.taxa};pattern=''.join(map(str,states.values()));rows.append(dict(seed=i,pattern=pattern));del h;gc.collect()
 if pattern=='1010':selected.append(i)
 (a.output/'rejection_extension.json').write_text(json.dumps(dict(rows=rows,selected=selected,stopping='target acceptances, maximum trials or history-time budget; not a fixed-n probability estimate'),indent=2))
(a.output/'cohort_plan.json').write_text(json.dumps(dict(requested=a.replicates,unconditional_seeds=list(range(a.replicates)),conditional_seeds=selected,history_seconds=time.perf_counter()-start,per_genealogy_seconds=a.seconds,no_failed_genealogy_replacement=True),indent=2))
for kind,seeds in [('A',list(range(a.replicates))),('B',selected)]:
 for i in seeds:
  j=dict(length=10000000,model='Moran_'+kind,duration=1000,kappa=.01,pattern='sampled',seed=i,history_seed=i,seconds=a.seconds)
  if kind=='B':j['condition']='1010'
  out=a.output/f'{kind}_seed{i}'
  if (out/'summary.json').exists():continue
  subprocess.run([sys.executable,str(root/'worker.py'),json.dumps(j),str(out)],timeout=a.seconds+30,check=True)
