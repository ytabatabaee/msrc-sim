"""Bounded fixed-frequency positive control; no ARG generator modifications."""
import sys,json,time,csv,signal,resource,traceback,os
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_biological_validation.theory_grid import case
from validation.moran_spatial.theory import make_tree,PrescribedHistory
from validation.chromosome_scalability.worker import validate,topology
from msrcsim.arg import simulate_msrc_arg
from msrcsim.arg.engine import simulate_arg,ARGConfig
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.arg.output.tree_sequence import build_tree_sequence

def run(job,out):
 out.mkdir(parents=True,exist_ok=True)
 try:fd=os.open(out/'worker.lock',os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(fd)
 except FileExistsError:return
 start=time.perf_counter();summary={'job':job,'status':'running'}
 def stop(*args):raise TimeoutError('120 seconds or 1.5 GiB limit')
 def tick(*args):
  if time.perf_counter()-start>120 or resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>1.5*2**30:stop()
 signal.signal(signal.SIGALRM,tick);signal.setitimer(signal.ITIMER_REAL,.5,.5)
 try:
  L=job['length'];c=case(job['duration'],job['kappa'],job.get('marker_fraction',.5)*L,job.get('rate',1e-8));c['suppressed_region']=[L*.3,L*.7];t=make_tree(c)
  h=PrescribedHistory(t,{bid:[(.5,b.older_age-b.younger_age)] for bid,b in t.branches.items()});states=dict(zip('ABCD',map(int,job['pattern'])))
  rec=Map(L,tuple(c['suppressed_region']),(c['map']['rates'][0],c['map']['rates'][0]*job['kappa'],c['map']['rates'][0])) if job['model']=='B' else Map.constant(L,c['map']['rates'][0])
  if job['model'] in ('A','B'):result=simulate_arg(t,ARGConfig(L,rec),np.random.default_rng(job['seed']))
  else:result=simulate_msrc_arg(t,rec,h,c['marker'],c['suppressed_region'],job['kappa'],'reservoir_v1',job['seed'],states)
  ts=build_tree_sequence(result).simplify();validate(ts,t);ts.dump(out/'chromosome.trees');rows=[]
  for gt in ts.trees():
   pairs=[gt.time(gt.mrca(i,j)) for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))]
   internal=sum(gt.branch_length(u) for u in gt.nodes() if gt.parent(u)!=-1 and gt.num_samples(u)==2)
   rows.append(dict(left=gt.interval.left,right=gt.interval.right,topology=topology(gt),root=gt.time(gt.root),internal=internal,mean_pair=np.mean(pairs),TAB=pairs[0],TAC=pairs[1],TAD=pairs[2],TBD=pairs[4]))
  with (out/'trees.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  from collections import Counter
  summary.update(status='complete',events=dict(Counter(e.event_type for e in result.events)),intervals=ts.num_trees,peak_carriers=result.peak_lineages)
 except Exception as e:summary.update(status='incomplete',error=repr(e),traceback=traceback.format_exc())
 finally:
  signal.setitimer(signal.ITIMER_REAL,0);summary.update(seconds=time.perf_counter()-start,rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);(out/'summary.json').write_text(json.dumps(summary,indent=2))
 print(out.name,summary['status'],round(summary['seconds'],2),flush=True)
if __name__=='__main__':
 if len(sys.argv)>1:
  j=json.loads(sys.argv[1]);run(j,Path(sys.argv[2]));sys.exit()
 import subprocess
 root=Path(__file__).parent/'results/runs'
 jobs=[]
 for L in (100000,1000000,10000000):
  for model,k in [('A',1.),('B',.01),('C',.01),('D',1.)]:
   for seed in range(100,120) if L==10000000 else (100,):jobs.append(dict(length=L,model=model,kappa=k,duration=1000,pattern='1010',seed=seed))
 for j in jobs:
  name=f"{j['model']}_L{j['length']}_s{j['seed']}";out=root/name
  if (out/'summary.json').exists():continue
  subprocess.run([sys.executable,__file__,json.dumps(j),str(out)],timeout=130,check=True)
