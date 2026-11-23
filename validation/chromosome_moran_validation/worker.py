"""Bounded genuine shared-history Moran experiment; no ARG generator modifications."""
import sys,json,time,csv,signal,resource,traceback,os
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_biological_validation.theory_grid import case
from validation.moran_spatial.theory import make_tree,PrescribedHistory
from validation.chromosome_scalability.worker import validate,topology
from msrcsim.arg import simulate_msrc_arg
from msrcsim.moran import simulate_frequency_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.arg.engine import simulate_arg,ARGConfig
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.arg.output.tree_sequence import build_tree_sequence

def run(job,out):
 out.mkdir(parents=True,exist_ok=True)
 try:fd=os.open(out/'worker.lock',os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(fd)
 except FileExistsError:return
 start=time.perf_counter();summary={'job':job,'status':'running'}
 def stop(*args):raise TimeoutError('configured time or 1.5 GiB limit')
 def tick(sig,frame):
  f=frame
  while f:
   if f.f_code.co_name=='simulate_arg':
    loc=f.f_locals;summary['partial_generation']=loc.get('t');summary['partial_carriers']=len(loc.get('lineages',[]));break
   f=f.f_back
  (out/'heartbeat.json').write_text(json.dumps(dict(seconds=time.perf_counter()-start,generation=summary.get('partial_generation'),carriers=summary.get('partial_carriers'))))
  if time.perf_counter()-start>job.get('seconds',180) or resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>1.5*2**30:stop()
 signal.signal(signal.SIGALRM,tick);signal.setitimer(signal.ITIMER_REAL,.5,.5)
 try:
  L=job['length'];c=case(job['duration'],job['kappa'],job.get('marker_fraction',.5)*L,job.get('rate',1e-8));c['suppressed_region']=[L*.3,L*.7];t=make_tree(c)
  s0=time.perf_counter();rng=np.random.default_rng(np.random.SeedSequence([20261123,job['history_seed'],0]));h=simulate_frequency_history(t,Rearrangement('MSRC','suppressed_recombination','ROOT',0,1000,job.get('selection',0)),rng)
  states={taxon:int(rng.random()<h.terminal_frequency(taxon)) for taxon in t.taxa};pattern=''.join(map(str,states.values()))
  if job.get('condition') is not None:assert pattern==job['condition'],'posterior selection/replay mismatch'
  summary.update(history_seconds=time.perf_counter()-s0,history_records=len(h.records),terminal_pattern=pattern,terminal_frequencies={taxon:h.terminal_frequency(taxon) for taxon in t.taxa},origin=h.structured_origin_age)
  traces=[]
  for bid,b in t.branches.items():
   for age in np.linspace(b.younger_age,b.older_age,101):traces.append(dict(population=bid,age=age,p=h.frequency_at(bid,age)))
  with (out/'history_trace.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(traces[0]));w.writeheader();w.writerows(traces)

  rec=Map(L,tuple(c['suppressed_region']),(c['map']['rates'][0],c['map']['rates'][0]*job['kappa'],c['map']['rates'][0])) if job['model']=='B' else Map.constant(L,c['map']['rates'][0])
  if job['model'] in ('A','B'):result=simulate_arg(t,ARGConfig(L,rec),np.random.default_rng(job['seed']))
  else:result=simulate_msrc_arg(t,rec,h,c['marker'],c['suppressed_region'],job['kappa'],'reservoir_v1',np.random.SeedSequence([20261123,job['history_seed'],1]),states)
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
 j=json.loads(sys.argv[1]);run(j,Path(sys.argv[2]))
