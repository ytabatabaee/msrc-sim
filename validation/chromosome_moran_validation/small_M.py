"""Independent finite-M multispecies target, followed by genuine prior draws."""
import sys,json,time,csv
from pathlib import Path
import numpy as np
from scipy.stats import beta
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.moran_spatial.theory import load_theory,make_tree,predict,rho_effective
from msrcsim.moran import simulate_frequency_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.arg import simulate_msrc_arg
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.arg.output.tree_sequence import build_tree_sequence
from validation.chromosome_scalability.worker import topology
if __name__=='__main__':
 out=Path(__file__).parent/'results';out.mkdir(exist_ok=True);mods,hashes=load_theory('/Users/ytabatabaee/Desktop/msrc-paper');pi=np.zeros(21);pi[10]=1
 c=dict(M=20,shape='balanced',newick='((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;',root_duration=10,pi_old=pi.tolist(),selection=0,marker=50,suppressed_region=[30,70],kappa=.01,map=dict(length=100,breakpoints=[],rates=[.001]))
 positions=[5,31,50,69,95];theory=[predict(mods,c,rho_effective(c,x)) for x in positions];t=make_tree(c);counts=np.zeros((2,5,3),int);patterns=np.zeros(16,int);failed=[];start=time.perf_counter()
 for i in range(4000):
  rng=np.random.default_rng(np.random.SeedSequence([20261123,20,i,0]));h=simulate_frequency_history(t,Rearrangement('audit','suppressed_recombination','ROOT',0,10,0),rng);states={x:int(rng.random()<h.terminal_frequency(x)) for x in t.taxa};pattern=''.join(map(str,states.values()));patterns[int(pattern,2)]+=1
  try:
   r=simulate_msrc_arg(t,Map.constant(100,.001),h,50,(30,70),.01,'reservoir_v1',np.random.SeedSequence([20261123,20,i,1]),states,record_events=False);ts=build_tree_sequence(r).simplify()
   for j,x in enumerate(positions):
    q=topology(ts.at(x));counts[0,j,q]+=1
    if pattern=='1010':counts[1,j,q]+=1
  except Exception as e:failed.append(dict(seed=i,error=repr(e)))
 result=dict(case=c,theory_hashes=hashes,n_prior=4000,pattern_counts=patterns.tolist(),theory_patterns=theory[0].pattern_probabilities.tolist(),counts=counts.tolist(),failures=failed,seconds=time.perf_counter()-start,positions=positions,theory_unconditional=[r.quartet_probabilities.tolist() for r in theory],theory_conditional1010=[r.conditional_quartet_probabilities[10].tolist() for r in theory])
 (out/'small_M.json').write_text(json.dumps(result,indent=2));print('small M complete',4000,len(failed),flush=True)
