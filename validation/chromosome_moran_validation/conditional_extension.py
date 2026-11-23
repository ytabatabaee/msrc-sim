"""Sequential rejection extension; initial fixed 64 trials alone estimate P(Z)."""
import sys,json,time,gc,subprocess
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_biological_validation.theory_grid import case
from validation.moran_spatial.theory import make_tree
from msrcsim.moran import simulate_frequency_history
from msrcsim.rearrangement import Rearrangement
root=Path(__file__).parent;out=root/'results';t=make_tree(case());rows=[];start=time.perf_counter()
for i in range(64,128):
 rng=np.random.default_rng(np.random.SeedSequence([20261123,i,0]));s=time.perf_counter();h=simulate_frequency_history(t,Rearrangement('MSRC','suppressed_recombination','ROOT',0,1000,0),rng);states={x:int(rng.random()<h.terminal_frequency(x)) for x in t.taxa};p=''.join(map(str,states.values()));rows.append(dict(seed=i,pattern=p,seconds=time.perf_counter()-s));del h;gc.collect();(out/'conditional_extension.json').write_text(json.dumps(dict(stopping_rule='first additional 1010 or 64 more prior attempts; no genealogy outcomes used',rows=rows,elapsed=time.perf_counter()-start),indent=2))
 if p=='1010':
  j=dict(length=10000000,model='Moran_B',duration=1000,kappa=.01,pattern='sampled',seed=i,history_seed=i,condition='1010',seconds=600)
  subprocess.run([sys.executable,str(root/'worker.py'),json.dumps(j),str(out/'runs'/f'B_seed{i}')],check=True,timeout=620);break
