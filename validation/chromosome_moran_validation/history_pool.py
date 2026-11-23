"""Fixed prior audit. Selection is based only on sampled terminal backgrounds."""
import sys,json,time,gc,csv
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_biological_validation.theory_grid import case
from validation.moran_spatial.theory import make_tree
from msrcsim.moran import simulate_frequency_history
from msrcsim.rearrangement import Rearrangement
if __name__=='__main__':
 out=Path(__file__).parent/'results';out.mkdir(exist_ok=True);rows=[];t=make_tree(case());start=time.perf_counter()
 for i in range(64):
  rng=np.random.default_rng(np.random.SeedSequence([20261123,i,0]));s=time.perf_counter();h=simulate_frequency_history(t,Rearrangement('MSRC','suppressed_recombination','ROOT',0,1000,0),rng);states={x:int(rng.random()<h.terminal_frequency(x)) for x in t.taxa}
  ends={x:h.terminal_frequency(x) for x in t.taxa};rows.append(dict(seed=i,pattern=''.join(map(str,states.values())),seconds=time.perf_counter()-s,records=len(h.records),terminal_frequencies=ends,lost_tips=sum(p==0 for p in ends.values()),fixed_tips=sum(p==1 for p in ends.values())))
  (out/'history_pool.json').write_text(json.dumps(dict(planned=64,completed=len(rows),rows=rows,elapsed=time.perf_counter()-start),indent=2));print(i,rows[-1]['pattern'],flush=True);del h;gc.collect()
