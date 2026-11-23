"""Exact-binomial uncertainty across independent chromosomes, never genomic windows."""
import json,csv
from pathlib import Path
import numpy as np
from scipy.stats import beta
root=Path(__file__).parent/'results';s=json.loads((root/'small_M.json').read_text());rows=[]
for group in (0,1):
 for j,x in enumerate(s['positions']):
  counts=np.array(s['counts'])[group,j];n=int(counts.sum());target=s['theory_unconditional' if group==0 else 'theory_conditional1010'][j]
  for qi,(k,q) in enumerate(zip(counts,target)):
   lo=0 if k==0 else beta.ppf(.025,k,n-k+1);hi=1 if k==n else beta.ppf(.975,k+1,n-k);status='equivalent_within_0.05' if lo-q>=-.05 and hi-q<=.05 else 'demonstrated_difference_beyond_0.05' if lo-q>.05 or hi-q<-.05 else 'unresolved_precision'
   rows.append(dict(condition='unconditional' if group==0 else '1010',position=x,metric=f'q{qi+1}',n=n,observed=k/n,theory=q,difference=k/n-q,ci_low=lo,ci_high=hi,status=status))
for name,rows2 in [('small_M_comparisons.csv',rows)]:
 with (root/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows2[0]));w.writeheader();w.writerows(rows2)
p=json.loads((root/'history_pool.json').read_text());n=p['completed'];k=sum(r['pattern']=='1010' for r in p['rows']);result=dict(fixed_prior_attempts=n,acceptances1010=k,probability=k/n,ci_low=0 if not k else beta.ppf(.025,k,n-k+1),ci_high=1 if k==n else beta.ppf(.975,k+1,n-k),accepted_seeds=[r['seed'] for r in p['rows'] if r['pattern']=='1010'],lost_tip_populations=sum(r['lost_tips'] for r in p['rows']),fixed_tip_populations=sum(r['fixed_tips'] for r in p['rows']),note='tip populations on the same chromosome are correlated; CI is for whole-chromosome terminal pattern',forward_history_seconds=sum(r['seconds'] for r in p['rows']))
(root/'conditioning_audit.json').write_text(json.dumps(result,indent=2));print(result)
