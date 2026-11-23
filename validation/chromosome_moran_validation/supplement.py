"""Small-M matched-theory comparisons and separate T200 sensitivity tracks."""
import sys,json,csv,gzip
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tskit
from scipy.stats import beta
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_moran_validation.windows import summarize
root=Path(__file__).parent/'results';s=json.loads((root/'small_M.json').read_text());f,axs=plt.subplots(2,3,figsize=(12,7),sharex=True)
for group in (0,1):
 a=np.array(s['counts'])[group];n=a.sum(axis=1);pred=np.array(s['theory_unconditional' if group==0 else 'theory_conditional1010'])
 for q in range(3):
  obs=a[:,q]/n;lo=np.array([0 if k==0 else beta.ppf(.025,k,nn-k+1) for k,nn in zip(a[:,q],n)]);hi=np.array([1 if k==nn else beta.ppf(.975,k+1,nn-k) for k,nn in zip(a[:,q],n)]);axs[group,q].errorbar(s['positions'],obs,yerr=[obs-lo,hi-obs],fmt='o',label='ARG');axs[group,q].plot(s['positions'],pred[:,q],'k.-',label='full Moran theory');axs[group,q].set_title(('Unconditional' if group==0 else 'Conditional1010')+f'; q{q+1}, n={n[0]}');axs[group,q].legend();axs[group,q].set_xlabel('Position (bp)');axs[group,q].axvspan(30,70,color='grey',alpha=.12)
f.tight_layout();f.savefig(root/'small_M_theory.png',dpi=180);f.savefig(root/'small_M_theory.svg')
p=root/'runs/sensitivity_T200_seed0/chromosome.trees'
if p.exists():
 ts=tskit.load(p);rows=[]
 from validation.chromosome_scalability.worker import topology
 from msrcsim.arg.output.genealogy_export import export_genealogies
 export_genealogies(ts,p.parent/'local_trees.nwk',p.parent/'local_trees.csv',chrom='chr1')
 for t in ts.trees():
  q=topology(t);rows.append([*t.interval,*[float(q==i) for i in range(3)],t.time(t.root),t.time(t.mrca(0,2)),t.time(t.mrca(0,1))])
 a=np.array(rows);f,axs=plt.subplots(2,3,figsize=(12,6),sharex=True);names=['q1','q2','q3','root','TAC','TAB']
 for w in [250000,50000,25000,10000]:
  x,v=summarize(a[:,0],a[:,1],a[:,2:],1e7,w)
  with gzip.open(root/f'T200_{w}.csv.gz','wt',newline='') as out:
   z=csv.writer(out);z.writerow(['position',*names]);z.writerows(zip(x,*v))
  if w==25000:
   for i,ax in enumerate(axs.flat):ax.plot(x/1e6,v[i]);ax.set_title(names[i]);ax.axvspan(3,7,color='grey',alpha=.12);ax.set_xlabel('Position (Mb)')
 f.suptitle('Genuine Moran T200 sensitivity: one unconditioned chromosome, no uncertainty estimate');f.tight_layout();f.savefig(root/'T200_sensitivity.png',dpi=180);f.savefig(root/'T200_sensitivity.svg')
