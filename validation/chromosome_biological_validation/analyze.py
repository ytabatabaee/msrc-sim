"""Independent chromosomes are the uncertainty unit; genealogy spans weight windows."""
import sys,json,csv
from functools import lru_cache
from pathlib import Path
import numpy as np
import tskit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import beta
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_biological_validation.theory_grid import case
from validation.moran_spatial.theory import load_theory,make_tree,predict,rho_effective,PrescribedHistory
root=Path(__file__).parent/'results';mods,_=load_theory('/Users/ytabatabaee/Desktop/msrc-paper')
positions=np.array([0,1e6,2.9e6,3.1e6,5e6,6.9e6,7.1e6,9e6,9999999.])
summary=[];cohorts={};pointrows=[]
for p in sorted((root/'runs').glob('*/summary.json')):
 s=json.loads(p.read_text());summary.append(dict(case=p.parent.name,**{k:v for k,v in s.items() if k not in ('job','traceback')},**s['job']))
 if s['status']!='complete' or s['job']['length']!=10000000:continue
 with (p.parent/'trees.csv').open() as f:r=list(csv.DictReader(f))
 rows={k:np.array([float(v[k]) for v in r]) for k in r[0]}
 ts=tskit.load(p.parent/'chromosome.trees');tracts=[];prev=None
 for gt in ts.trees():
  sig=tuple(gt.time(gt.mrca(i,j)) for i in range(4) for j in range(i))
  if sig==prev:tracts[-1]+=gt.span
  else:tracts.append(gt.span)
  prev=sig
 rows['_tracts']=np.array(tracts);cohorts.setdefault(s['job']['model'],[]).append(rows)

@lru_cache(None)
def target(model,x):
 if model in ('A','B'):return np.array([.3967750546426936,.3016124726786532,.3016124726786532])
 c=case(1000,.01 if model=='C' else 1.);t=make_tree(c);h=PrescribedHistory(t,{bid:[(.5,b.older_age-b.younger_age)] for bid,b in t.branches.items()})
 return predict(mods,c,rho_effective(c,x),h).conditional_quartet_probabilities[10]

fig,axes=plt.subplots(4,4,figsize=(15,10),sharex=True)
bins=np.linspace(0,1e7,41);centers=(bins[:-1]+bins[1:])/2
fig2,ax2=plt.subplots(6,4,figsize=(15,14),sharex=True)
fig3,ax3=plt.subplots(1,3,figsize=(14,4))
for col,(model,runs) in enumerate(sorted(cohorts.items())):
 if col>=4:continue
 n=len(runs);axes[0,col].set_title(f'{model}; {n} chromosomes');ax2[0,col].set_title(model)
 for qi in range(4):
  a=axes[qi,col];vals=[]
  for r in runs:
   weights=np.maximum(0,np.minimum(r['right'][:,None],bins[None,1:])-np.maximum(r['left'][:,None],bins[None,:-1]))/np.diff(bins)
   value=(r['topology']==qi) if qi<3 else (r['topology']==1).astype(float)-(r['topology']==0)
   vals.append(value@weights)
  v=np.array(vals);a.plot(centers/1e6,v.mean(0),label='ARG window mean');a.fill_between(centers/1e6,v.min(0),v.max(0),alpha=.15)
  curve=np.unique(np.r_[np.linspace(0,1e7,101),5e6+np.array([-2e5,-1e5,-5e4,-1e4,1e4,5e4,1e5,2e5])]);q=np.array([target(model,x) for x in curve]);a.plot(curve/1e6,q[:,qi] if qi<3 else q[:,1]-q[:,0],'k.-',label='point theory');a.axvspan(3,7,color='grey',alpha=.12)
  pointq=np.array([[int(r['topology'][np.searchsorted(r['right'],x,side='right')]) for x in positions] for r in runs]);pointvalue=np.mean(pointq==qi,axis=0) if qi<3 else np.mean((pointq==1).astype(float)-(pointq==0),axis=0);a.plot(positions/1e6,pointvalue,'s',color='tab:red',markersize=3,label='ARG at selected points')
  if col==0:a.set_ylabel(('q1','q2','q3','q2-q1')[qi])
  if qi==3:a.set_xlabel('Position (Mb)')
 for row,key in enumerate(('root','internal','mean_pair','TAC','TBD','TAB')):
  v=[]
  for r in runs:
   weights=np.maximum(0,np.minimum(r['right'][:,None],bins[None,1:])-np.maximum(r['left'][:,None],bins[None,:-1]))/np.diff(bins);v.append(r[key]@weights)
  ax2[row,col].plot(centers/1e6,np.mean(v,0));ax2[row,col].fill_between(centers/1e6,np.min(v,0),np.max(v,0),alpha=.15);ax2[row,col].axvspan(3,7,color='grey',alpha=.12)
  if col==0:ax2[row,col].set_ylabel(key+' (generations)')
 for x in positions:
  draws=[int(r['topology'][np.searchsorted(r['right'],x,side='right')]) for r in runs];q=target(model,x)
  for qi in range(3):
   k=draws.count(qi);lo=0 if k==0 else beta.ppf(.025,k,n-k+1);hi=1 if k==n else beta.ppf(.975,k+1,n-k)
   pointrows.append(dict(model=model,position=x,metric=f'q{qi+1}',n=n,count=k,observed=k/n,theory=q[qi],ci_low=lo,ci_high=hi,difference=k/n-q[qi],status='incompatible_95CI' if q[qi]<lo or q[qi]>hi else 'not_rejected_not_equivalence'))
 grid=np.arange(0,1e7,10000);lags=np.array([1,5,10,50,100,200]);persist=[];corr=[];tracts=[]
 for r in runs:
  ix=np.searchsorted(r['right'],grid,side='right');tops=r['topology'][ix];times=r['mean_pair'][ix]
  persist.append([np.mean(tops[:-d]==tops[d:]) for d in lags]);corr.append([np.corrcoef(times[:-d],times[d:])[0,1] for d in lags]);tracts.extend(r['_tracts'][1:-1])
 ax3[0].plot(lags*.01,np.mean(persist,0),label=model);ax3[1].plot(lags*.01,np.mean(corr,0),label=model);v=np.sort(tracts);ax3[2].plot(v/1e6,np.arange(1,len(v)+1)/len(v),label=model)
ax3[2].set_xscale('log')
for col in range(4):ax2[-1,col].set_xlabel('Position (Mb)')
axes[0,0].legend(fontsize=7);fig.suptitle('Conditional fixed p=0.5, Z=1010; bands: chromosome ranges, not CIs')
fig2.suptitle('Simulation branch lengths; fixed p=0.5, Z=1010; observed chromosome ranges')
for a in ax3:a.legend();a.set_xlabel('Distance / tract (Mb)')
ax3[0].set_ylabel('Topology agreement');ax3[1].set_ylabel('Pair-time correlation');ax3[2].set_ylabel('Genealogy tract ECDF (censored ends excluded)')
for f,name in ((fig,'quartets'),(fig2,'branches'),(fig3,'spatial')):
 f.tight_layout();f.savefig(root/(name+'.png'),dpi=180);f.savefig(root/(name+'.svg'))
for name,rows in [('point_estimates',pointrows),('run_summary',summary)]:
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (root/(name+'.csv')).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
# Independent theory sensitivity, existing precomputed grid.
with (root/'theory_grid.csv').open() as f:rs=list(csv.DictReader(f))
f,axs=plt.subplots(1,3,figsize=(13,4))
for duration in (50,200,500,1000,2000):
 a=[r for r in rs if r['pattern']=='1010' and float(r['marker'])==5e6 and float(r['kappa'])==.01 and int(r['duration'])==duration];axs[0].plot([float(r['position'])/1e6 for r in a],[float(r['q2']) for r in a],label=str(duration))
for k in (0,.01,.1,1):
 a=[r for r in rs if r['pattern']=='1010' and float(r['marker'])==5e6 and float(r['kappa'])==k and int(r['duration'])==1000];axs[1].plot([float(r['position'])/1e6 for r in a],[float(r['q2']) for r in a],label=str(k))
for pattern in ('1010','0111','0100','0000'):
 a=[r for r in rs if r['pattern']==pattern and float(r['marker'])==0 and float(r['kappa'])==.01 and int(r['duration'])==50];axs[2].plot([float(r['position'])/1e6 for r in a],[float(r['q2']) for r in a],label=pattern)
for a,title in zip(axs,('Structured root duration','Cross suppression kappa','Original geometry / terminal patterns')):a.set_title(title);a.legend();a.set_xlabel('Position (Mb)');a.set_ylabel('Theoretical q2');a.axvspan(3,7,color='grey',alpha=.12)
f.tight_layout();f.savefig(root/'sensitivity.png',dpi=180);f.savefig(root/'sensitivity.svg')
