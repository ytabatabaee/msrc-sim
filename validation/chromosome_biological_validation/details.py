"""Span-weighted branch distributions and independent sensitivity replicates."""
import csv,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_biological_validation.theory_grid import case
from validation.moran_spatial.theory import load_theory,make_tree,PrescribedHistory,predict,rho_effective
root=Path(__file__).parent/'results';mods,_=load_theory('/Users/ytabatabaee/Desktop/msrc-paper');rows=[];groups={}
for p in sorted((root/'runs').glob('*/summary.json')):
 s=json.loads(p.read_text());j=s['job']
 if s['status']!='complete' or j['length']!=10000000:continue
 with (p.parent/'trees.csv').open() as f:r=list(csv.DictReader(f))
 d={k:np.array([float(v[k]) for v in r]) for k in r[0]};groups.setdefault(j['model'],[]).append(d)
 for region in ('inside','outside'):
  w=np.maximum(0,np.minimum(d['right'],7e6)-np.maximum(d['left'],3e6))
  if region=='outside':w=d['right']-d['left']-w
  rows.append(dict(model=j['model'],seed=j['seed'],region=region,**{k:float(d[k]@w/w.sum()) for k in ('root','internal','mean_pair','TAC','TBD','TAB','TAD')}))
f,axs=plt.subplots(1,3,figsize=(13,4))
for model,ds in groups.items():
 for region,style in (('inside','-'),('outside','--')):
  for a,key in zip(axs,('root','internal','mean_pair')):
   v=[];weights=[]
   for d in ds:
    w=np.maximum(0,np.minimum(d['right'],7e6)-np.maximum(d['left'],3e6))
    if region=='outside':w=d['right']-d['left']-w
    v.extend(d[key]);weights.extend(w/w.sum()/len(ds))
   order=np.argsort(v);a.plot(np.array(v)[order],np.cumsum(np.array(weights)[order]),style,label=model+' '+region);a.set_xlabel(key+' (generations)');a.set_ylabel('Span-weighted ECDF')
axs[0].legend(fontsize=7);f.tight_layout();f.savefig(root/'branch_distributions.png',dpi=180)
with (root/'branch_by_chromosome.csv').open('w',newline='') as out:
 w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
f,axs=plt.subplots(2,4,figsize=(15,8));points=np.array([0,2.9e6,3.1e6,5e6,6.9e6,7.1e6,9999999.]);allrows=[]
for i,a in enumerate(axs.flat):
 paths=list((root/'sensitivity').glob(f'case{i}_*/summary.json'))
 if not paths:a.set_visible(False);continue
 ss=[json.loads(p.read_text()) for p in paths];j=ss[0]['job'];c=case(j['duration'],j['kappa'],j.get('marker_fraction',.5)*1e7,j.get('rate',1e-8));t=make_tree(c);h=PrescribedHistory(t,{bid:[(.5,b.older_age-b.younger_age)] for bid,b in t.branches.items()});pred=np.array([predict(mods,c,rho_effective(c,x),h).conditional_quartet_probabilities[int(j['pattern'],2),1] for x in points]);observations=[]
 for p,s in zip(paths,ss):
  if s['status']!='complete':continue
  with (p.parent/'trees.csv').open() as z:rs=list(csv.DictReader(z))
  right=np.array([float(r['right']) for r in rs]);top=np.array([int(r['topology']) for r in rs]);observations.append((top[np.searchsorted(right,points,side='right')]==1).astype(int))
 if observations:
  obs=np.mean(observations,0);a.plot(points/1e6,obs,'o--',label=f'ARG n={len(observations)}')
  for x,y,q in zip(points,obs,pred):allrows.append(dict(case=i,position=x,n=len(observations),observed_q2=y,theory_q2=q,job=json.dumps(j)))
 a.plot(points/1e6,pred,'k.-',label='theory');a.axvspan(3,7,color='grey',alpha=.12);a.set_ylim(0,1);a.set_title(f"T={j['duration']}, k={j['kappa']}, z={c['marker']/1e6:g}\nZ={j['pattern']}, r={j.get('rate',1e-8):g}",fontsize=9);a.legend(fontsize=7);a.set_xlabel('Position (Mb)');a.set_ylabel('q2')
f.suptitle('Sensitivity: three chromosomes per setting; descriptive point frequencies')
f.tight_layout();f.savefig(root/'sensitivity_observed.png',dpi=180)
if allrows:
 with (root/'sensitivity_estimates.csv').open('w',newline='') as out:
  w=csv.DictWriter(out,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
