"""One chromosome is the replication unit; windows intersect native ARG intervals."""
import sys,csv,json,gzip,collections
from pathlib import Path
import numpy as np
import tskit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_moran_validation.windows import summarize
from validation.chromosome_scalability.worker import topology
root=Path(__file__).parent/'results';root.mkdir(exist_ok=True)

def load(path,pattern=None):
 ts=tskit.load(path);rows=[];signatures=[]
 if path.is_relative_to(root):
  from msrcsim.arg.output.genealogy_export import export_genealogies
  export_genealogies(ts,path.parent/'local_trees.nwk',path.parent/'local_trees.csv',chrom='chr1')
 for t in ts.trees():
  pairs=np.array([t.time(t.mrca(i,j)) for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))]);q=topology(t)
  rows.append([*t.interval,*[float(q==i) for i in range(3)],float(q==1)-float(q==0),t.time(t.root),pairs.mean(),sum(t.branch_length(u) for u in t.nodes() if t.parent(u)!=-1 and t.num_samples(u)==2),pairs[1],pairs[4],pairs[0],pairs[2],np.mean([v for v,(i,j) in zip(pairs,((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))) if pattern and pattern[i]==pattern[j]]) if pattern and any(pattern[i]==pattern[j] for i in range(4) for j in range(i)) else np.nan,np.mean([v for v,(i,j) in zip(pairs,((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))) if pattern and pattern[i]!=pattern[j]]) if pattern and len(set(pattern))>1 else np.nan])
  signatures.append(tuple(pairs))
 a=np.array(rows);return a,signatures

def write(path,rows):
 if not rows:return
 opener=gzip.open if str(path).endswith('.gz') else open
 with opener(path,'wt',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

groups=collections.defaultdict(list);attempts=[]
for p in sorted((root/'runs').glob('*/summary.json')):
 s=json.loads(p.read_text());attempts.append(s)
 if s['status']=='complete' and s['job']['duration']==1000:
  key='conditional1010' if s['job'].get('condition')=='1010' else 'unconditional';a,sig=load(p.parent/'chromosome.trees',s['terminal_pattern']);groups[key].append((p.parent.name,a,sig,s))
# Replay duplicates are alternate measurements of the same chromosome, never extra replication.
for key in list(groups):
 seen=set();unique=[]
 for entry in sorted(groups[key],key=lambda z:z[3]['job'].get('seconds',180),reverse=True):
  seed=entry[3]['job']['history_seed']
  if seed not in seen:unique.append(entry);seen.add(seed)
 groups[key]=unique
old=Path(__file__).resolve().parents[1]/'chromosome_biological_validation/results/runs'
for model,key in [('A','neutral'),('C','fixed_half')]:
 for p in sorted(old.glob(f'{model}_L10000000_*/chromosome.trees')):
  a,sig=load(p);groups[key].append((p.parent.name,a,sig,{}))
keys=['neutral','fixed_half','unconditional','conditional1010'];widths=[250000,50000,25000,10000];tracks={};metrics=['q1','q2','q3','D_topology','root','pair_mean','internal','TAC','TBD','TAB','TAD','same_background','opposite_background'];rng=np.random.default_rng(20261123);summary=[];spatial=[];branch=[];tractrows=[];within_topology=[]
for key in keys:
 for width in widths:
  rows=[];vs=[]
  for name,a,sig,s in groups[key]:
   x,v=summarize(a[:,0],a[:,1],a[:,2:],1e7,width);vs.append(v)
   for j,xx in enumerate(x):rows.append(dict(chromosome=name,left=xx-width/2,right=xx+width/2,**dict(zip(metrics,v[:,j]))))
  tracks[key,width]=(x if vs else np.array([]),np.array(vs));write(root/f'{key}_{width}.csv.gz',rows)
 for name,a,sig,s in groups[key]:
  for kind,labels in [('genealogy',sig),('topology',list(a[:,2:5].argmax(1)))]:
   merged=[]
   for j,label in enumerate(labels):
    if j and label==labels[j-1]:merged[-1][1]=a[j,1]
    else:merged.append([a[j,0],a[j,1]])
   for j,(left,right) in enumerate(merged):tractrows.append(dict(model=key,chromosome=name,kind=kind,left=left,right=right,length=right-left,censored=j in (0,len(merged)-1)))
  changes=[a[j,0] for j in range(1,len(sig)) if sig[j]!=sig[j-1]];tops=a[:,2:5].argmax(1);topchanges=a[1:,0][tops[1:]!=tops[:-1]]
  summary.append(dict(model=key,chromosome=name,windows_primary=400,marginal_intervals=len(a),distinct_dated_genealogies=len(set(sig)),distinct_topologies=len(set(tops)),genealogy_changes=len(changes),topology_changes=len(topchanges)))
  for region in ['inside','outside']:
   w=np.maximum(0,np.minimum(a[:,1],7e6)-np.maximum(a[:,0],3e6));length=4e6
   if region=='outside':w=a[:,1]-a[:,0]-w;length=6e6
   for q in range(3):
    wt=w*(a[:,2:5].argmax(1)==q)
    if wt.sum()>0:within_topology.append(dict(model=key,chromosome=name,region=region,topology=q,span=wt.sum(),**dict(zip(metrics,np.average(a[:,2:],axis=0,weights=wt)))))
   branch.append(dict(model=key,chromosome=name,region=region,**dict(zip(metrics,np.average(a[:,2:],axis=0,weights=w))),genealogy_changes_per_Mb=sum((3e6<=v<7e6)==(region=='inside') for v in changes)/(length/1e6),topology_changes_per_Mb=sum((3e6<=v<7e6)==(region=='inside') for v in topchanges)/(length/1e6)))
  grid=np.arange(0,1e7,10000);ix=np.searchsorted(a[:,1],grid,side='right');t=tops[ix]
  for lag in [1,5,10,50,100,200]:
   for region in ['inside','outside','all']:
    inside=(grid>=3e6)&(grid<7e6);mask=np.ones(len(grid)-lag,bool) if region=='all' else (inside[:-lag]&inside[lag:] if region=='inside' else ~inside[:-lag]&~inside[lag:])
    d=dict(model=key,chromosome=name,region=region,distance=lag*10000,topology_agreement=float(np.mean((t[:-lag]==t[lag:])[mask])))
    left=np.eye(3)[t[:-lag][mask]];right=np.eye(3)[t[lag:][mask]];left-=left.mean(axis=0);right-=right.mean(axis=0);den=np.sqrt((left*left).sum()*(right*right).sum());d['topology_correlation']=float((left*right).sum()/den) if den else float('nan')
    for name2,col in [('pair_mean',7),('root',6),('internal',8)]:
     v=a[ix,col];d[name2+'_correlation']=float(np.corrcoef(v[:-lag][mask],v[lag:][mask])[0,1])
    spatial.append(d)
write(root/'within_topology_branches.csv',within_topology);write(root/'tracts.csv',tractrows);write(root/'chromosome_summary.csv',summary);write(root/'inside_outside.csv',branch);write(root/'spatial.csv',spatial)

bounds=[]
for key in ['unconditional','conditional1010']:
 x,v=tracks[key,25000];n=len(v)
 if n:
  for metric in ['q1','q2','q3','D_topology']:
   col=metrics.index(metric);total=v[:,col,:].sum(0);low=-1 if metric=='D_topology' else 0
   for j,xx in enumerate(x):bounds.append(dict(model=key,position=xx,metric=metric,planned=2,completed=n,cohort_lower=(total[j]+(2-n)*low)/2,cohort_upper=(total[j]+2-n)/2,interpretation='missing-chromosome bounds, not confidence intervals'))
write(root/'missing_chromosome_bounds.csv',bounds)

def band(ax,key,width,metric):
 x,v=tracks[key,width];n=len(v)
 if not n:ax.text(.5,.5,'No completed chromosomes',ha='center',transform=ax.transAxes);return
 k=metrics.index(metric);a=v[:,k,:];ax.plot(x/1e6,a.mean(0),label=f'{key}, completed n={n}'+(' of 2; descriptive' if key in ('unconditional','conditional1010') else ''))
 if n>1 and key not in ('unconditional','conditional1010'):
  boot=a[rng.integers(0,n,(500,n))].mean(1);lo,hi=np.quantile(boot,[.025,.975],axis=0);ax.fill_between(x/1e6,lo,hi,alpha=.15)
 ax.axvspan(3,7,color='grey',alpha=.12);ax.set_xlabel('Position (Mb)')

def save(f,name):
 f.tight_layout();f.savefig(root/(name+'.png'),dpi=180);f.savefig(root/(name+'.svg'));plt.close(f)
f,axs=plt.subplots(4,4,figsize=(15,11),sharex=True)
for j,key in enumerate(keys):
 axs[0,j].set_title(key)
 for i,m in enumerate(metrics[:4]):band(axs[i,j],key,25000,m);axs[i,j].set_ylabel(m);axs[i,j].set_ylim((-1,1) if i==3 else (0,1))
f.suptitle('25-kb tracks; Moran panels: completed examples, not validated ensemble expectations');save(f,'figure1_quartets')
f,axs=plt.subplots(2,1,figsize=(12,7))
for key in keys:
 band(axs[0],key,25000,'q2');band(axs[1],key,25000,'D_topology')
for a in axs:a.legend(fontsize=8)
axs[0].set_ylim(0,1);axs[0].set_ylabel('q2')
axs[1].set_ylim(-1,1);axs[1].set_ylabel('D_topology')
save(f,'figure2_conditioning')
f,axs=plt.subplots(5,4,figsize=(15,12),sharex=True)
for j,key in enumerate(keys):
 axs[0,j].set_title(key)
 for i,m in enumerate(['root','pair_mean','internal','TAC','TAB']):band(axs[i,j],key,25000,m);axs[i,j].set_ylabel(m+' (generations)')
save(f,'figure3_branches')
f,axs=plt.subplots(4,4,figsize=(15,10),sharex=True)
for j,key in enumerate(keys):
 axs[0,j].set_title(key)
 for i,w in enumerate(widths):band(axs[i,j],key,w,'q2');axs[i,j].set_ylabel(f'q2; {w//1000} kb');axs[i,j].set_ylim(0,1)
save(f,'figure4_resolution')
f,gridaxs=plt.subplots(2,4,figsize=(16,8));axs=gridaxs.flat
for key in keys:
 for i,m in enumerate(['topology_agreement','topology_correlation','pair_mean_correlation','root_correlation','internal_correlation']):
  rs=[r for r in spatial if r['model']==key and r['region']=='all'];ds=sorted({r['distance'] for r in rs})
  if ds:axs[i].plot(np.array(ds)/1e6,[np.mean([r[m] for r in rs if r['distance']==d]) for d in ds],label=key)
  axs[i].set_title(m);axs[i].set_xlabel('Distance (Mb)')
for key in keys:
 for j,kind in enumerate(['genealogy','topology']):
  rs=[r for r in tractrows if r['model']==key and r['kind']==kind and not r['censored']];names=sorted({r['chromosome'] for r in rs});values=[];weights=[]
  for name in names:
   ds=[r['length'] for r in rs if r['chromosome']==name];values.extend(ds);weights.extend([1/len(names)/len(ds)]*len(ds))
  if values:
   order=np.argsort(values);axs[5+j].plot(np.array(values)[order]/1e6,np.cumsum(np.array(weights)[order]),label=key)
  axs[5+j].set_xscale('log');axs[5+j].set_xlabel('Completed tract (Mb)');axs[5+j].set_ylabel(kind+' ECDF; censored ends excluded')
axs[7].set_visible(False);axs[0].legend(fontsize=7);save(f,'figure5_dependence')
pool=json.loads((root/'history_pool.json').read_text());f,axs=plt.subplots(1,3,figsize=(15,4))
seen=set()
complete_names={name for key in ['unconditional','conditional1010'] for name,a,sig,s in groups[key]}
for path in sorted((root/'runs').glob('*/history_trace.csv'),key=lambda p:(p.parent.name not in complete_names,p.parent.name)):
 label='conditional1010' if path.parent.name.startswith('B_') else 'unconditional'
 if label in seen:continue
 seen.add(label);rs=list(csv.DictReader(path.open()));origin=json.loads((path.parent/'summary.json').read_text())['origin'];ps=[r for r in rs if r['population']=='ROOT' and float(r['age'])<origin-1e-8];axs[0].plot([float(r['age']) for r in ps]+[origin],[float(r['p']) for r in ps]+[.5],label=path.parent.name);axs[0].axvline(origin,color='grey',linestyle=':')
counts=collections.Counter(r['pattern'] for r in pool['rows']);labels=[format(i,'04b') for i in range(16)];axs[1].bar(labels,[counts[k] for k in labels]);axs[1].tick_params(axis='x',rotation=90);axs[1].set_title(f"Terminal patterns; {pool['completed']} prior histories")
axs[2].bar(['lost','fixed','segregating'],[sum(r['lost_tips'] for r in pool['rows']),sum(r['fixed_tips'] for r in pool['rows']),4*pool['completed']-sum(r['lost_tips']+r['fixed_tips'] for r in pool['rows'])]);axs[2].set_title('Tip populations (correlated within chromosome)');axs[0].set_title('Representative root paths only');axs[0].set_xlabel('Backward age (generations)');axs[0].set_ylabel('p(t)');axs[0].legend(fontsize=7);save(f,'figure6_histories')
(root/'analysis_summary.json').write_text(json.dumps(dict(completed_unique={k:len(groups[k]) for k in keys},attempted=[dict(job=s['job'],status=s['status'],seconds=s['seconds'],rss=s['rss']) for s in attempts],pool_completed=pool['completed'],conditional_prior_acceptance=counts['1010']/pool['completed'],large_M_theory='unavailable; no approximation substituted'),indent=2))
