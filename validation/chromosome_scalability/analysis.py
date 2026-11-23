"""Chromosome-level diagnostics; genomic windows are never independent replicates."""
import argparse,collections,csv,hashlib,json,math,shutil,sys
from pathlib import Path
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from msrcsim.arg.output.genealogy_export import export_genealogies

COLORS={'neutral_uniform':'#2374ab','neutral_low':'#4d923c','MSRC_ARG':'#bd631d'}
LABELS={'neutral_uniform':'Neutral uniform','neutral_low':'Neutral low interval','MSRC_ARG':'MSRC-ARG'}


def write_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def load_runs(roots):
    runs=[]
    for root in roots:
        for path in sorted(root.glob('*/summary.json')):
            r=json.loads(path.read_text());r['directory']=path.parent;r['id']=path.parent.name
            if r['job']['kind']!='history':runs.append(r)
    return runs


def key(r):
    j=r['job'];return j['kind'],j['taxa'],j['length'],j['seed'],j['rate'],j['kappa']


def prefer(runs):
    chosen={}
    for r in runs:
        k=key(r)
        # Cohorts are the same fixed seeds rerun with increased resources,
        # never a new seed selected because it happened to finish.
        if k not in chosen or r['job']['seconds']>chosen[k]['job']['seconds']:chosen[k]=r
    return list(chosen.values())


def enrich(r):
    import tskit
    rows=list(csv.DictReader((r['directory']/'marginals.csv').open()))
    for row in rows:
        for k in row:row[k]=int(row[k]) if k=='topology' else float(row[k])
    ts=tskit.load(r['directory']/'chromosome.trees')
    assert len(rows)==ts.num_trees
    for row,gt in zip(rows,ts.trees()):
        pairs=[gt.time(gt.mrca(i,j)) for i in range(ts.num_samples) for j in range(i)]
        row['genealogy_id']=hashlib.sha256(np.asarray(pairs,dtype='<f8').tobytes()).hexdigest()
        row['quartet_internal_generations']=sum(gt.branch_length(u) for u in gt.nodes() if gt.parent(u)!=-1 and gt.num_samples(u)==2) if ts.num_samples==4 else float('nan')
    return rows,ts


def windows(rows,L,width):
    starts=np.arange(0,L,width);ends=np.minimum(starts+width,L);vals=[];i=0
    for a,b in zip(starts,ends):
        sums=np.zeros(6);cursor=i
        while cursor<len(rows) and rows[cursor]['left']<b:
            r=rows[cursor];overlap=max(0.,min(b,r['right'])-max(a,r['left']))
            if overlap:
                sums[:3]+=overlap*np.eye(3)[int(r['topology'])]
                sums[3:]+=overlap*np.array([r['root_age'],r['pair_mean'],r['quartet_internal_generations']])
            if r['right']<=b:cursor+=1
            else:break
        i=cursor;vals.append(sums/(b-a))
    return starts,ends,np.asarray(vals)


def binomial_interval(k,n,alpha=.05):
    if n==0:return float('nan'),float('nan')
    return (0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1)),
            1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k)))


def point_estimates(draws,neutral):
    rows=[]
    for model,x in sorted({(r['model'],r['position']) for r in draws}):
        values=[r['quartet'] for r in draws if r['model']==model and r['position']==x];n=len(values);q=np.bincount(values,minlength=3)/n
        for k,name in enumerate(('q1','q2','q3','D_topology')):
            observed=q[k] if k<3 else q[1]-q[0]
            if k<3:lo,hi=binomial_interval(sum(v==k for v in values),n)
            else:
                l0,u0=binomial_interval(sum(v==0 for v in values),n,.025);l1,u1=binomial_interval(sum(v==1 for v in values),n,.025);lo,hi=max(-1,l1-u0),min(1,u1-l0)
            target=None if model=='MSRC_ARG' or neutral is None else float(neutral[k] if k<3 else neutral[1]-neutral[0])
            status='target_unavailable' if target is None else 'PASS' if lo-target>=-.05 and hi-target<=.05 else 'FAIL_mismatch' if lo-target>.05 or hi-target<-.05 else 'UNRESOLVED'
            rows.append(dict(model=model,position=x,metric=name,n_chromosomes=n,observed=observed,ci_low=lo,ci_high=hi,theory=target,status=status,tolerance=.05,difference=None if target is None else observed-target,absolute_error=None if target is None else abs(observed-target),difference_ci_low=None if target is None else lo-target,difference_ci_high=None if target is None else hi-target,method='Clopper-Pearson q; Bonferroni difference of 97.5% binomial intervals for D'))
    return rows


def tracts(rows,field):
    out=[]
    for r in rows:
        if out and out[-1]['value']==r[field]:out[-1]['right']=r['right']
        else:out.append(dict(left=r['left'],right=r['right'],value=r[field]))
    L=rows[-1]['right']
    return [dict(**r,length=r['right']-r['left'],censored=r['left']==0 or r['right']==L) for r in out]


def span_ecdf(rows,field,inside,a,b):
    pairs=[]
    for r in rows:
        w=max(0.,min(b,r['right'])-max(a,r['left']))
        if not inside:w=(r['right']-r['left'])-w
        if w:pairs.append((r[field],w))
    pairs.sort();x=np.array([v for v,w in pairs]);weights=np.array([w for v,w in pairs]);return x,np.cumsum(weights)/sum(weights)


def savefig(fig,out,name):
    for ext in ('png','svg'):
        p=out/(name+'.'+ext);fig.savefig(p,dpi=200,bbox_inches='tight')
        if ext=='svg':p.write_text('\n'.join(line.rstrip() for line in p.read_text().splitlines())+'\n')
    plt.close(fig)


def analyze(roots,out,theory_check=None):
    out.mkdir(parents=True,exist_ok=True);runs=load_runs(roots);selected=prefer(runs)
    theoretical=json.loads(theory_check.read_text()) if theory_check else {}
    neutral=np.asarray(theoretical['neutral_MSC_quartets']) if 'neutral_MSC_quartets' in theoretical else None
    table=[]
    for r in runs:
        j=r['job'];events=r.get('event_counts',{});state=r.get('partial_state',{}) or {};event=lambda name:events.get(name,0) if r['status']=='complete' else None
        table.append(dict(id=r['id'],model=j['kind'],taxa=j['taxa'],length=j['length'],seed=j['seed'],workflow='inline' if j.get('history_mode')=='inline' else 'cached_history' if j['kind']=='MSRC_ARG' else 'neutral',ne=j['ne'],rate=j['rate'],kappa=j['kappa'],status=r['status'],budget_seconds=j['seconds'],process_wall_seconds=r.get('process_wall_seconds'),pipeline_wall_seconds=r.get('wall_seconds'),simulation_seconds=r.get('simulation_seconds'),peak_rss_bytes=r.get('peak_rss_bytes'),recombination_events=event('recombination'),background_transfer_events=event('background_transfer'),coalescence_events=event('coalescence'),marginal_intervals=r.get('marginal_trees'),distinct_genealogies=r.get('distinct_dated_genealogies'),nodes=r.get('ts_nodes'),edges=r.get('ts_edges'),tree_file_bytes=r.get('tree_file_bytes'),event_log_bytes=r.get('event_log_bytes'),export_seconds=r.get('export_seconds'),peak_carriers=r.get('peak_carriers'),peak_stored_segments=r.get('peak_segments',r.get('peak_stored_segments')),runtime_per_tree=r.get('runtime_per_tree'),failed_stage=r.get('failed_stage'),partial_generation=state.get('time_generations')))
    write_csv(out/'benchmark_metrics.csv',table)
    plt.rcParams.update({'font.size':9,'svg.fonttype':'none'})
    L=10000000;width=250000;cohorts={};spatial=[];wtable=[];ptr=[];profile=[];point=[]
    for model in COLORS:
        rs=[r for r in selected if r['job']['kind']==model and r['job']['taxa']==4 and r['job']['length']==L and r['job']['rate']==1e-8 and r['job']['kappa']==.01]
        byseed={r['job']['seed']:r for r in rs}
        if set(byseed)!={101,202,303} or any(r['status']!='complete' for r in rs):continue
        data=[];series=[];cors=[];pers=[];changes=[]
        for seed,r in sorted(byseed.items()):
            rows,ts=enrich(r);data.append(rows);a,b,v=windows(rows,L,width);series.append(v)
            sample_x=np.arange(5000,L,10000);index=np.searchsorted([row['right'] for row in rows],sample_x,side='right');times=np.array([rows[i]['pair_mean'] for i in index]);tops=np.array([rows[i]['topology'] for i in index]);cs=[];ps=[]
            for lag in (1,5,10,25,50,100,200,300):
                corr=np.corrcoef(times[:-lag],times[lag:])[0,1] if times[:-lag].std()>0 and times[lag:].std()>0 else np.nan
                agreement=float(np.mean(tops[:-lag]==tops[lag:]));cs.append(corr);ps.append(agreement)
                spatial.append(dict(model=model,seed=seed,distance_bp=lag*10000,pair_time_correlation=corr,quartet_agreement=agreement))
            cors.append(cs);pers.append(ps)
            counts=[]
            for field in ('genealogy_id','topology'):
                tr=tracts(rows,field)
                for tract in tr:ptr.append(dict(model=model,seed=seed,kind=field,**tract))
                positions=[tr[i]['left'] for i in range(1,len(tr))];ninside=sum(3000000<=x<7000000 for x in positions);counts.extend((ninside/4,(len(positions)-ninside)/6))
            changes.append(counts)
            for x in (500000,2900000,3100000,5000000,6900000,7100000,9500000):
                gt=ts.at(x);q=next(row['topology'] for row in rows if row['left']<=x<row['right']);point.append(dict(model=model,seed=seed,position=x,quartet=q))
            if seed==101:
                artifact=out/'examples'/model;artifact.mkdir(parents=True,exist_ok=True);tables=ts.dump_tables();tables.time_units='generations';annotated=tables.tree_sequence();annotated.dump(artifact/'chromosome.trees');export_genealogies(annotated,artifact/'local_trees.nwk',artifact/'local_trees.csv',chrom='chr1')
                write_csv(artifact/'marginals.csv',rows);(artifact/'input.json').write_text(json.dumps(r['job'],indent=2)+'\n')
            for i,(aa,bb) in enumerate(zip(a,b)):wtable.append(dict(model=model,seed=seed,left=aa,right=bb,q1=v[i,0],q2=v[i,1],q3=v[i,2],D_topology=v[i,1]-v[i,0],root_generations=v[i,3],pair_generations=v[i,4],quartet_internal_generations=v[i,5]))
        cohorts[model]=dict(data=data,series=np.asarray(series),cors=np.asarray(cors),pers=np.asarray(pers),changes=np.asarray(changes),x=(a+b)/2/1e6)
    write_csv(out/'change_rates_by_chromosome.csv',[dict(model=m,seed=seed,genealogy_inside_per_Mb=v[0],genealogy_outside_per_Mb=v[1],quartet_inside_per_Mb=v[2],quartet_outside_per_Mb=v[3]) for m,c in cohorts.items() for seed,v in zip((101,202,303),c['changes'])]);write_csv(out/'window_by_chromosome.csv',wtable);write_csv(out/'spatial_by_chromosome.csv',spatial);write_csv(out/'tracts.csv',ptr);write_csv(out/'selected_position_draws.csv',point);write_csv(out/'selected_position_estimates.csv',point_estimates(point,neutral))
    # 1. Window-average quartets; bands show observed chromosome ranges, not CIs.
    fig,axes=plt.subplots(4,3,figsize=(13,10),sharex=True,constrained_layout=True)
    for col,model in enumerate(COLORS):
        axes[0,col].set_title(LABELS[model]+'; n=3 chromosomes' if model in cohorts else LABELS[model]+'; incomplete cohort')
        for row in range(4):
            ax=axes[row,col];ax.axvspan(3,7,color='.8',alpha=.3)
            if model in cohorts:
                c=cohorts[model];v=c['series'][:,:,row] if row<3 else c['series'][:,:,1]-c['series'][:,:,0];mean=v.mean(axis=0);lo,hi=v.min(axis=0),v.max(axis=0);ax.plot(c['x'],mean,color=COLORS[model]);ax.fill_between(c['x'],lo,hi,color=COLORS[model],alpha=.18)
            if model!='MSRC_ARG' and neutral is not None:ax.axhline(neutral[row] if row<3 else neutral[1]-neutral[0],color='k',ls=':',lw=1)
            ax.set_ylabel(f'q{row+1}' if row<3 else 'D_topology');ax.grid(alpha=.15)
        axes[-1,col].set_xlabel('Position (Mb)')
    savefig(fig,out,'quartet_windows')
    # 2. Branch summaries and span-weighted inside/outside distributions.
    fig,axes=plt.subplots(3,2,figsize=(12,10),constrained_layout=True)
    for k,field in enumerate(('root_age','pair_mean','quartet_internal_generations')):
        for model,c in cohorts.items():
            v=c['series'][:,:,k+3];mean=v.mean(axis=0);lo,hi=v.min(axis=0),v.max(axis=0);axes[k,0].plot(c['x'],mean,color=COLORS[model],label=LABELS[model]);axes[k,0].fill_between(c['x'],lo,hi,color=COLORS[model],alpha=.12)
            for inside in (True,False):
                allrows=[r for rows in c['data'] for r in rows];x,y=span_ecdf(allrows,field,inside,3000000,7000000);axes[k,1].plot(x,y,color=COLORS[model],ls='-' if inside else '--',label=LABELS[model]+(' inside' if inside else ' outside'))
        axes[k,0].axvspan(3,7,color='.8',alpha=.3);axes[k,0].set_ylabel(field+' (generations)');axes[k,0].set_xlabel('Position (Mb)');axes[k,1].set_xlabel(field+' (generations)');axes[k,1].set_ylabel('Span-weighted ECDF')
    axes[0,0].legend(fontsize=7);axes[0,1].legend(fontsize=7);savefig(fig,out,'branch_lengths')
    # 3. Descriptive dependence statistics: no tract/window independence assumed.
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True);dist=np.array([1,5,10,25,50,100,200,300])/100
    for model,c in cohorts.items():
        lengths=sorted(r['length']/1e6 for r in ptr if r['model']==model and r['kind']=='genealogy_id' and not r['censored'])
        if lengths:axes[0,0].plot(lengths,np.arange(1,len(lengths)+1)/len(lengths),color=COLORS[model],label=LABELS[model])
        axes[0,1].plot(dist,c['pers'].mean(axis=0),'o-',color=COLORS[model],label=LABELS[model]);axes[0,1].axhline(np.mean(np.sum(c['series'][:,:,:3].mean(axis=1)**2,axis=1)),color=COLORS[model],ls=':',lw=1);axes[1,0].plot(dist,np.nanmean(c['cors'],axis=0),'o-',color=COLORS[model],label=LABELS[model])
    axes[0,0].set_xscale('log');axes[0,0].set_xlabel('Completed genealogy tract (Mb)');axes[0,0].set_ylabel('Empirical CDF; censored tracts excluded');axes[0,0].legend(fontsize=8)
    axes[0,1].set_xlabel('Distance (Mb)');axes[0,1].set_ylabel('Quartet topology agreement');axes[0,1].set_title('Dotted: same-chromosome random-position baseline',fontsize=8);axes[1,0].set_xlabel('Distance (Mb)');axes[1,0].set_ylabel('Pair-time correlation')
    offset=0;labels=[]
    for model,c in cohorts.items():
        for k,kind in enumerate(('Genealogy','Quartet topology')):
            v=c['changes'][:,2*k:2*k+2];axes[1,1].bar([offset-.15,offset+.15],v.mean(axis=0),width=.28,color=COLORS[model],alpha=.8);axes[1,1].errorbar([offset-.15,offset+.15],v.mean(axis=0),yerr=[v.mean(axis=0)-v.min(axis=0),v.max(axis=0)-v.mean(axis=0)],fmt='none',ecolor='k',capsize=2);labels.append(LABELS[model]+'\n'+kind);offset+=1
    axes[1,1].set_xticks(range(offset),labels,rotation=35,ha='right',fontsize=7);axes[1,1].set_ylabel('Changes / Mb; inside then outside');savefig(fig,out,'spatial_dependence')
    # 4. Measured complete values and censored resource-limited attempts.
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for model in COLORS:
        for r in selected:
            j=r['job']
            if j['kind']!=model or j['rate']!=1e-8 or j['kappa']!=.01:continue
            complete=r['status']=='complete';marker='o' if complete else '^'
            if j['taxa']==4:
                x=j['length']/1e6;y=r.get('process_wall_seconds',r.get('wall_seconds'))
                if y is not None:axes[0,0].scatter(x,y,color=COLORS[model],marker=marker,s=25,facecolors=COLORS[model] if complete else 'none')
                if r.get('peak_rss_bytes'):axes[0,1].scatter(x,r['peak_rss_bytes']/1024**3,color=COLORS[model],marker=marker,s=25)
                if complete:axes[1,0].scatter(x,r['marginal_trees'],color=COLORS[model],s=25)
            if j['length']==L:
                y=r.get('process_wall_seconds',r.get('wall_seconds'))
                if y is not None:axes[1,1].scatter(j['taxa'],y,color=COLORS[model],marker=marker,s=25,facecolors=COLORS[model] if complete else 'none')
        axes[0,0].plot([],[],color=COLORS[model],marker='o',ls='none',label=LABELS[model])
    for r in runs:
        j=r['job']
        if j['kind']=='MSRC_ARG' and j['taxa']==4 and j.get('history_mode')=='inline' and r['status']=='complete':
            axes[0,0].scatter(j['length']/1e6,r['process_wall_seconds'],marker='*',s=65,color=COLORS['MSRC_ARG'])
            axes[0,1].scatter(j['length']/1e6,r['peak_rss_bytes']/1024**3,marker='*',s=65,color=COLORS['MSRC_ARG'])
    axes[0,0].set_title('Circles: controlled cached grid; stars: fresh MSRC')
    axes[0,1].set_title('History caching changes RSS, not the model')
    axes[0,0].legend(fontsize=8);axes[0,0].set_ylabel('Process wall time (s); △ limited');axes[0,1].set_ylabel('Peak RSS (GiB)');axes[1,0].set_ylabel('Marginal genealogy intervals');axes[1,1].set_ylabel('10-Mb process wall time (s); △ limited')
    for ax in (axes[0,0],axes[0,1],axes[1,0]):ax.set_xlabel('Chromosome length (Mb)');ax.set_xscale('log')
    for ax in (axes[0,0],axes[1,0],axes[1,1]):ax.set_yscale('log')
    axes[1,1].set_xlabel('Sampled taxa');savefig(fig,out,'scalability')
    summary=dict(total_attempts=len(runs),statuses=dict(collections.Counter(r['status'] for r in runs)),complete_10Mb_4taxon_cohorts=list(cohorts),cohort_seeds=[101,202,303],replication='one independent chromosome per seed; paired settings reuse histories/streams; no genomic windows pooled',interval='window bands are observed min-max across three chromosomes, not CIs; selected positions use exact binomial/Bonferroni intervals',exact_theory_status=theoretical.get('status','not_requested'),branch_units='generations',statistical_limit='three chromosomes are descriptive, not sufficient for tight marginal equivalence tests',phase8='deferred while 30-taxon MSRC scalability and matched-theory validation remain unresolved')
    (out/'analysis_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--benchmarks',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--theory-check',type=Path);a=p.parse_args();analyze(a.benchmarks,a.output,a.theory_check)
