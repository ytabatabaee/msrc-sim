"""Compact summaries from primary grid; chromosome bootstrap uncertainty."""
import argparse,csv,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from msrcsim.species_tree import SpeciesTree
from msrcsim.smc.engine import _nodes,_desc


def read(path):
    with path.open() as f:return list(csv.DictReader(f))


def write(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def clades_gene(nwk):
    root=SpeciesTree(nwk,20,200,ultrametric_tolerance=1e-7).root
    taxa=sorted(_desc(root));pairs=[]
    for i,a in enumerate(taxa):
        for b in taxa[:i]:pairs.append(min(n.age for n in _nodes(root) if {a,b}<=_desc(n)))
    return matrix_clades(pairs,len(taxa)),tuple(pairs)


def matrix_clades(pairs,nt):
    a=np.zeros((nt,nt));j=0
    for i in range(nt):
        for k in range(i):a[i,k]=a[k,i]=pairs[j];j+=1
    clades=set()
    for i in range(nt):
        for k in range(i):
            group=tuple(x for x in range(nt) if a[i,x]<=a[i,k]+1e-9)
            if 1<len(group)<nt:clades.add(group)
    return clades


def main(src,out):
    out.mkdir(parents=True,exist_ok=True);tests=[];spatial=[];cases=[];wait=[];tract=[];paths=[]
    rng=np.random.default_rng(9971)
    for directory in sorted(src.iterdir()):
        if not (directory/'summary.json').exists():continue
        summary=json.loads((directory/'summary.json').read_text());case=summary['case'];cases.append(summary)
        for row in read(directory/'marginal_tests.csv'):
            row={k:(v if k in ('case','method','metric') else float(v)) for k,v in row.items()};tests.append(row)
        data=read(directory/'spatial.csv');nt=int(case[0]);replicates={}
        cached={}
        for r in data:
            rep=int(r['replicate']);x=int(r['position'])
            sm=(r['smc_newick'])
            if sm not in cached:cached[sm]=clades_gene(sm)
            cls,sig=cached[sm]
            asig=tuple(float(x) for x in r['arg_signature'].split(';'))
            replicates.setdefault(rep,{})[x]=(r,cls,sig,matrix_clades(asig,nt),asig)
        for x,y in ((5,25),(25,50),(25,75),(5,95)):
            smpairs=[];arpairs=[];smroot=[];arroot=[];smrf=[];arrf=[];smeq=[];areq=[];smtop=[];artop=[]
            for rep,rows in replicates.items():
                a,b=rows[x],rows[y];ra,rb=a[0],b[0]
                smpairs.append((float(ra['smc_pair']),float(rb['smc_pair'])));arpairs.append((float(ra['arg_pair']),float(rb['arg_pair'])))
                smroot.append((float(ra['smc_root']),float(rb['smc_root'])));arroot.append((float(ra['arg_root']),float(rb['arg_root'])))
                smrf.append(len(a[1]^b[1]));arrf.append(len(a[3]^b[3]));smeq.append(a[2]==b[2]);areq.append(a[4]==b[4])
                smtop.append(ra['smc_top']==rb['smc_top']);artop.append(ra['arg_top']==rb['arg_top'])
            for metric,sm,ar in [('tree_agreement',smeq,areq),('topology_agreement',smtop,artop),('rooted_RF',smrf,arrf),('pair_time_correlation',smpairs,arpairs),('root_time_correlation',smroot,arroot)]:
                a=np.asarray(sm,float);b=np.asarray(ar,float);iscorr=a.ndim==2
                f=lambda z:np.corrcoef(z.T)[0,1] if iscorr else np.mean(z)
                v=float(f(a));ref=float(f(b));boot=[]
                for _ in range(250):
                    ia=rng.integers(len(a),size=len(a));ib=rng.integers(len(b),size=len(b))
                    boot.append(f(a[ia])-f(b[ib]))
                lo,hi=np.quantile(boot,[.025,.975])
                spatial.append(dict(case=case,start=x,end=y,distance=y-x,metric=metric,smc=v,reservoir=ref,difference=v-ref,ci_low=lo,ci_high=hi,n_chromosomes=len(a)))
        # Current local-path dwell means and switching/exposure diagnostics.
        pp=read(directory/'path_summaries.csv')
        for key in ('switches','dwell0','dwell1'):
            a=np.array([float(r[key]) for r in pp]);se=a.std(ddof=1)/np.sqrt(len(a))
            paths.append(dict(case=case,metric=key,mean=a.mean(),se=se,ci_low=a.mean()-1.96*se,ci_high=a.mean()+1.96*se))
        evfile=directory/'events.csv';events=read(evfile) if evfile.exists() else []
        byrep={i:[] for i in replicates}
        for e in events:byrep[int(e['replicate'])].append(e)
        for rep,evs in byrep.items():
            for kind,field in [('recombination_proposal',None),('genealogy_change','genealogy_changed'),('topology_change','topology_changed')]:
                positions=[float(e['position']) for e in evs if field is None or e[field]=='True'];last=0.
                for pos in positions:
                    wait.append(dict(case=case,replicate=rep,kind=kind,length=pos-last,right_censored=False));last=pos
                wait.append(dict(case=case,replicate=rep,kind=kind,length=100-last,right_censored=True))
                if kind=='topology_change':
                    for a,b in zip([0]+positions,positions+[100]):tract.append(dict(case=case,replicate=rep,length=b-a,boundary_censored=a==0 or b==100))
    order=sorted(range(len(tests)),key=lambda i:tests[i]['pvalue']);prev=0
    for rank,i in enumerate(order):
        prev=max(prev,min(1.,(len(tests)-rank)*tests[i]['pvalue']));tests[i]['holm_p']=prev
        tests[i]['practical_equivalence']=tests[i]['ci_low']>=-tests[i]['tolerance'] and tests[i]['ci_high']<=tests[i]['tolerance']
    write(out/'marginal_tests.csv',tests);write(out/'spatial_summary.csv',spatial);write(out/'case_summary.csv',cases);write(out/'path_summary.csv',paths)
    # Keep compact waiting summaries, preserving censoring rather than treating
    # chromosome-end intervals as completed waits.
    waits=[]
    for case in [r['case'] for r in cases]:
        for kind in ('recombination_proposal','genealogy_change','topology_change'):
            rr=[r for r in wait if r['case']==case and r['kind']==kind]
            complete=[r['length'] for r in rr if not r['right_censored']]
            waits.append(dict(case=case,kind=kind,n_completed=len(complete),n_right_censored=sum(r['right_censored'] for r in rr),mean_completed=float(np.mean(complete)) if complete else None,median_completed=float(np.median(complete)) if complete else None))
    write(out/'waiting_summary.csv',waits)
    tractsummary=[]
    for case in [r['case'] for r in cases]:
        rr=[r for r in tract if r['case']==case];interior=[r['length'] for r in rr if not r['boundary_censored']]
        tractsummary.append(dict(case=case,n_tracts=len(rr),n_interior=len(interior),mean_interior_length=float(np.mean(interior)) if interior else None))
    write(out/'tract_summary.csv',tractsummary)
    result=dict(cases=len(cases),seed=20261121,primary_replicates_per_case=1000,
                scientific_scope='frozen reference marginal at 50, uniform map and whole-region cross suppression',
                targets={'quartet_and_D_absolute_difference':.05,'time_mean_relative_difference':.1,'spatial_primary_tolerance':'quantified rather than calibrated: no acceptance asserted'},
                comparisons={m:dict(n=sum(r['method']==m for r in tests),holm_significant=sum(r['method']==m and r['holm_p']<.05 for r in tests),practical_CI_pass=sum(r['method']==m and r['practical_equivalence'] for r in tests)) for m in ('SMC','reservoir','SMC_analytic')},
                max_abs_log_acceptance=max([r['max_abs_log_acceptance'] for r in cases],default=0))
    (out/'validation_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,3,figsize=(14,4))
    rr=[r for r in tests if r['method']=='SMC' and r['metric']=='q2']
    ax[0].errorbar(range(len(rr)),[r['difference'] for r in rr],yerr=[1.96*r['se'] for r in rr],fmt='o');ax[0].axhline(0,color='k',lw=.5)
    ax[0].axhline(.05,color='r',ls='--');ax[0].axhline(-.05,color='r',ls='--');ax[0].set(title='Reference-position q2 differences',xlabel='Four-taxon case',ylabel='SMC − legacy (95% CI)')
    rr=[r for r in tests if r['method']=='SMC' and r['metric']=='root_age'];ax[1].plot([r['reference'] for r in rr],[r['observed'] for r in rr],'o');ax[1].plot([40,160],[40,160],'k--');ax[1].set(title='Reference root-age means',xlabel='Legacy',ylabel='SMC')
    for m,label in [('smc','SMC'),('reservoir','Reservoir ARG')]:
        rr=[r for r in spatial if r['case']=='4tax_p0.5_k1.0' and r['metric']=='pair_time_correlation']
        ax[2].plot([r['distance'] for r in rr],[r[m] for r in rr],'o-',label=label)
    ax[2].legend();ax[2].set(title='Pair-time dependence, p=.5 κ=1',xlabel='Separation (bp)',ylabel='Correlation')
    fig.tight_layout();fig.savefig(out/'validation.png',dpi=160);plt.close(fig)
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--input',type=Path,required=True);a.add_argument('--output',type=Path,required=True)
    args=a.parse_args();main(args.input,args.output)
