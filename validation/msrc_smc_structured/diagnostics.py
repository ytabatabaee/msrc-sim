"""Off-reference marginal error, topology baseline and background-rate audit."""
import argparse,csv,json,re,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from scipy.stats import norm


def read(p):
    with p.open() as f:return list(csv.DictReader(f))


def write(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def run(src,bg,out):
    out.mkdir(parents=True,exist_ok=True);off=[];excess=[];rates=[]
    for directory in sorted(src.iterdir()):
        if not (directory/'spatial.csv').exists():continue
        rows=read(directory/'spatial.csv');case=directory.name;nt=int(case[0])
        for x in (5,25,50,75,95):
            rr=[r for r in rows if float(r['position'])==x]
            for metric,sa,sb in [('root_age','smc_root','arg_root'),('pair_time','smc_pair','arg_pair')]+([('q2','smc_top','arg_top')] if nt==4 else []):
                a=np.array([float(r[sa])==1 if metric=='q2' else float(r[sa]) for r in rr],float)
                b=np.array([float(r[sb])==1 if metric=='q2' else float(r[sb]) for r in rr],float)
                diff=a.mean()-b.mean();se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))
                off.append(dict(case=case,position=x,metric=metric,smc=a.mean(),reservoir=b.mean(),difference=diff,se=se,ci_low=diff-1.96*se,ci_high=diff+1.96*se))
        if nt==4:
            for x,y in ((5,25),(25,50),(25,75),(5,95)):
                a=[r for r in rows if float(r['position'])==x];b=[r for r in rows if float(r['position'])==y]
                for method,key in [('SMC','smc_top'),('reservoir','arg_top')]:
                    tx=np.array([int(float(r[key])) for r in a]);ty=np.array([int(float(r[key])) for r in b])
                    baseline=sum(np.mean(tx==q)*np.mean(ty==q) for q in range(3));observed=np.mean(tx==ty)
                    excess.append(dict(case=case,start=x,end=y,distance=y-x,method=method,observed_agreement=observed,independent_marginal_baseline=baseline,excess_agreement=observed-baseline))
    bdata=read(bg/'background_draws.csv')
    for case in sorted({r['case'] for r in bdata}):
        p=float(case.split('_p')[1].split('_k')[0]);k=float(case.split('_k')[1]);rho=.0002*50*k
        for method in ('SMC','legacy','reservoir'):
            rr=[r for r in bdata if r['case']==case and r['method']==method]
            n=np.array([float(r['switches']) for r in rr]);target=rho*np.array([p*float(r['dwell0'])+(1-p)*float(r['dwell1']) for r in rr]);delta=n-target
            se=delta.std(ddof=1)/np.sqrt(len(delta))
            rates.append(dict(case=case,method=method,mean_switches=n.mean(),mean_integrated_target=target.mean(),difference=delta.mean(),se=se,ci_low=delta.mean()-1.96*se,ci_high=delta.mean()+1.96*se,pvalue=float(2*norm.sf(abs(delta.mean()/se))) if se else 1.))
    order=sorted(range(len(rates)),key=lambda i:rates[i]['pvalue']);prev=0.
    for rank,i in enumerate(order):
        prev=max(prev,min(1.,(len(rates)-rank)*rates[i]['pvalue']));rates[i]['holm_p']=prev
    write(out/'off_reference_marginals.csv',off);write(out/'topology_excess_agreement.csv',excess);write(out/'background_rate_tests.csv',rates)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(14,4))
    for metric,ax in [('root_age',axes[0]),('q2',axes[1])]:
        rr=[r for r in off if r['case']=='4tax_p0.5_k1.0' and r['metric']==metric]
        ax.errorbar([r['position'] for r in rr],[r['difference'] for r in rr],yerr=[1.96*r['se'] for r in rr],fmt='o-');ax.axhline(0,color='k',lw=.5);ax.axvline(50,color='grey',ls='--');ax.set(title=f'Off-reference {metric} error',xlabel='Coordinate',ylabel='SMC − reservoir')
    rr=read(src/'4tax_p0.5_k1.0'/'draws.csv')
    for m in ('SMC','legacy','reservoir'):
        vals=np.sort([float(r['root_age']) for r in rr if r['method']==m]);axes[2].step(vals,np.arange(1,len(vals)+1)/len(vals),where='post',label=m)
    axes[2].legend();axes[2].set(title='Reference root-time CDF, p=.5 κ=1',xlabel='Generations',ylabel='CDF')
    fig.tight_layout();fig.savefig(out/'off_reference.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--input',type=Path,required=True);a.add_argument('--backgrounds',type=Path,required=True);a.add_argument('--output',type=Path,required=True)
    c=a.parse_args();run(c.input,c.backgrounds,c.output)
