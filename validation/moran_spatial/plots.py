"""Publication-oriented quartet profiles; uncertainty uses independent chromosomes."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def profiles(root,cases,ids,name):
    plt.rcParams.update({'font.size':9,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(5,len(ids),figsize=(6.6*len(ids),12),squeeze=False,sharex='col',constrained_layout=True)
    colors={'ARG':'#c06014','SMC':'#196caa','single_locus':'#438c39'}
    for col,id in enumerate(ids):
        c=next(c for c in cases if c['id']==id)
        rows=list(csv.DictReader((root/id/'comparisons.csv').open()))
        for j,metric in enumerate(('q1','q2','q3','D')):
            ax=axes[j,col]
            arg=sorted([r for r in rows if r['method']=='ARG' and r['pattern']=='all' and r['metric']==metric],key=lambda r:float(r['position']))
            xs=np.array([float(r['position']) for r in arg]);target=np.array([float(r['theory']) for r in arg]);ax.plot(xs,target,'k-',lw=1.5,label='Exact multispecies theory')
            frozen=sorted([r for r in rows if r['method']=='SMC_frozen_target' and r['pattern']=='all' and r['metric']==metric],key=lambda r:float(r['position']))
            ax.plot(xs,[float(r['theory']) for r in frozen],color='.5',ls=':',lw=1,label='Frozen reference target')
            for method in ('ARG','SMC'):
                rs=sorted([r for r in rows if r['method']==method and r['pattern']=='all' and r['metric']==metric],key=lambda r:float(r['position']))
                good=[r for r in rs if r['observed']!='']
                if good:
                    x=[float(r['position'])+(-.5 if method=='ARG' else .5) for r in good];y=np.array([float(r['observed']) for r in good]);q=np.array([float(r['theory']) for r in good]);lo=q+np.array([float(r['ci_low']) for r in good]);hi=q+np.array([float(r['ci_high']) for r in good])
                    ax.errorbar(x,y,yerr=[y-lo,hi-y],fmt='o' if method=='ARG' else 's',ms=3,capsize=2,color=colors[method],label="SMC′" if method=='SMC' else method)
            ax.axvspan(*c['suppressed_region'],color='#859477',alpha=.12);ax.axvline(c['reference'],color='.4',ls='--',lw=.7);ax.set_ylabel(metric if metric=='D' else '$q_'+metric[-1]+'$');ax.grid(alpha=.15)
            if j==0:
                ax.set_title(c['shape'].capitalize()+' tree; '+('unconditional samples' if c['sampling']=='unconditional' else 'Z='+c['sampling'])+'\n'+f"M(root)={c['M']}; s={c.get('selection',0)}; κ={c['kappa']}; reference=50 bp")
                ax.legend(fontsize=7,loc='best')
        ax=axes[4,col];m=c['map'];x=sorted({0,m['length'],*m['breakpoints'],*c['suppressed_region']});y=[]
        for a,b in zip(x,x[1:]):
            mid=(a+b)/2;r=m['rates'][int(np.searchsorted(m['breakpoints'],mid,side='right'))];y.append(r*(c['kappa'] if c['suppressed_region'][0]<=mid<c['suppressed_region'][1] else 1))
        ax.stairs(y,x,color='k');ax.set_xlabel('Genomic position (bp)');ax.set_ylabel('Cross-background\nmap rate / bp / generation');ax.axvspan(*c['suppressed_region'],color='#859477',alpha=.12)
    for ext in ('png','svg'):
        path=root/(name+'.'+ext);fig.savefig(path,dpi=220)
        if ext=='svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--config',type=Path,required=True);a=p.parse_args();c=json.loads(a.config.read_text())['cases']
    profiles(a.output,c,['balanced_M40_cond1010','pectinate_M40_cond0101'],'quartet_profiles')
    profiles(a.output,c,['balanced_plateau_cond1010','pectinate_plateau_cond0101'],'suppressed_region_profiles')

    if all((a.output/case/'comparisons.csv').exists() for case in ('balanced_stochastic_M10','pectinate_stochastic_M20')):
        profiles(a.output,c,['balanced_stochastic_M10','pectinate_stochastic_M20'],'stochastic_moran_profiles')
