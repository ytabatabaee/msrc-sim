"""Plot the runnable window example without treating its intervals as replicates."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def read(p):
    with p.open() as f:return list(csv.DictReader(f))


def main(src,dest):
    meta=json.loads((src/'metadata.json').read_text());cfg=meta['configuration']['structured_smc']
    tracks=read(src/'quartet_tracks.csv');events=read(src/'recombination_events.csv');origin=meta['origin_age']
    fig,ax=plt.subplots(4,1,figsize=(11,10),gridspec_kw={'height_ratios':[1,1,1,1.5]})
    p=cfg['frequency_history']['p'];ax[0].plot([0,origin],[p,p],color='purple');ax[0].axvline(origin,color='k',ls='--');ax[0].set(xlabel='Generations backward',ylabel='p(t)',ylim=(0,1),title='One shared frequency history; background partition erased at origin')
    L=cfg['chromosome_length'];a,b=cfg['suppressed_region'];start,end=meta['scan_start'],meta['scan_end']
    ax[1].plot([0,L],[1,1],color='blue',label='Baseline r₀ (constant)');ax[1].plot([0,a,a,b,b,L],[1,1,meta['cross_retention'],meta['cross_retention'],1,1],color='orange',label='Cross retention κ')
    ax[1].add_patch(Rectangle((start,0),end-start,1.05,fill=False,edgecolor='green',lw=2,label='Simulated window'))
    ax[1].axvline(meta['marker'],color='purple',ls=':');ax[1].text(meta['marker']+10,.5,'Marker');ax[1].set(xlim=(-20,L),ylim=(-.05,1.15),xlabel='Map coordinate',ylabel='Relative rate / κ');ax[1].legend(loc='upper right',fontsize=8)
    colors=['#4c78a8','#f58518','#54a24b']
    for row in tracks:
        lo,hi=float(row['start']),float(row['end']);q=int(row['quartet_index'])
        ax[2].broken_barh([(lo,hi-lo)],(.1,.6),facecolors=colors[q] if q>=0 else 'grey')
        ax[2].text((lo+hi)/2,.4,f'q{q+1}',ha='center',va='center',fontsize=8)
    for e in events:
        x=float(e['position']);ax[2].axvline(x,color='k',lw=.6)
        if e['background_changed']=='True':ax[2].plot(x,.85,'rv',ms=6)
    ax[2].axvline(meta['reference_position'],color='purple',ls='--');ax[2].set(xlim=(start,end),ylim=(0,1),yticks=[],xlabel='Window coordinate',title='Local quartet tracks: candidate cuts (black), background-mark changes (red)')
    xx=[float(r['start']) for r in tracks]+[float(tracks[-1]['end'])];yy=[float(r['root_age']) for r in tracks];ax[3].step(xx,yy+[yy[-1]],where='post');ax[3].set(xlim=(start,end),xlabel='Window coordinate',ylabel='Root age (generations)',title='Dependent dated genealogies (reference marginal frozen at u*)')
    fig.suptitle('Structured SMC′ window example: an explicit approximation',fontsize=14);fig.tight_layout();dest.parent.mkdir(parents=True,exist_ok=True);fig.savefig(dest,dpi=160);plt.close(fig)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--input',type=Path,required=True);a.add_argument('--output',type=Path,required=True);c=a.parse_args();main(c.input,c.output)
