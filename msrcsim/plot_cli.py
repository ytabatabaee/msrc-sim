from __future__ import annotations
import argparse,csv,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import matplotlib.pyplot as plt


def _rows(path):
    with open(path,newline='') as h: return list(csv.DictReader(h))

def _f(row,k,default=np.nan):
    try:return float(row.get(k,default))
    except:return default

def _pattern_class(p):
    p=str(p).zfill(4)
    if p in {'0000','1111'}: return '4:0'
    if p.count('1') in {1,3}: return '3:1'
    if p in {'0011','1100'}: return 'concordant 2:2'
    return 'discordant 2:2'

def main():
    ap=argparse.ArgumentParser(description='Create automated MSRC summary figures')
    ap.add_argument('--input',required=True)
    ap.add_argument('--output',default='figures')
    ap.add_argument('--format',choices=['png','pdf','svg'],default='png')
    args=ap.parse_args(); rows=[r for r in _rows(args.input) if str(r.get('accepted','True')).lower() not in {'false','0'}]
    out=Path(args.output); out.mkdir(parents=True,exist_ok=True); ext=args.format
    # simplex
    fig,ax=plt.subplots(figsize=(7,6)); verts=np.array([[0,0],[1,0],[.5,np.sqrt(3)/2],[0,0]]); ax.plot(verts[:,0],verts[:,1])
    for cls in ['4:0','3:1','concordant 2:2','discordant 2:2']:
        rr=[r for r in rows if _pattern_class(r.get('terminal_pattern',''))==cls]
        if rr: ax.scatter([_f(r,'simplex_x') for r in rr],[_f(r,'simplex_y') for r in rr],s=18,alpha=.7,label=cls)
    # MSC arms centroid to vertices
    c=np.array([.5,np.sqrt(3)/6])
    for v in verts[:3]: ax.plot([c[0],v[0]],[c[1],v[1]],linestyle='--',linewidth=1)
    ax.set_aspect('equal'); ax.axis('off'); ax.legend(); ax.set_title('Quartet simplex and MSC arms'); fig.tight_layout(); fig.savefig(out/f'quartet_simplex.{ext}',dpi=250); plt.close(fig)
    # distance by class
    classes=defaultdict(list)
    for r in rows: classes[_pattern_class(r.get('terminal_pattern',''))].append(_f(r,'distance_to_nearest_msc_arm'))
    order=['4:0','3:1','concordant 2:2','discordant 2:2']; vals=[classes[x] for x in order]
    fig,ax=plt.subplots(figsize=(8,5)); ax.boxplot(vals,tick_labels=order,showfliers=False); ax.set_ylabel('Distance to nearest MSC arm'); ax.tick_params(axis='x',rotation=20); fig.tight_layout(); fig.savefig(out/f'off_arm_distance.{ext}',dpi=250); plt.close(fig)
    # terminal pattern prevalence
    counts=Counter(str(r.get('terminal_pattern','')).zfill(4) for r in rows); pats=sorted(counts)
    fig,ax=plt.subplots(figsize=(10,5)); ax.bar(pats,[counts[x] for x in pats]); ax.set_ylabel('Replicates'); ax.set_xlabel('Terminal arrangement pattern'); ax.tick_params(axis='x',rotation=45); fig.tight_layout(); fig.savefig(out/f'terminal_pattern_prevalence.{ext}',dpi=250); plt.close(fig)
    # AIC evidence
    ds=[_f(r,'delta_aic_network_vs_msc') for r in rows if np.isfinite(_f(r,'delta_aic_network_vs_msc'))]
    fig,ax=plt.subplots(figsize=(7,5)); ax.hist(ds,bins=30); ax.axvline(0,linestyle='--'); ax.axvline(-4,linestyle=':'); ax.set_xlabel('Delta AIC (network - MSC)'); ax.set_ylabel('Replicates'); fig.tight_layout(); fig.savefig(out/f'model_comparison.{ext}',dpi=250); plt.close(fig)
    manifest={'input':args.input,'num_rows':len(rows),'figures':[f'quartet_simplex.{ext}',f'off_arm_distance.{ext}',f'terminal_pattern_prevalence.{ext}',f'model_comparison.{ext}']}
    (out/'figure_manifest.json').write_text(json.dumps(manifest,indent=2)); print(f'Wrote {out}')
if __name__=='__main__': main()
