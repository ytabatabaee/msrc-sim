"""Reproducible four-taxon SMC' map comparison and diagnostic figure."""
from pathlib import Path
import argparse
import numpy as np
import matplotlib.pyplot as plt
from msrcsim.species_tree import SpeciesTree
from msrcsim.smc import PiecewiseConstantMap, SMCConfig, simulate_msc_smc
from msrcsim.smc.engine import write_result

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output",default="smc_demo",help="directory for reproducible trees and figure")
OUT=Path(parser.parse_args().output)
OUT.mkdir(exist_ok=True)
tree=SpeciesTree("((A:10000,B:10000)AB:5000,(C:10000,D:10000)CD:5000)ROOT;",10000,100000)
L=200000
maps={
    "constant_low":PiecewiseConstantMap.constant(L,5e-9),
    "constant":PiecewiseConstantMap.constant(L,1e-8),
    "constant_high":PiecewiseConstantMap.constant(L,2e-8),
    "suppressed":PiecewiseConstantMap(L,(50000.,150000.),(1e-8,1e-9,1e-8)),
}
colors={"constant_low":"#65a765","constant":"#315a9a","constant_high":"#8064a2","suppressed":"#d36b3f"}
results={}
for i,(name,recmap) in enumerate(maps.items()):
    results[name]=simulate_msc_smc(tree,SMCConfig(L,recmap,20261116+i),np.random.default_rng(20261116+i))
    write_result(results[name],OUT/name)

fig,axes=plt.subplots(6,1,figsize=(12,12),sharex=True,gridspec_kw={"height_ratios":[1,1,1,1,0.7,1.5]})
names=("constant_low","constant","constant_high","suppressed")
for row,name in enumerate(names):
    r=results[name]
    # Display a categorical topology code over each local interval.
    cats=sorted({t.topology for t in r.trees})
    cat={c:i for i,c in enumerate(cats)}
    ax=axes[row]
    for t in r.trees:
        ax.plot([t.start,t.end],[cat[t.topology]]*2,color=colors[name],lw=5,solid_capstyle="butt")
    for e in r.events:
        if e.event_type!="no_genealogy_change": ax.axvline(e.position,color="black",lw=.55,alpha=.4)
    ax.set_ylabel(name+"\ntopology"); ax.set_yticks([])
for name,recmap in maps.items():
    xx=[0,*recmap.breakpoints,L]; yy=[*recmap.rates,recmap.rates[-1]]
    axes[4].step(xx,yy,where="post",color=colors[name],label=name,alpha=.9)
axes[4].set_ylabel("r(x)"); axes[4].ticklabel_format(axis="y",style="sci",scilimits=(0,0))
axes[4].axvspan(50000,150000,color="#e9c8b8",alpha=.25)
# Autocorrelation: fraction of pairs at each genomic lag with identical local topology.
grid=np.linspace(0,L,1001)
for name,r in results.items():
    intervals=r.trees
    def at(x):
        for t in intervals:
            if t.start<=x<t.end: return t.topology
        return intervals[-1].topology
    lag_bp=np.linspace(0,50000,60); ac=[]
    for d in lag_bp:
        xs=np.linspace(0,L-d-1,300)
        ac.append(np.mean([at(x)==at(x+d) for x in xs]))
    axes[5].plot(lag_bp,ac,label=name,color=colors[name])
axes[5].set_ylabel("P(same topology)"); axes[5].set_xlabel("Genomic position / distance (bp)")
axes[5].legend(frameon=False); axes[4].legend(frameon=False,ncol=4); fig.tight_layout()
fig.savefig(OUT/"smc_map_comparison.png",dpi=180)
fig.savefig(OUT/"smc_map_comparison.pdf")
for name,r in results.items():
    print(name,"trees",len(r.trees),"recombination events",len(r.events),"changed",sum(e.event_type!="no_genealogy_change" for e in r.events))
