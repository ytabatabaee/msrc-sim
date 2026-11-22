"""Explicit terminal-pattern conditioning with a fixed number of prior draws.

Draw the unmodified forward Moran/history and terminal-sampling law. Simulate
chromosomes only for the declared patterns. This is posterior conditioning on
Z, not filtering on simulator success, frequency survival, or genealogy.
The attempt count is fixed in advance; failed simulations are never replaced.
"""
import argparse,json,csv,time
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.moran_spatial.run import (load_theory,make_tree,moran_history,rho_effective,predict,AMap,SMap,simulate_msrc_arg,build_tree_sequence,ts_topology,StructuredConfig,simulate_structured_smc,SpeciesTree,tree_topology,simulate_genealogy_reference,assess,write_csv)


def run(c,mods,out,attempts):
    t=make_tree(c);m=c['map'];amap=AMap(m['length'],tuple(m['breakpoints']),tuple(m['rates']));smap=SMap(m['length'],tuple(m['breakpoints']),tuple(m['rates']));x=c['reference'];rho=rho_effective(c,x);theory=predict(mods,c,rho)
    patterns=('1010','0101');methods=('ARG','SMC','single_locus');counts={(p,m):np.zeros(3,int) for p in patterns for m in methods};failures={p:{m:[] for m in methods} for p in patterns};selected={p:0 for p in patterns};start=time.perf_counter()
    for i in range(attempts):
        seeds=np.random.SeedSequence([20261123,c['case_number'],i]).generate_state(4);rng=np.random.default_rng(seeds[0]);h=moran_history(c,t,rng);states={taxon:int(rng.random()<h.terminal_frequency(taxon)) for taxon in t.taxa};p=''.join(map(str,states.values()))
        if p not in patterns:continue
        selected[p]+=1
        for j,method in enumerate(methods):
            try:
                if method=='ARG':
                    r=simulate_msrc_arg(t,amap,h,c['marker'],tuple(c['suppressed_region']),c['kappa'],'reservoir_v1',int(seeds[1]),states,record_events=False);q=ts_topology(build_tree_sequence(r).simplify().at(x))
                elif method=='SMC':
                    cfg=StructuredConfig(smap,c['marker'],tuple(c['suppressed_region']),c['kappa'],int(seeds[2]),reference_position=x);r=simulate_structured_smc(t,h,cfg,states);row=next(row for row in r.trees if row['start']<=x<row['end']);q=tree_topology(SpeciesTree(row['newick'],20,1,ultrametric_tolerance=1e-6))
                else:q=simulate_genealogy_reference(i,t,h,states,rho,1,np.random.default_rng(seeds[3]),False).topology_index
                counts[p,method][q]+=1
            except (ValueError,RuntimeError) as e:failures[p][method].append(dict(attempt=i,error=str(e)))
        if sum(selected.values())%1000==0:print(c['id'],'posterior draws',selected,'prior attempts',i+1,flush=True)
    rows=[]
    for p in patterns:
        for method in methods:
            for row in assess(counts[p,method],theory.conditional_quartet_probabilities[int(p,2)],.05,len(failures[p][method])):
                rows.append(dict(case=c['id'],method=method,position=x,pattern=p,**row))
    out.mkdir(parents=True,exist_ok=True);write_csv(out/'comparisons.csv',rows)
    summary=dict(case=c['id'],prior_attempts=attempts,selected=selected,expected_pattern_probabilities={p:float(theory.pattern_probabilities[int(p,2)]) for p in patterns},failures=failures,seconds=time.perf_counter()-start,seed=20261123,selection_rule='only sampled terminal pattern, before any genealogy; fixed attempt budget; no replacement on simulator failure')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--theory-source',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--attempts',type=int,default=160000);a=p.parse_args();mods,_=load_theory(a.theory_source);cases=json.loads(a.config.read_text())['cases']
    for c in cases:
        if c['id'] in ('balanced_stochastic_M10','pectinate_stochastic_M20'):run(c,mods,a.output/c['id'],a.attempts)
