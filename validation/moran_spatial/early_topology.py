"""Check pectinate topology mass fixed before all four lineages share a population."""
import argparse,json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.moran_spatial.run import (load_theory,make_tree,prescribed,predict,rho_effective,AMap,SMap,simulate_msrc_arg,build_tree_sequence,ts_topology,StructuredConfig,simulate_structured_smc,SpeciesTree,tree_topology,simulate_genealogy_reference,write_csv)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--theory-source',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--replicates',type=int,default=4000);a=p.parse_args();mods,_=load_theory(a.theory_source)
    c=next(c for c in json.loads(a.config.read_text())['cases'] if c['id']=='pectinate_M20_unconditional');t=make_tree(c);h=prescribed(c,t);x=c['reference'];rho=rho_effective(c,x);theory=predict(mods,c,rho,h);target=theory.early_topology_contributions;m=c['map'];amap=AMap(m['length'],tuple(m['breakpoints']),tuple(m['rates']));smap=SMap(m['length'],tuple(m['breakpoints']),tuple(m['rates']));counts={m:np.zeros(4,int) for m in ('ARG','SMC','single_locus')}
    for i in range(a.replicates):
        s=np.random.SeedSequence([20261124,i]).generate_state(4);rng=np.random.default_rng(s[0]);states={taxon:int(rng.random()<h.terminal_frequency(taxon)) for taxon in t.taxa}
        arg=simulate_msrc_arg(t,amap,h,c['marker'],tuple(c['suppressed_region']),c['kappa'],'reservoir_v1',int(s[1]),states,record_events=False);ts=build_tree_sequence(arg).simplify().at(x);age=min(ts.time(ts.mrca(i,j)) for i in range(4) for j in range(i));counts['ARG'][ts_topology(ts) if age<t.root.age else 3]+=1
        cfg=StructuredConfig(smap,c['marker'],tuple(c['suppressed_region']),c['kappa'],int(s[2]),reference_position=x);r=simulate_structured_smc(t,h,cfg,states);row=next(row for row in r.trees if row['start']<=x<row['end']);g=SpeciesTree(row['newick'],20,1,ultrametric_tolerance=1e-6);age=min(n.age for n in g.nodes() if not n.is_tip());counts['SMC'][tree_topology(g) if age<t.root.age else 3]+=1
        g=simulate_genealogy_reference(i,t,h,states,rho,1,np.random.default_rng(s[3]),False);counts['single_locus'][g.topology_index if g.coalescences[0].time<t.root.age else 3]+=1
    rows=[]
    for method,cnt in counts.items():
        for k,name in enumerate(('early_q1','early_q2','early_q3')):
            q=cnt[k]/a.replicates;se=np.sqrt(q*(1-q)/a.replicates);diff=q-target[k];lo=diff-1.96*se;hi=diff+1.96*se
            rows.append(dict(method=method,metric=name,theory=target[k],observed=q,difference=diff,mc_se=se,ci_low=lo,ci_high=hi,tolerance=.05,status='PASS' if lo>=-.05 and hi<=.05 else 'FAIL',n=a.replicates))
    a.output.mkdir(parents=True,exist_ok=True);write_csv(a.output/'comparisons.csv',rows);(a.output/'provenance.json').write_text(json.dumps(dict(seed=20261124,n=a.replicates,case=c,position=x,criterion='95% early topology probability error interval within +/-0.05; no histories resampled'),indent=2)+'\n')
    print(rows,flush=True)
