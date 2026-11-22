"""Independent-theory validation; one chromosome is the replication unit."""
import argparse
import csv
import json
from pathlib import Path
import sys
import time
import platform
import subprocess
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.moran_spatial.theory import load_theory,make_tree,predict,rho_effective,PrescribedHistory,freeze_history
from msrcsim.arg import simulate_msrc_arg
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as AMap
from msrcsim.arg.output.tree_sequence import build_tree_sequence
from msrcsim.smc.recombination_map import PiecewiseConstantMap as SMap
from msrcsim.smc.structured import StructuredConfig,simulate_structured_smc
from msrcsim.species_tree import SpeciesTree
from msrcsim.moran import simulate_frequency_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.structured_coalescent import simulate_genealogy_reference


def moran_history(case,tree,rng):
    k=int(rng.choice(case['M']+1,p=case['pi_old']))
    rr=Rearrangement('theory_validation','suppressed_recombination','ROOT',0,k,case.get('selection',0.))
    return simulate_frequency_history(tree,rr,rng)


def prescribed(case,tree):
    if case['history']=='realized_moran':
        h=moran_history(case,tree,np.random.default_rng(case['history_seed']))
        return freeze_history(tree,h)
    segments={bid:[(case['p'],b.older_age-b.younger_age)] for bid,b in tree.branches.items()}
    return PrescribedHistory(tree,segments)


def tree_topology(tree):
    # The first coalescing pair fixes the unrooted quartet; includes pre-root coalescence.
    taxa=sorted(t.name for t in tree.nodes() if t.is_tip())
    def descendants(n):
        return {n.name} if n.is_tip() else set().union(*(descendants(c) for c in n.children))
    pair=descendants(min((n for n in tree.nodes() if not n.is_tip()),key=lambda n:n.age))
    return next(i for i,p in enumerate(({taxa[0],taxa[1]},{taxa[0],taxa[2]},{taxa[0],taxa[3]})) if pair==p or set(taxa)-pair==p)


def ts_topology(t):
    _,a,b=min((t.time(t.mrca(a,b)),a,b) for a in range(4) for b in range(a))
    pair={a,b}
    return next(i for i,p in enumerate(({0,1},{0,2},{0,3})) if pair==p or set(range(4))-pair==p)


def write_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def assess(counts,theory,tolerance,failed):
    n=int(sum(counts));q=np.array(counts)/n if n else np.full(3,np.nan)
    rows=[]
    for i,name in enumerate(('q1','q2','q3','D')):
        obs=q[i] if i<3 else q[1]-q[0];target=theory[i] if i<3 else theory[1]-theory[0]
        if n:
            variance=q[i]*(1-q[i]) if i<3 else q[0]+q[1]-(q[1]-q[0])**2
            se=float(np.sqrt(max(0,variance)/n))
            if i<3: # Wilson interval avoids zero-width intervals at q=0 or q=1.
                z=1.96;center=(obs+z*z/(2*n))/(1+z*z/n)
                half=z*np.sqrt(obs*(1-obs)/n+z*z/(4*n*n))/(1+z*z/n)
                lo,hi=center-half-target,center+half-target
            else:lo,hi=obs-target-1.96*se,obs-target+1.96*se
        else:se=lo=hi=float('nan')
        status=('NOT_DEFINED_zero_probability_pattern' if not np.isfinite(target) and n==0 else 'FAIL_incompatible' if failed else 'PASS' if lo>=-tolerance and hi<=tolerance
                else 'FAIL_mismatch' if lo>tolerance or hi<-tolerance else 'FAIL_uncertainty')
        rows.append(dict(metric=name,theory=float(target),observed=float(obs) if not failed else None,
                         difference=float(obs-target) if not failed else None,absolute_error=float(abs(obs-target)) if not failed else None,
                         mc_se=se if not failed else None,ci_low=lo if not failed else None,ci_high=hi if not failed else None,
                         tolerance=tolerance,n=n,failed_chromosomes=failed,status=status))
    return rows


def run_case(case,modules,out,n,seed,tolerance):
    start=time.perf_counter();tree=make_tree(case);fixed=case['history']!='unconditional_moran'
    hist=prescribed(case,tree) if fixed else None
    # Honest Moran histories can include absorbing counts: assess actual occupied
    # states and failures rather than discarding paths at fixation/loss.
    if fixed and case['sampling']!='unconditional':
        mass=np.prod([
                      hist.terminal_frequency(t) if bit=='1' else 1-hist.terminal_frequency(t)
                      for t,bit in zip(tree.taxa,case['sampling'])])
        if mass==0:raise ValueError('Conditioned pattern has probability zero: '+case['id'])
    positions=case['positions'];m=case['map'];amap=AMap(m['length'],tuple(m['breakpoints']),tuple(m['rates']));smap=SMap(m['length'],tuple(m['breakpoints']),tuple(m['rates']))
    theories={x:predict(modules,case,rho_effective(case,x),hist) for x in positions}
    frozen=predict(modules,case,rho_effective(case,case['reference']),hist)
    theory_rows=[]
    for x,t in theories.items():
        assert np.all(t.joint>=-1e-10) and abs(t.joint.sum()-1)<1e-9
        assert np.allclose(t.joint.sum(axis=0),t.quartet_probabilities,atol=1e-10)
        theory_rows.append(dict(position=x,rho=rho_effective(case,x),q=t.quartet_probabilities.tolist(),joint=t.joint.tolist(),pattern=t.pattern_probabilities.tolist(),diagnostics=t.diagnostics))
    methods=('ARG','SMC','single_locus');patterns=('all','1010','0101')
    counts={(method,x,z):np.zeros(3,dtype=int) for method in methods for x in positions for z in patterns}
    failures={method:[] for method in methods};raw=[]
    for i in range(n):
        seeds=np.random.SeedSequence([seed,case['case_number'],i]).generate_state(5)
        rng=np.random.default_rng(seeds[0]);h=hist if fixed else moran_history(case,tree,rng)
        if case['sampling']=='unconditional':states={taxon:int(rng.random()<h.terminal_frequency(taxon)) for taxon in tree.taxa}
        else:states=dict(zip(tree.taxa,map(int,case['sampling'])))
        pattern=''.join(map(str,states.values()))
        for method in methods:
            tops={}
            try:
                if method=='ARG':
                    r=simulate_msrc_arg(tree,amap,h,case['marker'],tuple(case['suppressed_region']),case['kappa'],'reservoir_v1',int(seeds[1]),states,record_events=False)
                    ts=build_tree_sequence(r).simplify();tops={x:ts_topology(ts.at(x)) for x in positions}
                elif method=='SMC':
                    cfg=StructuredConfig(smap,case['marker'],tuple(case['suppressed_region']),case['kappa'],int(seeds[2]),reference_position=case['reference'])
                    r=simulate_structured_smc(tree,h,cfg,states)
                    for x in positions:
                        row=next(row for row in r.trees if row['start']<=x<row['end'])
                        tops[x]=tree_topology(SpeciesTree(row['newick'],20,1,ultrametric_tolerance=1e-6))
                else:
                    for j,x in enumerate(positions):
                        g=simulate_genealogy_reference(i,tree,h,states,rho_effective(case,x),1.,np.random.default_rng(np.random.SeedSequence([int(seeds[3]),j])),False)
                        tops[x]=g.topology_index
            except (ValueError,RuntimeError) as e:
                failures[method].append(dict(replicate=i,pattern=pattern,error=str(e)))
                tops={} # Whole-chromosome comparison fails; partial output never used.
            for x,top in tops.items():
                counts[method,x,'all'][top]+=1
                if pattern in patterns:counts[method,x,pattern][top]+=1
                raw.append(dict(replicate=i,method=method,position=x,pattern=pattern,topology=int(top)))
        if (i+1)%1000==0:print(case['id'],i+1,'/',n,flush=True)
    comparisons=[]
    for method in methods:
        for x in positions:
            for z in patterns:
                if case['sampling']!='unconditional' and z not in ('all',case['sampling']):continue
                pat=case['sampling'] if z=='all' and case['sampling']!='unconditional' else z
                index=int(pat,2) if pat!='all' else None
                q=theories[x].quartet_probabilities if index is None else theories[x].conditional_quartet_probabilities[index]
                for row in assess(counts[method,x,z],q,tolerance,len(failures[method])):
                    comparisons.append(dict(case=case['id'],method=method,position=x,pattern=z,rho=rho_effective(case,x),**row))
                if method=='SMC':
                    qf=frozen.quartet_probabilities if index is None else frozen.conditional_quartet_probabilities[index]
                    for row in assess(counts[method,x,z],qf,tolerance,len(failures[method])):
                        comparisons.append(dict(case=case['id'],method='SMC_frozen_target',position=x,pattern=z,rho=rho_effective(case,case['reference']),**row))
    dest=out/case['id'];dest.mkdir(parents=True,exist_ok=True)
    (dest/'case.json').write_text(json.dumps(case,indent=2)+'\n')
    write_csv(dest/'comparisons.csv',comparisons);write_csv(dest/'draws.csv',raw)
    (dest/'theory.json').write_text(json.dumps(theory_rows,indent=2)+'\n')
    (dest/'failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    (dest/'history.json').write_text(json.dumps(hist.segments if fixed else {'law':'unconditional_moran','pi_old':case['pi_old']},indent=2)+'\n')
    summary=dict(case=case['id'],n=n,seconds=time.perf_counter()-start,failures={k:len(v) for k,v in failures.items()},statuses={s:sum(r['status']==s for r in comparisons) for s in sorted({r['status'] for r in comparisons})})
    (dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary,flush=True)
    return summary


def main():
    p=argparse.ArgumentParser();p.add_argument('--theory-source',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--replicates',type=int,default=4000);p.add_argument('--case',action='append');a=p.parse_args()
    cfg=json.loads(a.config.read_text());mods,hashes=load_theory(a.theory_source);a.output.mkdir(parents=True,exist_ok=True)
    provenance=dict(theory_commit=cfg['theory_commit'],theory_source_hashes=hashes,simulator_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[2],text=True).strip(),numpy=np.__version__,python=platform.python_version(),seed=cfg['seed'],replicates=a.replicates,quartet_tolerance=cfg['quartet_tolerance'],D_tolerance=cfg['quartet_tolerance'],interval='95% pointwise Wilson q / multinomial normal D',replication='independent chromosomes; conditional histories fixed when declared',configuration=cfg)
    (a.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    summaries=[]
    for case in cfg['cases']:
        if a.case and case['id'] not in a.case:continue
        summaries.append(run_case(case,mods,a.output,a.replicates,cfg['seed'],cfg['quartet_tolerance']))
        (a.output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')

if __name__=='__main__':main()
