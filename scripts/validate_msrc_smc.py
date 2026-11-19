"""Independent statistical validation of the neutral multispecies SMC'.

The script uses independent MSC draws plus msprime SMCPrime and Hudson ARG
simulations. Chromosomes, not adjacent local trees, are the independent units
for confidence intervals. The default suite is intentionally moderate; use
--replicates and --reference-trees to increase precision.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass, field
from itertools import combinations
import json
from pathlib import Path
import time
import tracemalloc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import msprime
import numpy as np
import scipy.stats as st

from msrcsim.species_tree import SpeciesTree, Node, _Parser
from msrcsim.smc import PiecewiseConstantMap, SMCConfig, simulate_msc_smc
from msrcsim.smc.engine import write_result


@dataclass
class SNode:
    name: str
    children: list["SNode"] = field(default_factory=list)
    age: float = 0.0


def make_species_tree(n_taxa: int, tau: float, ne: float, variable_ne=False):
    """Build an ultrametric approximately balanced tree with internals tau apart."""
    counter=0
    def build(names):
        nonlocal counter
        if len(names)==1: return SNode(names[0])
        mid=len(names)//2
        kids=[build(names[:mid]),build(names[mid:])]
        age=max(k.age for k in kids)+tau
        counter+=1
        return SNode(f"I{counter}",kids,age)
    root=build([f"T{i+1}" for i in range(n_taxa)])
    root.name="ROOT"
    def render(n,parent_age=None):
        body=n.name if not n.children else "("+",".join(render(c,n.age) for c in n.children)+")"+n.name
        return body if parent_age is None else body+":"+format(parent_age-n.age,".12g")
    newick=render(root)+";"
    branches={}
    def assign(n,depth=0):
        size=ne*(0.5 if variable_ne and depth%2 else (2.0 if variable_ne else 1.0))
        branches[n.name]={"effective_population_size":size}
        for c in n.children: assign(c,depth+1)
    assign(root)
    tree=SpeciesTree(newick,int(ne),max(1e7,root.age+100*ne),branches)
    return tree


def msprime_demography(tree):
    d=msprime.Demography()
    # With one haploid sample per species (ploidy=1), msprime's pairwise rate
    # is 1/N. Doubling sizes matches this project's 1/(2 Ne) convention.
    for b in tree.branches.values(): d.add_population(name=b.branch_id,initial_size=2*b.effective_population_size)
    for n in tree.nodes():
        if n.children:
            d.add_population_split(time=n.age,derived=[c.name for c in n.children],ancestral=n.name)
    d.sort_events()
    return d


def _metric_newick(newick):
    root=_Parser(newick).parse()
    tips=[]; clades=[]; pairs={}; internal=[]; total=0.0
    def rec(n):
        nonlocal total
        if n.is_tip():
            n.age=0.0; desc=frozenset((n.name,)); tips.append(n.name)
            return desc
        childsets=[]
        child_ages=[]
        for c in n.children:
            d=rec(c); childsets.append(d); child_ages.append(c.age+c.length)
            total+=c.length
            if not c.is_tip() and not n.is_tip(): internal.append(c.length)
        if max(child_ages)-min(child_ages)>1e-5: raise ValueError("non-ultrametric output gene tree")
        n.age=float(np.mean(child_ages)); desc=frozenset().union(*childsets)
        clades.append(desc)
        for i,a in enumerate(childsets):
            for b in childsets[i+1:]:
                for x in a:
                    for y in b: pairs[tuple(sorted((x,y)))]=n.age
        return desc
    rec(root)
    return _pack_metric(tips,root.age,clades,pairs,internal,total)


def _metric_tskit(tree, labels):
    tips=[]; clades=[]; pairs={}; internal=[]; total=0.0
    roots=tree.roots
    if len(roots)!=1: raise ValueError(f"local tree has {len(roots)} roots")
    def rec(u):
        nonlocal total
        children=list(tree.children(u))
        if not children:
            label=labels[u]; tips.append(label); return frozenset((label,))
        sets=[]
        for v in children:
            d=rec(v); sets.append(d)
            bl=tree.time(u)-tree.time(v); total+=bl
            if list(tree.children(v)): internal.append(bl)
        desc=frozenset().union(*sets); clades.append(desc)
        for i,a in enumerate(sets):
            for b in sets[i+1:]:
                for x in a:
                    for y in b: pairs[tuple(sorted((x,y)))]=tree.time(u)
        return desc
    rec(roots[0])
    return _pack_metric(tips,tree.time(roots[0]),clades,pairs,internal,total)


def _pack_metric(tips,root_age,clades,pairs,internal,total):
    tips=tuple(sorted(tips)); nontrivial=tuple(sorted(tuple(sorted(c)) for c in set(clades) if 1<len(c)<len(tips)))
    pair_names=tuple(combinations(tips,2)); pair_values=tuple(float(pairs[x]) for x in pair_names)
    qfreq={}
    for q in combinations(tips,4):
        qset=set(q); pair_clades=[]
        for c in nontrivial:
            side=tuple(sorted(qset.intersection(c)))
            if len(side)==2:
                other=tuple(sorted(qset-set(side)))
                pair_clades.append(tuple(sorted((side,other))))
        if pair_clades:
            key=";".join(",".join(x) for x in min(pair_clades))
            qfreq["/".join(q)+"="+key]=1.0
    internal_mean=float(np.mean(internal)) if internal else 0.0
    return {"tips":tips,"topology":nontrivial,"root_age":float(root_age),"pair_mean":float(np.mean(pair_values)) if pair_values else float(root_age),"pair_values":pair_values,"pair_names":pair_names,"internal_mean":internal_mean,"internal_lengths":tuple(internal),"total_length":float(total),"quartets":qfreq}


def _weighted_mean(values,weights):
    return float(np.average(values,weights=weights))

def _weighted_quantile(values,weights,q):
    idx=np.argsort(values); v=np.asarray(values)[idx]; w=np.asarray(weights)[idx]
    c=(np.cumsum(w)-0.5*w)/w.sum()
    return float(np.interp(q,c,v,left=v[0],right=v[-1]))


def _smc_chrom(tree,L,rate,seed):
    start=time.perf_counter()
    res=simulate_msc_smc(tree,SMCConfig(L,PiecewiseConstantMap.constant(L,rate),seed),np.random.default_rng(seed))
    elapsed=time.perf_counter()-start
    metrics=[_metric_newick(x.newick) for x in res.trees]
    widths=[x.end-x.start for x in res.trees]
    return metrics,widths,len(res.events),[e.event_type for e in res.events],elapsed,res


def _msprime_chrom(tree,L,rate,seed,model,full_arg=False):
    dem=msprime_demography(tree)
    samples=[msprime.SampleSet(1,population=x,ploidy=1) for x in tree.taxa]
    start=time.perf_counter()
    ts=msprime.sim_ancestry(samples=samples,demography=dem,sequence_length=L,recombination_rate=rate,
                            model=model,random_seed=seed,record_full_arg=full_arg,ploidy=1)
    elapsed=time.perf_counter()-start
    # Full-ARG recording inserts unary recombination nodes in marginal trees.
    # They are not genealogy branches; simplify before measuring tree times.
    metric_ts=ts.simplify()
    labels={int(u):str(tree.taxa[i]) for i,u in enumerate(metric_ts.samples())}
    metrics=[]; widths=[]; last_top=None; tree_changes=topo_changes=0
    for t in metric_ts.trees():
        metric=_metric_tskit(t,labels); metrics.append(metric); widths.append(t.interval.right-t.interval.left)
        if last_top is not None:
            prev=metrics[-2]
            if metric["topology"]!=last_top: topo_changes+=1
            if metric["topology"]!=last_top or metric["pair_values"]!=prev["pair_values"]: tree_changes+=1
        last_top=metric["topology"]
    if full_arg:
        re_nodes=int(np.count_nonzero(ts.tables.nodes.flags & msprime.NODE_IS_RE_EVENT))
        if re_nodes%2: raise ValueError("msprime full ARG has an unmatched recombination-node flag")
        # msprime records the two complementary ancestral products as two
        # NODE_IS_RE_EVENT nodes for each biological crossover.
        re_count=re_nodes//2
    else: re_count=None
    return metrics,widths,re_count,tree_changes,topo_changes,elapsed,ts


def _chrom_summary(metrics,widths):
    fields=("root_age","pair_mean","internal_mean","total_length")
    result={k:_weighted_mean([m[k] for m in metrics],widths) for k in fields}
    for k in ("root_age","pair_mean","internal_mean","total_length"):
        result[k+"_q10"]=_weighted_quantile([m[k] for m in metrics],widths,.1)
        result[k+"_q50"]=_weighted_quantile([m[k] for m in metrics],widths,.5)
        result[k+"_q90"]=_weighted_quantile([m[k] for m in metrics],widths,.9)
    keys=set().union(*(m["quartets"].keys() for m in metrics))
    for k in keys:
        result["qfreq:"+k]=_weighted_mean([m["quartets"].get(k,0.0) for m in metrics],widths)
    return result


def _waiting_summary(events,L):
    out={}
    for name,positions in events.items():
        positions=sorted(float(x) for x in positions)
        # The trailing interval is right-censored at the chromosome end and is
        # excluded from empirical waiting-time means.
        gaps=np.diff([0.0,*positions])
        out[name+"_count"]=len(positions)
        out[name+"_mean_wait_bp"]=float(np.mean(gaps)) if len(gaps) else np.nan
        out[name+"_median_wait_bp"]=float(np.median(gaps)) if len(gaps) else np.nan
    return out

def _waiting_samples(events):
    rows=[]
    for name,positions in events.items():
        xs=sorted(float(x) for x in positions)
        gaps=np.diff([0.0,*xs])
        rows.extend({"event_type":name,"wait_bp":float(x)} for x in gaps)
    return rows


def _transition_positions(metrics,widths):
    pos=[]; topo=[]
    x=0.0
    for i in range(1,len(metrics)):
        x+=widths[i-1]
        a,b=metrics[i-1],metrics[i]
        if a["topology"]!=b["topology"] or a["pair_values"]!=b["pair_values"]: pos.append(x)
        if a["topology"]!=b["topology"]: topo.append(x)
    return pos,topo


def _smc_hazard_audit(result,metrics,widths,recmap):
    audit=[]; prev=0.0
    for i,e in enumerate(result.events):
        left=prev; right=e.position; g=metrics[i]["total_length"]
        cum=0.0; cursor=left
        while cursor<right:
            end=min(right,recmap.next_boundary(cursor))
            cum+=(end-cursor)*recmap.rate_at(cursor)*g
            cursor=end
        audit.append(cum); prev=right
    return audit


def _summarize_method(summaries,ref_values,case_name,method):
    rows=[]
    keys=set().union(*(s.keys() for s in summaries))
    keys.update("qfreq:"+q for m in ref_values for q in m["quartets"])
    keys=sorted(keys)
    for key in keys:
        # A quartet topology absent from an entire chromosome contributes a
        # true zero. Omitting that chromosome conditions on seeing the topology
        # and inflates linked-tree frequencies, especially for rare splits.
        default=0.0 if key.startswith("qfreq:") else np.nan
        vals=np.asarray([s.get(key,default) for s in summaries],float)
        vals=vals[np.isfinite(vals)]
        if key.startswith("qfreq:"):
            # Independent MSC trees are the sampling unit for the reference.
            qkey=key[6:]
            ref=np.asarray([m["quartets"].get(qkey,0.0) for m in ref_values],float)
        else:
            if key.endswith(("_q10","_q50","_q90")): continue
            ref=np.asarray([m[key] for m in ref_values],float)
        if len(vals)<2 or len(ref)<2: continue
        test=st.ttest_ind(vals,ref,equal_var=False)
        diff=float(np.mean(vals)-np.mean(ref))
        se=float(np.sqrt(np.var(vals,ddof=1)/len(vals)+np.var(ref,ddof=1)/len(ref)))
        va=np.var(vals,ddof=1)/len(vals); vb=np.var(ref,ddof=1)/len(ref)
        den=va**2/(len(vals)-1)+vb**2/(len(ref)-1)
        df=(va+vb)**2/den if den>0 else float("inf")
        if se==0:
            pvalue=1.0 if diff==0 else 0.0
        else: pvalue=float(test.pvalue)
        crit=float(st.t.ppf(.975,df)) if np.isfinite(df) else 1.959963984540054
        rows.append({"case":case_name,"method":method,"metric":key,"reference_mean":float(np.mean(ref),),
                     "method_mean":float(np.mean(vals)),"difference":diff,"difference_ci95_low":diff-crit*se,
                     "difference_ci95_high":diff+crit*se,"standard_error":se,"welch_t":float(test.statistic) if np.isfinite(test.statistic) else (0.0 if diff==0 else np.inf),
                     "welch_df":float(df),"p_value":pvalue,"n_chromosomes":len(vals),"n_independent_trees":len(ref)})
    return rows


def _holm(rows):
    order=sorted(range(len(rows)),key=lambda i:rows[i]["p_value"]); prev=0.0; m=len(order)
    for rank,i in enumerate(order):
        adj=min(1.0,max(prev,(m-rank)*rows[i]["p_value"])); rows[i]["p_holm_all_metrics"]=adj; prev=adj


def _case_rows():
    return [
        (2,.05,100,False,1,"2taxa_short_small_low"),
        (3,.05,100,False,1,"3taxa_short_small_low"),
        (3,2.0,10000,False,5,"3taxa_long_large_moderate"),
        (4,.05,10000,True,5,"4taxa_short_large_variable_moderate"),
        (4,2.0,100,False,20,"4taxa_long_small_high"),
        (8,.05,100,True,5,"8taxa_short_small_variable_moderate"),
        (8,2.0,10000,False,20,"8taxa_long_large_high"),
    ]


def run_marginal(args,out):
    rows=[]; testrows=[]; summary=[]; waiting_rows=[]; hazard_rows=[]; waiting_samples=[]
    for idx,(n,tau,ne,var,qrate,name) in enumerate(_case_rows()):
        tree=make_species_tree(n,tau*2*ne,ne,var); L=args.length; rate=qrate/(4*ne*L)
        dem=msprime_demography(tree)
        sample_sets=[msprime.SampleSet(1,population=x,ploidy=1) for x in tree.taxa]
        ref=[]
        reference_seed=100000+idx*10000
        ref_iterator=msprime.sim_ancestry(samples=sample_sets,demography=dem,sequence_length=1,
                recombination_rate=0,ploidy=1,model="hudson",random_seed=reference_seed,
                num_replicates=args.reference_trees)
        for ts in ref_iterator:
            labels={int(u):str(tree.taxa[i]) for i,u in enumerate(ts.samples())}
            metric=_metric_tskit(ts.first(),labels); ref.append(metric)
            rows.append({"case":name,"method":"independent_MSC","replicate":len(ref)-1,"start":0.0,"end":1.0,
                         "span":1.0,"topology":metric["topology"],"root_age":metric["root_age"],"pair_mean":metric["pair_mean"],
                         "internal_mean":metric["internal_mean"],"total_length":metric["total_length"],**{"qfreq:"+k:v for k,v in metric["quartets"].items()}})
        method_summaries={"msrc_smc":[],"msprime_smc_prime":[],"msprime_hudson":[],"msrc_arg":[]}
        recomb_counts={k:[] for k in method_summaries}; change_counts={k:[] for k in method_summaries}; topo_counts={k:[] for k in method_summaries}
        tree_counts={k:[] for k in method_summaries}
        for j in range(args.replicates):
            seed=500000+idx*10000+j
            mm,ww,nrec,types,elapsed,res=_smc_chrom(tree,L,rate,seed)
            method_summaries["msrc_smc"].append(_chrom_summary(mm,ww)); recomb_counts["msrc_smc"].append(nrec)
            change_counts["msrc_smc"].append(sum(x!="no_genealogy_change" for x in types)); topo_counts["msrc_smc"].append(sum(x=="topology_change" for x in types))
            tree_counts["msrc_smc"].append(len(mm))
            smc_treepos=[e.position for e in res.events if e.event_type!="no_genealogy_change"]
            smc_toppos=[e.position for e in res.events if e.event_type=="topology_change"]
            waiting_rows.append({"case":name,"method":"msrc_smc","replicate":j,"recombination_events":nrec,
                **_waiting_summary({"recombination": [e.position for e in res.events],"tree_change":smc_treepos,"topology_change":smc_toppos},L)})
            waiting_samples.extend({"case":name,"method":"msrc_smc","replicate":j,**x} for x in _waiting_samples(
                {"recombination":[e.position for e in res.events],"tree_change":smc_treepos,"topology_change":smc_toppos}))
            hvals=_smc_hazard_audit(res,mm,ww,PiecewiseConstantMap.constant(L,rate))
            for hidx,hval in enumerate(hvals): hazard_rows.append({"case":name,"replicate":j,"event_index":hidx,"integrated_hazard":hval})
            for method,key in (("msprime_smc_prime","smc_prime"),("msprime_hudson","hudson")):
                m,w,nr,ng,nt,el,ts=_msprime_chrom(tree,L,rate,seed+1000000,key,True)
                method_summaries[method].append(_chrom_summary(m,w)); recomb_counts[method].append(nr or 0); change_counts[method].append(ng); topo_counts[method].append(nt)
                tree_counts[method].append(len(m))
                bp,tp=_transition_positions(m,w)
                waiting_rows.append({"case":name,"method":method,"replicate":j,"recombination_events":nr,
                    "recombination_count_only":nr,"tree_change_count":ng,"topology_change_count":nt,
                    **_waiting_summary({"tree_change":bp,"topology_change":tp},L)})
                waiting_samples.extend({"case":name,"method":method,"replicate":j,**x} for x in _waiting_samples({"tree_change":bp,"topology_change":tp}))
                local_start=0.0
                for mm,span in zip(m,w):
                    rows.append({"case":name,"method":method,"replicate":j,"start":local_start,"end":local_start+span,"span":span,
                                 "topology":mm["topology"],"root_age":mm["root_age"],"pair_mean":mm["pair_mean"],
                                 "internal_mean":mm["internal_mean"],"total_length":mm["total_length"],
                                 **{"qfreq:"+k:v for k,v in mm["quartets"].items()}})
                    local_start+=span
            from msrcsim.arg.engine import ARGConfig, simulate_arg
            from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as ARGMap
            from msrcsim.arg.output.tree_sequence import build_tree_sequence
            start_arg=time.perf_counter()
            arg=simulate_arg(tree,ARGConfig(L,ARGMap.constant(L,rate)),np.random.default_rng(seed+3000000))
            arg_ts=build_tree_sequence(arg).simplify()
            from msrcsim.arg.output.genealogy_export import _labels as arg_labels
            arg_metrics=[]; arg_widths=[]; arg_tree_changes=arg_top_changes=0; last=None
            labs=arg_labels(arg_ts)
            for at in arg_ts.trees():
                am=_metric_tskit(at,labs); arg_metrics.append(am); arg_widths.append(at.interval.right-at.interval.left)
                if last is not None:
                    prev=arg_metrics[-2]
                    if am["topology"]!=last: arg_top_changes+=1
                    if am["topology"]!=last or am["pair_values"]!=prev["pair_values"]: arg_tree_changes+=1
                last=am["topology"]
            method_summaries["msrc_arg"].append(_chrom_summary(arg_metrics,arg_widths))
            recomb_counts["msrc_arg"].append(sum(e.event_type=="recombination" for e in arg.events))
            change_counts["msrc_arg"].append(arg_tree_changes); topo_counts["msrc_arg"].append(arg_top_changes)
            tree_counts["msrc_arg"].append(len(arg_metrics))
            arg_recs=[e.position for e in arg.events if e.event_type=="recombination" and e.position is not None]
            arg_treepos,arg_toppos=_transition_positions(arg_metrics,arg_widths)
            waiting_rows.append({"case":name,"method":"msrc_arg","replicate":j,"recombination_events":len(arg_recs),
                **_waiting_summary({"recombination":arg_recs,"tree_change":arg_treepos,"topology_change":arg_toppos},L)})
            waiting_samples.extend({"case":name,"method":"msrc_arg","replicate":j,**x} for x in _waiting_samples(
                {"recombination":arg_recs,"tree_change":arg_treepos,"topology_change":arg_toppos}))
            local_start=0.0
            for mm,span in zip(arg_metrics,arg_widths):
                rows.append({"case":name,"method":"msrc_arg","replicate":j,"start":local_start,"end":local_start+span,"span":span,
                             "topology":mm["topology"],"root_age":mm["root_age"],"pair_mean":mm["pair_mean"],
                             "internal_mean":mm["internal_mean"],"total_length":mm["total_length"],
                             **{"qfreq:"+k:v for k,v in mm["quartets"].items()}})
                local_start+=span
            for t in res.trees:
                metric=_metric_newick(t.newick); span=t.end-t.start
                rows.append({"case":name,"method":"msrc_smc","replicate":j,"start":t.start,"end":t.end,"span":span,"topology":t.topology,
                             "root_age":metric["root_age"],"pair_mean":metric["pair_mean"],"internal_mean":metric["internal_mean"],
                             "total_length":metric["total_length"],**{"qfreq:"+k:v for k,v in metric["quartets"].items()}})
        all_tests=[]
        for method,sums in method_summaries.items():
            trows=_summarize_method(sums,ref,name,method); testrows.extend(trows); all_tests.extend(trows)
            rows.append({"case":name,"method":method,"n_chromosomes":len(sums),"mean_recombination_events":float(np.mean(recomb_counts[method])),
                         "mean_tree_changes":float(np.mean(change_counts[method])),"mean_topology_changes":float(np.mean(topo_counts[method])),
                         "mean_local_trees":float(np.mean(tree_counts[method]))})
        summary.append({"case":name,"n_taxa":n,"internal_branch_tau_over_2Ne":tau,"default_Ne":ne,"variable_branch_Ne":var,
                        "rate_proxy_q":qrate,"rate_per_bp_per_generation":rate,"chromosome_length":L,
                        "reference_independent_trees":len(ref),"chromosome_replicates":args.replicates,
                        "method_event_means":{m:{"recombinations":float(np.mean(recomb_counts[m])),"tree_changes":float(np.mean(change_counts[m])),"topology_changes":float(np.mean(topo_counts[m])),"local_trees":float(np.mean(tree_counts[m]))} for m in method_summaries},
                        "n_mean_tests_ci_excludes_zero":sum(x["difference_ci95_low"]*x["difference_ci95_high"]>0 for x in all_tests if not x["metric"].startswith("qfreq:"))})
        print(name,"done",{m:round(float(np.mean(recomb_counts[m])),1) for m in method_summaries})
    _holm(testrows)
    return rows,testrows,summary,waiting_rows,hazard_rows,waiting_samples


def _tree_at(metrics,widths,x):
    c=np.cumsum(widths)
    i=min(int(np.searchsorted(c,x,side="right")),len(metrics)-1)
    return metrics[i]

def spatial_stats(metric_chroms,L,seed,points=250):
    metrics,widths=metric_chroms
    rng=np.random.default_rng(seed); starts=rng.random(points)
    result=[]
    for frac in (.01,.05,.1,.25,.5):
        d=frac*L; xs=starts*(L-d)
        pairs=[(_tree_at(metrics,widths,x),_tree_at(metrics,widths,x+d)) for x in xs]
        same_gene=[a["topology"]==b["topology"] and np.allclose(a["pair_values"],b["pair_values"],rtol=0,atol=1e-8) for a,b in pairs]
        same_top=[a["topology"]==b["topology"] for a,b in pairs]
        rf=[len(set(a["topology"])^set(b["topology"])) for a,b in pairs]
        pta=[np.mean(a["pair_values"]) for a,b in pairs]; ptb=[np.mean(b["pair_values"]) for a,b in pairs]
        la=[a["total_length"] for a,b in pairs]; lb=[b["total_length"] for a,b in pairs]
        def corr(a,b): return float(np.corrcoef(a,b)[0,1]) if np.std(a)>0 and np.std(b)>0 else np.nan
        result.append({"distance_bp":d,"same_genealogy":float(np.mean(same_gene)),"same_topology":float(np.mean(same_top)),
                       "rooted_rf":float(np.mean(rf)),"pair_time_corr":corr(pta,ptb),"total_length_corr":corr(la,lb)})
    return result


def run_spatial(args,out):
    tree=make_species_tree(4,.05*2*100,100,False); L=args.spatial_length
    rates=[.5/(4*100*L),2/(4*100*L),8/(4*100*L)]
    labels=["low","moderate","high"]; observations=[]
    for k,(label,rate) in enumerate(zip(labels,rates)):
        for rep in range(args.spatial_replicates):
            seed=700000+k*10000+rep
            m,w,nr,types,elapsed,res=_smc_chrom(tree,L,rate,seed)
            ss=spatial_stats((m,w),L,seed+1)
            for row in ss: observations.append({"method":"msrc_smc","rate":label,"rate_per_bp_per_generation":rate,"replicate":rep,**row})
            for method,key in (("msprime_smc_prime","smc_prime"),("msprime_hudson","hudson")):
                rm,rw,re,ng,nt,el,ts=_msprime_chrom(tree,L,rate,seed+2000000,key,True)
                rs=spatial_stats((rm,rw),L,seed+2)
                for row in rs: observations.append({"method":method,"rate":label,"rate_per_bp_per_generation":rate,"replicate":rep,**row})
        print("spatial",label,"done")
    # Rate-by-method estimates with chromosome-level t intervals.
    aggregate=[]
    for method in sorted(set(x["method"] for x in observations)):
        for label in labels:
            for d in sorted(set(x["distance_bp"] for x in observations)):
                group=[x for x in observations if x["method"]==method and x["rate"]==label and x["distance_bp"]==d]
                for metric in ("same_genealogy","same_topology","rooted_rf","pair_time_corr","total_length_corr"):
                    vals=np.asarray([x[metric] for x in group],float); vals=vals[np.isfinite(vals)]
                    if not len(vals): continue
                    se=float(np.std(vals,ddof=1)/np.sqrt(len(vals))) if len(vals)>1 else np.nan
                    crit=float(st.t.ppf(.975,len(vals)-1)) if len(vals)>1 else np.nan
                    aggregate.append({"method":method,"rate":label,"distance_bp":d,"metric":metric,"mean":float(np.mean(vals)),"se_chromosome":se,
                                      "ci95_low":float(np.mean(vals)-crit*se) if len(vals)>1 else np.nan,"ci95_high":float(np.mean(vals)+crit*se) if len(vals)>1 else np.nan,"n_chromosomes":len(vals)})
    fig,axs=plt.subplots(2,2,figsize=(12,8))
    for ax,metric,title in zip(axs.flat,("same_genealogy","same_topology","rooted_rf","pair_time_corr"),
                               ("Identical genealogy","Identical topology","Rooted RF distance","Pairwise coalescence correlation")):
        for method in ("msrc_smc","msprime_smc_prime","msprime_hudson"):
            for label in labels:
                z=[x for x in aggregate if x["method"]==method and x["rate"]==label and x["metric"]==metric]
                z.sort(key=lambda x:x["distance_bp"])
                if not z: continue
                style={"msrc_smc":"-","msprime_smc_prime":"--","msprime_hudson":":"}[method]
                ax.plot([x["distance_bp"] for x in z],[x["mean"] for x in z],linestyle=style,marker="o",label=f"{method}/{label}")
                ax.fill_between([x["distance_bp"] for x in z],[x["ci95_low"] for x in z],[x["ci95_high"] for x in z],alpha=.08)
        ax.set_title(title); ax.set_xlabel("distance (bp)"); ax.set_ylabel(metric)
    axs[0,0].legend(fontsize=7,ncol=2); fig.tight_layout(); fig.savefig(out/"spatial_dependence.png",dpi=180); plt.close(fig)
    return observations,aggregate


def _region_spatial(metrics,widths,L,ranges,seed,points=100):
    rng=np.random.default_rng(seed); out=[]
    for region,intervals in ranges.items():
        for frac in (.01,.05,.1):
            d=frac*L
            valid=[(lo,hi) for lo,hi in intervals if hi-lo>d]
            if not valid: continue
            lens=np.array([b-a-d for a,b in valid],float); total=lens.sum()
            if total<=0: continue
            xs=[]
            for _ in range(points):
                k=int(rng.choice(len(valid),p=lens/total)); a,b=valid[k]
                xs.append(a+rng.random()*(b-a-d))
            pairs=[(_tree_at(metrics,widths,x),_tree_at(metrics,widths,x+d)) for x in xs]
            same=[a["topology"]==b["topology"] and np.allclose(a["pair_values"],b["pair_values"],atol=1e-8,rtol=0) for a,b in pairs]
            same_top=[a["topology"]==b["topology"] for a,b in pairs]
            rf=[len(set(a["topology"])^set(b["topology"])) for a,b in pairs]
            pa=[np.mean(a["pair_values"]) for a,b in pairs]; pb=[np.mean(b["pair_values"]) for a,b in pairs]
            la=[a["total_length"] for a,b in pairs]; lb=[b["total_length"] for a,b in pairs]
            corr=lambda x,y: float(np.corrcoef(x,y)[0,1]) if np.std(x)>0 and np.std(y)>0 else np.nan
            out.append({"region":region,"distance_bp":d,"same_genealogy":float(np.mean(same)),"same_topology":float(np.mean(same_top)),
                        "rooted_rf":float(np.mean(rf)),"pair_time_corr":corr(pa,pb),"total_length_corr":corr(la,lb)})
    return out


def run_suppression(args,out):
    tree=make_species_tree(4,10,100,False); L=args.spatial_length; a=.25*L; b=.75*L
    r0=8/(4*100*L); kappas=(0,.01,.1,.5,1.0); per_chrom=[]; ac_rows=[]
    dem=msprime_demography(tree); samples=[msprime.SampleSet(1,population=x,ploidy=1) for x in tree.taxa]
    refs=[]
    for ts in msprime.sim_ancestry(samples=samples,demography=dem,sequence_length=1,recombination_rate=0,ploidy=1,model="hudson",
                                   random_seed=61001,num_replicates=args.reference_trees):
        refs.append(_metric_tskit(ts.first(),{int(u):str(tree.taxa[i]) for i,u in enumerate(ts.samples())}))
    for ki,kappa in enumerate(kappas):
        recmap=PiecewiseConstantMap(L,(a,b),(r0,kappa*r0,r0))
        for rep in range(args.spatial_replicates):
            seed=900000+ki*10000+rep
            result=simulate_msc_smc(tree,SMCConfig(L,recmap,seed),np.random.default_rng(seed))
            metrics=[_metric_newick(t.newick) for t in result.trees]; widths=[t.end-t.start for t in result.trees]
            areas={"inside":0.0,"outside":0.0}
            weighted={key:[] for key in areas}; weightv={key:[] for key in areas}
            quart={key:{} for key in areas}
            for loc,met in zip(result.trees,metrics):
                ins=max(0.0,min(loc.end,b)-max(loc.start,a)); outlen=(loc.end-loc.start)-ins
                for region,w in (("inside",ins),("outside",outlen)):
                    if w<=0: continue
                    areas[region]+=w; weighted[region].append(met); weightv[region].append(w)
                    for q,v in met["quartets"].items(): quart[region][q]=quart[region].get(q,0.0)+v*w
            records={}
            for region in ("inside","outside"):
                mets=weighted[region]; ws=weightv[region]
                s=_chrom_summary(mets,ws)
                for key,val in s.items():
                    if key.startswith("qfreq:"): continue
                    records[region+"_"+key]=val
                for q,v in quart[region].items(): records[region+"_qfreq:"+q]=v/areas[region]
            inside_events=[e for e in result.events if a<=e.position<b]
            outside_events=[e for e in result.events if not a<=e.position<b]
            records.update({"kappa":kappa,"replicate":rep,"r0":r0,"inside_recomb_density":len(inside_events)/(b-a),
                "outside_recomb_density":len(outside_events)/(L-(b-a)),
                "inside_tree_change_density":sum(e.event_type!="no_genealogy_change" for e in inside_events)/(b-a),
                "outside_tree_change_density":sum(e.event_type!="no_genealogy_change" for e in outside_events)/(L-(b-a)),
                "inside_topology_change_density":sum(e.event_type=="topology_change" for e in inside_events)/(b-a),
                "outside_topology_change_density":sum(e.event_type=="topology_change" for e in outside_events)/(L-(b-a))})
            per_chrom.append(records)
            ranges={"inside":[(a,b)],"outside":[(0,a),(b,L)]}
            for s in _region_spatial(metrics,widths,L,ranges,seed+1): ac_rows.append({"kappa":kappa,"replicate":rep,**s})
        print("suppression kappa",kappa,"done")
    # Replicate-level summaries and uncertainty. Regions remain clustered by chromosome.
    summary=[]
    measures=["recomb_density","tree_change_density","topology_change_density","root_age","pair_mean","total_length"]
    for kappa in kappas:
        grp=[x for x in per_chrom if x["kappa"]==kappa]
        for region in ("inside","outside"):
            for measure in measures:
                key=f"{region}_{measure}" if measure.endswith("_density") else f"{region}_{measure}"
                vals=np.asarray([x[key] for x in grp if key in x],float); vals=vals[np.isfinite(vals)]
                if len(vals)<2: continue
                mean=float(np.mean(vals)); se=float(np.std(vals,ddof=1)/np.sqrt(len(vals))); c=float(st.t.ppf(.975,len(vals)-1))
                summary.append({"kappa":kappa,"region":region,"measure":measure,"mean":mean,"se_chromosome":se,
                                "ci95_low":mean-c*se,"ci95_high":mean+c*se,"n_chromosomes":len(vals)})
    ac_summary=[]
    for kappa in kappas:
        for region in ("inside","outside"):
            for metric in ("same_genealogy","same_topology","rooted_rf","pair_time_corr","total_length_corr"):
                for d in sorted(set(x["distance_bp"] for x in ac_rows)):
                    vals=np.asarray([x[metric] for x in ac_rows if x["kappa"]==kappa and x["region"]==region and x["distance_bp"]==d],float); vals=vals[np.isfinite(vals)]
                    if len(vals)<2: continue
                    mean=float(np.mean(vals)); se=float(np.std(vals,ddof=1)/np.sqrt(len(vals))); c=float(st.t.ppf(.975,len(vals)-1))
                    ac_summary.append({"kappa":kappa,"region":region,"metric":metric,"distance_bp":d,"mean":mean,"se_chromosome":se,
                                       "ci95_low":mean-c*se,"ci95_high":mean+c*se,"n_chromosomes":len(vals)})
    # Region-specific MSC marginal checks use chromosome means and cluster uncertainty.
    marginal_tests=[]
    for kappa in kappas:
        grp=[x for x in per_chrom if x["kappa"]==kappa]
        for region in ("inside","outside"):
            for metric in ("root_age","pair_mean","total_length"):
                vals=np.asarray([x[f"{region}_{metric}"] for x in grp],float); ref=np.asarray([x[metric] for x in refs],float)
                test=st.ttest_ind(vals,ref,equal_var=False); diff=float(vals.mean()-ref.mean())
                se=float(np.sqrt(vals.var(ddof=1)/len(vals)+ref.var(ddof=1)/len(ref)))
                marginal_tests.append({"kappa":kappa,"region":region,"metric":metric,"msc_reference_mean":float(ref.mean()),"smc_mean":float(vals.mean()),
                    "difference":diff,"difference_se":se,"p_value":float(test.pvalue) if np.isfinite(test.pvalue) else 1.0,
                    "n_chromosomes":len(vals),"n_independent_MSC_trees":len(ref)})
            qkeys=sorted({key.split(":",1)[1] for row in grp for key in row if key.startswith(region+"_qfreq:")})
            for qkey in qkeys:
                vals=np.asarray([x.get(region+"_qfreq:"+qkey,0.0) for x in grp],float)
                ref=np.asarray([m["quartets"].get(qkey,0.0) for m in refs],float)
                test=st.ttest_ind(vals,ref,equal_var=False); diff=float(vals.mean()-ref.mean())
                se=float(np.sqrt(vals.var(ddof=1)/len(vals)+ref.var(ddof=1)/len(ref)))
                marginal_tests.append({"kappa":kappa,"region":region,"metric":"qfreq:"+qkey,"msc_reference_mean":float(ref.mean()),"smc_mean":float(vals.mean()),
                    "difference":diff,"difference_se":se,"p_value":float(test.pvalue) if np.isfinite(test.pvalue) else 1.0,
                    "n_chromosomes":len(vals),"n_independent_MSC_trees":len(ref)})
    fig,axs=plt.subplots(2,2,figsize=(11,8))
    for ax,measure,title in zip(axs.flat,("recomb_density","tree_change_density","topology_change_density","same_topology"),
                                ("Recombination density","Tree-change density","Topology-change density","Same topology by separation")):
        for region,color in (("inside","#d36b3f"),("outside","#315a9a")):
            if measure=="same_topology":
                rows=[x for x in ac_summary if x["metric"]==measure and x["region"]==region and np.isclose(x["distance_bp"],.05*L)]
            else: rows=[x for x in summary if x["measure"]==measure and x["region"]==region]
            rows.sort(key=lambda x:x["kappa"])
            ax.errorbar([x["kappa"] for x in rows],[x["mean"] for x in rows],yerr=[1.96*x["se_chromosome"] for x in rows],marker="o",capsize=2,label=region,color=color)
        ax.set_title(title); ax.set_xlabel("kappa"); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(out/"suppressed_region_validation.png",dpi=180); plt.close(fig)
    return per_chrom,summary,ac_rows,ac_summary,marginal_tests


def run_benchmarks(out):
    data=[]
    for n in (4,8,16,32):
        tree=make_species_tree(n,0,1000,False)
        for L in (5000,20000):
            for q in (.5,2.0):
                rate=q/(4*1000*L)
                tracemalloc.start(); start=time.perf_counter()
                res=simulate_msc_smc(tree,SMCConfig(L,PiecewiseConstantMap.constant(L,rate),800000+n+L),np.random.default_rng(800000+n+L))
                runtime=time.perf_counter()-start; _,peak=tracemalloc.get_traced_memory(); tracemalloc.stop()
                data.append({"taxa":n,"chromosome_length_bp":L,"rate_per_bp_per_generation":rate,"q_proxy":q,"recombination_events":len(res.events),
                             "local_genealogies":len(res.trees),"genealogy_changes":sum(e.event_type!="no_genealogy_change" for e in res.events),
                             "topology_changes":sum(e.event_type=="topology_change" for e in res.events),"runtime_seconds":runtime,"peak_python_bytes":peak,
                             "seconds_per_event":runtime/max(1,len(res.events)),"seconds_per_transition":runtime/max(1,len(res.events)),"seconds_per_local_tree":runtime/max(1,len(res.trees)),"method":"msrc_smc_corrected"})
                from msrcsim.arg.engine import ARGConfig, simulate_arg
                from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as ARGMap
                from msrcsim.arg.output.tree_sequence import build_tree_sequence
                tracemalloc.start(); start=time.perf_counter()
                arg=simulate_arg(tree,ARGConfig(L,ARGMap.constant(L,rate)),np.random.default_rng(900000+n+L))
                ts=build_tree_sequence(arg).simplify()
                arg_runtime=time.perf_counter()-start; _,arg_peak=tracemalloc.get_traced_memory(); tracemalloc.stop()
                from msrcsim.arg.output.genealogy_export import _labels as arg_labels
                labels=arg_labels(ts)
                metrics=[]; widths=[]
                for t in ts.trees(): metrics.append(_metric_tskit(t,labels)); widths.append(t.interval.right-t.interval.left)
                bp,tp=_transition_positions(metrics,widths)
                arg_recombinations=sum(e.event_type=="recombination" for e in arg.events)
                data.append({"taxa":n,"chromosome_length_bp":L,"rate_per_bp_per_generation":rate,"q_proxy":q,
                             "recombination_events":arg_recombinations,"local_genealogies":len(metrics),"genealogy_changes":len(bp),"topology_changes":len(tp),
                             "runtime_seconds":arg_runtime,"peak_python_bytes":arg_peak,"seconds_per_event":arg_runtime/max(1,arg_recombinations),
                             "seconds_per_transition":arg_runtime/max(1,arg_recombinations),"seconds_per_local_tree":arg_runtime/max(1,len(metrics)),"method":"msrc_arg"})
                print("benchmark",n,L,q,"smc",len(res.events),round(runtime,3),"arg",arg_recombinations,round(arg_runtime,3))
    return data


def write_csv(path,rows):
    import csv
    keys=sorted(set().union(*(x.keys() for x in rows))) if rows else []
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def resummarize_marginals(out):
    """Recompute marginal tests from raw span rows, including absent quartets as zero."""
    import csv
    with (out/"marginal_genealogies.csv").open() as f:
        data=list(csv.DictReader(f))
    cases=list(dict.fromkeys(x["case"] for x in data))
    rows=[]; continuous=("root_age","pair_mean","internal_mean","total_length")
    for case in cases:
        ref=[x for x in data if x["case"]==case and x["method"]=="independent_MSC"]
        qkeys=sorted({key[6:] for x in ref for key,value in x.items() if key.startswith("qfreq:") and value not in ("",None)})
        for method in ("msrc_smc","msprime_smc_prime","msprime_hudson","msrc_arg"):
            linked=[x for x in data if x["case"]==case and x["method"]==method and x.get("span")]
            reps=sorted({x["replicate"] for x in linked})
            values={metric:[] for metric in (*continuous,*("qfreq:"+q for q in qkeys))}
            for rep in reps:
                trees=[x for x in linked if x["replicate"]==rep]
                weights=np.asarray([float(x["span"]) for x in trees]); weights=weights/weights.sum()
                for metric in continuous:
                    values[metric].append(float(np.average([float(x[metric]) for x in trees],weights=weights)))
                for qkey in qkeys:
                    col="qfreq:"+qkey
                    values[col].append(float(np.average([float(x[col]) if x.get(col) else 0.0 for x in trees],weights=weights)))
            for metric,linked_values in values.items():
                refs=np.asarray([float(x[metric[6:]]) if metric.startswith("qfreq:") and x.get(metric[6:]) else
                                 (float(x["qfreq:"+metric[6:]]) if metric.startswith("qfreq:") and x.get("qfreq:"+metric[6:]) else 0.0)
                                 if metric.startswith("qfreq:") else float(x[metric]) for x in ref])
                vals=np.asarray(linked_values,float)
                if len(vals)<2 or len(refs)<2: continue
                diff=float(vals.mean()-refs.mean()); va=float(vals.var(ddof=1)/len(vals)); vb=float(refs.var(ddof=1)/len(refs)); se=float(np.sqrt(va+vb))
                den=va**2/(len(vals)-1)+vb**2/(len(refs)-1); df=(va+vb)**2/den if den>0 else float("inf")
                stat=0.0 if se==0 and diff==0 else (np.inf if se==0 else diff/se)
                pvalue=1.0 if se==0 and diff==0 else (0.0 if se==0 else float(2*st.t.sf(abs(stat),df)))
                crit=float(st.t.ppf(.975,df)) if np.isfinite(df) else 1.959963984540054
                rows.append({"case":case,"method":method,"metric":metric,"reference_mean":float(refs.mean()),"method_mean":float(vals.mean()),
                             "difference":diff,"difference_ci95_low":diff-crit*se,"difference_ci95_high":diff+crit*se,
                             "standard_error":se,"welch_t":stat,"welch_df":float(df),"p_value":pvalue,
                             "n_chromosomes":len(vals),"n_independent_trees":len(refs)})
    _holm(rows)
    write_csv(out/"marginal_metric_tests.csv",rows)
    summary_path=out/"validation_summary.json"
    if summary_path.exists():
        report=json.loads(summary_path.read_text())
        report["marginal_test_results"]={m:{"comparisons":sum(x["method"]==m for x in rows),
            "holm_significant_at_0.05":sum(x["method"]==m and x["p_holm_all_metrics"]<.05 for x in rows),
            "95pct_CI_excludes_zero":sum(x["method"]==m and x["difference_ci95_low"]*x["difference_ci95_high"]>0 for x in rows)}
            for m in ("msrc_smc","msprime_smc_prime","msprime_hudson","msrc_arg")}
        summary_path.write_text(json.dumps(report,indent=2)+"\n")
    return rows

def plot_diagnostics(out):
    import csv
    def read(name):
        path=out/name
        return list(csv.DictReader(path.open())) if path.exists() else []
    tests=read("marginal_metric_tests.csv")
    if tests:
        cases=list(dict.fromkeys(x["case"] for x in tests if x["method"]=="msrc_smc" and x["metric"]=="root_age"))
        fig,ax=plt.subplots(figsize=(11,5))
        ax.axhline(1,color="black",lw=1)
        for method,marker in (("msrc_smc","o"),("msprime_smc_prime","s"),("msprime_hudson","^"),("msrc_arg","D")):
            chosen=[x for x in tests if x["method"]==method and x["metric"]=="root_age"]
            chosen.sort(key=lambda x:cases.index(x["case"]))
            if not chosen: continue
            x=[cases.index(z["case"]) for z in chosen]
            ref=np.asarray([float(z["reference_mean"]) for z in chosen]); ratio=np.asarray([float(z["method_mean"]) for z in chosen])/ref
            low=np.asarray([float(z["difference_ci95_low"]) for z in chosen])/ref
            high=np.asarray([float(z["difference_ci95_high"]) for z in chosen])/ref
            ax.errorbar(x,ratio,yerr=np.vstack([ratio-(1+low),1+high-ratio]),marker=marker,capsize=2,label=method)
        ax.set_xticks(range(len(cases)),cases,rotation=35,ha="right");ax.set_ylabel("Root-age mean / independent MSC mean")
        ax.legend(frameon=False);fig.tight_layout();fig.savefig(out/"root_age_comparison.png",dpi=180);plt.close(fig)
        quartet=[x for x in tests if x["metric"].startswith("qfreq:")]
        if quartet:
            fig,ax=plt.subplots(figsize=(6,6));ax.plot([0,1],[0,1],color="black",lw=1)
            for method,marker in (("msrc_smc","o"),("msprime_smc_prime","s"),("msprime_hudson","^"),("msrc_arg","D")):
                z=[x for x in quartet if x["method"]==method]
                if z: ax.scatter([float(a["reference_mean"]) for a in z],[float(a["method_mean"]) for a in z],s=15,alpha=.55,label=method,marker=marker)
            ax.set(xlim=(-.03,1.03),ylim=(-.03,1.03),xlabel="Independent MSC quartet frequency",ylabel="Linked chromosome span-weighted frequency")
            ax.legend(frameon=False);fig.tight_layout();fig.savefig(out/"quartet_frequency_comparison.png",dpi=180);plt.close(fig)
    tree_rows=read("marginal_genealogies.csv")
    if tree_rows:
        cases=["4taxa_long_small_high","8taxa_long_large_high"]
        fig,axs=plt.subplots(1,2,figsize=(11,4))
        methods=("independent_MSC","msrc_smc","msprime_smc_prime","msprime_hudson","msrc_arg")
        colors=dict(zip(methods,("black","#d36b3f","#315a9a","#65a765","#8064a2")))
        for ax,case in zip(axs,cases):
            all_values=[float(r["root_age"]) for r in tree_rows if r["case"]==case and r.get("root_age")]
            bins=np.linspace(min(all_values),max(all_values),36) if all_values else 35
            for method in methods:
                z=[r for r in tree_rows if r["case"]==case and r["method"]==method and r.get("root_age")]
                if not z: continue
                vals=np.asarray([float(r["root_age"]) for r in z]); weights=np.asarray([float(r.get("span") or 1) for r in z])
                ax.hist(vals,bins=bins,weights=weights/weights.sum(),density=False,histtype="step",lw=1.5,color=colors[method],label=method)
            ax.set_title(case);ax.set_xlabel("Root age (generations)");ax.set_ylabel("Span-weighted frequency")
        axs[0].legend(fontsize=7,frameon=False);fig.tight_layout();fig.savefig(out/"root_age_distributions.png",dpi=180);plt.close(fig)
    waits=read("waiting_distances.csv")
    if waits:
        fig,axs=plt.subplots(1,3,figsize=(13,4),sharey=True)
        methods=("msrc_smc","msprime_smc_prime","msprime_hudson","msrc_arg")
        for ax,event in zip(axs,("recombination","tree_change","topology_change")):
            for method in methods:
                vals=[float(r["wait_bp"]) for r in waits if r["method"]==method and r["event_type"]==event]
                if vals: ax.hist(vals,bins=35,density=True,histtype="step",lw=1.3,label=method)
            ax.set_title(event.replace("_"," "));ax.set_xlabel("Waiting distance (bp)")
        axs[0].set_ylabel("Density");axs[0].legend(fontsize=7,frameon=False);fig.tight_layout();fig.savefig(out/"waiting_distance_distributions.png",dpi=180);plt.close(fig)
    bench=read("benchmark.csv")
    if bench:
        fig,axs=plt.subplots(1,2,figsize=(11,4))
        for method in sorted(set(x["method"] for x in bench)):
            for L in sorted(set(int(x["chromosome_length_bp"]) for x in bench)):
                z=[x for x in bench if x["method"]==method and int(x["chromosome_length_bp"])==L]
                z.sort(key=lambda a:int(a["taxa"]))
                axs[0].plot([int(x["taxa"]) for x in z],[float(x["runtime_seconds"]) for x in z],marker="o",label=f"{method}, L={L}")
        axs[0].set(xlabel="Taxa",ylabel="Runtime (s)",yscale="log");axs[0].legend(fontsize=7,frameon=False)
        taxa=sorted(set(int(x["taxa"]) for x in bench))
        for method in sorted(set(x["method"] for x in bench)):
            vals=[np.mean([float(x["seconds_per_transition"]) for x in bench if x["method"]==method and int(x["taxa"])==n]) for n in taxa]
            axs[1].plot(taxa,vals,marker="o",label=method)
        axs[1].set(xlabel="Taxa",ylabel="Seconds per recombination event");axs[1].legend(fontsize=7,frameon=False)
        fig.tight_layout();fig.savefig(out/"runtime_scaling.png",dpi=180);plt.close(fig)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="validation_msrc_smc")
    p.add_argument("--replicates",type=int,default=20)
    p.add_argument("--reference-trees",type=int,default=1000)
    p.add_argument("--length",type=int,default=10000)
    p.add_argument("--spatial-replicates",type=int,default=20)
    p.add_argument("--spatial-length",type=int,default=20000)
    p.add_argument("--skip-benchmark",action="store_true")
    p.add_argument("--suppression-only",action="store_true",help="run only the region-of-suppressed-recombination suite")
    p.add_argument("--benchmark-only",action="store_true",help="run only the corrected SMC and neutral ARG performance comparison")
    p.add_argument("--resummarize-only",action="store_true",help="recompute marginal tests from existing span-level CSV output")
    args=p.parse_args(); out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    if args.resummarize_only:
        resummarize_marginals(out); plot_diagnostics(out); return
    if args.benchmark_only:
        write_csv(out/"benchmark.csv",run_benchmarks(out)); plot_diagnostics(out); return
    if args.suppression_only:
        suppressed,suppressed_summary,suppressed_ac,suppressed_ac_summary,suppressed_tests=run_suppression(args,out)
        write_csv(out/"suppressed_region_chromosomes.csv",suppressed)
        write_csv(out/"suppressed_region_summary.csv",suppressed_summary)
        write_csv(out/"suppressed_region_autocorrelation.csv",suppressed_ac)
        write_csv(out/"suppressed_region_autocorrelation_summary.csv",suppressed_ac_summary)
        write_csv(out/"suppressed_region_marginal_tests.csv",suppressed_tests)
        return
    marginal,tests,summaries,waiting,hazards,waiting_samples=run_marginal(args,out)
    write_csv(out/"marginal_genealogies.csv",marginal); write_csv(out/"marginal_metric_tests.csv",tests)
    write_csv(out/"waiting_distance_summary.csv",waiting); write_csv(out/"waiting_distances.csv",waiting_samples); write_csv(out/"integrated_hazards.csv",hazards)
    hazard_values=np.asarray([x["integrated_hazard"] for x in hazards],float)
    hazard_values=hazard_values[np.isfinite(hazard_values)]
    hazard_summary={"n":int(len(hazard_values)),"mean":float(hazard_values.mean()) if len(hazard_values) else None,
                    "expected_mean":1.0,"ks_exponential_pvalue":float(st.kstest(hazard_values,"expon").pvalue) if len(hazard_values)>0 else None}
    spatial,spatial_agg=run_spatial(args,out)
    write_csv(out/"spatial_chromosome_metrics.csv",spatial); write_csv(out/"spatial_summary.csv",spatial_agg)
    suppressed,suppressed_summary,suppressed_ac,suppressed_ac_summary,suppressed_tests=run_suppression(args,out)
    write_csv(out/"suppressed_region_chromosomes.csv",suppressed)
    write_csv(out/"suppressed_region_summary.csv",suppressed_summary)
    write_csv(out/"suppressed_region_autocorrelation.csv",suppressed_ac)
    write_csv(out/"suppressed_region_autocorrelation_summary.csv",suppressed_ac_summary)
    write_csv(out/"suppressed_region_marginal_tests.csv",suppressed_tests)
    if not args.skip_benchmark:
        bench=run_benchmarks(out); write_csv(out/"benchmark.csv",bench)
    plot_diagnostics(out)
    report={"versions":{"msprime":msprime.__version__,"numpy":np.__version__,"scipy":__import__('scipy').__version__},
            "parameters":{"chromosome_replicates":args.replicates,"independent_MSC_gene_trees_per_case":args.reference_trees,
                          "marginal_length_bp":args.length,"spatial_replicates_per_rate_and_method":args.spatial_replicates,
                          "spatial_length_bp":args.spatial_length},"marginal_cases":summaries,
            "integrated_hazard_exponential_audit":hazard_summary,
            "limitations":["ipcoal unavailable in this environment; published conditional waiting-distance calculator not called",
                           "msprime was the independent SMC-prime and exact Hudson ARG reference",
                           "Monte Carlo confidence intervals use independent chromosome replicates for linked results"]}
    (out/"validation_summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print("wrote validation to",out)
if __name__=="__main__": main()
