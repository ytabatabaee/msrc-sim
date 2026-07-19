from __future__ import annotations
import numpy as np
from .species_tree import SpeciesTree
from .wright_fisher import FrequencyHistory
from .lineages import Lineage,GenealogyEvent,CoalescenceRecord,GenealogyResult

def quartet_index(desc):
    pairs=[frozenset(x) for x in desc if len(x)==2]
    taxa=sorted(set().union(*desc))
    mapping=[({taxa[0],taxa[1]},{taxa[2],taxa[3]}),({taxa[0],taxa[2]},{taxa[1],taxa[3]}),({taxa[0],taxa[3]},{taxa[1],taxa[2]})]
    for i,(a,b) in enumerate(mapping):
        if frozenset(a) in pairs or frozenset(b) in pairs: return i
    return 0

def simulate_genealogy(locus_id,tree:SpeciesTree,history:FrequencyHistory,sampled_arrangements,recombination_rate,effective_fraction,rng,record_events=True):
    rho=recombination_rate*effective_fraction
    lineages=[]; next_id=0
    for taxon in tree.taxa:
        lineages.append(Lineage(next_id,frozenset([taxon]),tree.branch_for_tip(taxon),int(sampled_arrangements[taxon]),taxon,0.0)); next_id+=1
    t=0.0; events=[]; coals=[]; coal_times=[]; eidx=0
    max_time=tree.root.age+tree.root_extension
    while len(lineages)>1:
        # move lineages at boundaries
        moved=False
        for L in lineages:
            b=tree.branches[L.population_id]
            if t>=b.older_age-1e-10 and b.parent_branch_id is not None:
                old=L.population_id; L.population_id=b.parent_branch_id; moved=True
                if record_events:
                    events.append(GenealogyEvent(locus_id,eidx,t,'species_boundary',old,L.population_id,str(L.lineage_id),'',','.join(sorted(L.descendants)),str(L.arrangement),str(L.arrangement),None,None,None,None,None,None,None)); eidx+=1
        if moved: continue
        possible=[]
        # calculate per-lineage switching
        for L in lineages:
            p=history.frequency_at(L.population_id,min(t,tree.branches[L.population_id].older_age))
            m01=rho*p; m10=rho*(1-p); rate=m01 if L.arrangement==0 else m10
            if rate>0: possible.append(('switch',(L,),rate,p,m01,m10,None,None))
        # coalescence pairs
        for i in range(len(lineages)):
            for j in range(i+1,len(lineages)):
                a,b=lineages[i],lineages[j]
                if a.population_id!=b.population_id or a.arrangement!=b.arrangement: continue
                p=history.frequency_at(a.population_id,min(t,tree.branches[a.population_id].older_age)); N=tree.branches[a.population_id].effective_population_size
                l0=1/(2*N*max(1-p,1/(2*N))); l1=1/(2*N*max(p,1/(2*N))); rate=l0 if a.arrangement==0 else l1
                possible.append(('coal',(a,b),rate,p,rho*p,rho*(1-p),l0,l1))
        total=sum(x[2] for x in possible)
        next_boundary=min([tree.branches[L.population_id].older_age for L in lineages if tree.branches[L.population_id].parent_branch_id is not None]+[float('inf')])
        if total<=0:
            if next_boundary<float('inf'): t=next_boundary; continue
            # root fallback: ordinary panmictic coalescence
            a,b=lineages[0],lineages[1]; wait=rng.exponential(2*tree.branches[a.population_id].effective_population_size); t+=wait; chosen=('coal',(a,b),1/wait,0.5,0,0,1,1)
        else:
            wait=rng.exponential(1/total)
            if t+wait>=next_boundary: t=next_boundary; continue
            t+=wait; u=rng.random()*total; c=0; chosen=possible[-1]
            for x in possible:
                c+=x[2]
                if u<=c: chosen=x; break
        kind,objs,rate,p,m01,m10,l0,l1=chosen
        if kind=='switch':
            L=objs[0]; old=L.arrangement; L.arrangement=1-old
            if record_events:
                events.append(GenealogyEvent(locus_id,eidx,t,'switch',L.population_id,L.population_id,str(L.lineage_id),'',','.join(sorted(L.descendants)),str(old),str(L.arrangement),rate,total,p,m01,m10,l0,l1)); eidx+=1
        else:
            a,b=objs; parent=Lineage(next_id,a.descendants|b.descendants,a.population_id,a.arrangement,'',t)
            bl1=t-a.birth_time; bl2=t-b.birth_time; parent.node_newick=f"({a.node_newick}:{bl1:.8f},{b.node_newick}:{bl2:.8f})"; next_id+=1
            lineages=[L for L in lineages if L not in (a,b)]+[parent]; coal_times.append(t)
            coals.append(CoalescenceRecord(locus_id,len(coals),t,a.population_id,a.lineage_id,b.lineage_id,parent.lineage_id,a.arrangement,','.join(sorted(a.descendants)),','.join(sorted(b.descendants))))
            if record_events:
                events.append(GenealogyEvent(locus_id,eidx,t,'coalescence',a.population_id,a.population_id,f"{a.lineage_id};{b.lineage_id}",str(parent.lineage_id),f"{','.join(sorted(a.descendants))}|{','.join(sorted(b.descendants))}",f"{a.arrangement};{b.arrangement}",str(parent.arrangement),rate,total,p,m01,m10,l0,l1)); eidx+=1
    nw=lineages[0].node_newick+';'
    # topology from first coalescence descendant sets
    first=coals[0]; pair=frozenset(first.descendant_set_1.split(',')+first.descendant_set_2.split(','))
    taxa=list(tree.taxa); sets=[frozenset([taxa[0],taxa[1]]),frozenset([taxa[0],taxa[2]]),frozenset([taxa[0],taxa[3]])]
    top=sets.index(pair) if pair in sets else (sets.index(frozenset(set(taxa)-set(pair))) if frozenset(set(taxa)-set(pair)) in sets else 0)
    return GenealogyResult(locus_id,nw,top,tuple(coal_times),tuple(events),tuple(coals))
