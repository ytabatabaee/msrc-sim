from __future__ import annotations
import numpy as np
from .species_tree import SpeciesTree
from .population_process import PopulationFrequencyHistory
from .lineages import Lineage,GenealogyEvent,CoalescenceRecord,GenealogyResult

def quartet_index(desc):
    pairs=[frozenset(x) for x in desc if len(x)==2]
    taxa=sorted(set().union(*desc))
    mapping=[({taxa[0],taxa[1]},{taxa[2],taxa[3]}),({taxa[0],taxa[2]},{taxa[1],taxa[3]}),({taxa[0],taxa[3]},{taxa[1],taxa[2]})]
    for i,(a,b) in enumerate(mapping):
        if frozenset(a) in pairs or frozenset(b) in pairs: return i
    return 0

def _topology_from_coalescences(coals, taxa):
    if len(taxa) != 4:
        return -1
    first=coals[0]; pair=frozenset(first.descendant_set_1.split(',')+first.descendant_set_2.split(','))
    taxa=list(taxa); sets=[frozenset([taxa[0],taxa[1]]),frozenset([taxa[0],taxa[2]]),frozenset([taxa[0],taxa[3]])]
    return sets.index(pair) if pair in sets else (sets.index(frozenset(set(taxa)-set(pair))) if frozenset(set(taxa)-set(pair)) in sets else 0)

def simulate_genealogy_reference(locus_id,tree:SpeciesTree,history:PopulationFrequencyHistory,sampled_arrangements,recombination_rate,effective_fraction,rng,record_events=True):
    """Simulate backward genealogy with an explicit structured-origin boundary.

    A generated frequency history carries ``structured_origin_age``. Before
    that age, the process is the finite-frequency structured coalescent. At
    and beyond it, arrangement labels cease to have biological meaning and
    the process becomes an ordinary unstructured MSC within each population
    branch. Older histories without this provenance use the historical
    ``tree.root.age + tree.root_extension`` boundary as a compatibility
    fallback.
    """
    rho=recombination_rate*effective_fraction
    lineages=[]; next_id=0
    for taxon in tree.taxa:
        lineages.append(Lineage(next_id,frozenset([taxon]),tree.branch_for_tip(taxon),int(sampled_arrangements[taxon]),taxon,0.0)); next_id+=1
    t=0.0; events=[]; coals=[]; coal_times=[]; eidx=0
    max_time=tree.root.age+tree.root_extension
    origin_age=getattr(history,"structured_origin_age",None)
    structured_origin_age=max_time if origin_age is None else float(origin_age)
    while len(lineages)>1:
        # Move lineages at species-tree boundaries before drawing events.
        moved=False
        for L in lineages:
            b=tree.branches[L.population_id]
            if t>=b.older_age-1e-10 and b.parent_branch_id is not None:
                old=L.population_id; L.population_id=b.parent_branch_id; moved=True
                if record_events:
                    events.append(GenealogyEvent(locus_id,eidx,t,'species_boundary',old,L.population_id,str(L.lineage_id),'',','.join(sorted(L.descendants)),str(L.arrangement),str(L.arrangement),None,None,None,None,None,None,None)); eidx+=1
        if moved: continue

        structured = t < structured_origin_age-1e-10
        possible=[]
        if structured:
            # Structured MSRC: switching and same-background coalescence use
            # the frequency history exactly as recorded on each branch.
            for L in lineages:
                p=history.frequency_at(L.population_id,min(t,tree.branches[L.population_id].older_age))
                m01=rho*p; m10=rho*(1-p); rate=m01 if L.arrangement==0 else m10
                if rate>0: possible.append(('switch',(L,),rate,p,m01,m10,None,None))
            for i in range(len(lineages)):
                for j in range(i+1,len(lineages)):
                    a,b=lineages[i],lineages[j]
                    if a.population_id!=b.population_id or a.arrangement!=b.arrangement: continue
                    p=history.frequency_at(a.population_id,min(t,tree.branches[a.population_id].older_age)); N=tree.branches[a.population_id].effective_population_size
                    l0=1/(2*N*max(1-p,1/(2*N))); l1=1/(2*N*max(p,1/(2*N))); rate=l0 if a.arrangement==0 else l1
                    possible.append(('coal',(a,b),rate,p,rho*p,rho*(1-p),l0,l1))
        else:
            # Pre-rearrangement ordinary ancestor: labels are ignored,
            # switching stops, and every same-population pair coalesces at
            # the ordinary MSC rate 1/(2 Ne).
            for i in range(len(lineages)):
                for j in range(i+1,len(lineages)):
                    a,b=lineages[i],lineages[j]
                    if a.population_id!=b.population_id: continue
                    rate=1/(2*tree.branches[a.population_id].effective_population_size)
                    possible.append(('coal',(a,b),rate,None,0.0,0.0,rate,rate))

        total=sum(x[2] for x in possible)
        next_species_boundary=min([tree.branches[L.population_id].older_age for L in lineages if tree.branches[L.population_id].parent_branch_id is not None]+[float('inf')])
        if structured:
            next_frequency_boundary=min([history.next_frequency_boundary(L.population_id, t) for L in lineages]+[float('inf')])
            next_boundary=min(next_species_boundary,next_frequency_boundary,structured_origin_age)
        else:
            next_boundary=next_species_boundary
        if total<=0:
            if next_boundary<float('inf') and next_boundary>t+1e-12:
                t=next_boundary; continue
            raise RuntimeError("genealogy has no event and no finite species boundary")
        wait=rng.exponential(1/total)
        if t+wait>=next_boundary:
            t=next_boundary; continue
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
            a,b=objs; parent=Lineage(next_id,a.descendants|b.descendants,a.population_id,0 if not structured else a.arrangement,'',t)
            bl1=t-a.birth_time; bl2=t-b.birth_time; parent.node_newick=f"({a.node_newick}:{bl1:.8f},{b.node_newick}:{bl2:.8f})"; next_id+=1
            lineages=[L for L in lineages if L not in (a,b)]+[parent]; coal_times.append(t)
            coals.append(CoalescenceRecord(locus_id,len(coals),t,a.population_id,a.lineage_id,b.lineage_id,parent.lineage_id,parent.arrangement,','.join(sorted(a.descendants)),','.join(sorted(b.descendants))))
            if record_events:
                events.append(GenealogyEvent(locus_id,eidx,t,'coalescence',a.population_id,a.population_id,f"{a.lineage_id};{b.lineage_id}",str(parent.lineage_id),f"{','.join(sorted(a.descendants))}|{','.join(sorted(b.descendants))}",f"{a.arrangement};{b.arrangement}",str(parent.arrangement),rate,total,p,m01,m10,l0,l1)); eidx+=1
    nw=lineages[0].node_newick+';'
    # topology from first coalescence descendant sets
    top=_topology_from_coalescences(coals, tree.taxa)
    return GenealogyResult(locus_id,nw,top,tuple(coal_times),tuple(events),tuple(coals))

def simulate_msc_genealogy(locus_id,tree:SpeciesTree,rng,record_events=True):
    """Simulate an ordinary unstructured MSC genealogy on the species tree."""
    lineages=[]; next_id=0
    for taxon in tree.taxa:
        lineages.append(Lineage(next_id,frozenset([taxon]),tree.branch_for_tip(taxon),0,taxon,0.0)); next_id+=1
    t=0.0; events=[]; coals=[]; coal_times=[]; eidx=0
    while len(lineages)>1:
        moved=False
        for L in lineages:
            b=tree.branches[L.population_id]
            if t>=b.older_age-1e-10 and b.parent_branch_id is not None:
                old=L.population_id; L.population_id=b.parent_branch_id; moved=True
                if record_events:
                    events.append(GenealogyEvent(locus_id,eidx,t,'species_boundary',old,L.population_id,str(L.lineage_id),'',','.join(sorted(L.descendants)),'0','0',None,None,None,None,None,None,None)); eidx+=1
        if moved: continue
        possible=[]
        for i in range(len(lineages)):
            for j in range(i+1,len(lineages)):
                a,b=lineages[i],lineages[j]
                if a.population_id!=b.population_id: continue
                N=tree.branches[a.population_id].effective_population_size
                possible.append((a,b,1/(2*N)))
        total=sum(x[2] for x in possible)
        next_species_boundary=min([tree.branches[L.population_id].older_age for L in lineages if tree.branches[L.population_id].parent_branch_id is not None]+[float('inf')])
        if total<=0:
            t=next_species_boundary
            continue
        wait=rng.exponential(1/total)
        if t+wait>=next_species_boundary:
            t=next_species_boundary
            continue
        t+=wait; u=rng.random()*total; c=0; chosen=possible[-1]
        for x in possible:
            c+=x[2]
            if u<=c:
                chosen=x; break
        a,b,rate=chosen
        parent=Lineage(next_id,a.descendants|b.descendants,a.population_id,0,'',t)
        bl1=t-a.birth_time; bl2=t-b.birth_time; parent.node_newick=f"({a.node_newick}:{bl1:.8f},{b.node_newick}:{bl2:.8f})"; next_id+=1
        lineages=[L for L in lineages if L not in (a,b)]+[parent]; coal_times.append(t)
        coals.append(CoalescenceRecord(locus_id,len(coals),t,a.population_id,a.lineage_id,b.lineage_id,parent.lineage_id,0,','.join(sorted(a.descendants)),','.join(sorted(b.descendants))))
        if record_events:
            events.append(GenealogyEvent(locus_id,eidx,t,'coalescence',a.population_id,a.population_id,f"{a.lineage_id};{b.lineage_id}",str(parent.lineage_id),f"{','.join(sorted(a.descendants))}|{','.join(sorted(b.descendants))}",'0;0','0',rate,total,None,None,None,rate,rate)); eidx+=1
    nw=lineages[0].node_newick+';'
    top=_topology_from_coalescences(coals, tree.taxa)
    return GenealogyResult(locus_id,nw,top,tuple(coal_times),tuple(events),tuple(coals))


def _remove_pair(group, rng):
    """Draw an unordered pair uniformly from a mutable lineage list."""
    n = len(group)
    i = int(rng.integers(n))
    j = int(rng.integers(n - 1))
    if j >= i:
        j += 1
    a, b = group[i], group[j]
    for index in sorted((i, j), reverse=True):
        group.pop(index)
    return a, b


def simulate_genealogy(locus_id, tree: SpeciesTree, history: PopulationFrequencyHistory,
                       sampled_arrangements, recombination_rate, effective_fraction,
                       rng, record_events=True):
    """Exact grouped implementation of the structured backward process.

    The original lineage/pair enumeration remains available as
    ``simulate_genealogy_reference``.  This implementation groups lineages by
    population and background, calculating aggregate hazards analytically and
    selecting a lineage or pair uniformly within the selected group.  It uses
    the same rates, boundaries, and transition rules; changing RNG consumption
    is intentional and does not imply a different stochastic model.
    """
    rho = recombination_rate * effective_fraction
    groups = {}
    next_id = 0
    for taxon in tree.taxa:
        lineage = Lineage(next_id, frozenset([taxon]), tree.branch_for_tip(taxon),
                          int(sampled_arrangements[taxon]), taxon, 0.0)
        groups.setdefault((lineage.population_id, lineage.arrangement), []).append(lineage)
        next_id += 1
    t = 0.0
    events, coals, coal_times = [], [], []
    eidx = 0
    max_time = tree.root.age + tree.root_extension
    origin_value = getattr(history, "structured_origin_age", None)
    structured_origin_age = max_time if origin_value is None else float(origin_value)

    def all_lineages():
        return [lineage for members in groups.values() for lineage in members]

    def rebuild_population_groups():
        rebuilt = {}
        for lineage in all_lineages():
            rebuilt.setdefault((lineage.population_id, lineage.arrangement), []).append(lineage)
        return rebuilt

    while sum(len(members) for members in groups.values()) > 1:
        moved = False
        for members in list(groups.values()):
            for lineage in members:
                branch = tree.branches[lineage.population_id]
                if t >= branch.older_age - 1e-10 and branch.parent_branch_id is not None:
                    old = lineage.population_id
                    lineage.population_id = branch.parent_branch_id
                    moved = True
                    if record_events:
                        events.append(GenealogyEvent(
                            locus_id, eidx, t, "species_boundary", old,
                            lineage.population_id, str(lineage.lineage_id), "",
                            ",".join(sorted(lineage.descendants)), str(lineage.arrangement),
                            str(lineage.arrangement), None, None, None, None, None, None, None))
                        eidx += 1
        if moved:
            groups = rebuild_population_groups()
            continue

        structured = t < structured_origin_age - 1e-10
        possible = []
        if structured:
            # One frequency lookup per occupied population, rather than one
            # lookup for every lineage and every same-background pair.
            population_frequency = {}
            for population_id, _arrangement in groups:
                branch = tree.branches[population_id]
                population_frequency[population_id] = history.frequency_at(
                    population_id, min(t, branch.older_age))
            for (population_id, arrangement), members in groups.items():
                p = population_frequency[population_id]
                rate = rho * p if arrangement == 0 else rho * (1 - p)
                if rate > 0:
                    possible.append(("switch", (population_id, arrangement),
                                    len(members) * rate, p, rho * p, rho * (1 - p), None, None))
                if len(members) >= 2:
                    branch = tree.branches[population_id]
                    c0 = 1 / (2 * branch.effective_population_size * max(1 - p, 1 / (2 * branch.effective_population_size)))
                    c1 = 1 / (2 * branch.effective_population_size * max(p, 1 / (2 * branch.effective_population_size)))
                    rate = c1 if arrangement else c0
                    possible.append(("coal", (population_id, arrangement),
                                    len(members) * (len(members) - 1) * rate / 2,
                                    p, rho * p, rho * (1 - p), c0, c1))
        else:
            by_population = {}
            for (population_id, _arrangement), members in groups.items():
                by_population.setdefault(population_id, []).extend(members)
            for population_id, members in by_population.items():
                if len(members) < 2:
                    continue
                rate = 1 / (2 * tree.branches[population_id].effective_population_size)
                possible.append(("coal", (population_id, None),
                                len(members) * (len(members) - 1) * rate / 2,
                                None, 0.0, 0.0, rate, rate))

        total = sum(option[2] for option in possible)
        boundaries = [tree.branches[lineage.population_id].older_age
                      for lineage in all_lineages()
                      if tree.branches[lineage.population_id].parent_branch_id is not None]
        next_species_boundary = min(boundaries + [float("inf")])
        if structured:
            frequency_boundaries = [history.next_frequency_boundary(population_id, t)
                                    for population_id, _arrangement in groups]
            next_boundary = min(next_species_boundary, min(frequency_boundaries + [float("inf")]),
                                structured_origin_age)
        else:
            next_boundary = next_species_boundary
        if total <= 0:
            if next_boundary < float("inf") and next_boundary > t + 1e-12:
                t = next_boundary
                continue
            raise RuntimeError("genealogy has no event and no finite species boundary")
        wait = rng.exponential(1 / total)
        if t + wait >= next_boundary:
            t = next_boundary
            continue
        t += wait
        draw = rng.random() * total
        cumulative = 0.0
        chosen = possible[-1]
        for option in possible:
            cumulative += option[2]
            if draw <= cumulative:
                chosen = option
                break
        kind, key, rate, p, m01, m10, l0, l1 = chosen
        population_id, arrangement = key
        if kind == "switch":
            members = groups[key]
            index = int(rng.integers(len(members)))
            lineage = members.pop(index)
            old = lineage.arrangement
            lineage.arrangement = 1 - old
            groups.setdefault((population_id, lineage.arrangement), []).append(lineage)
            if not members:
                del groups[key]
            if record_events:
                events.append(GenealogyEvent(
                    locus_id, eidx, t, "switch", population_id, population_id,
                    str(lineage.lineage_id), "", ",".join(sorted(lineage.descendants)),
                    str(old), str(lineage.arrangement), rate, total, p, m01, m10, l0, l1))
                eidx += 1
        else:
            if arrangement is None:
                members = []
                for group_key in list(groups):
                    if group_key[0] == population_id:
                        members.extend(groups.pop(group_key))
            else:
                members = groups[key]
            a, b = _remove_pair(members, rng)
            if members:
                groups[key] = members
            parent_arrangement = 0 if not structured else a.arrangement
            parent = Lineage(next_id, a.descendants | b.descendants, population_id,
                             parent_arrangement, "", t)
            bl1 = t - a.birth_time
            bl2 = t - b.birth_time
            parent.node_newick = f"({a.node_newick}:{bl1:.8f},{b.node_newick}:{bl2:.8f})"
            next_id += 1
            groups.setdefault((population_id, parent_arrangement), []).append(parent)
            coal_times.append(t)
            coals.append(CoalescenceRecord(
                locus_id, len(coals), t, population_id, a.lineage_id, b.lineage_id,
                parent.lineage_id, parent.arrangement, ",".join(sorted(a.descendants)),
                ",".join(sorted(b.descendants))))
            if record_events:
                events.append(GenealogyEvent(
                    locus_id, eidx, t, "coalescence", population_id, population_id,
                    f"{a.lineage_id};{b.lineage_id}", str(parent.lineage_id),
                    f"{','.join(sorted(a.descendants))}|{','.join(sorted(b.descendants))}",
                    f"{a.arrangement};{b.arrangement}", str(parent.arrangement),
                    rate, total, p, m01, m10, l0, l1))
                eidx += 1
    final = all_lineages()[0]
    return GenealogyResult(locus_id, final.node_newick + ";",
                           _topology_from_coalescences(coals, tree.taxa),
                           tuple(coal_times), tuple(events), tuple(coals))
