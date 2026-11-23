"""One resource-bounded benchmark process; production RNG/event logic is unchanged."""
import argparse,collections,cProfile,csv,dataclasses,gzip,hashlib,json,math,os,pickle,pstats,resource,signal,sys,time,traceback
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from msrcsim.species_tree import SpeciesTree
from msrcsim.arg import simulate_msrc_arg
from msrcsim.arg.engine import ARGConfig,simulate_arg,write_event_log
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.arg.output.tree_sequence import build_tree_sequence
from msrcsim.arg.output.genealogy_export import export_genealogies
from msrcsim.moran import simulate_frequency_history
from msrcsim.rearrangement import Rearrangement


def species(n,ne,step,extension):
    names=list('ABCD') if n==4 else [f'T{i:02d}' for i in range(n)]
    def rec(ns):
        if len(ns)==1:return ns[0],0.
        j=len(ns)//2;l,la=rec(ns[:j]);r,ra=rec(ns[j:]);age=max(la,ra)+step
        label='AB' if ns==['A','B'] else 'CD' if ns==['C','D'] else ''
        return f'({l}:{age-la},{r}:{age-ra})'+label,age
    newick,_=rec(names)
    return SpeciesTree(newick+'ROOT;',ne,extension)


def rss():return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)


def topology(tree,n=4):
    if n!=4:return -1
    _,a,b=min((tree.time(tree.mrca(i,j)),i,j) for i in range(n) for j in range(i))
    pair={a,b};return next(i for i,p in enumerate(({0,1},{0,2},{0,3})) if pair==p or set(range(4))-pair==p)


def validate(ts,tree):
    taxa=tuple(tree.taxa);samples=list(ts.samples());assert len(samples)==len(taxa)
    assert tuple(json.loads(ts.node(u).metadata)['taxon'] for u in samples)==taxa
    species_clades=sorted([(frozenset(tree.descendants(node.name)),node.age) for node in tree.nodes()],key=lambda x:x[1])
    species_mrca={}
    covered=0.;minimum=float('inf');eligible=0
    # Validate the complete marginal interval set, not just selected positions.
    for gt in ts.trees():
        assert gt.num_roots==1 and gt.num_samples()==len(taxa)
        left,right=gt.interval;assert left==covered and right>left;covered=right
        descendants={u:{taxa[i]} for i,u in enumerate(samples)}
        for u in gt.nodes(order='postorder'):
            if u not in descendants:descendants[u]=set().union(*(descendants[c] for c in gt.children(u)))
            p=gt.parent(u)
            if p!=-1:
                dt=gt.time(p)-gt.time(u);assert dt>=0;minimum=min(minimum,dt)
            if not gt.is_sample(u):
                # A local ancestor cannot predate the species MRCA of its descendants.
                names=frozenset(descendants[u])
                if names not in species_mrca:species_mrca[names]=next(age for desc,age in species_clades if names<=desc)
                assert gt.time(u)>=species_mrca[names]-1e-8;eligible+=1
        assert gt.num_edges==2*len(taxa)-2
    assert covered==ts.sequence_length
    return dict(valid=True,covered_length=covered,sample_taxa=list(taxa),minimum_branch_generations=minimum,validated_internal_nodes=eligible)


def validate_background_events(result,tree,history,states):
    live=dict(enumerate(states[t] for t in tree.taxa));young=old=0
    origin=history.structured_origin_age
    for e in result.events:
        ids=e.lineage_ids
        if e.event_type=='recombination':
            live.pop(ids[0]);left,right=map(int,e.background_after.split(','));live[ids[1]]=left;live[ids[2]]=right
        elif e.event_type=='background_transfer':live[ids[0]]=int(e.background_after)
        elif e.event_type=='coalescence':
            a,b=live.pop(ids[0]),live.pop(ids[1]);pop=e.population_before[0];ne=tree.branches[pop].effective_population_size
            if e.time<origin-1e-10:
                assert a==b
                p=history.frequency_at(pop,e.time);ph=p if a else 1-p;assert ph>0
                expected=1/(2*ne*ph);young+=1
            else:expected=1/(2*ne);old+=1
            assert np.isclose(e.rate,expected,rtol=1e-12,atol=0);live[ids[2]]=a
        if e.event_type in ('recombination','background_transfer'):assert 0<=e.position<=result.sequence_length
    return dict(valid=True,structured_coalescences=young,ordinary_ancestral_coalescences=old)


class Limit(Exception):pass


def run(job,out):
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();stage='setup';trace=[];profile=cProfile.Profile() if job.get('profile') else None
    summary=dict(job=job,status='running',wall_seconds=None,peak_rss_bytes=None,validation_algorithm='cached_species_MRCA_v1');last_events=0;event_counts=collections.Counter();last_snapshot={};result=None
    deadline=job.get('seconds',30);memory=job.get('rss_gib',1.5)*1024**3
    def sample(sig,frame):
        nonlocal last_events,last_snapshot
        snap=dict(wall=time.perf_counter()-start,stage=stage,peak_rss_bytes=rss())
        f=frame
        while f:
            if f.f_code.co_name=='simulate_arg' and f.f_code.co_filename.endswith('/arg/engine.py'):
                loc=f.f_locals;ls=loc.get('lineages',[]);ev=loc.get('events',[])
                event_counts.update(e.event_type for e in ev[last_events:]);last_events=len(ev)
                snap.update(time_generations=loc.get('t'),active_carriers=len(ls),stored_segments=sum(len(l.segments) for l in ls),active_intervals=sum(len(l.active) for l in ls),nodes=len(loc.get('nodes',[])),edges=len(loc.get('edges',[])),logged_events=len(ev),event_counts=dict(event_counts),candidate_options=len(loc.get('options',[])))
                break
            f=f.f_back
        trace.append(snap);last_snapshot=snap
        (out/'heartbeat.json').write_text(json.dumps(snap)+'\n')
        if snap['peak_rss_bytes']>memory:raise Limit('rss_limit')
        if snap['wall']>deadline:raise Limit('time_limit')
    signal.signal(signal.SIGALRM,sample);signal.setitimer(signal.ITIMER_REAL,.25,.25)
    try:
        t=species(job['taxa'],job['ne'],job['divergence_step'],job['root_extension'])
        if job['kind']=='history':
            stage='history_generation';rng=np.random.default_rng(np.random.SeedSequence([job['seed'],job['taxa'],0]));h=simulate_frequency_history(t,Rearrangement('MSRC','suppressed_recombination','ROOT',0,job['ne'],job.get('selection',0.)),rng)
            states={taxon:int(rng.random()<h.terminal_frequency(taxon)) for taxon in t.taxa}
            summary.update(generation_seconds=time.perf_counter()-start,history_records=len(h.records),origin=h.structured_origin_age,states=states)
            stage='history_write'
            with Path(job['history_path']).open('wb') as f:pickle.dump((h,states),f,protocol=5)
            summary['history_bytes']=Path(job['history_path']).stat().st_size;summary['status']='complete'
        else:
            h=states=None;length=job['length'];region=tuple(f*length for f in job.get('region_fraction',[.3,.7]));rate=job['rate'];kind=job['kind']
            if kind=='MSRC_ARG':
                s=time.perf_counter()
                if job.get('history_mode')=='inline':
                    stage='history_generation';rng=np.random.default_rng(np.random.SeedSequence([job['seed'],job['taxa'],0]))
                    h=simulate_frequency_history(t,Rearrangement('MSRC','suppressed_recombination','ROOT',0,job['ne'],job.get('selection',0.)),rng)
                    states={taxon:int(rng.random()<h.terminal_frequency(taxon)) for taxon in t.taxa}
                    summary['history_generation_seconds']=time.perf_counter()-s
                else:
                    stage='history_load'
                    with Path(job['history_path']).open('rb') as f:h,states=pickle.load(f)
                    summary['history_load_seconds']=time.perf_counter()-s
                summary.update(history_records=len(h.records),origin=h.structured_origin_age,states=states)
            rec=Map(length,region,(rate,rate*job['kappa'],rate)) if kind=='neutral_low' else Map.constant(length,rate)
            stage='ARG_simulation';s=time.perf_counter();trace.append(dict(wall=s-start,stage=stage,time_generations=0.,active_carriers=job['taxa'],stored_segments=job['taxa'],active_intervals=job['taxa'],logged_events=0))
            if profile:profile.enable()
            if kind=='MSRC_ARG':result=simulate_msrc_arg(t,rec,h,job.get('marker_fraction',0.)*length,region,job['kappa'],'reservoir_v1',np.random.SeedSequence([job['seed'],job['taxa'],1]),states)
            else:result=simulate_arg(t,ARGConfig(length,rec),np.random.default_rng(np.random.SeedSequence([job['seed'],job['taxa'],1])))
            if profile:profile.disable()
            summary.update(simulation_seconds=time.perf_counter()-s,raw_nodes=len(result.nodes),raw_edges=len(result.edges),peak_carriers=result.peak_lineages,peak_stored_segments=result.peak_segments,event_counts=dict(collections.Counter(e.event_type for e in result.events)),metadata=result.metadata)
            stage='tree_sequence_build';s=time.perf_counter();raw=build_tree_sequence(result);summary['build_seconds']=time.perf_counter()-s
            stage='simplification';s=time.perf_counter();ts=raw.simplify();summary['simplify_seconds']=time.perf_counter()-s
            stage='validity';s=time.perf_counter();summary['validity']=validate(ts,t);summary['validity_seconds']=time.perf_counter()-s
            if h is not None:summary['background_validity']=validate_background_events(result,t,h,states)
            stage='tree_sequence_write';s=time.perf_counter();ts.dump(out/'chromosome.trees');summary['tree_file_bytes']=(out/'chromosome.trees').stat().st_size;summary['dump_seconds']=time.perf_counter()-s
            stage='event_export';s=time.perf_counter();write_event_log(result.events,out/'events.csv');summary['event_export_seconds']=time.perf_counter()-s;summary['event_log_bytes']=(out/'events.csv').stat().st_size
            stage='newick_export';s=time.perf_counter();export_genealogies(ts,out/'local_trees.nwk',out/'local_trees.csv',chrom=job.get('chrom','chr1'));summary['export_seconds']=time.perf_counter()-s
            summary.update(marginal_trees=ts.num_trees,ts_nodes=ts.num_nodes,ts_edges=ts.num_edges,runtime_per_tree=summary['simulation_seconds']/ts.num_trees,newick_bytes=(out/'local_trees.nwk').stat().st_size,coordinate_tree_bytes=(out/'local_trees.csv').stat().st_size)
            stage='diagnostics';s=time.perf_counter();rows=[];unique=set();previous=None;changes=top_changes=0
            for gt in ts.trees():
                pair_times=[gt.time(gt.mrca(i,j)) for i in range(job['taxa']) for j in range(i)]
                # Internal branch total in generations, not substitutions.
                internal=sum(gt.branch_length(u) for u in gt.nodes() if not gt.is_sample(u) and gt.parent(u)!=-1)
                signature=hashlib.sha256(np.asarray(pair_times,dtype='<f8').tobytes()).hexdigest();unique.add(signature)
                top=topology(gt,job['taxa'])
                if previous is not None:
                    changes+=int(previous[0]!=signature);top_changes+=int(previous[1]!=top)
                previous=signature,top
                rows.append(dict(left=gt.interval.left,right=gt.interval.right,topology=topology(gt,job['taxa']),root_age=gt.time(gt.root),pair_mean=float(np.mean(pair_times)),internal_generations=internal))
            with (out/'marginals.csv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
            summary['diagnostic_seconds']=time.perf_counter()-s
            summary.update(distinct_dated_genealogies=len(unique),genealogy_changes=changes,quartet_topology_changes=top_changes if job['taxa']==4 else None)
            summary['event_positions_inside']={kind:sum(e.event_type==kind and e.position is not None and region[0]<=e.position<region[1] for e in result.events) for kind in ('recombination','background_transfer')}
            summary['status']='complete'
    except Limit as e:summary.update(status=str(e),failed_stage=stage,last_snapshot=last_snapshot)
    except Exception as e:summary.update(status='error',failed_stage=stage,error=repr(e),traceback=traceback.format_exc())
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        if profile:
            profile.disable();profile.dump_stats(out/'profile.pstats')
            with (out/'profile.txt').open('w') as f:pstats.Stats(profile,stream=f).strip_dirs().sort_stats('cumulative').print_stats(35)
        summary.update(wall_seconds=time.perf_counter()-start,peak_rss_bytes=rss(),last_stage=stage)
        if summary['status']!='complete':summary.update(partial_event_counts=dict(event_counts),partial_state=last_snapshot)
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        if trace:
            fields=sorted(set().union(*(row.keys() for row in trace)))
            with (out/'ancestry_trace.csv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(trace)
    print(json.dumps({k:summary[k] for k in ('status','wall_seconds','peak_rss_bytes')}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--job',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(json.loads(a.job.read_text()),a.output)
