"""Run sequential, isolated, resource-limited jobs with immutable checkpoints."""
import argparse,csv,hashlib,json,os,subprocess,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_scalability.worker import species

ROOT=Path(__file__).resolve().parents[2]


def run_job(job,path):
    path.mkdir(parents=True,exist_ok=True)
    if (path/'summary.json').exists():
        prior=json.loads((path/'summary.json').read_text())
        def normalized(j):
            j=dict(j);j.setdefault('region_fraction',[.3,.7]);j.setdefault('marker_fraction',0.)
            return j
        if normalized(prior['job'])!=normalized(job):raise ValueError('checkpoint configuration differs; use a fresh output: '+str(path))
        return prior
    (path/'job.json').write_text(json.dumps(job,indent=2)+'\n');start=time.perf_counter()
    with (path/'stdout.log').open('w') as f:
        try:
            p=subprocess.run([sys.executable,str(Path(__file__).with_name('worker.py')),'--job',str(path/'job.json'),'--output',str(path)],stdout=f,stderr=subprocess.STDOUT,timeout=job['seconds']+15,cwd=ROOT)
            if not (path/'summary.json').exists():
                result=dict(job=job,status='worker_exit',returncode=p.returncode,last_snapshot=json.loads((path/'heartbeat.json').read_text()) if (path/'heartbeat.json').exists() else None)
                (path/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
        except subprocess.TimeoutExpired:
            result=dict(job=job,status='parent_timeout',last_snapshot=json.loads((path/'heartbeat.json').read_text()) if (path/'heartbeat.json').exists() else None)
            (path/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    result=json.loads((path/'summary.json').read_text());result['process_wall_seconds']=time.perf_counter()-start
    (path/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(path.name,result['status'],round(result['process_wall_seconds'],2),flush=True);return result


def flatten(results,out):
    rows=[]
    for id,r in results:
        j=r['job'];events=r.get('event_counts',{});part=r.get('partial_state',r.get('last_snapshot',{})) or {}
        rows.append(dict(id=id,model=j['kind'],taxa=j['taxa'],length=j.get('length'),seed=j['seed'],ne=j['ne'],rate=j.get('rate'),kappa=j.get('kappa'),status=r['status'],wall_seconds=r.get('wall_seconds'),process_wall_seconds=r.get('process_wall_seconds'),simulation_seconds=r.get('simulation_seconds'),history_load_seconds=r.get('history_load_seconds'),peak_rss_bytes=r.get('peak_rss_bytes'),recombinations=events.get('recombination'),background_marks=events.get('background_transfer'),coalescences=events.get('coalescence'),marginal_intervals=r.get('marginal_trees'),distinct_genealogies=r.get('distinct_dated_genealogies'),nodes=r.get('ts_nodes'),edges=r.get('ts_edges'),tree_file_bytes=r.get('tree_file_bytes'),event_log_bytes=r.get('event_log_bytes'),export_seconds=r.get('export_seconds'),peak_carriers=r.get('peak_carriers'),peak_stored_segments=r.get('peak_stored_segments'),build_seconds=r.get('build_seconds'),simplify_seconds=r.get('simplify_seconds'),validity_seconds=r.get('validity_seconds'),runtime_per_tree=r.get('runtime_per_tree'),failed_stage=r.get('failed_stage'),partial_simulation_time=part.get('time_generations'),partial_logged_events=part.get('logged_events')))
    if rows:
        with (out/'metrics.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    (out/'summary.json').write_text(json.dumps([dict(id=id,**r) for id,r in results],indent=2)+'\n')


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--suite',choices=['primary','multi','replicates','sensitivity','all'],default='all');a=p.parse_args();cfg=json.loads(a.config.read_text());out=a.output.resolve();out.mkdir(parents=True,exist_ok=True);cache=out/'histories';cache.mkdir(exist_ok=True)
    source={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'msrcsim/arg/engine.py',ROOT/'msrcsim/arg/coalescence.py',ROOT/'msrcsim/arg/rates.py']}
    if (out/'provenance.json').exists():
        prior=json.loads((out/'provenance.json').read_text())
        if prior['config']!=cfg or prior['engine_source_sha256']!=source:
            raise ValueError('configuration or engine changed; use a fresh benchmark output')
    (out/'provenance.json').write_text(json.dumps(dict(config=cfg,base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),engine_source_sha256=source,python=sys.version,rng_streams='SeedSequence([replicate_seed,taxa,0]) for Moran and terminal sampling; independent stream 1 for ARG',telemetry='SIGALRM wall-clock samples; no RNG draws; baseline timings unprofiled',limits='each job in a fresh process; wall and RSS limits; parent timeout backup'),indent=2)+'\n')
    tasks=[]
    def add(n,L,seed,kind,rate=None,kappa=None):
        tasks.append(dict(kind=kind,taxa=n,length=L,seed=seed,rate=cfg['rate'] if rate is None else rate,kappa=cfg['kappa'] if kappa is None else kappa))
    models=('neutral_uniform','neutral_low','MSRC_ARG')
    if a.suite in ('primary','all'):
        for L in cfg['lengths']:
            for model in models:add(4,L,cfg['seeds'][0],model)
    if a.suite in ('multi','all'):
        for n in cfg['taxa'][1:]:
            for L in (cfg['lengths'][0],cfg['lengths'][-1]):
                for model in models:add(n,L,cfg['seeds'][0],model)
    if a.suite in ('replicates','all'):
        for seed in cfg['seeds'][1:]:
            for L in cfg['replicated_lengths']:
                for model in models:add(4,L,seed,model)
    if a.suite in ('sensitivity','all'):
        for k in (0.,.1,1.):add(4,10000000,101,'MSRC_ARG',kappa=k)
        for r in (5e-9,2e-8):
            for kind in ('neutral_uniform','MSRC_ARG'):add(4,10000000,101,kind,rate=r)
    results=[];histories={};preflight=[]
    for task in tasks:
        job=dict(ne=cfg['ne'],divergence_step=cfg['divergence_step'],root_extension=cfg['root_extension'],selection=cfg['selection'],region_fraction=cfg['region_fraction'],marker_fraction=cfg['marker_fraction'],rss_gib=cfg['rss_gib'],seconds=cfg['large_seconds'] if task['length']>=5000000 else cfg['small_seconds'],**task)
        id=f"{job['kind']}_n{job['taxa']}_L{job['length']}_s{job['seed']}_r{job['rate']}_k{job['kappa']}"
        tree=species(job['taxa'],job['ne'],job['divergence_step'],job['root_extension']);total=sum(b.older_age-b.younger_age for b in tree.branches.values())
        integrated=job['rate']*job['length']*(1-(cfg['region_fraction'][1]-cfg['region_fraction'][0])*(1-job['kappa']) if job['kind']=='neutral_low' else 1.)
        preflight.append(dict(id=id,chromosomal_rate_per_lineage=integrated,population_scaled_4Ne_map=4*job['ne']*integrated,origin_age=tree.branches['ROOT'].older_age,initial_lineage_max_structured_marks=job['taxa']*integrated*tree.branches['ROOT'].older_age,moran_expected_changes_upper_bound=total*job['ne']/2))
        if job['kind']=='MSRC_ARG':
            key=job['taxa'],job['seed'];hfile=cache/f'n{key[0]}_s{key[1]}.pkl'
            if key not in histories:
                hjob={k:v for k,v in job.items() if k not in ('length','rate','kappa')};hjob.update(kind='history',seconds=cfg['history_seconds'],history_path=str(hfile))
                histories[key]=run_job(hjob,out/f'history_n{key[0]}_s{key[1]}')
                results.append((f'history_n{key[0]}_s{key[1]}',histories[key]))
            job['history_path']=str(hfile)
            if histories[key]['status']!='complete':
                result=dict(job=job,status='blocked_history',reason=histories[key]['status']);(out/id).mkdir(exist_ok=True);(out/id/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
            else:result=run_job(job,out/id)
        else:result=run_job(job,out/id)
        results.append((id,result));flatten(results,out)
        (out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')
    print('suite complete',len(results),'jobs including history generation',flush=True)

if __name__=='__main__':main()
