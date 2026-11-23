"""Escalate the same prespecified seeds, and bypass benchmark-only history serialization."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_scalability.benchmark import run_job,flatten

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--benchmarks',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--seconds',type=int,default=300);a=p.parse_args();source=a.benchmarks.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((source/'provenance.json').read_text());(out/'provenance.json').write_text(json.dumps(dict(parent=cfg,extended_seconds=a.seconds,design='same seeds and demographics; all three original 4-taxon seeds replayed, no replacement of failed chromosomes; inline 30-taxon history avoids cache serialization'),indent=2)+'\n')
    tasks=[]
    for seed in cfg['config']['seeds']:
        path=next(p for p in source.glob('MSRC_ARG_n4_L10000000_s*_r1e-08_k0.01/summary.json') if json.loads(p.read_text())['job']['seed']==seed)
        job=json.loads(path.read_text())['job'];job['seconds']=a.seconds;tasks.append((path.parent.name+'_extended',job))
    for kind in ('neutral_uniform','neutral_low','MSRC_ARG'):
        path=next(source.glob(kind+'_n30_L10000000_s101_r1e-08_k0.01/summary.json'));job=json.loads(path.read_text())['job'];job['seconds']=a.seconds
        if kind=='MSRC_ARG':job.pop('history_path',None);job['history_mode']='inline'
        tasks.append((path.parent.name+'_extended_inline' if kind=='MSRC_ARG' else path.parent.name+'_extended',job))
    results=[]
    for id,job in tasks:
        result=run_job(job,out/id);results.append((id,result));flatten(results,out)
    print('extended suite complete',flush=True)
