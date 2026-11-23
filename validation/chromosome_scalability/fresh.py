"""Measure full generation-to-export cost for the original three seeds, without a history cache."""
import argparse,hashlib,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.chromosome_scalability.benchmark import run_job,flatten

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--benchmarks',type=Path,required=True);p.add_argument('--cached-extended',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--seconds',type=int,default=300);a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True);source=a.benchmarks.resolve();cfg=json.loads((source/'provenance.json').read_text())
    (out/'provenance.json').write_text(json.dumps(dict(parent=cfg,design='all original 4-taxon seeds replayed inline; same chromosomes, not additional independent replicates',seconds=a.seconds),indent=2)+'\n');results=[];checks=[]
    for seed in cfg['config']['seeds']:
        old=next(p for p in source.glob('MSRC_ARG_n4_L10000000_s*_r1e-08_k0.01/summary.json') if json.loads(p.read_text())['job']['seed']==seed)
        j=json.loads(old.read_text())['job'];j.pop('history_path',None);j.update(history_mode='inline',seconds=a.seconds)
        id=f'MSRC_ARG_n4_L10000000_s{seed}_fresh';r=run_job(j,out/id);results.append((id,r));flatten(results,out)
        cached=next(p.parent for p in a.cached_extended.glob('MSRC_ARG_n4_L10000000_s*_r1e-08_k0.01_extended/summary.json') if json.loads(p.read_text())['job']['seed']==seed)
        if r['status']=='complete':
            hashes={file:hashlib.sha256((out/id/file).read_bytes()).hexdigest()==hashlib.sha256((cached/file).read_bytes()).hexdigest() for file in ('events.csv','local_trees.csv','local_trees.nwk','marginals.csv')}
            checks.append(dict(seed=seed,identical_scientific_files=hashes));assert all(hashes.values()),'inline history replay changed scientific data'
    (out/'cache_vs_fresh_regression.json').write_text(json.dumps(checks,indent=2)+'\n');print('fresh cohort complete',flush=True)
