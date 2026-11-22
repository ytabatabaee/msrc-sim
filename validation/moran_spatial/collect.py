"""Collect completed case checkpoints without pooling their chromosome counts."""
import argparse,json,shutil,csv,math
from collections import Counter
from pathlib import Path


def collect(inputs,out):
    out.mkdir(parents=True,exist_ok=True)
    for i,source in enumerate(inputs):
        provenance=json.loads((source/'provenance.json').read_text())
        (out/('provenance.json' if i==0 else f'source_{i}_provenance.json')).write_text(json.dumps(provenance,indent=2)+'\n')
        cases={c['id']:c for c in provenance['configuration']['cases']}
        for path in sorted(source.iterdir()):
            if not path.is_dir() or not (path/'summary.json').exists():continue
            dest=out/path.name
            if dest.exists():raise FileExistsError('use a fresh collection output directory: '+str(dest))
            shutil.copytree(path,dest)
            (dest/'case.json').write_text(json.dumps(cases[path.name],indent=2)+'\n')
            # Early checkpoints used an uncertainty label for undefined
            # zero-probability patterns; normalize the label, not any criterion.
            with (dest/'comparisons.csv').open() as f:rows=list(csv.DictReader(f))
            for row in rows:
                if int(row['n'])==0 and not math.isfinite(float(row['theory'])):
                    row['status']='NOT_DEFINED_zero_probability_pattern'
            with (dest/'comparisons.csv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
            summary=json.loads((dest/'summary.json').read_text());summary['statuses']=dict(Counter(r['status'] for r in rows));(dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();collect(a.inputs,a.output)
