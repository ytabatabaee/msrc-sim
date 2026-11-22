"""Compact joint/conditional tables, with failures preventing survivor-only claims."""
import argparse,csv,json
from pathlib import Path
from collections import Counter
import numpy as np


def aggregate(root):
    summaries=[];joint=[];allrows=[]
    for dest in sorted(p for p in root.iterdir() if p.is_dir() and (p/'draws.csv').exists()):
        summary=json.loads((dest/'summary.json').read_text());summaries.append(summary)
        case=json.loads((dest/'case.json').read_text());unconditional=case['sampling']=='unconditional'
        with (dest/'comparisons.csv').open() as f:allrows.extend(csv.DictReader(f))
        counts=Counter()
        with (dest/'draws.csv').open() as f:
            for row in csv.DictReader(f):counts[row['method'],float(row['position']),row['pattern'],int(row['topology'])]+=1
        theory={float(t['position']):t for t in json.loads((dest/'theory.json').read_text())}
        for method,x in sorted({(m,x) for m,x,z,q in counts}):
            failed=summary['failures'][method]
            for z in range(16):
                pattern=f'{z:04b}';n=sum(counts[method,x,pattern,q] for q in range(3));J=theory[x]['joint'][z]
                for q in range(3):
                    count=counts[method,x,pattern,q]
                    joint.append(dict(case=dest.name,method=method,position=x,pattern=pattern,topology=q+1,attempted=summary['n'],pattern_n=n,failed_chromosomes=failed,joint_theory=J[q],joint_empirical=count/summary['n'] if not failed and unconditional else None,pattern_theory=sum(J),pattern_empirical=n/summary['n'] if not failed and unconditional else None,conditional_theory=J[q]/sum(J) if sum(J) else None,conditional_empirical=count/n if n and not failed else None,joint_mc_se=float(np.sqrt((count/summary['n'])*(1-count/summary['n'])/summary['n'])) if not failed and unconditional else None,conditional_mc_se=float(np.sqrt((count/n)*(1-count/n)/n)) if n and not failed else None,sampling=case['sampling']))
    for name,rows in [('all_comparisons.csv',allrows),('joint_pattern_comparisons.csv',joint)]:
        with (root/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    (root/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    verdict={}
    for method in ('ARG','SMC','SMC_frozen_target','single_locus'):
        rows=[r for r in allrows if r['method']==method and r['pattern']=='all']
        verdict[method]=dict(Counter(r['status'] for r in rows))
    (root/'verdict_counts.json').write_text(json.dumps(verdict,indent=2)+'\n');print(json.dumps(verdict,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);aggregate(p.parse_args().output)
