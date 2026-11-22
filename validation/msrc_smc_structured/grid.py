"""Run independently reproducible p-specific blocks of the primary grid."""
import argparse
from pathlib import Path
from validate import run

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);a.add_argument('--p',type=float,required=True);a.add_argument('--replicates',type=int,default=1000)
    cfg=a.parse_args()
    for nt in (2,4):
        for k in (0.,.01,.1,.5,1.):
            run(cfg.output/f'{nt}tax_p{cfg.p}_k{k}',cfg.replicates,nt,cfg.p,k,20261121)
