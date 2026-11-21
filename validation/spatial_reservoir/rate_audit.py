"""Independent focal mark-rate and missing-marker audit, not tree validation."""
import csv,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from msrcsim.spatial_reservoir import ReservoirBackground
from msrcsim.arg.backgrounds.msrc_two_state import MSRCTwoBackground
from msrcsim.arg.recombination_map import PiecewiseConstantRecombinationMap as Map
from msrcsim.arg.segments import Segment
from validate import Frozen


def run(out):
    rec=Map.constant(100,.01);active=((80.,90.),);segments=(Segment(80,90,0),)
    rows=[];rng=np.random.default_rng(11991);n=12000
    for p in (.2,.5,.8):
        for h in (0,1):
            for k in (0.,.1,1.):
                model=ReservoirBackground(Frozen(p),0,(0,100),k)
                rate=model.recombination_hazard('ROOT',h,active,segments,0,rec)
                switches=0
                for _ in range(n):
                    e=model.sample_recombination('ROOT',h,active,segments,0,rec,rng)
                    switches += e['position']<85 and e['parent_background']!=h
                q=rate*switches/n;se=rate*np.sqrt((switches/n)*(1-switches/n)/n)
                cross=p if h==0 else 1-p;target=.85*k*cross
                rows.append(dict(p=p,background=h,kappa=k,proposals=n,observed_switch_rate=q,target_switch_rate=target,se=se,difference=q-target,ci_low=q-target-1.96*se,ci_high=q-target+1.96*se,current_ARG_missing_marker_rate=.05*k*cross))
    with Path(out).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


if __name__=='__main__':run(sys.argv[1])
