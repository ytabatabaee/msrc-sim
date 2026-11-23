"""Prespecified sensitivity replays; three seeds per setting are descriptive."""
import json,subprocess,sys
from pathlib import Path
root=Path(__file__).parent
changes=[dict(duration=200),dict(duration=2000),dict(kappa=0.),dict(kappa=.1),dict(marker_fraction=0.),dict(pattern='0101'),dict(rate=2e-8)]
for i,change in enumerate(changes):
 for seed in (100,101,102):
  j=dict(length=10000000,model='C',kappa=.01,duration=1000,pattern='1010',seed=seed,**{})
  j.update(change);out=root/'results/sensitivity'/f'case{i}_s{seed}'
  if not (out/'summary.json').exists():subprocess.run([sys.executable,str(root/'run.py'),json.dumps(j),str(out)],check=True,timeout=130)
