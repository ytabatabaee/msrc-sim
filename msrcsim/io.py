from __future__ import annotations
import csv,json
from pathlib import Path
import yaml
from .analytic import TOPOLOGY_NAMES
from .wright_fisher import FrequencyHistory


def write_frequency_history(path: Path,history: FrequencyHistory):
    with path.open("w",newline="") as h:
        w=csv.writer(h); w.writerow(["branch","forward_generation","age","count","frequency"])
        for bid,b in history.branches.items():
            for g,k in enumerate(b.counts_forward):
                age=b.older_age-g
                w.writerow([bid,g,age,int(k),float(k/(2*b.effective_size))])

def write_mechanistic_output(outdir: Path,config: dict,history: FrequencyHistory,
                             arrangements: dict,genealogies):
    outdir.mkdir(parents=True,exist_ok=True)
    with (outdir/"config.resolved.yaml").open("w") as h: yaml.safe_dump(config,h,sort_keys=False)
    write_frequency_history(outdir/"frequency_history.csv",history)
    with (outdir/"sampled_arrangements.csv").open("w",newline="") as h:
        w=csv.writer(h); w.writerow(["taxon","terminal_frequency","sampled_arrangement"])
        for t in sorted(arrangements): w.writerow([t,history.terminal_frequencies[t],arrangements[t]])
    with (outdir/"true_gene_trees.nwk").open("w") as h:
        for g in genealogies: h.write(g.newick+"\n")
    counts=[0,0,0]
    for g in genealogies: counts[g.topology_index]+=1
    n=len(genealogies)
    with (outdir/"quartet_counts.csv").open("w",newline="") as h:
        w=csv.writer(h); w.writerow(["topology","count","frequency"])
        for i,name in enumerate(TOPOLOGY_NAMES): w.writerow([name,counts[i],counts[i]/n])
    summary={"num_loci":n,"counts":counts,"frequencies":[x/n for x in counts],
             "terminal_frequencies":history.terminal_frequencies,"sampled_arrangements":arrangements}
    with (outdir/"summary.json").open("w") as h: json.dump(summary,h,indent=2)
