from __future__ import annotations
from pathlib import Path
import csv,json,yaml
from dataclasses import asdict

def write_outputs(config,tree,rearrangement,history,sampled,results):
    out=Path(config['output']['directory']); out.mkdir(parents=True,exist_ok=True)
    if config['output']['record_resolved_config']:
        with open(out/'config.resolved.yaml','w') as f: yaml.safe_dump(config,f,sort_keys=False)
    if config['output']['record_frequency_history']:
        rows=[asdict(r) for r in history.records]
        with open(out/'frequency_history.csv','w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    if config['output']['record_sampled_arrangements']:
        with open(out/'sampled_arrangements.csv','w',newline='') as f:
            w=csv.writer(f); w.writerow(['taxon','arrangement','terminal_frequency_A1']);
            for t,z in sampled.items(): w.writerow([t,z,history.terminal_frequency(t)])
    if config['output']['record_gene_trees']:
        with open(out/'true_gene_trees.nwk','w') as f:
            for r in results: f.write(r.newick+'\n')
    with open(out/'coalescence_times.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['locus_id','coalescence_index','time'])
        for r in results:
            for i,t in enumerate(r.coalescence_times): w.writerow([r.locus_id,i,t])
    with open(out/'coalescence_events.csv','w',newline='') as f:
        rows=[asdict(x) for r in results for x in r.coalescences]
        if rows: w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    if config['output']['record_backward_events']:
        rows=[asdict(x) for r in results for x in r.events]
        if rows:
            with open(out/'genealogy_events.csv','w',newline='') as f: w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    counts=[sum(r.topology_index==i for r in results) for i in range(3)]
    with open(out/'summary.json','w') as f: json.dump({'version':'0.3.0','taxa':tree.taxa,'sampled_arrangements':sampled,'topology_counts':counts,'topology_frequencies':[c/len(results) for c in counts]},f,indent=2)
    return out
