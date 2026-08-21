from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from .history_io import load_frozen_history
from .history_summary import summarize_frequency_history
from .experiments import _tree_from_config,_recombination
from .structured_coalescent import simulate_genealogy
from .model_fitting import compare_models


def main() -> None:
    ap=argparse.ArgumentParser(description="Replay genealogies conditional on a frozen rearrangement history")
    ap.add_argument("--history", required=True)
    ap.add_argument("--num-loci", type=int, required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--output", default="replay_output")
    args=ap.parse_args()
    config,hist,sampled,metadata=load_frozen_history(args.history)
    tree=_tree_from_config(config); base,frac=_recombination(config); rng=np.random.default_rng(args.seed)
    results=[simulate_genealogy(i,tree,hist,sampled,base,frac,rng,record_events=False) for i in range(args.num_loci)]
    counts=np.bincount([r.topology_index for r in results],minlength=3); q=counts/args.num_loci
    model=compare_models(counts)
    out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    summary={"history_file":str(args.history),"history_metadata":metadata,"genealogy_seed":args.seed,"num_loci":args.num_loci,"n1":int(counts[0]),"n2":int(counts[1]),"n3":int(counts[2]),"q1":float(q[0]),"q2":float(q[1]),"q3":float(q[2]),**model}
    (out/'replay_summary.json').write_text(json.dumps(summary,indent=2))
    rows=summarize_frequency_history(hist)
    with (out/'branch_history_summary.csv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with (out/'true_gene_trees.nwk').open('w') as h:
        for r in results: h.write(r.newick+'\n')
    print(f"Wrote {out}")

if __name__=="__main__": main()
