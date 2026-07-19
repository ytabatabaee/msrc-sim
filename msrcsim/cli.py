from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from .config import load_config
from .structured_coalescent import simulate_conditional_quartets,simulate_genealogy
from .species_tree import BalancedQuartetTree
from .rearrangement import Rearrangement
from .wright_fisher import simulate_frequency_history,sample_terminal_arrangements
from .io import write_mechanistic_output


def _conditional(raw):
    s=raw["structured_interval"]; m=s["migration"]; c=s["coalescence"]
    r=simulate_conditional_quartets(s["configuration"],float(s["duration"]),float(m["m01"]),float(m["m10"]),
                                    float(c["lambda0"]),float(c["lambda1"]),int(raw["num_loci"]),int(raw.get("seed",1)))
    out=Path(raw.get("output_dir","conditional_output")); out.mkdir(parents=True,exist_ok=True)
    d={"counts":r.counts.tolist(),"empirical":r.empirical_probabilities.tolist(),
       "exact":r.exact_probabilities.tolist(),"absolute_error":r.absolute_error.tolist()}
    (out/"summary.json").write_text(json.dumps(d,indent=2)); print(json.dumps(d,indent=2))


def _mechanistic(raw):
    st=raw["species_tree"]
    tree=BalancedQuartetTree.create(int(st["cherry_age"]),int(st["root_age"]),int(st["origin_age"]),
                                    int(st["effective_population_size"]),int(st.get("ancestral_population_size",st["effective_population_size"])))
    rr=raw["rearrangement"]
    rearr=Rearrangement(initial_copies=int(rr.get("initial_copies",1)),selection_coefficient=float(rr.get("selection_coefficient",0)),
                         recombination_rate=float(rr["recombination_rate"]),suppression_factor=float(rr.get("suppression_factor",1)))
    seed=int(raw.get("seed",1)); rng=np.random.default_rng(seed)
    history=simulate_frequency_history(tree,rearr,rng)
    arrangements=sample_terminal_arrangements(history,tree,rng)
    genealogies=[simulate_genealogy(tree,history,arrangements,rearr.effective_recombination_rate,rng)
                 for _ in range(int(raw["num_loci"]))]
    out=Path(raw.get("output_dir","mechanistic_output")); write_mechanistic_output(out,raw,history,arrangements,genealogies)
    print(f"sampled arrangement pattern: {''.join(str(arrangements[str(i)]) for i in range(1,5))}")
    print(f"output directory: {out}")


def main():
    p=argparse.ArgumentParser(); p.add_argument("--config",required=True); a=p.parse_args()
    raw=load_config(a.config).raw
    _conditional(raw) if raw["mode"]=="conditional" else _mechanistic(raw)
if __name__=="__main__": main()
