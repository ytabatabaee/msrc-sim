from __future__ import annotations
import argparse
import numpy as np
import yaml
from .experiments import _tree_from_config, _rearrangement_from_config
from .wright_fisher import simulate_frequency_history
from .conditioning import evaluate_conditioning, terminal_pattern
from .history_io import save_frozen_history


def main() -> None:
    ap=argparse.ArgumentParser(description="Simulate and save one accepted rearrangement history")
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", default="frozen_history.yaml")
    ap.add_argument("--max-attempts", type=int, default=None)
    args=ap.parse_args()
    with open(args.config) as h: config=yaml.safe_load(h)
    tree=_tree_from_config(config); rearr=_rearrangement_from_config(config)
    conditioning=config.get("conditioning", {"mode":"none"})
    max_attempts=args.max_attempts or int(conditioning.get("max_attempts",100000))
    master=np.random.SeedSequence(int(config.get("seed",1)))
    for attempt,ss in enumerate(master.spawn(max_attempts)):
        rng=np.random.default_rng(ss)
        hist=simulate_frequency_history(tree,rearr,rng)
        sampled={t:int(rng.random()<hist.terminal_frequency(t)) for t in tree.taxa}
        result=evaluate_conditioning(conditioning,hist,sampled,tree.taxa)
        if result.accepted:
            save_frozen_history(args.output,config,hist,sampled,{"attempt_id":attempt,"history_seed":int(ss.generate_state(1)[0]),"terminal_pattern":terminal_pattern(sampled,tree.taxa),"conditioning_reason":result.reason})
            print(f"Wrote {args.output}")
            return
    raise RuntimeError(f"No accepted history after {max_attempts} attempts")

if __name__=="__main__": main()
