import argparse, yaml
import numpy as np
from ..species_tree import SpeciesTree
from .engine import SMCConfig, simulate_msc_smc, write_result
from .recombination_map import PiecewiseConstantMap

def main(argv=None):
    p=argparse.ArgumentParser(description="Simulate linked neutral MSC genealogies under multispecies SMC'")
    p.add_argument("--config",required=True); p.add_argument("--output",default=None)
    a=p.parse_args(argv)
    cfg=yaml.safe_load(open(a.config))
    t=cfg["species_tree"]
    tree=SpeciesTree(t["newick"],t.get("default_ne",t.get("default_effective_population_size",10000)),t.get("root_extension",1e7),t.get("branches"))
    c=cfg["smc"]
    m=c["recombination_map"]
    rec=PiecewiseConstantMap(float(c["chromosome_length"]),tuple(map(float,m.get("breakpoints",[]))),tuple(map(float,m["rates"])))
    sampling=cfg.get("sampling",{"samples_per_species":1,"time_generations":0.0})
    result=simulate_msc_smc(tree,SMCConfig(float(c["chromosome_length"]),rec,int(c.get("seed",cfg.get("seed",1))),sampling),np.random.default_rng(int(c.get("seed",cfg.get("seed",1)))))
    result.metadata["species_tree_newick"]=t["newick"]
    write_result(result,a.output or cfg.get("output",{}).get("directory","smc_output"))
    return 0
if __name__=="__main__": raise SystemExit(main())
