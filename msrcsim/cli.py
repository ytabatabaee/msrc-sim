import argparse
from .config import load_config
from .io import write_result
from .structured_coalescent import simulate_conditional_quartets

def main():
    p=argparse.ArgumentParser(description="Simulate conditional quartet genealogies under the MSRC")
    p.add_argument("--config",required=True); args=p.parse_args(); c=load_config(args.config)
    r=simulate_conditional_quartets(c.configuration,c.duration,c.m01,c.m10,c.lambda0,c.lambda1,c.num_loci,c.seed)
    out=write_result(c,r)
    print(f"Configuration: {c.configuration}")
    print(f"Counts: {r.counts.tolist()}")
    print(f"Empirical: {r.empirical_probabilities.tolist()}")
    print(f"Exact: {r.exact_probabilities.tolist()}")
    print(f"Absolute error: {r.absolute_error.tolist()}")
    print(f"Output directory: {out}")
if __name__=="__main__": main()
