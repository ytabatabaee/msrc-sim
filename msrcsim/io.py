import csv,json
from .analytic import TOPOLOGY_NAMES

def write_result(config,result):
    out=config.output_dir; out.mkdir(parents=True,exist_ok=True)
    summary={"mode":config.mode,"seed":config.seed,"num_loci":config.num_loci,"configuration":"".join(map(str,result.configuration)),"duration":config.duration,"m01":config.m01,"m10":config.m10,"lambda0":config.lambda0,"lambda1":config.lambda1,"topologies":list(TOPOLOGY_NAMES),"counts":result.counts.tolist(),"empirical_probabilities":result.empirical_probabilities.tolist(),"exact_probabilities":result.exact_probabilities.tolist(),"absolute_error":result.absolute_error.tolist()}
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    with (out/'quartet_probabilities.csv').open('w',newline='') as h:
        w=csv.writer(h); w.writerow(["topology","count","empirical_probability","exact_probability","absolute_error"])
        for i,t in enumerate(TOPOLOGY_NAMES): w.writerow([t,int(result.counts[i]),float(result.empirical_probabilities[i]),float(result.exact_probabilities[i]),float(result.absolute_error[i])])
    return out
