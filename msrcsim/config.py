from dataclasses import dataclass
from pathlib import Path
import yaml

@dataclass(frozen=True)
class ConditionalConfig:
    mode:str; seed:int; num_loci:int; configuration:str; duration:float
    m01:float; m10:float; lambda0:float; lambda1:float; output_dir:Path

def load_config(path):
    with Path(path).open() as h: data=yaml.safe_load(h)
    if data.get("mode")!="conditional": raise ValueError("Milestone 1 supports only mode='conditional'")
    x=data["structured_interval"]; mig=x["migration"]; coal=x["coalescence"]
    return ConditionalConfig("conditional",int(data.get("seed",1)),int(data["num_loci"]),str(x["configuration"]),float(x["duration"]),float(mig["m01"]),float(mig["m10"]),float(coal["lambda0"]),float(coal["lambda1"]),Path(data.get("output_dir","msrc_output")))
