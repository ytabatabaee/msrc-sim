from __future__ import annotations
from pathlib import Path
import yaml

def load_config(path):
    with open(path) as f: d=yaml.safe_load(f)
    if d.get('mode')!='mechanistic': raise ValueError("This release config supports mode: mechanistic")
    d.setdefault('seed',1); d.setdefault('num_loci',1000)
    d.setdefault('output',{}); d['output'].setdefault('directory','mechanistic_output')
    d['output'].setdefault('record_resolved_config',True); d['output'].setdefault('record_frequency_history',True)
    d['output'].setdefault('record_sampled_arrangements',True); d['output'].setdefault('record_gene_trees',True)
    d['output'].setdefault('record_coalescence_times',True); d['output'].setdefault('record_backward_events',False)
    d['output'].setdefault('event_log_loci',{'first_n':100})
    if d.get('sampling',{}).get('samples_per_species',1)!=1: raise ValueError("Only one sample per species is supported")
    return d
