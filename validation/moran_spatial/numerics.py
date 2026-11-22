"""Numerical cross-checks on original operators, including full-tree recursions."""
import argparse
from contextlib import contextmanager
import json
from pathlib import Path
import sys
import time
import tracemalloc
import resource
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.moran_spatial.theory import load_theory,predict
from validation.moran_spatial.run import prescribed,make_tree


@contextmanager
def dense_exponentials(modules):
    saved={}
    for name,m in modules.items():
        if hasattr(m,'expm_multiply'):
            saved[name]=m.expm_multiply
            m.expm_multiply=lambda A,B,**kw:expm(A.toarray() if hasattr(A,'toarray') else np.asarray(A))@B
    try:yield
    finally:
        for name,f in saved.items():modules[name].expm_multiply=f


def check(mods,cases):
    base=mods['moran_quartet_operator'];rows=[]
    for shape in ('balanced','pectinate'):
        c=next(c.copy() for c in cases if c['shape']==shape)
        c.update(M=4,branch_M={},root_duration=.7,pi_old=[0,0,1,0,0],selection=.2)
        start=time.perf_counter();sparse=predict(mods,c,.17)
        with dense_exponentials(mods):dense=predict(mods,c,.17)
        error=float(np.max(abs(sparse.joint-dense.joint)))
        assert error<1e-10 and abs(sparse.joint.sum()-1)<1e-10
        rows.append(dict(check='full_recursion_sparse_vs_dense',shape=shape,M=4,max_error=error,seconds=time.perf_counter()-start))
    for M in (8,20,80):
        start=time.perf_counter();tracemalloc.start();pi=np.zeros(M+1);pi[M//2]=1
        a=base.integrated_quartet_probabilities(M,.17,.7,pi,selection=.2)
        b=mods['matrix_free_quartet'].matrix_free_shared_population_quartet_probabilities(M,.17,.7,pi,selection=.2)
        _,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
        error=float(np.max(abs(a.joint-b.joint)));assert error<1e-10
        rows.append(dict(check='shared_sparse_vs_matrix_free',M=M,dimension=19*(M+1),max_error=error,normalization_error=float(abs(a.joint.sum()-1)),seconds=time.perf_counter()-start,peak_python_traced_bytes=peak,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)))
    for n,d in ((1,2),(2,6),(3,11),(4,19)):
        A=base.joint_operator_n(8,.17,n,.2);v=np.linspace(0,1,A.shape[0])
        a=expm_multiply(A*.7,v);b=expm_multiply(A*.35,expm_multiply(A*.35,v))
        assert np.max(abs(a-b))<1e-10
        for k in (0,3,8):assert np.max(abs(base.genealogy_generator(n,k,8,.17).sum(axis=1)))<1e-12
        rows.append(dict(check='semigroup_and_generator',n=n,dimension=d*9,max_error=float(np.max(abs(a-b)))))
    c=next(c.copy() for c in cases if c['shape']=='balanced');c.update(p=.5,branch_M={})
    h=prescribed(c,make_tree(c));a=predict(mods,c,.017,h)
    b=mods['full_species_tree'].fixed_p_half_species_tree_quartet_probabilities(M=c['M'],rho=.017,root_length=c['root_duration'],sister_lengths=(2.,2.),terminal_lengths=(2.,2.,2.,2.))
    assert np.max(abs(a.quartet_probabilities-b))<1e-10
    for pattern in ('1010','0101'):
        b=mods['full_species_tree'].fixed_p_half_species_tree_conditional(pattern,M=c['M'],rho=.017,root_length=c['root_duration'],sister_lengths=(2.,2.),terminal_lengths=(2.,2.,2.,2.))
        assert np.max(abs(a.conditional_quartet_probabilities[int(pattern,2)]-b))<1e-10
    rows.append(dict(check='full_species_fixed_half_regression',max_error=float(abs(a.joint.sum()-1))))
    return rows

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--theory-source',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();mods,_=load_theory(a.theory_source);rows=check(mods,json.loads(a.config.read_text())['cases']);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
