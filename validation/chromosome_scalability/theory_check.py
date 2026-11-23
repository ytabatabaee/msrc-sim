"""Attempt the unmodified finite-M multispecies theory at benchmark parameters.

A failed numerical dependency is reported, never replaced by a constant-p
prediction for the stochastic Moran model. Neutral MSC limits are separate.
"""
import argparse,json,resource,signal,sys,time,traceback
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.moran_spatial.theory import load_theory,predict,rho_effective,make_tree,PrescribedHistory

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--theory-source',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--seconds',type=int,default=120);a=p.parse_args();cfg=json.loads(a.config.read_text());mods,hashes=load_theory(a.theory_source);start=time.perf_counter();M=2*cfg['ne'];step=cfg['divergence_step'];L=max(cfg['lengths']);pi=np.zeros(M+1);pi[cfg['ne']]=1
    case=dict(shape='balanced',newick=f'((A:{step},B:{step})AB:{step},(C:{step},D:{step})CD:{step})ROOT;',M=M,root_duration=cfg['root_extension'],selection=cfg['selection'],pi_old=pi.tolist(),marker=cfg['marker_fraction']*L,suppressed_region=[v*L for v in cfg['region_fraction']],kappa=cfg['kappa'],map=dict(length=L,breakpoints=[],rates=[cfg['rate']]))
    result=dict(case=case,theory_source_sha256=hashes,positions=[],status='running',method='original full finite-M Moran recursion, unmodified',numerical_scope='same Ne, selection, branch durations, binomial inheritance and root count law as MSRC-ARG')
    # A cheap exact neutral control, explicitly not the stochastic Moran target.
    t=make_tree(case);h=PrescribedHistory(t,{bid:[(0.,b.older_age-b.younger_age)] for bid,b in t.branches.items()});neutral=predict(mods,case,0.,h);result['neutral_MSC_quartets']=neutral.quartet_probabilities.tolist()
    def alarm(sig,frame):raise TimeoutError('exact-theory time budget exceeded')
    signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,a.seconds)
    try:
        result['stage']='binomial_inheritance_kernel_preflight'
        H=mods['full_species_tree'].inheritance_kernel(M,M)
        assert np.allclose(H.sum(axis=1),1,atol=1e-10);del H
        for x in [v*L for v in (.05,.29,.31,.5,.69,.71,.95)]:
            result['stage']=f'full_Moran_at_{x}';rho=rho_effective(case,x);theory=predict(mods,case,rho)
            assert abs(theory.quartet_probabilities.sum()-1)<1e-9
            result['positions'].append(dict(position=x,rho_effective=rho,quartets=theory.quartet_probabilities.tolist()))
        result['status']='complete'
    except Exception as e:
        result.update(status='unavailable',error=repr(e),traceback=traceback.format_exc(),interpretation='the requested exact target was not evaluated; no approximation substituted')
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);result['seconds']=time.perf_counter()-start;result['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
