"""Independent full multispecies prescribed-frequency predictions, before ARG runs."""
import sys,json,csv,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from validation.moran_spatial.theory import load_theory,make_tree,PrescribedHistory,predict,rho_effective

def case(duration=1000,kappa=.01,marker=5e6,rate=1e-8):
 return dict(M=2000,shape='balanced',newick='((A:100,B:100)AB:100,(C:100,D:100)CD:100)ROOT;',root_duration=duration,marker=marker,suppressed_region=[3e6,7e6],kappa=kappa,map={'breakpoints':[],'rates':[rate]})
if __name__=='__main__':
 start=time.perf_counter();mods,hashes=load_theory('/Users/ytabatabaee/Desktop/msrc-paper');rows=[]
 for duration in (50,200,500,1000,2000):
  for k in (0.,.01,.1,1.):
   for marker in (0.,5e6):
    c=case(duration,k,marker);t=make_tree(c);h=PrescribedHistory(t,{bid:[(.5,b.older_age-b.younger_age)] for bid,b in t.branches.items()})
    for x in (0.,1e6,2.9e6,3.1e6,5e6,6.9e6,7.1e6,9e6,1e7):
     result=predict(mods,c,rho_effective(c,x),h)
     for pattern in ('1010','0101','0111','0100','0000'):
      q=result.conditional_quartet_probabilities[int(pattern,2)]
      rows.append(dict(duration=duration,kappa=k,marker=marker,position=x,pattern=pattern,rho=rho_effective(c,x),q1=q[0],q2=q[1],q3=q[2]))
 out=Path(__file__).parent/'results';out.mkdir(exist_ok=True)
 with (out/'theory_grid.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 (out/'theory_provenance.json').write_text(json.dumps(dict(hashes=hashes,seconds=time.perf_counter()-start,model='full balanced multispecies prescribed p=0.5; no Moran approximation'),indent=2))
 print('theory complete',len(rows),flush=True)
