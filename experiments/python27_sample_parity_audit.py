"""Mirror CPython2.7 float sample semantics, compare actual2.7 archived vectors."""
import json,random,math,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sample27(r,n,k):
 result=[];size=21
 if k>5:size+=4**math.ceil(math.log(k*3,4))
 if n<=size:
  pool=list(range(n))
  for i in range(k):
   j=int(r.random()*(n-i));result.append(pool[j]);pool[j]=pool[n-i-1]
 else:
  seen=set()
  for i in range(k):
   j=int(r.random()*n)
   while j in seen:j=int(r.random()*n)
   seen.add(j);result.append(j)
 return result
if __name__=='__main__':
 import argparse
 a=argparse.ArgumentParser();a.add_argument('--make-input',action='store_true');v=a.parse_args();K=67088;delta=.001;c=.025;S=c*math.log(K/delta)*math.sqrt(K);pivot=int(math.floor(K/S));tau=[S/K/d for d in range(1,pivot)]+[S/K*math.log(S/delta)]+[0]*(K-pivot);rho=[1./K]+[1./(d*(d-1)) for d in range(2,K+1)];Z=sum(rho)+sum(tau);mu=[(rho[d]+tau[d])/Z for d in range(K)];cdf=[];acc=0
 for u in mu:acc+=u;cdf.append(acc)
 inp={'K':K,'cdf':cdf,'seeds':sorted({r['seed'] for r in json.loads((ROOT/'results/dna_fountain_receiver_screen_audit.json').read_text())['records']})}
 if v.make_input:(ROOT/'results/python27_sample_input.json').write_text(json.dumps(inp,separators=(',',':')));raise SystemExit
 vectors=json.loads((ROOT/'results/python27_sample_vectors.json').read_text());errors=[];modern_diff=0
 for x in vectors['vectors']:
  r=random.Random(x['seed']);u=r.random();d=next(i+1 for i,c in enumerate(cdf) if c>u);ids=sample27(r,K,d)
  if u!=x['first_uniform'] or d!=x['degree'] or ids!=x['indices']:errors.append(x['seed'])
  r=random.Random(x['seed']);r.random();modern_diff+=r.sample(range(K),d)!=x['indices']
 j={'plan_sha256':hashlib.sha256((ROOT/'research/python27_sample_parity_plan.json').read_bytes()).hexdigest(),'input_sha256':hashlib.sha256((ROOT/'results/python27_sample_input.json').read_bytes()).hexdigest(),'vectors_sha256':hashlib.sha256((ROOT/'results/python27_sample_vectors.json').read_bytes()).hexdigest(),'vector_count':len(vectors['vectors']),'mismatch_seeds':errors,'modern_python3_sample_different_vector_count':modern_diff,'python27_version':vectors['python_version'],'scope':'Compiled2.7.18nativebranch only,not originalruntime/npflag/payloadidentity;vectorsgeneratednow notoriginalauthororacle'}
 (ROOT/'results/python27_sample_parity_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j)
