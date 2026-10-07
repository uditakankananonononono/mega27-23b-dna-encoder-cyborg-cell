"""Perfectly shared substitutions across copies at unchanged emitted budget."""
import hashlib,json,sys
from collections import Counter
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import bound_nonce_file_audit as F

def compute():
 raw=(ROOT/'experiments/bound_nonce_shared_noise_plan.json').read_bytes();p=json.loads(raw);old=json.loads((ROOT/'results/bound_nonce_file_audit.json').read_text());rng=np.random.default_rng(74119);messages=[rng.bytes(n) for n in [128,2048] for _ in range(2)];rows=[];wrong=[]
 lookup={(r['message'],r['stripes'],r['rate'],r['replicate']):r for r in old['trials']}
 for mi,msg in enumerate(messages):
  for stripes in [1,4]:
   seq,ledger=F.G.encode_geometry(msg,stripes)
   if seq is None:continue
   for rate in p['rates']:
    for rep in range(p['replicates']):
     seed=p['seed_base']+mi*1000+rep*10;noisy,mut=F.strict(seq,rate,seed);prior=lookup[(mi,stripes,rate,rep)];digest=hashlib.sha256(noisy.encode()).hexdigest();assert prior['noisy_sha256'][0]==digest
     result,error,out=F.outcome([noisy]*3,stripes,msg);rows.append({'message':mi,'bytes':len(msg),'stripes':stripes,'rate':rate,'replicate':rep,'seed':seed,'total_nt':3*len(seq),'changed_nt_over_copies':3*mut,'shared_noisy_sha256':digest,'independent_outcome':prior['outcome'],'shared_outcome':result,'error':error})
     if result=='silent_wrong':wrong.append({'row':rows[-1],'expected_hex':msg.hex(),'actual_hex':out.hex(),'noisy_strand':noisy})
 summaries=[]
 for size in [128,2048]:
  for stripes in [1,4]:
   for rate in p['rates']:
    rr=[r for r in rows if (r['bytes'],r['stripes'],r['rate'])==(size,stripes,rate)]
    summaries.append({'bytes':size,'stripes':stripes,'rate':rate,'n':len(rr),'independent_counts':dict(Counter(r['independent_outcome'] for r in rr)),'shared_counts':dict(Counter(r['shared_outcome'] for r in rr))})
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'prior_ledger_sha256':hashlib.sha256((ROOT/'results/bound_nonce_file_audit.json').read_bytes()).hexdigest(),'rows':rows,'summaries':summaries,'silent_wrong':wrong,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/bound_nonce_shared_noise_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j['summaries'],indent=2))
