"""Per-position shared-generator mixture, not a physical correlation estimate."""
import hashlib,json,sys
from collections import Counter
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import bound_nonce_file_audit as F

def mix(independent,gate,q):
 if len(independent)!=3 or len({len(s) for s in independent})!=1 or len(gate)!=len(independent[0]):raise ValueError('mixture geometry')
 if not 0<=q<=1:raise ValueError('mixture q')
 return [''.join(independent[0][i] if gate[i]<q else s[i] for i in range(len(s))) for s in independent]

def compute():
 raw=(ROOT/'experiments/bound_nonce_mixture_plan.json').read_bytes();p=json.loads(raw);rng=np.random.default_rng(p['message_seed']);messages=[rng.bytes(n) for n in p['sizes'] for _ in range(p['identities_per_size'])]
 old=json.loads((ROOT/'results/bound_nonce_file_audit.json').read_text());lookup={(r['message'],r['stripes'],r['rate'],r['replicate']):r for r in old['trials']};shared=json.loads((ROOT/'results/bound_nonce_shared_noise_audit.json').read_text());shared_lookup={(r['message'],r['stripes'],r['rate'],r['replicate']):r for r in shared['rows']};rows=[];wrong=[]
 for mi,msg in enumerate(messages):
  for stripes in p['stripes']:
   seq,ledger=F.G.encode_geometry(msg,stripes)
   if seq is None:continue
   for rate in p['rates']:
    for rep in range(p['replicates']):
     seeds=[p['noise_seed_base']+mi*1000+rep*10+c for c in range(3)];independent=[F.strict(seq,rate,seed)[0] for seed in seeds];prior=lookup[(mi,stripes,rate,rep)];assert [hashlib.sha256(s.encode()).hexdigest() for s in independent]==prior['noisy_sha256']
     gate_seed=p['mixture_seed_base']+mi*1000+rep*10;gate=np.random.default_rng(gate_seed).random(len(seq))
     for q in p['shared_fractions']:
      reads=mix(independent,gate,q);result,error,out=F.outcome(reads,stripes,msg)
      if q==0:assert result==prior['outcome']
      if q==1:assert result==shared_lookup[(mi,stripes,rate,rep)]['shared_outcome']
      rows.append({'message':mi,'bytes':len(msg),'stripes':stripes,'rate':rate,'replicate':rep,'q':q,'noise_seeds':seeds,'gate_seed':gate_seed,'shared_positions':int(np.sum(gate<q)),'strand_nt':len(seq),'total_nt':3*len(seq),'changes_per_copy':[sum(a!=b for a,b in zip(seq,s)) for s in reads],'noisy_sha256':[hashlib.sha256(s.encode()).hexdigest() for s in reads],'source_sha256':hashlib.sha256(msg).hexdigest(),'clean_sha256':hashlib.sha256(seq.encode()).hexdigest(),'outcome':result,'error':error})
      if result=='silent_wrong':wrong.append({'trial':rows[-1],'expected_hex':msg.hex(),'actual_hex':out.hex(),'strands':reads})
   print('finished',mi,stripes,flush=True)
 summaries=[]
 for size in p['sizes']:
  for stripes in p['stripes']:
   for rate in p['rates']:
    for q in p['shared_fractions']:
     rr=[r for r in rows if (r['bytes'],r['stripes'],r['rate'],r['q'])==(size,stripes,rate,q)];summaries.append({'bytes':size,'stripes':stripes,'rate':rate,'q':q,'n':len(rr),'counts':dict(Counter(r['outcome'] for r in rr))})
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'prior_sha256':hashlib.sha256((ROOT/'results/bound_nonce_file_audit.json').read_bytes()).hexdigest(),'rows':rows,'summaries':summaries,'silent_wrong':wrong,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/bound_nonce_mixture_audit.json').write_text(json.dumps(j,indent=2)+'\n')
 print({q:dict(Counter(r['outcome'] for r in j['rows'] if r['q']==q)) for q in [0,.25,.5,.75,1]})
