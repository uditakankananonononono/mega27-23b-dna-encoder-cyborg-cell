"""Audit-local whole-file integration. Production encoder is unchanged."""
import hashlib,json,sys
from collections import Counter
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'))
import bound_nonce_geometry_audit as G
from strict_substitution_benchmark import strict
B=G.B;E=B.E

def decode(strands,stripes):
 if isinstance(stripes,bool) or not isinstance(stripes,int) or stripes<1:raise ValueError('stripes')
 if not isinstance(strands,(list,tuple)) or not strands or any(not isinstance(s,str) for s in strands):raise ValueError('geometry')
 if len(set(map(len,strands)))!=1 or not strands[0] or len(strands[0])%64 or len(strands[0])//64<=stripes:raise ValueError('geometry')
 if any(set(s)-set(E.BASES) for s in strands):raise ValueError('alphabet')
 voted=''.join(max(E.BASES,key=lambda b:sum(s[i]==b for s in strands)) for i in range(len(strands[0])))
 blocks=[]
 for start in range(0,len(voted),64):
  candidates=set()
  for source in [voted]+list(strands):
   try:candidates.add(tuple(B.decode(source[start:start+64])))
   except ValueError:pass
  if len(candidates)>1:raise ValueError('valid copy conflict')
  blocks.append(list(next(iter(candidates))) if candidates else None)
 data_count=len(blocks)-stripes
 for g in range(stripes):
  ids=list(range(g,data_count,stripes));bad=[i for i in ids if blocks[i] is None];parity=blocks[data_count+g]
  if len(bad)==1 and parity is not None:blocks[bad[0]]=[(parity[k]-sum(blocks[i][k] for i in ids if blocks[i] is not None))%3 for k in range(32)]
  elif bad:raise ValueError('unrecoverable stripe')
  if parity is not None and any(sum(blocks[i][k] for i in ids)%3!=parity[k] for k in range(32)):raise ValueError('stripe parity mismatch')
 padded=sum(blocks[:data_count],[])
 t=E.descramble(padded)
 if len(t)<24:raise ValueError('header')
 n=sum(x*3**(11-i) for i,x in enumerate(t[:12]))
 if n%6 or n>len(t)-24:raise ValueError('header')
 if any(padded[24+n:]):raise ValueError('padding')
 body=t[12:12+n]
 if t[12+n:24+n]!=E.checksum_trits(body):raise ValueError('checksum')
 if any(sum(body[i+j]*3**(5-j) for j in range(6))>255 for i in range(0,n,6)):raise ValueError('byte')
 return E.trits_to_bytes(body)

def outcome(strands,stripes,msg):
 try:
  out=decode(strands,stripes)
  return ('exact_success' if out==msg else 'silent_wrong'),None,out
 except ValueError as e:return 'loud_failure',str(e),None

def compute():
 raw=(ROOT/'experiments/bound_nonce_file_plan.json').read_bytes();p=json.loads(raw);rng=np.random.default_rng(p['seed'])
 messages=[rng.bytes(n) for n in p['sizes'] for _ in range(p['identities_per_size'])]
 geometry=[];trials=[];controlled=[];wrong=[]
 for mi,msg in enumerate(messages):
  for stripes in p['stripes']:
   seq,ledger=G.encode_geometry(msg,stripes)
   g={'message':mi,'bytes':len(msg),'stripes':stripes,'source_sha256':hashlib.sha256(msg).hexdigest(),'included':seq is not None,'search_ledger':ledger};geometry.append(g)
   if seq is None:continue
   g.update({'strand_nt':len(seq),'total_nt':len(seq)*p['copies'],'bits_per_nt':8*len(msg)/(len(seq)*p['copies']),'strand_sha256':hashlib.sha256(seq.encode()).hexdigest()})
   cases={'single_data_erasure':['A'*64+seq[64:]]*3,'same_stripe_double_erasure':[('A'*64+seq[64:64*stripes]+'A'*64+seq[64*(stripes+1):])]*3,'one_copy_truncated':[seq,seq,seq[:-1]],'invalid_base':[seq,seq,'N'+seq[1:]]}
   for name,reads in cases.items():
    result,error,_=outcome(reads,stripes,msg);controlled.append({'message':mi,'stripes':stripes,'case':name,'outcome':result,'error':error})
   for rate in p['noise_rates']:
    for rep in range(p['noise_replicates']):
     reads=[];changes=0;seeds=[]
     for copy in range(p['copies']):
      seed=p['noise_seed_base']+mi*1000+rep*10+copy;seeds.append(seed);s,k=strict(seq,rate,seed);reads.append(s);changes+=k
     result,error,out=outcome(reads,stripes,msg)
     trials.append({'message':mi,'bytes':len(msg),'stripes':stripes,'rate':rate,'replicate':rep,'seeds':seeds,'changes':changes,'total_nt':g['total_nt'],'outcome':result,'error':error,'noisy_sha256':[hashlib.sha256(s.encode()).hexdigest() for s in reads]})
     if result=='silent_wrong':wrong.append({'trial':trials[-1],'expected_hex':msg.hex(),'actual_hex':out.hex(),'noisy_strands':reads})
   print('finished',mi,stripes,len(seq),flush=True)
 summaries=[]
 for size in p['sizes']:
  for stripes in p['stripes']:
   for rate in p['noise_rates']:
    rows=[r for r in trials if (r['bytes'],r['stripes'],r['rate'])==(size,stripes,rate)]
    summaries.append({'bytes':size,'stripes':stripes,'rate':rate,'n_trials':len(rows),'counts':dict(Counter(r['outcome'] for r in rows))})
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'geometry':geometry,'trials':trials,'controlled':controlled,'summaries':summaries,'silent_wrong':wrong,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/bound_nonce_file_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j['summaries'],indent=2))
