"""Exhaust every first-block nonce substitution, separating core/file acceptance."""
import hashlib,json,sys
from pathlib import Path
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import local_constrained_audit as L
def compute():
 raw=(ROOT/'experiments/local_nonce_audit_plan.json').read_bytes();p=json.loads(raw);pool=json.loads((ROOT/'experiments/local_constrained_plan.json').read_text());rng=np.random.default_rng(pool['seed']);rows=[];wrong=[]
 for i,size in enumerate([n for n in pool['sizes'] for _ in range(4)]):
  msg=rng.bytes(size);seq,ledger=L.encode(msg);assert seq is not None
  for d,s in seq.items():
   original_t=L.E.decode_never_same(s[4:52],s[3]);original_nonce=sum(x*3**(3-k) for k,x in enumerate(L.E.decode_never_same(s[:4],'A')));original_payload=[(x-y)%3 for x,y in zip(original_t[:32],L.mask(original_nonce,0))]
   for pos in range(4):
    for base in L.E.BASES:
     if base==s[pos]:continue
     corrupt=s[:pos]+base+s[pos+1:];core_outcome=None
     try:
      nonce=sum(x*3**(3-k) for k,x in enumerate(L.E.decode_never_same(corrupt[:4],'A')))
      t=L.E.decode_never_same(corrupt[4:52],corrupt[3]);canon=L.E.encode_trits_never_same(t,L.E.BNS_SEEDS[0]);payload=L.M.block_decode(canon,0,d)
      if payload is None:core_outcome='core_reject'
      else:core_outcome='core_accept_correctpayload' if [(x-y)%3 for x,y in zip(payload,L.mask(nonce,0))]==original_payload else 'core_accept_wrongpayload'
     except ValueError:core_outcome='prefix_or_core_spacing_reject'
     try:out=L.decode([corrupt]*3,d);file_outcome='wholefile_exact' if out==msg else 'wholefile_silentwrong'
     except ValueError:file_outcome='wholefile_loud'
     rows.append({'message':i,'bytes':size,'detector':d,'position':pos,'replacement':base,'core_outcome':core_outcome,'file_outcome':file_outcome})
     if file_outcome=='wholefile_silentwrong':wrong.append({'message':i,'detector':d,'position':pos,'expected_hex':msg.hex(),'actual_hex':out.hex(),'corrupted_strand':corrupt})
 summary=[]
 for d in p['detectors']:
  rr=[r for r in rows if r['detector']==d];summary.append({'detector':d,'n_events':len(rr),'core_counts':dict(Counter(r['core_outcome'] for r in rr)),'file_counts':dict(Counter(r['file_outcome'] for r in rr))})
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'pool_plan_sha256':hashlib.sha256((ROOT/'experiments/local_constrained_plan.json').read_bytes()).hexdigest(),'rows':rows,'summary':summary,'wrong':wrong,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/local_nonce_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j['summary'])
