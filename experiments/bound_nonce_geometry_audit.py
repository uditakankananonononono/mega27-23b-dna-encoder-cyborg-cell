"""Prospective feasibility and charged geometry for nonce-bound frames."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import bound_nonce_frame_audit as B
P=B.T.F.P;E=B.E
def encode_geometry(msg,stripes):
 body=E.bytes_to_trits(msg);n=len(body);header=[(n//3**p)%3 for p in range(11,-1,-1)];t=E.scramble(header+body+E.checksum_trits(body));data=[(t[i:i+32]+[0]*32)[:32] for i in range(0,len(t),32)];blocks=data+[[sum(data[b][k] for b in range(g,len(data),stripes))%3 for k in range(32)] for g in range(stripes)];seq='';ledger=[]
 # Frame implementation uses a fixed block-zero mask, audited here without file decoder.
 for bi,payload in enumerate(blocks):
  for nonce in range(81):
   cand=B.frame(payload,nonce);joined=seq+cand;passed=E.max_homopolymer(joined)<=2 and .35<=E.gc_content(cand)<=.65 and (len(joined)<64 or (lambda r:r[0]>=.35 and r[1]<=.65)(P.S.gc_range(joined,64)))
   if passed:break
  else:return None,ledger+[{'block':bi,'attempts':81,'nonce':None}]
  assert B.decode(cand)==payload;ledger.append({'block':bi,'attempts':nonce+1,'nonce':nonce});seq=joined
 if not .45<=E.gc_content(seq)<=.55:return None,ledger
 return seq,ledger
def compute():
 raw=(ROOT/'experiments/bound_nonce_geometry_plan.json').read_bytes();p=json.loads(raw);rng=np.random.default_rng(p['seed']);rows=[]
 for mi,msg in enumerate([rng.bytes(n) for n in p['sizes'] for _ in range(p['identities_per_size'])]):
  for stripes in p['stripes']:
   seq,ledger=encode_geometry(msg,stripes);r={'message':mi,'bytes':len(msg),'stripes':stripes,'source_sha256':hashlib.sha256(msg).hexdigest(),'included':seq is not None,'search_ledger':ledger,'candidate_attempts':sum(x['attempts'] for x in ledger)}
   if seq is not None:r.update({'strand_nt':len(seq),'total_nt_three_copies':3*len(seq),'bits_per_nt':8*len(msg)/(3*len(seq)),'whole_gc':E.gc_content(seq),'max_homopolymer':E.max_homopolymer(seq),'sliding64_gc_range':P.S.gc_range(seq,64),'frame_gc_range':[min(E.gc_content(seq[i:i+64]) for i in range(0,len(seq),64)),max(E.gc_content(seq[i:i+64]) for i in range(0,len(seq),64))]})
   rows.append(r);print(mi,stripes,r['included'],r.get('total_nt_three_copies'),r['candidate_attempts'],flush=True)
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'rows':rows,'limits':p['limits'],'mask_scope':'fixed block-zero mask as isolated candidate;notdrop-infilecodec'}
if __name__=='__main__':
 j=compute();(ROOT/'results/bound_nonce_geometry_audit.json').write_text(json.dumps(j,indent=2)+'\n')
