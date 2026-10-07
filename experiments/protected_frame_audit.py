"""Aligned distance-five frame audit. Does not replace production or file decoder."""
import hashlib,itertools,json,sys
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import protected_nonce_audit as P
def frame(payload,nonce,bi=0):
 t=[(nonce//3**p)%3 for p in range(3,-1,-1)];masked=[(x+y)%3 for x,y in zip(payload,P.mask(nonce,bi))];return P.E.encode_trits_never_same(t,'A')*3+P.core(masked,bi,'distance5','A')
def decode_frame(s,bi=0):
 if len(s)!=60:raise ValueError('geometry')
 z=[P.E.decode_never_same(s[k:k+4],'A') for k in [0,4,8]]
 if not z[0]==z[1]==z[2]:raise ValueError('nonce disagreement')
 nonce=sum(x*3**(3-i) for i,x in enumerate(z[0]));t=P.E.decode_never_same(s[12:],'A');canon=P.E.encode_trits_never_same(t,P.E.BNS_SEEDS[bi%4]);p=P.M.block_decode(canon,bi,'distance5')
 if p is None:raise ValueError('syndrome')
 return [(x-y)%3 for x,y in zip(p,P.mask(nonce,bi))]
def compute():
 raw=(ROOT/'experiments/protected_frame_plan.json').read_bytes();plan=json.loads(raw);rows=[]
 for payload_id,payload in enumerate([[0]*32,[i%3 for i in range(32)]]):
  for nonce in plan['nonces']:
   s=frame(payload,nonce);assert decode_frame(s)==payload;counts=Counter();partition=Counter()
   for order in [1,2]:
    for positions in itertools.combinations(range(60),order):
     for letters in itertools.product(*[[b for b in P.E.BASES if b!=s[i]] for i in positions]):
      changed=list(s)
      for i,b in zip(positions,letters):changed[i]=b
      try:out=decode_frame(''.join(changed));outcome='accepted_unchanged' if out==payload else 'accepted_wrong'
      except ValueError:outcome='rejected'
      counts[f'{order}:{outcome}']+=1;region='prefix' if max(positions)<12 else ('core' if min(positions)>=12 else 'mixed');partition[f'{order}:{region}:{outcome}']+=1
   rows.append({'payload_id':payload_id,'nonce':nonce,'frame_sha256':hashlib.sha256(s.encode()).hexdigest(),'counts':dict(counts),'partitions':dict(partition)});print(rows[-1],flush=True)
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'rows':rows,'limits':plan['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/protected_frame_audit.json').write_text(json.dumps(j,indent=2)+'\n')
