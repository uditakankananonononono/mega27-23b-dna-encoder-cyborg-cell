"""Audit-local nonce-bound [13,9,5] core with explicitly larger frame."""
import hashlib,itertools,json,sys
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import protected_frame_three_audit as T
D=T.F.P.M.D;E=T.F.P.E
H=[[D.power(i,r) for i in range(13)] for r in range(4)]
def syndrome(symbols):
 out=[]
 for row in H:
  s=0
  for h,x in zip(row,symbols):s=D.ADD[s][D.MUL[h][x]]
  out.append(s)
 return out
def nonce_trits(n):return [(n//3**p)%3 for p in range(3,-1,-1)]
def frame(payload,nonce):
 data=nonce_trits(nonce)+[(x+y)%3 for x,y in zip(payload,T.F.P.mask(nonce,0))];u=D.symbols(data);parity=D.solve([r[9:] for r in H],[D.neg(x) for x in syndrome(u+[0]*4)]);core=E.encode_trits_never_same(data+sum([D.digits(x) for x in parity],[]),'A');return T.prefix(nonce)+core
def decode(s):
 if len(s)!=64:raise ValueError('geometry')
 z=[E.decode_never_same(s[k:k+4],'A') for k in [0,4,8]]
 if not z[0]==z[1]==z[2]:raise ValueError('replicas')
 t=E.decode_never_same(s[12:],'A')
 if any(syndrome(D.symbols(t))):raise ValueError('syndrome')
 if t[:4]!=z[0]:raise ValueError('nonce association')
 n=sum(x*3**(3-i) for i,x in enumerate(t[:4]));return [(x-y)%3 for x,y in zip(t[4:36],T.F.P.mask(n,0))]
def compute():
 raw=(ROOT/'experiments/bound_nonce_frame_plan.json').read_bytes();rows=[];payloads=[[0]*32,[i%3 for i in range(32)]]
 for pi,p in enumerate(payloads):
  for nonce in [0,80]:
   s=frame(p,nonce);assert decode(s)==p;counts=Counter()
   for order in [1,2]:
    for pos in itertools.combinations(range(64),order):
     for letters in itertools.product(*[[b for b in E.BASES if b!=s[i]] for i in pos]):
      changed=list(s)
      for i,b in zip(pos,letters):changed[i]=b
      try:out=decode(''.join(changed));outcome='accepted_same' if out==p else 'accepted_wrong'
      except ValueError:outcome='rejected'
      counts[f'{order}:{outcome}']+=1
   rows.append({'payload_id':pi,'nonce':nonce,'frame_sha256':hashlib.sha256(s.encode()).hexdigest(),'counts':dict(counts)})
 prefix_counts=Counter()
 for pi,p in enumerate(payloads):
  for a in range(81):
   s=frame(p,a)
   for b in range(81):
    if sum(x!=y for x,y in zip(T.prefix(a),T.prefix(b)))!=3:continue
    try:out=decode(T.prefix(b)+s[12:]);outcome='accepted_same' if out==p else 'accepted_wrong'
    except ValueError:outcome='rejected'
    prefix_counts[outcome]+=1
 # Every four columns have distinct-coordinate Vandermonde factors.
 independent=sum(all(D.ADD[b][D.neg(a)]!=0 for a,b in itertools.combinations(cols,2)) for cols in itertools.combinations(range(13),4))
 coefficients=D.solve([[r[i] for i in range(4)] for r in H],[D.neg(r[4]) for r in H])+[1];w=coefficients+[0]*8;assert all(coefficients) and not any(syndrome(w))
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'frame_nt':64,'payload_trits':32,'checked_nonce_trits':4,'parity_nt':16,'field_code':[13,9,5],'four_column_sets':independent,'weight5_witness':w,'rows':rows,'prefix_three_counts':dict(prefix_counts),'limits':json.loads(raw)['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/bound_nonce_frame_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j)
