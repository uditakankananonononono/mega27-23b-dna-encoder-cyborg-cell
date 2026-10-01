"""Prospective joint finite-nonce framing with all emitted nucleotide costs charged."""
import hashlib,json,sys
from pathlib import Path
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import block_stripe_frontier_audit as S
M=S.M;E=S.E
def mask(nonce,bi):
 h=hashlib.sha256(f'local-nonce:{bi}:{nonce}'.encode()).digest();return [b%3 for b in h]
def core(payload,bi,d,seed):
 t=E.decode_never_same(M.pack(payload,bi,d),E.BNS_SEEDS[bi%4]);return E.encode_trits_never_same(t,seed)
def local_ok(seq,block):
 return E.max_homopolymer(seq)<=2 and .35<=E.gc_content(block)<=.65 and (len(seq)<64 or (lambda z:z[0]>=.35 and z[1]<=.65)(S.gc_range(seq,64)))
def encode(msg):
 body=E.bytes_to_trits(msg);n=len(body);header=[(n//3**p)%3 for p in range(11,-1,-1)];t=E.scramble(header+body+E.checksum_trits(body));blocks=[(t[i:i+32]+[0]*32)[:32] for i in range(0,len(t),32)];blocks.append([sum(b[k] for b in blocks)%3 for k in range(32)])
 seq={d:'' for d in ['low','full','distance5']};ledger=[]
 for bi,payload in enumerate(blocks):
  for nonce in range(81):
   trits=[(nonce//3**p)%3 for p in range(3,-1,-1)];prefix=E.encode_trits_never_same(trits,'A');masked=[(x+y)%3 for x,y in zip(payload,mask(nonce,bi))];cand={d:prefix+core(masked,bi,d,prefix[-1]) for d in seq}
   if all(local_ok(seq[d]+cand[d],cand[d]) for d in seq):break
  else:return None,ledger+[{'block':bi,'attempts':81,'accepted_nonce':None}]
  ledger.append({'block':bi,'attempts':nonce+1,'accepted_nonce':nonce})
  for d in seq:seq[d]+=cand[d]
 if not all(.45<=E.gc_content(s)<=.55 for s in seq.values()):return None,ledger
 return seq,ledger
def decode(strands,d):
 if not strands or len(set(map(len,strands)))!=1 or len(strands[0])%52:raise ValueError('geometry')
 voted=''.join(max(E.BASES,key=lambda b:sum(s[i]==b for s in strands)) for i in range(len(strands[0])));blocks=[]
 for bi in range(len(voted)//52):
  payload=None
  for source in [voted]+strands:
   block=source[bi*52:(bi+1)*52]
   try:
    nonce=sum(x*3**(3-i) for i,x in enumerate(E.decode_never_same(block[:4],'A')));raw=E.decode_never_same(block[4:],block[3]);canon=E.encode_trits_never_same(raw,E.BNS_SEEDS[bi%4]);p=M.block_decode(canon,bi,d)
    if p is not None:payload=[(x-y)%3 for x,y in zip(p,mask(nonce,bi))];break
   except ValueError:pass
  blocks.append(payload)
 parity=blocks.pop();bad=[i for i,b in enumerate(blocks) if b is None]
 if len(bad)==1 and parity is not None:blocks[bad[0]]=[(parity[k]-sum(b[k] for b in blocks if b is not None))%3 for k in range(32)]
 elif bad:raise ValueError('erasure')
 t=E.descramble(sum(blocks,[]));n=sum(x*3**(11-i) for i,x in enumerate(t[:12]))
 if n%6 or n>len(t)-24:raise ValueError('header')
 body=t[12:12+n]
 if t[12+n:24+n]!=E.checksum_trits(body):raise ValueError('checksum')
 if any(sum(body[i+j]*3**(5-j) for j in range(6))>255 for i in range(0,n,6)):raise ValueError('byte')
 return E.trits_to_bytes(body)
def main():
 raw=(ROOT/'experiments/local_constrained_plan.json').read_bytes();p=json.loads(raw);rng=np.random.default_rng(p['seed']);geometry=[];rows=[];wrong=[]
 for i,size in enumerate([n for n in p['sizes'] for _ in range(p['identities_per_size'])]):
  msg=rng.bytes(size);seq,ledger=encode(msg);g={'message':i,'bytes':size,'source_sha256':hashlib.sha256(msg).hexdigest(),'included':seq is not None,'search_ledger':ledger,'candidate_attempts':sum(r['attempts'] for r in ledger),'candidate_rejections':sum(r['attempts']-int(r['accepted_nonce'] is not None) for r in ledger)};geometry.append(g)
  if seq is None:print('EXCLUDED',i,g['candidate_attempts'],flush=True);continue
  g['metrics']={d:{'nt_per_copy':len(s),'total_nt':len(s)*p['copies'],'bits_per_nt':size*8/(len(s)*p['copies']),'hp':E.max_homopolymer(s),'whole_gc':E.gc_content(s),'sliding64_gc_range':S.gc_range(s,64),'block52_gc_range':[min(E.gc_content(s[k:k+52]) for k in range(0,len(s),52)),max(E.gc_content(s[k:k+52]) for k in range(0,len(s),52))]} for d,s in seq.items()}
  for d,s in seq.items():
   assert decode([s]*p['copies'],d)==msg
   for rate in p['strict_p']:
    for rep in range(p['noise_repeats']):
     noisy=[S.F.strict(s,rate,990000+i*100+rep*10+c)[0] for c in range(p['copies'])]
     try:out=decode(noisy,d);outcome='exact_success' if out==msg else 'silent_wrong'
     except ValueError:outcome='loud_failure'
     rows.append({'message':i,'bytes':size,'detector':d,'p':rate,'repeat':rep,'outcome':outcome})
     if outcome=='silent_wrong':wrong.append({'message':i,'detector':d,'p':rate,'repeat':rep,'expected_hex':msg.hex(),'actual_hex':out.hex(),'noisy_strands':noisy})
  print('INCLUDED',i,g['candidate_attempts'],g['metrics']['distance5']['total_nt'],flush=True)
 summaries=[]
 for size in p['sizes']:
  for d in p['detectors']:
   for rate in p['strict_p']:
    rr=[r for r in rows if r['bytes']==size and r['detector']==d and r['p']==rate];summaries.append({'bytes':size,'detector':d,'p':rate,'n_trials':len(rr),'counts':dict(Counter(r['outcome'] for r in rr))})
 j={'plan_sha256':hashlib.sha256(raw).hexdigest(),'geometry':geometry,'rows':rows,'summaries':summaries,'limits':p['limits']};(ROOT/'results/local_constrained_audit.json').write_text(json.dumps(j,indent=2)+'\n');(ROOT/'results/local_constrained_wrong.json').write_text(json.dumps(wrong,indent=2)+'\n');print(summaries)
if __name__=='__main__':main()
