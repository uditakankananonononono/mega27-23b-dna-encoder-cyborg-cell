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
def encode(msg,stripes):
 body=E.bytes_to_trits(msg);n=len(body);header=[(n//3**p)%3 for p in range(11,-1,-1)];t=E.scramble(header+body+E.checksum_trits(body));blocks=[(t[i:i+32]+[0]*32)[:32] for i in range(0,len(t),32)];data=blocks[:];blocks.extend([[sum(data[b][k] for b in range(g,len(data),stripes))%3 for k in range(32)] for g in range(stripes)])
 seq={d:'' for d in ['low','full','distance5']};ledger=[]
 for bi,payload in enumerate(blocks):
  for nonce in range(81):
   trits=[(nonce//3**p)%3 for p in range(3,-1,-1)];prefix=E.encode_trits_never_same(trits,'A')*3;masked=[(x+y)%3 for x,y in zip(payload,mask(nonce,bi))];cand={d:prefix+core(masked,bi,d,'A') for d in seq}
   if all(local_ok(seq[d]+cand[d],cand[d]) for d in seq):break
  else:return None,ledger+[{'block':bi,'attempts':81,'accepted_nonce':None}]
  ledger.append({'block':bi,'attempts':nonce+1,'accepted_nonce':nonce})
  for d in seq:seq[d]+=cand[d]
 if not all(.45<=E.gc_content(s)<=.55 for s in seq.values()):return None,ledger
 return seq,ledger
def decode(strands,d,stripes):
 if d not in ['low','full','distance5']:raise ValueError('detector')
 if isinstance(stripes,bool) or not isinstance(stripes,int) or stripes<1:raise ValueError('stripes')
 if not isinstance(strands,(list,tuple)) or not strands or any(not isinstance(s,str) for s in strands):raise ValueError('geometry')
 if len(set(map(len,strands)))!=1 or not strands[0] or len(strands[0])%60 or len(strands[0])//60<=stripes:raise ValueError('geometry')
 if any(set(s)-set(E.BASES) for s in strands):raise ValueError('alphabet')
 voted=''.join(max(E.BASES,key=lambda b:sum(s[i]==b for s in strands)) for i in range(len(strands[0])));blocks=[]
 for bi in range(len(voted)//60):
  payload=None
  for source in [voted]+strands:
   block=source[bi*60:(bi+1)*60]
   try:
    prefixes=[E.decode_never_same(block[k:k+4],'A') for k in [0,4,8]];
    if not prefixes[0]==prefixes[1]==prefixes[2]:continue
    nonce=sum(x*3**(3-i) for i,x in enumerate(prefixes[0]));raw=E.decode_never_same(block[12:],'A');canon=E.encode_trits_never_same(raw,E.BNS_SEEDS[bi%4]);p=M.block_decode(canon,bi,d)
    if p is not None:payload=[(x-y)%3 for x,y in zip(p,mask(nonce,bi))];break
   except ValueError:pass
  blocks.append(payload)
 data_count=len(blocks)-stripes
 for g in range(stripes):
  ids=list(range(g,data_count,stripes));bad=[i for i in ids if blocks[i] is None];parity=blocks[data_count+g]
  if len(bad)==1 and parity is not None:blocks[bad[0]]=[(parity[k]-sum(blocks[i][k] for i in ids if blocks[i] is not None))%3 for k in range(32)]
  elif bad:raise ValueError('erasure')
 blocks=blocks[:data_count]
 t=E.descramble(sum(blocks,[]));n=sum(x*3**(11-i) for i,x in enumerate(t[:12]))
 if n%6 or n>len(t)-24:raise ValueError('header')
 body=t[12:12+n]
 if t[12+n:24+n]!=E.checksum_trits(body):raise ValueError('checksum')
 if any(sum(body[i+j]*3**(5-j) for j in range(6))>255 for i in range(0,n,6)):raise ValueError('byte')
 return E.trits_to_bytes(body)
def main():
 raw=(ROOT/'experiments/protected_nonce_plan.json').read_bytes();p=json.loads(raw);rng=np.random.default_rng(p['seed']);messages=[rng.bytes(n) for n in p['sizes'] for _ in range(p['identities_per_size'])];geometry=[];rows=[];wrong=[]
 for i,msg in enumerate(messages):
  for stripes in p['stripes']:
   seq,ledger=encode(msg,stripes);g={'message':i,'bytes':len(msg),'stripes':stripes,'source_sha256':hashlib.sha256(msg).hexdigest(),'included':seq is not None,'search_ledger':ledger,'candidate_attempts':sum(r['attempts'] for r in ledger),'candidate_rejections':sum(r['attempts']-int(r['accepted_nonce'] is not None) for r in ledger)};geometry.append(g)
   if seq is None:print('EXCLUDED',i,stripes,flush=True);continue
   g['metrics']={d:{'total_nt':len(s)*3,'bits_per_nt':len(msg)*8/(len(s)*3),'whole_gc':E.gc_content(s),'hp':E.max_homopolymer(s),'sliding64_gc_range':S.gc_range(s,64)} for d,s in seq.items()}
   for d,s in seq.items():
    assert decode([s]*3,d,stripes)==msg
    for rate in p['strict_p']:
     for rep in range(2):
      noisy=[S.F.strict(s,rate,1100000+i*100+rep*10+c)[0] for c in range(3)]
      try:out=decode(noisy,d,stripes);outcome='exact_success' if out==msg else 'silent_wrong'
      except ValueError:outcome='loud_failure'
      rows.append({'message':i,'bytes':len(msg),'stripes':stripes,'detector':d,'p':rate,'repeat':rep,'outcome':outcome})
      if outcome=='silent_wrong':wrong.append({'message':i,'stripes':stripes,'detector':d,'expected_hex':msg.hex(),'actual_hex':out.hex(),'noisy_strands':noisy})
   print('INCLUDED',i,stripes,g['candidate_attempts'],g['metrics']['distance5']['total_nt'],flush=True)
 summaries=[]
 for size in p['sizes']:
  for stripes in p['stripes']:
   for d in p['detectors']:
    for rate in p['strict_p']:
     rr=[r for r in rows if (r['bytes'],r['stripes'],r['detector'],r['p'])==(size,stripes,d,rate)];summaries.append({'bytes':size,'stripes':stripes,'detector':d,'p':rate,'n_trials':len(rr),'counts':dict(Counter(r['outcome'] for r in rr))})
 j={'plan_sha256':hashlib.sha256(raw).hexdigest(),'geometry':geometry,'rows':rows,'summaries':summaries,'wrong':wrong,'limits':p['limits'],'proof_scope':p['proof_scope']};(ROOT/'results/protected_nonce_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(summaries)
if __name__=='__main__':main()
