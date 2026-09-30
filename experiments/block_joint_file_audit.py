"""Whole-file clean/strict-noise detector comparison at charged, unequal budgets.
Audit candidate framing and deterministic tie policy; NOT production replacement.
"""
import hashlib,json,sys
from collections import Counter
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
from dnacell import encoder as E
import block_distance_code_audit as D
from strict_substitution_benchmark import strict
LENGTH={'low':38,'full':44,'distance5':48}
def pack(payload,bi,detector):
 if detector=='distance5':return D.encode(payload,bi)
 c=E.checksum_trits(payload)
 return E.encode_trits_never_same(payload+(c[-6:] if detector=='low' else c),E.BNS_SEEDS[bi%4])
def block_decode(seq,bi,detector):
 try:t=E.decode_never_same(seq,E.BNS_SEEDS[bi%4])
 except ValueError:return None
 if len(t)!=LENGTH[detector]:return None
 if detector=='distance5':return t[:32] if not any(D.syndrome(D.symbols(t))) else None
 c=E.checksum_trits(t[:32]);return t[:32] if t[32:]==(c[-6:] if detector=='low' else c) else None

def encode(message,copies,detector):
 body=E.bytes_to_trits(message);n=len(body);header=[(n//3**p)%3 for p in range(11,-1,-1)]
 t=E.scramble(header+body+E.checksum_trits(body));blocks=[(t[i:i+32]+[0]*32)[:32] for i in range(0,len(t),32)]
 parity=[sum(b[i] for b in blocks)%3 for i in range(32)];blocks.append(parity)
 dna=''.join(pack(b,i,detector) for i,b in enumerate(blocks));return [dna]*copies

def decode(strands,detector):
 if not strands or len({len(s) for s in strands})!=1:raise ValueError('strand geometry')
 blen=LENGTH[detector];n=len(strands[0])
 if n%blen or n//blen<2:raise ValueError('block geometry')
 # Fixed A,C,G,T tie ordering makes behavior reproducible.
 voted=''.join(max(E.BASES,key=lambda b:sum(s[i]==b for s in strands)) for i in range(n))
 blocks=[]
 for bi in range(n//blen):
  p=block_decode(voted[bi*blen:(bi+1)*blen],bi,detector)
  if p is None:
   for s in strands:
    p=block_decode(s[bi*blen:(bi+1)*blen],bi,detector)
    if p is not None:break
  blocks.append(p)
 parity=blocks.pop();bad=[i for i,b in enumerate(blocks) if b is None]
 if len(bad)==1 and parity is not None:
  blocks[bad[0]]=[(parity[k]-sum(b[k] for b in blocks if b is not None))%3 for k in range(32)]
 elif bad:raise ValueError('unrecoverable blocks')
 t=E.descramble(sum(blocks,[]));n=sum(x*3**(11-i) for i,x in enumerate(t[:12]))
 if n%6 or n>len(t)-24:raise ValueError('invalid header')
 body=t[12:12+n];chk=t[12+n:24+n]
 if chk!=E.checksum_trits(body):raise ValueError('global checksum')
 if any(sum(body[i+j]*3**(5-j) for j in range(6))>255 for i in range(0,n,6)):raise ValueError('invalid byte')
 return E.trits_to_bytes(body)

def main():
 rng=np.random.default_rng(930481);messages=[rng.bytes(n) for n in [128]*16+[2048]*4];rows=[];summaries=[];wrong=[]
 for detector in LENGTH:
  for copies in [1,3]:
   packed=[encode(m,copies,detector) for m in messages]
   assert all(decode(s,detector)==m for m,s in zip(messages,packed))
   for p in [0,.005,.01,.02,.03]:
    tally=Counter();bases=changes=0
    for i,(msg,strands) in enumerate(zip(messages,packed)):
     for replicate in range(3):
      noisy=[];mut=0
      for c,s in enumerate(strands):
       z,k=strict(s,p,100000+i*1000+replicate*20+c);noisy.append(z);mut+=k
      try:
       out=decode(noisy,detector);outcome='exact_success' if out==msg else 'silent_wrong';error=None
       if outcome=='silent_wrong':
        production=None
        if detector=='low':
         try:production=E.decode_message_v2(noisy).hex()
         except Exception as e:production=type(e).__name__+': '+str(e)
        wrong.append({'detector':detector,'message':i,'replicate':replicate,'p':p,'copies':copies,'expected_hex':msg.hex(),'actual_hex':out.hex(),'noisy_strands':noisy,'production_replay':production})
      except ValueError as e:outcome='loud_failure';error=str(e)
      tally[outcome]+=1;nt=sum(map(len,strands));bases+=nt;changes+=mut
      seq=strands[0];rows.append({'detector':detector,'copies':copies,'p':p,'message':i,'bytes':len(msg),'replicate':replicate,'outcome':outcome,'error':error,'changes':mut,'total_nt':nt,'strand_nt':len(seq),'blocks':len(seq)//LENGTH[detector],'rate':len(msg)*8/nt,'max_homopolymer':E.max_homopolymer(seq),'gc_fraction':E.gc_content(seq),'message_sha256':hashlib.sha256(msg).hexdigest(),'clean_dna_sha256':hashlib.sha256(seq.encode()).hexdigest()})
    summaries.append({'detector':detector,'copies':copies,'p':p,'n_trials':60,'n_distinct_messages':20,'counts':dict(tally),'changed_fraction':changes/bases,'total_exposed_nt':bases});print(detector,copies,p,dict(tally),flush=True)
 (ROOT/'results/block_joint_file_silent_wrong.json').write_text(json.dumps(wrong,indent=2)+'\n')
 raw=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in rows);(ROOT/'results/block_joint_file_trials.jsonl').write_text(raw)
 out={'message_seed':930481,'message_sizes':[128]*16+[2048]*4,'rows':summaries,'ledger_sha256':hashlib.sha256(raw.encode()).hexdigest(),'detector_block_nt':LENGTH,'limits':['Audit-local whole-file candidates use same one-block outer parity, deterministic plurality and copy fallback','Low/full candidate decoder adds explicit length/header/byte validation; not an unchanged production benchmark','Unequal detector lengths and copy counts mean unequal budgets and channel exposure; not a matched frontier','Twenty message identities; three noise realizations per condition are repeated measures','No oligo segmentation, primers, indexing, molecular copy model, indels, dropout or physical data','Whole-file seed resets can allow seam homopolymers; GC and homopolymer recorded, not screened']}
 (ROOT/'results/block_joint_file_audit.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
