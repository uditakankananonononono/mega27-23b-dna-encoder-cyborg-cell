"""Preregistered fixed-geometry detector/copy allocation comparison.
Low/full checks receive charged ignored padding, not free protection.
"""
import hashlib,json,sys
from collections import Counter
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import importlib.util
spec=importlib.util.spec_from_file_location('matched_file_wrapper',ROOT/'experiments/block_joint_file_audit.py')
F=importlib.util.module_from_spec(spec);spec.loader.exec_module(F)
E=F.E;D=F.D
ORIG_PACK=F.pack;ORIG_BLOCK=F.block_decode

def pack(payload,bi,detector):
 seq=ORIG_PACK(payload,bi,detector);t=E.decode_never_same(seq,E.BNS_SEEDS[bi%4])
 return E.encode_trits_never_same(t+[0]*(48-len(t)),E.BNS_SEEDS[bi%4])
def block_decode(seq,bi,detector):
 if len(seq)!=48:return None
 # Core decode includes transition from payload into check. Ignored tail is
 # known padding and charged exposure, not a detector or recovery safeguard.
 core={'low':38,'full':44,'distance5':48}[detector]
 try:t=E.decode_never_same(seq[:core],E.BNS_SEEDS[bi%4])
 except ValueError:return None
 if detector=='distance5':return t[:32] if not any(D.syndrome(D.symbols(t))) else None
 chk=E.checksum_trits(t[:32]);return t[:32] if t[32:]==(chk[-6:] if detector=='low' else chk) else None
F.pack=pack;F.block_decode=block_decode;F.LENGTH={k:48 for k in F.LENGTH}
def main():
 rawplan=(ROOT/'experiments/block_matched_plan.json').read_bytes();plan=json.loads(rawplan)
 rng=np.random.default_rng(plan['message_seed']);messages=[rng.bytes(n) for n in plan['sizes'] for _ in range(plan['identities_per_size'])]
 packed={};exclusions=[];geometry=[];included=[]
 for i,msg in enumerate(messages):
  candidate={d:F.encode(msg,1,d) for d in plan['detectors']}
  metrics={d:{'max_homopolymer':E.max_homopolymer(s[0]),'gc':E.gc_content(s[0]),'strand_nt':len(s[0])} for d,s in candidate.items()}
  passed=all(r['max_homopolymer']<=2 and .45<=r['gc']<=.55 for r in metrics.values())
  geometry.append({'message':i,'bytes':len(msg),'included':passed,'metrics':metrics,'input_sha256':hashlib.sha256(msg).hexdigest()})
  if not passed:exclusions.append(i);continue
  included.append(i)
  assert len({len(s[0]) for s in candidate.values()})==1
  for d,s in candidate.items():
   for copies in plan['copies']:
    ss=s*copies;assert F.decode(ss,d)==msg;packed[i,d,copies]=ss
 rows=[];summaries=[];wrong=[]
 for size in plan['sizes']:
  identities=[i for i in included if len(messages[i])==size]
  for copies in plan['copies']:
   for p in plan['strict_p']:
    for d in plan['detectors']:
     counts=Counter();changes=total=0
     for i in identities:
      for rep in range(plan['replicates_per_identity']):
       noisy=[];mut=0;strands=packed[i,d,copies]
       for c,s in enumerate(strands):
        z,k=F.strict(s,p,200000+i*1000+rep*20+c);noisy.append(z);mut+=k
       try:out=F.decode(noisy,d);outcome='exact_success' if out==messages[i] else 'silent_wrong';error=None
       except ValueError as e:outcome='loud_failure';error=str(e)
       if outcome=='silent_wrong':wrong.append({'message':i,'detector':d,'copies':copies,'p':p,'replicate':rep,'expected_hex':messages[i].hex(),'actual_hex':out.hex(),'noisy_strands':noisy})
       nt=sum(map(len,strands));changes+=mut;total+=nt;counts[outcome]+=1
       rows.append({'bytes':size,'message':i,'detector':d,'copies':copies,'p':p,'replicate':rep,'outcome':outcome,'error':error,'changed_nt':mut,'total_nt':nt,'rate':size*8/nt})
     summaries.append({'bytes':size,'detector':d,'copies':copies,'p':p,'n_messages':len(identities),'n_trials':len(identities)*3,'counts':dict(counts),'changed_fraction':changes/total if total else None,'total_exposed_nt':total});print(size,copies,p,d,dict(counts),flush=True)
 raw=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in rows)
 (ROOT/'results/block_matched_geometry_trials.jsonl').write_text(raw)
 (ROOT/'results/block_matched_geometry_wrong.json').write_text(json.dumps(wrong,indent=2)+'\n')
 out={'preregistered_plan_sha256':hashlib.sha256(rawplan).hexdigest(),'ledger_sha256':hashlib.sha256(raw.encode()).hexdigest(),'geometry':geometry,'excluded_message_ids':exclusions,'summaries':summaries,'limits':['Exact same payload identities, 48nt block geometry, copy count, constraints and total nt per detector condition','Noise masks and alternative-letter indices paired by equal length/seed, but different sequences yield different changed letters','Padding is charged but ignored; this isolates detector allocation, not optimized baseline performance','Whole-strand GC only, no sliding-window gate, motif gate or physical oligo geometry','New sixteen identities; three repeats are not independent messages','Local software comparison, no external codec and no optimized frontier claim']}
 (ROOT/'results/block_matched_geometry_audit.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
