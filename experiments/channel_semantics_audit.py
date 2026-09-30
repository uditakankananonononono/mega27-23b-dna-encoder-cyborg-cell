"""Audit actual substitution probability and differential-code error locality.
No implicit downstream error-propagation claim; exhaustive single-base changes.
"""
import json,sys,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from dnacell import encoder as E
def raw_decode(s):
 prev='A';out=[]
 for ch in s:out.append((E.B2I[ch]-E.B2I[prev]-1)%4);prev=ch
 return out

def main():
 rng=np.random.default_rng(941);events=[];invalid=0;changed_counts=[];outside=0
 for t in range(8):
  tr=rng.integers(0,3,120).tolist();s=E.encode_trits_never_same(tr)
  for pos in range(len(s)):
   for b in E.BASES:
    if b==s[pos]:continue
    altered=s[:pos]+b+s[pos+1:];z=raw_decode(altered);ix=[i for i,(a,c) in enumerate(zip(tr,z)) if a!=c];bad=[i for i,c in enumerate(z) if c==3]
    changed_counts.append(len(ix));invalid+=bool(bad);outside+=any(i not in [pos,pos+1] for i in ix)
    if len(events)<12:events.append({'position':pos,'changed_trits':ix,'out_of_alphabet_trits':bad})
 assert outside==0
 empirical=[]
 s='ACGT'*25000
 for p in [.01,.02,.03,.1]:
  rates=[]
  for seed in range(5):z=E.introduce_errors(s,p,seed);rates.append(sum(a!=b for a,b in zip(s,z))/len(s))
  empirical.append({'attempt_probability':p,'expected_realized_substitution_probability':.75*p,'realized_rates':rates,'mean':float(np.mean(rates))})
 out={'single_base_events':len(changed_counts),'changed_trit_min':min(changed_counts),'changed_trit_max':max(changed_counts),'events_changing_trits_outside_i_i_plus_1':outside,'events_with_invalid_trit_3':invalid,'examples':events,'noise_generator':empirical,'proof':'Each decoded trit uses only current and previous DNA base. Changing base i changes decoded trit i and (if present) i+1. Full-strand decode can reject on invalid spacing/checksum, but this is detection, not downstream symbol desynchronization. Insertions/deletions are different.','source_code_sha256':hashlib.sha256((ROOT/'src/dnacell/encoder.py').read_bytes()).hexdigest(),'limits':['Audit applies to implemented differential mapping, not every published rotating code','Historical nominal-rate benchmark outputs are unchanged; a strict-different-base rerun is required for realized-rate comparisons']}
 (ROOT/'results/channel_semantics_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['examples']},indent=2))
if __name__=='__main__':main()
