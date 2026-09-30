"""Fresh smaller benchmark with guaranteed-different-base substitutions.
Reimplemented baselines, not original published codecs. No historical overwrite.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(Path(__file__).parent))
import benchmark_codecs as B
from rs_inner import rsns_encode,rsns_decode
from dnacell import encoder as E
def strict(s,p,seed):
 rng=np.random.default_rng(seed);out=[];changes=0
 for ch in s:
  if rng.random()<p:
   choices=[x for x in E.BASES if x!=ch];out.append(choices[int(rng.integers(3))]);changes+=1
  else:out.append(ch)
 return ''.join(out),changes

def main():
 rng=np.random.default_rng(30930);messages=[rng.bytes(128) for _ in range(32)];defs={'v1':(lambda m:(E.encode_message(m),None),lambda s,p:E.decode_message(s)),'v2':(lambda m:(E.encode_message_v2(m),None),lambda s,p:E.decode_message_v2(s)),'goldman_style':(lambda m:(B.goldman_encode(m),None),lambda s,p:B.goldman_decode(s)),'fountain_crc8':(B.fountain_encode,lambda s,p:B.fountain_decode(s,p)[:128]),'rsns':(lambda m:([rsns_encode(m)],None),lambda s,p:rsns_decode(s[0])[:128])};rows=[]
 for name,(enc,dec) in defs.items():
  packed=[enc(m) for m in messages]
  for rate in [0,.01,.02,.03]:
   trials=[];actual=0;bases=0
   for i,(m,(strands,p)) in enumerate(zip(messages,packed)):
    noisy=[]
    for c,s in enumerate(strands):
     a,n=strict(s,rate,90000+97*i+c);actual+=n;bases+=len(s);noisy.append(a)
    try:ok=dec(noisy,p)==m;error=None
    except Exception as e:ok=False;error=type(e).__name__+': '+str(e)[:100]
    trials.append({'message':i,'recovered':bool(ok),'error':error})
   ok=sum(t['recovered'] for t in trials);z=1.96;nn=len(trials);ph=ok/nn;den=1+z*z/nn;center=(ph+z*z/(2*nn))/den;half=z*np.sqrt(ph*(1-ph)/nn+z*z/(4*nn*nn))/den
   row={'codec':name,'requested_realized_substitution_p':rate,'observed_changed_base_fraction':actual/bases,'recovered':ok,'n':nn,'recovery':ph,'wilson95':[float(center-half),float(center+half)],'mean_payload_bits_per_synthesized_base':float(np.mean([8*128/sum(map(len,s)) for s,p in packed])),'trials':trials};rows.append(row);print(name,rate,ok,actual/bases,flush=True)
 out={'payload_bytes':128,'message_seed':30930,'rows':rows,'channel':'independent Bernoulli attempt and uniform alternative among exactly three different bases','limits':['Only 32 messages per condition at one payload size','Same message identities, but unequal strand lengths consume RNG differently, not paired identical corrupted bases','Fountain uses CRC8 detection, uniform low-degree sampling and side metadata, NOT the published implementation or equivalent RS inner correction','RSNS implemented erasure handling needs its own audit; this records actual behavior without endorsing mechanism','No synthesis/sequencing measurements']}
 (ROOT/'results/strict_substitution_benchmark.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
