"""Geometry control using identical payload stream to fountain_port_audit.py.
One read per indexed oligo, ideal known order, no outer redundancy or primers.
"""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'))
from rs_inner import rsns_encode,rsns_decode
from strict_substitution_benchmark import strict

def run():
 rng=np.random.default_rng(511);rows=[]
 for i in range(4):
  msg=rng.bytes(2048);parts=[msg[k:k+32] for k in range(0,2048,32)];dna=[rsns_encode(p) for p in parts]
  for rate in [0,.005,.01,.02,.03]:
   for seed in range(3):
    correct=[];changes=0;errors=[]
    for c,(p,s) in enumerate(zip(parts,dna)):
     z,n=strict(s,rate,91000+i*1000+seed*300+c);changes+=n
     try:ok=rsns_decode(z)==p;err=None
     except Exception as e:ok=False;err=type(e).__name__+': '+str(e)
     correct.append(ok);errors.append(err)
    rows.append({'message':i,'noise_seed':seed,'p':rate,'whole_file_recovered':all(correct),'chunks_recovered':sum(correct),'chunks':64,'oligo_nt':len(dna[0]),'total_encoded_nt':sum(map(len,dna)),'payload_bits_per_encoded_nt':8*2048/sum(map(len,dna)),'realized_substitution_fraction':changes/sum(map(len,dna)),'input_sha256':hashlib.sha256(msg).hexdigest(),'chunk_success':correct,'errors':errors})
 out={'design':'same 4x2KB payloads split into 64 independent 32-byte RSNS oligos; no outer code; one read each','rows':rows,'limits':['Oligo length 495 nt, still longer than 152-nt Fountain oligos','Ideal known order and chunk identity, no index charged','No RSNS outer redundancy; neither budget optimized','Whole-file success requires all oligos to decode; chunk successes are not independent whole-message replicates','No physical synthesis validation or primers']}
 (ROOT/'results/rsns_segmented_audit.json').write_text(json.dumps(out,indent=2)+'\n')
 for rate in [0,.005,.01,.02,.03]:
  a=[x for x in rows if x['p']==rate];print(rate,sum(x['whole_file_recovered'] for x in a),'/12',sum(x['chunks_recovered'] for x in a),'/768')
if __name__=='__main__':run()
