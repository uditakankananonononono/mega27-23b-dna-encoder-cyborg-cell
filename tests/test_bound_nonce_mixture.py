import json,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import bound_nonce_mixture_audit as M

def test_mixture_reproduces_endpoint_and_partial_records():
 j=json.loads((ROOT/'results/bound_nonce_mixture_audit.json').read_text());assert M.compute()==j
 assert len(j['rows'])==360 and not j['silent_wrong']
 assert [sum(r['outcome']=='exact_success' for r in j['rows'] if r['q']==q) for q in [0,.25,.5,.75,1]]==[49,8,1,0,0]
 for mi in range(4):
  for stripe in [1,4]:
   for rate in [.005,.01,.02]:
    for rep in range(3):
     rows=[r for r in j['rows'] if (r['message'],r['stripes'],r['rate'],r['replicate'])==(mi,stripe,rate,rep)]
     assert len({r['noisy_sha256'][0] for r in rows})==1
     assert [r['shared_positions'] for r in rows]==sorted(r['shared_positions'] for r in rows)
     assert rows[0]['shared_positions']==0 and rows[-1]['shared_positions']==rows[-1]['strand_nt']
     assert len(set(rows[-1]['noisy_sha256']))==1

def test_mixture_fixed_gate_geometry():
 s=['ACGT','TGCA','GTAC'];gate=[.1,.3,.6,.9]
 assert M.mix(s,gate,0)==s
 assert M.mix(s,gate,1)==[s[0]]*3
 assert M.mix(s,gate,.5)==['ACGT','ACCA','ACAC']
 with pytest.raises(ValueError):M.mix(s,[.1],.5)
 with pytest.raises(ValueError):M.mix(s,gate,2)
