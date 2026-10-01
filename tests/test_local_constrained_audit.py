import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import local_constrained_audit as L
def test_clean_roundtrip_and_charged_constraints():
 p=json.loads((ROOT/'experiments/local_constrained_plan.json').read_text());rng=np.random.default_rng(p['seed']);j=json.loads((ROOT/'results/local_constrained_audit.json').read_text())
 for i,size in enumerate([n for n in p['sizes'] for _ in range(4)]):
  msg=rng.bytes(size);seq,ledger=L.encode(msg);g=j['geometry'][i];assert seq is not None and ledger==g['search_ledger']
  assert sum(x['attempts'] for x in ledger)==g['candidate_attempts']
  for d,s in seq.items():
   assert L.decode([s]*3,d)==msg and len(s)*3==g['metrics'][d]['total_nt']
   assert L.E.max_homopolymer(s)<=2 and .45<=L.E.gc_content(s)<=.55
   assert all(.35<=L.E.gc_content(s[k:k+52])<=.65 for k in range(0,len(s),52))
   assert L.S.gc_range(s,64)[0]>=.35 and L.S.gc_range(s,64)[1]<=.65
def test_outcome_denominators_and_negative():
 j=json.loads((ROOT/'results/local_constrained_audit.json').read_text());assert len(j['rows'])==96 and len(j['geometry'])==8
 assert all(r['included'] for r in j['geometry'])
 for r in j['summaries']:
  assert r['n_trials']==8 and sum(r['counts'].values())==8
  if r['p']==0:assert r['counts']=={'exact_success':8}
  elif r['bytes']==2048:assert r['counts']=={'loud_failure':8}
 assert json.loads((ROOT/'results/local_constrained_wrong.json').read_text())==[]
