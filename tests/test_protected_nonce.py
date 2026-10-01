import itertools,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import protected_nonce_audit as P
def test_prefix_two_substitution_detection_exhaustive():
 events=0
 for nonce in range(81):
  t=[(nonce//3**p)%3 for p in range(3,-1,-1)];prefix=P.E.encode_trits_never_same(t,'A')*3
  for order in [1,2]:
   for positions in itertools.combinations(range(12),order):
    for letters in itertools.product(*[[b for b in P.E.BASES if b!=prefix[i]] for i in positions]):
     s=list(prefix)
     for i,b in zip(positions,letters):s[i]=b
     s=''.join(s);events+=1
     try:z=[P.E.decode_never_same(s[k:k+4],'A') for k in [0,4,8]]
     except ValueError:continue
     assert not z[0]==z[1]==z[2] or z[0]==t
 assert events==51030
def test_prospective_geometry_and_clean_recovery():
 j=json.loads((ROOT/'results/protected_nonce_audit.json').read_text());p=json.loads((ROOT/'experiments/protected_nonce_plan.json').read_text());rng=np.random.default_rng(p['seed'])
 for i,msg in enumerate([rng.bytes(n) for n in p['sizes'] for _ in range(2)]):
  for stripes in [1,4]:
   seq,ledger=P.encode(msg,stripes);g=next(g for g in j['geometry'] if (g['message'],g['stripes'])==(i,stripes));assert g['included'] and ledger==g['search_ledger']
   for d,s in seq.items():
    assert P.decode([s]*3,d,stripes)==msg and len(s)*3==g['metrics'][d]['total_nt']
    assert all(.35<=P.E.gc_content(s[k:k+60])<=.65 for k in range(0,len(s),60))
    assert P.S.gc_range(s,64)[0]>=.35 and P.S.gc_range(s,64)[1]<=.65
 assert len(j['rows'])==96 and j['wrong']==[]
 assert all(r['counts']=={'loud_failure':4} for r in j['summaries'] if r['bytes']==2048 and r['p']==.03)
