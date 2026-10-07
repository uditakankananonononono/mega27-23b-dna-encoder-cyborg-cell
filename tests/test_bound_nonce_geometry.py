import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import bound_nonce_geometry_audit as G
def test_prospective_geometry_complete_and_charged():
 j=json.loads((ROOT/'results/bound_nonce_geometry_audit.json').read_text());assert G.compute()==j and len(j['rows'])==8
 for r in j['rows']:
  assert r['included'] and r['total_nt_three_copies']==3*r['strand_nt']==192*len(r['search_ledger'])
  assert r['max_homopolymer']<=2 and .45<=r['whole_gc']<=.55
  assert r['sliding64_gc_range'][0]>=.35 and r['sliding64_gc_range'][1]<=.65
  assert r['frame_gc_range'][0]>=.35 and r['frame_gc_range'][1]<=.65
  assert all(0<=x['nonce']<81 and x['attempts']==x['nonce']+1 for x in r['search_ledger'])
