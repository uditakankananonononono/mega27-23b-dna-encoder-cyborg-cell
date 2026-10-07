import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import bound_nonce_frame_audit as B
def test_bound_nonce_complete_audit():
 j=json.loads((ROOT/'results/bound_nonce_frame_audit.json').read_text());assert B.compute()==j
 assert j['frame_nt']==64 and j['field_code']==[13,9,5] and j['four_column_sets']==715
 assert j['prefix_three_counts']=={'rejected':972}
 assert all(r['counts']=={'1:rejected':192,'2:rejected':18144} for r in j['rows'])
def test_clean_all_nonces_and_explicit_nonce_binding():
 for n in range(81):
  p=[i%3 for i in range(32)];s=B.frame(p,n);assert len(s)==64 and B.decode(s)==p
