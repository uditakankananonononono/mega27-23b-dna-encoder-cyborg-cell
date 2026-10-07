import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import bound_nonce_shared_noise_audit as S

def test_shared_copy_grid_reproduces_and_preserves_negative():
 j=json.loads((ROOT/'results/bound_nonce_shared_noise_audit.json').read_text());assert S.compute()==j
 assert len(j['rows'])==72 and not j['silent_wrong']
 assert all(r['shared_outcome']=='loud_failure' and r['error']=='unrecoverable stripe' for r in j['rows'])
 assert sum(r['independent_outcome']=='exact_success' for r in j['rows'])==49
 assert all(r['changed_nt_over_copies']%3==0 and r['total_nt']%192==0 for r in j['rows'])
