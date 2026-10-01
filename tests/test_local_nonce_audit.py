import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from local_nonce_audit import compute
def test_exhaustive_first_nonce_grid():
 j=json.loads((ROOT/'results/local_nonce_audit.json').read_text());assert compute()==j
 assert len(j['rows'])==288 and j['wrong']==[]
 assert len({(r['message'],r['detector'],r['position'],r['replacement']) for r in j['rows']})==288
 for r in j['summary']:
  assert r['n_events']==96 and r['core_counts']=={'prefix_or_core_spacing_reject':56,'core_accept_wrongpayload':24,'core_reject':16}
  assert r['file_counts']=={'wholefile_exact':72,'wholefile_loud':24}
