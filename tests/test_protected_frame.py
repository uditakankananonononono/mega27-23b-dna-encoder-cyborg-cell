import hashlib,json,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import protected_frame_audit as F
def test_full_frame_saved_enumeration():
 j=json.loads((ROOT/'results/protected_frame_audit.json').read_text());assert F.compute()==j and len(j['rows'])==4
 assert j['plan_sha256']==hashlib.sha256((ROOT/'experiments/protected_frame_plan.json').read_bytes()).hexdigest()
 for r in j['rows']:
  assert r['counts']=={'1:rejected':180,'2:rejected':15930}
  assert sum(r['partitions'].values())==16110 and r['partitions']['2:mixed:rejected']==5184
def test_clean_frame_and_geometry():
 for nonce in range(81):
  p=[i%3 for i in range(32)];s=F.frame(p,nonce);assert F.decode_frame(s)==p
  with pytest.raises(ValueError):F.decode_frame(s[:-1])
