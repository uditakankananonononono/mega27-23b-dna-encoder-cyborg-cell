import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import protected_frame_three_audit as T
def test_three_change_sharpness_complete_ledger():
 j=json.loads((ROOT/'results/protected_frame_three_audit.json').read_text());new,raw=T.compute();assert new==j and raw==(ROOT/'results/protected_frame_three.jsonl').read_text()
 assert j['nonce_pair_count']==486 and j['event_count']==972 and j['counts']=={'accepted_wrong':972}
 for r in map(json.loads,raw.splitlines()):
  p=[[0]*32,[i%3 for i in range(32)]][r['payload_id']];s=T.F.frame(p,r['original_nonce']);changed=T.prefix(r['changed_nonce'])+s[12:]
  assert sum(x!=y for x,y in zip(s,changed))==3 and all(i<12 for i in r['positions'])
  assert T.F.decode_frame(changed)==r['decoded_payload']!=p
