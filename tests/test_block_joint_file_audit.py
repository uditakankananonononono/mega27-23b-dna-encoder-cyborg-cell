import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import block_joint_file_audit as F

def test_clean_whole_file_geometry_and_seams():
 for det,blen in F.LENGTH.items():
  for copies in [1,3]:
   for n in [0,1,128,2048]:
    msg=bytes((17*i+31)%256 for i in range(n));s=F.encode(msg,copies,det)
    assert F.decode(s,det)==msg
    assert len(s)==copies and len(s[0])==((6*n+24+31)//32+1)*blen

def test_whole_file_trials_and_negative_replays():
 j=json.loads((ROOT/'results/block_joint_file_audit.json').read_text());raw=(ROOT/'results/block_joint_file_trials.jsonl').read_bytes()
 assert hashlib.sha256(raw).hexdigest()==j['ledger_sha256']
 assert len(raw.splitlines())==1800
 assert all(r['counts']=={'exact_success':60} for r in j['rows'] if r['p']==0)
 wrong=json.loads((ROOT/'results/block_joint_file_silent_wrong.json').read_text())
 assert len(wrong)==2
 for w in wrong:
  out=F.decode(w['noisy_strands'],w['detector'])
  assert out.hex()==w['actual_hex']!=w['expected_hex']
 assert sum(r['counts'].get('silent_wrong',0) for r in j['rows'])==2

def test_unmodified_production_has_a_seeded_silent_wrong_counterexample():
 import os,subprocess
 env=dict(os.environ,PYTHONHASHSEED='0',PYTHONPATH=str(ROOT/'src'))
 code="import json;from dnacell.encoder import decode_message_v2;w=json.load(open('results/block_joint_file_silent_wrong.json'))[0];out=decode_message_v2(w['noisy_strands']);assert out.hex()==w['actual_hex']!=w['expected_hex']"
 subprocess.run([sys.executable,'-c',code],cwd=ROOT,env=env,check=True)
