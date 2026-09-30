import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import block_matched_geometry_audit as M

def test_matched_clean_geometry_and_no_test_wrapper_pollution():
 import block_joint_file_audit as baseline
 assert baseline.LENGTH=={'low':38,'full':44,'distance5':48}
 for d in ['low','full','distance5']:
  for size in [128,2048]:
   msg=bytes((17*i+31)%256 for i in range(size));s=M.F.encode(msg,3,d)
   assert M.F.decode(s,d)==msg
   assert len(s[0])==((6*size+24+31)//32+1)*48

def test_matched_ledger_plan_and_gates():
 j=json.loads((ROOT/'results/block_matched_geometry_audit.json').read_text())
 plan=(ROOT/'experiments/block_matched_plan.json').read_bytes();raw=(ROOT/'results/block_matched_geometry_trials.jsonl').read_bytes()
 assert hashlib.sha256(plan).hexdigest()==j['preregistered_plan_sha256']
 assert hashlib.sha256(raw).hexdigest()==j['ledger_sha256']
 rows=[json.loads(s) for s in raw.splitlines()]
 assert len(rows)==1440 and j['excluded_message_ids']==[]
 assert all(r['counts']=={'exact_success':24} for r in j['summaries'] if r['p']==0)
 # All candidates expose identical masks/counts, lengths and copy budget.
 groups={}
 for r in rows:groups.setdefault((r['message'],r['copies'],r['p'],r['replicate']),[]).append(r)
 assert all(len(g)==3 and len({(r['changed_nt'],r['total_nt']) for r in g})==1 for g in groups.values())
 wrong=json.loads((ROOT/'results/block_matched_geometry_wrong.json').read_text())
 assert len(wrong)==1
 for w in wrong:assert M.F.decode(w['noisy_strands'],w['detector']).hex()==w['actual_hex']!=w['expected_hex']
