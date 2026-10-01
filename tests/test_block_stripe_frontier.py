import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import block_stripe_frontier_audit as S

def test_clean_and_single_stripe_erasure_recovery():
 msg=bytes(range(128))
 for detector in ['low','full','distance5']:
  for stripes in [1,4,16]:
   strand=S.encode(msg,1,detector,stripes)[0];assert S.decode([strand],detector,stripes)==msg
   # An invalid block is a detected erasure; one in each stripe remains repairable.
   damaged='A'*48*stripes+strand[48*stripes:]
   assert S.decode([damaged],detector,stripes)==msg
   damaged='A'*48*(stripes+1)+strand[48*(stripes+1):]
   import pytest
   with pytest.raises(ValueError):S.decode([damaged],detector,stripes)

def test_saved_grid_counts_costs_and_wrong_replay():
 j=json.loads((ROOT/'results/block_stripe_frontier_audit.json').read_text());p=(ROOT/'experiments/block_stripe_frontier_plan.json').read_bytes()
 assert j['plan_sha256']==hashlib.sha256(p).hexdigest()
 assert len(j['rows'])==216 and sum(r['n_trials'] for r in j['rows'])==3456 and not j['excluded_geometry']
 assert all(sum(r['counts'].values())==r['n_trials']==16 for r in j['rows'])
 assert all(r['counts']=={'exact_success':16} for r in j['rows'] if r['p']==0)
 for r in j['rows']:assert r['payload_bits_per_nt']==r['bytes']*8/r['total_nt_per_file']
 assert sum(not m['sliding_gc_gate'] for g in j['geometry'] for m in g['detectors'].values())==70
 wrong=json.loads((ROOT/'results/block_stripe_frontier_wrong.json').read_text());assert len(wrong)==2
 for r in wrong:assert S.decode(r['noisy_strands'],r['detector'],r['stripes']).hex()==r['actual_hex']!=r['expected_hex']

def test_chunk_merge_exact_reproduction():
 from merge_block_stripe_frontier import merge
 assert merge()==json.loads((ROOT/'results/block_stripe_frontier_audit.json').read_text())
