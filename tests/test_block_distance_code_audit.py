import hashlib,itertools,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import block_distance_code_audit as C

def test_field_and_all_four_column_vandermonde_minors():
 assert C.irreducible() and all(C.INV.values())
 for cols in itertools.combinations(range(12),4):
  # Nonzero determinant via product of pairwise coordinate differences.
  det=1
  for i,j in itertools.combinations(cols,2):det=C.MUL[det][C.ADD[j][C.neg(i)]]
  assert det!=0

def test_explicit_distance_five_witness_and_clean_roundtrip():
 j=json.loads((ROOT/'results/block_distance_code_audit.json').read_text())
 w=j['distance_upper_bound_witness']
 assert sum(x!=0 for x in w)==5 and C.syndrome(w)==[0]*4
 for bi in range(4):
  payload=[(i+bi)%3 for i in range(32)];dna=C.encode(payload,bi)
  assert len(dna)==48 and C.classify(dna,bi,payload)=='accepted_same'

def test_enumerated_whole_block_mutation_ledger():
 j=json.loads((ROOT/'results/block_distance_code_audit.json').read_text())
 raw=(ROOT/'results/block_distance_code_trials.jsonl').read_bytes()
 assert len(raw.splitlines())==j['total_events']==82368
 assert hashlib.sha256(raw).hexdigest()==j['ledger_sha256']
 assert sum(v for k,v in j['counts'].items() if k.startswith('1:'))==1152
 assert sum(v for k,v in j['counts'].items() if k.startswith('2:'))==81216
 assert not any('accepted' in k for k in j['counts'])
 assert j['max_valid_trit_changes']==j['max_valid_symbol_changes']==4
 assert j['block_nt']==48 and j['parity_nt']==16
