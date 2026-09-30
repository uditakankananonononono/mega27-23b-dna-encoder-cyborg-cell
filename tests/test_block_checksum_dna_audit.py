import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
from dnacell import encoder as E
from block_checksum_dna_audit import classify

def test_single_dna_substitution_reaches_compensating_trit_collision():
 j=json.loads((ROOT/'results/block_checksum_dna_audit.json').read_text());r=j['examples']['1:payload_only']
 assert sum(a!=b for a,b in zip(r['original_dna'],r['changed_dna']))==1
 kind,payload=classify(r['changed_dna'],r['block'],r['original_trits'])
 assert kind=='accepted_wrong_payload'
 assert E._bns_try_block(r['changed_dna'],r['block'])==payload
 assert E.checksum_trits(payload)[-6:]==E.checksum_trits(r['original_trits'])[-6:]
 assert E.checksum_trits(payload)!=E.checksum_trits(r['original_trits'])

def test_complete_retained_enumeration_integrity():
 j=json.loads((ROOT/'results/block_checksum_dna_audit.json').read_text());raw=(ROOT/'results/block_checksum_dna_trials.jsonl').read_bytes()
 assert hashlib.sha256(raw).hexdigest()==j['trial_file_sha256']
 assert len(raw.splitlines())==j['total_trials']==51528
 for row in j['rows']:assert sum(row['counts'].values())==(114 if row['order']==1 else 6327)
 assert j['aggregate']['1:accepted_wrong_payload']==285
 assert j['aggregate']['2:accepted_wrong_payload']==5149
