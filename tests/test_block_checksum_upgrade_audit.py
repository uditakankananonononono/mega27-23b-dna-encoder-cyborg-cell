import json,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
from dnacell import encoder as E
from block_checksum_upgrade_audit import classify

def test_double_dna_error_can_preserve_full_fletcher():
 j=json.loads((ROOT/'results/block_checksum_upgrade_audit.json').read_text());r=j['full_fletcher_misses']['2']
 assert sum(a!=b for a,b in zip(r['original_dna'],r['changed_dna']))==2
 assert r['original_payload']!=r['changed_payload']
 assert E.checksum_trits(r['original_payload'])==E.checksum_trits(r['changed_payload'])
 assert classify(r['changed_dna'],r['block'],r['original_payload'],True)=='accepted_wrong'

def test_paired_upgrade_trial_integrity_and_accounting():
 j=json.loads((ROOT/'results/block_checksum_upgrade_audit.json').read_text());raw=(ROOT/'results/block_checksum_upgrade_trials.jsonl').read_bytes()
 assert len(raw.splitlines())==j['total_paired_events']==36480
 assert hashlib.sha256(raw).hexdigest()==j['trial_sha256']
 singles=sum(n for k,n in j['counts'].items() if k.startswith('1:'))
 assert singles==768
 assert not any(k.startswith('1:') and k.endswith('|accepted_wrong') for k in j['counts'])
 assert j['counts']['2:accepted_wrong|accepted_wrong']==1653
 assert j['per_block_nt']=={'low_byte':38,'full_fletcher':44}
