import json,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from dna_fountain_receiver_screen_audit import screen
def test_predicate_boundaries():
 assert screen('ACGT'*36)['pass']
 assert not screen('AAAA'+'CG'*70)['pass']
 assert not screen('AT'*72)['pass']
def test_snapshot_decision_conservation():
 j=json.loads((ROOT/'results/dna_fountain_receiver_screen_audit.json').read_text());r=j['records'];c=j['counts']
 assert j['plan_sha256']==hashlib.sha256((ROOT/'research/dna_fountain_receiver_screen_plan.json').read_bytes()).hexdigest()
 assert len(r)==j['RS_exact_records']==4873
 assert c['unique_seed_screen144_pass']+c['unique_seed_screen144_fail']+c['duplicate_seed_pre_screen']==4873
 assert c['screen144_pass']==sum(x['screen144']['pass'] for x in r)
 assert c['predicate_disagreement']==sum(x['screen144']['pass']!=x['screen152']['pass'] for x in r)
 assert not j['decoding_performed']
