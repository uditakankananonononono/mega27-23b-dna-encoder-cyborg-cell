"""Mapping-free equation rank ceiling, not actual rank, coverage, or decoding."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def compute():
 pth=ROOT/'research/dna_fountain_rank_bound_plan.json';plan=json.loads(pth.read_text());hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in plan['inputs']};r=json.loads((ROOT/plan['inputs'][0]).read_text());s=json.loads((ROOT/plan['inputs'][1]).read_text());K=s['receiver']['readme_parameters']['chunks'];rows=r['records'];seen=set();accepted=set();all_pass=set()
 for x in rows:
  if x['screen144']['pass']:all_pass.add(x['seed'])
  if x['seed'] not in seen and x['screen144']['pass']:accepted.add(x['seed'])
  seen.add(x['seed'])
 assert len(accepted)==r['counts']['unique_seed_screen144_pass'];assert r['seed_payload_conflict_count']==0
 scenarios=[]
 for name,m in [('first_seen_144nt_predicate',len(accepted)),('any_144nt_passing_unique_seed',len(all_pass)),('unscreened_RS_exact_unique_seed',len(seen))]:
  scenarios.append({'scenario':name,'unique_equations_upper_count':m,'coefficient_rank_upper_bound':min(K,m),'coefficient_nullity_lower_bound':max(0,K-m),'rank_fraction_upper_bound':min(K,m)/K,'additional_independent_equations_necessary_at_least':max(0,K-m),'payload_bits_unresolved_lower_bound':max(0,K-m)*32*8})
 return {'plan_sha256':hashlib.sha256(pth.read_bytes()).hexdigest(),'input_sha256':hashes,'conditional_source_chunks':K,'chunk_bytes':32,'scenarios':scenarios,'coverage_measured':False,'coordinate_coverage_upper_bound':K,'actual_rank_measured':False,'peeling_performed':False,'decoding_performed':False,'limits':plan['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/dna_fountain_rank_bound_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j)
