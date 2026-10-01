"""Merge six bounded chunks; compute only empirical tested-choice nondominance."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results'
def merge():
 parts=[json.loads((R/f'block_stripe_frontier_audit_{n}_{g}.json').read_text()) for n in [128,2048] for g in [1,4,16]]
 assert len({x['plan_sha256'] for x in parts})==1
 rows=sum([x['rows'] for x in parts],[]);geo=parts[0]['geometry'];assert all(x['geometry']==geo for x in parts)
 front=[]
 for r in rows:
  choices=[z for z in rows if z['bytes']==r['bytes'] and z['p']==r['p'] and z['n_trials']==r['n_trials'] and z['n_trials']]
  if not any(z['total_nt_per_file']<=r['total_nt_per_file'] and z['counts'].get('exact_success',0)>=r['counts'].get('exact_success',0) and z['counts'].get('silent_wrong',0)<=r['counts'].get('silent_wrong',0) and (z['total_nt_per_file']<r['total_nt_per_file'] or z['counts'].get('exact_success',0)>r['counts'].get('exact_success',0) or z['counts'].get('silent_wrong',0)<r['counts'].get('silent_wrong',0)) for z in choices):front.append(r)
 wrong=sum([json.loads((R/f'block_stripe_frontier_wrong_{n}_{g}.json').read_text()) for n in [128,2048] for g in [1,4,16]],[])
 out={k:v for k,v in parts[0].items() if k not in ['rows','empirical_nondominated_tested_choices']};out.update(rows=rows,empirical_nondominated_tested_choices=front,complete_chunk_count=6)
 (R/'block_stripe_frontier_audit.json').write_text(json.dumps(out,indent=2)+'\n');(R/'block_stripe_frontier_wrong.json').write_text(json.dumps(wrong,indent=2)+'\n')
 return out
if __name__=='__main__':merge()
