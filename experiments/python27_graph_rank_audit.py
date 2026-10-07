"""Exact coefficient rank under conditional archived vectors; no payload bytes."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def rank(rows):
 piv={}
 for r in rows:
  bits=0
  for i in r['indices']:bits^=1<<i
  while bits:
   p=bits.bit_length()-1
   if p not in piv:piv[p]=bits;break
   bits^=piv[p]
 return len(piv)
def compute():
 p=ROOT/'research/python27_graph_rank_plan.json';plan=json.loads(p.read_text());v=json.loads((ROOT/plan['inputs'][0]).read_text())['vectors'];s=json.loads((ROOT/plan['inputs'][1]).read_text());passing={r['seed'] for r in s['records'] if not r['duplicate_seed_pre_screen'] and r['screen144']['pass']};scenarios={}
 for name,rows in [('first_seen144passing',[r for r in v if r['seed'] in passing]),('all_RS_exact_unique',v)]:
  n=rank(rows);scenarios[name]={'equations':len(rows),'exact_conditional_coefficient_rank':n,'dependent_equations':len(rows)-n,'coefficient_nullity':67088-n}
 return {'plan_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'input_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in plan['inputs']},'K':67088,'scenarios':scenarios,'payload_bytes_read':False,'payload_decoding_performed':False,'limits':plan['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/python27_graph_rank_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j)
