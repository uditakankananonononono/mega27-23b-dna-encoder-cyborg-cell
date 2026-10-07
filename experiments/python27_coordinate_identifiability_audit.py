"""Coefficient-only coordinate identifiability. Never assigns payload values."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def identifiable(rows):
 piv={}
 for r in rows:
  bits=0
  for i in r['indices']:bits^=1<<i
  while bits:
   p=bits.bit_length()-1
   if p not in piv:piv[p]=bits;break
   bits^=piv[p]
 order=sorted(piv)
 for n,p in enumerate(order):
  mask=1<<p
  for q in order[n+1:]:
   if piv[q]&mask:piv[q]^=piv[p]
 return [p for p in order if piv[p]==1<<p]
def compute():
 p=ROOT/'research/python27_coordinate_identifiability_plan.json';plan=json.loads(p.read_text());vectors=json.loads((ROOT/plan['inputs'][0]).read_text())['vectors'];s=json.loads((ROOT/plan['inputs'][1]).read_text());g=json.loads((ROOT/plan['inputs'][2]).read_text());passing={r['seed'] for r in s['records'] if not r['duplicate_seed_pre_screen'] and r['screen144']['pass']};out={}
 for name,rs in [('first_seen144passing',[r for r in vectors if r['seed'] in passing]),('all_RS_exact_unique',vectors)]:
  ids=identifiable(rs);peel=g['scenarios'][name]['structurally_peeled_coordinates'];out[name]={'equations':len(rs),'identifiable_coordinates_count':len(ids),'identifiable_coordinate_indices':ids,'structurally_peeled_coordinates':peel,'additional_identifiable_count_beyond_peeling':len(ids)-peel}
 return {'plan_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'input_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in plan['inputs']},'K':67088,'scenarios':out,'payload_bytes_read':False,'chunk_values_assigned':False,'limits':plan['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/python27_coordinate_identifiability_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j)
