"""Conditional graph structure only. Never reads or solves payload bytes."""
import json,hashlib
from pathlib import Path
from collections import Counter,defaultdict,deque
ROOT=Path(__file__).resolve().parents[1]
def graph(rows,K):
 edges=[set(r['indices']) for r in rows];covered=set().union(*edges);adj=defaultdict(set);parent=list(range(K))
 def find(x):
  while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
  return x
 for e,s in enumerate(edges):
  assert len(s)==rows[e]['degree'] and all(0<=v<K for v in s)
  first=next(iter(s))
  for v in s:adj[v].add(e);parent[find(v)]=find(first)
 comp=Counter(find(v) for v in range(K));queue=deque(i for i,s in enumerate(edges) if len(s)==1);resolved=set()
 while queue:
  i=queue.popleft()
  if len(edges[i])!=1:continue
  v=next(iter(edges[i]))
  if v in resolved:continue
  resolved.add(v)
  for k in list(adj[v]):
   edges[k].discard(v)
   if len(edges[k])==1:queue.append(k)
 return {'equations':len(rows),'covered_coordinates':len(covered),'uncovered_coordinates':K-len(covered),'degree_histogram':dict(sorted(Counter(r['degree'] for r in rows).items())),'initial_singletons':sum(r['degree']==1 for r in rows),'component_count_including_isolates':len(comp),'largest_component_size':max(comp.values()),'structurally_peeled_coordinates':len(resolved),'residual_nonempty_equations':sum(bool(s) for s in edges)}
def compute():
 p=ROOT/'research/python27_graph_coverage_plan.json';plan=json.loads(p.read_text());hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in plan['inputs']};v=json.loads((ROOT/plan['inputs'][0]).read_text())['vectors'];screen=json.loads((ROOT/plan['inputs'][1]).read_text());K=json.loads((ROOT/plan['inputs'][2]).read_text())['K'];passing={r['seed'] for r in screen['records'] if not r['duplicate_seed_pre_screen'] and r['screen144']['pass']}
 return {'plan_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'input_sha256':hashes,'K':K,'scenarios':{'first_seen144passing':graph([r for r in v if r['seed'] in passing],K),'all_RS_exact_unique':graph(v,K)},'payload_bytes_read':False,'payload_decoding_performed':False,'GF2_rank_measured':False,'limits':plan['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/python27_graph_coverage_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j)
