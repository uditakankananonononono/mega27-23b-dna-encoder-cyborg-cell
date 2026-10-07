import importlib.util,json,random
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];s=importlib.util.spec_from_file_location('coord',ROOT/'experiments/python27_coordinate_identifiability_audit.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def test_no_singleton_linear_control():
 assert m.identifiable([{'indices':[0,1]},{'indices':[1,2]},{'indices':[0,1,2]}])==[0,1,2]
def test_small_bruteforce_rowspace():
 r=random.Random(99)
 for _ in range(20):
  rows=[{'indices':[i for i in range(6) if r.random()<.5]} for _ in range(5)];space={0}
  for row in rows:
   v=sum(1<<i for i in row['indices']);space|={x^v for x in list(space)}
  assert m.identifiable(rows)==[i for i in range(6) if 1<<i in space]
def test_recorded_identity():
 j=json.loads((ROOT/'results/python27_coordinate_identifiability_audit.json').read_text());assert not j['chunk_values_assigned'];assert all(r['identifiable_coordinates_count']==10 for r in j['scenarios'].values())
