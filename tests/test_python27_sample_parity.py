import importlib.util,json,random
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('sample27',ROOT/'experiments/python27_sample_parity_audit.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def test_all_archived_cpython27_sample_vectors():
 j=json.loads((ROOT/'results/python27_sample_input.json').read_text());vs=json.loads((ROOT/'results/python27_sample_vectors.json').read_text())['vectors'];assert len(vs)==4677
 for x in vs:
  r=random.Random(x['seed']);assert r.random()==x['first_uniform'];assert m.sample27(r,j['K'],x['degree'])==x['indices']
