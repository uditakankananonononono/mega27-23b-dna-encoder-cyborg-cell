import importlib.util,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('g',ROOT/'experiments/python27_graph_coverage_audit.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def test_conditional_graph_replay():
 j=m.compute();assert json.loads(json.dumps(j))==json.loads((ROOT/'results/python27_graph_coverage_audit.json').read_text());assert not j['payload_decoding_performed'];assert not j['GF2_rank_measured']
 for r in j['scenarios'].values():assert r['covered_coordinates']+r['uncovered_coordinates']==j['K'];assert r['structurally_peeled_coordinates']<=r['equations']
