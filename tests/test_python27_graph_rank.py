import importlib.util,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('rank',ROOT/'experiments/python27_graph_rank_audit.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def test_gf2_rank_and_replay():
 assert m.rank([{'indices':[0,1]},{'indices':[1,2]},{'indices':[0,2]}])==2
 j=m.compute();assert j==json.loads((ROOT/'results/python27_graph_rank_audit.json').read_text());assert not j['payload_decoding_performed']
