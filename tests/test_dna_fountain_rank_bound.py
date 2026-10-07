import importlib.util,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('bound',ROOT/'experiments/dna_fountain_rank_bound_audit.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def test_mapping_free_rank_bounds():
 j=m.compute();assert j==json.loads((ROOT/'results/dna_fountain_rank_bound_audit.json').read_text())
 assert [r['coefficient_rank_upper_bound'] for r in j['scenarios']]==[4595,4595,4677]
 assert j['scenarios'][0]['coefficient_nullity_lower_bound']==62493
 assert not j['actual_rank_measured'] and not j['decoding_performed']
