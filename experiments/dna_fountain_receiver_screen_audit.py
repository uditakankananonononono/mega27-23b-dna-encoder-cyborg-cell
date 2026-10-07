"""Pinned receiver predicate geometry, not complete receiver or payload decode."""
import argparse,hashlib,json,subprocess
from pathlib import Path
from collections import Counter,defaultdict
import dna_fountain_bounded_read_audit as B
ROOT=Path(__file__).resolve().parents[1]
def screen(s):
 gc=(s.count('G')+s.count('C'))/len(s);homo=any(c*4 in s for c in 'ACGT');return {'pass':not homo and .45<=gc<=.55,'gc':gc,'homopolymer_fail':homo,'gc_fail':not .45<=gc<=.55}
def compute(prefix,upstream):
 p=ROOT/'research/dna_fountain_receiver_screen_plan.json';prior=json.loads((ROOT/'results/dna_fountain_bounded_read_audit.json').read_text());fresh=B.compute(prefix)
 if json.loads(json.dumps(fresh))!=prior:raise ValueError('prefix replay')
 if subprocess.check_output(['git','-C',upstream,'rev-parse','HEAD']).decode().strip()!='8ee2777aa5e9e101e5d756f7e2449f6672b66f1f':raise ValueError('source pin')
 aa=[B.read_prefix(Path(prefix)/f'mate{k}.gz')[0] for k in [1,2]];counts=Counter();rows=[];seen=set();payloads=defaultdict(set)
 for r in prior['RS_exact_records']:
  i=r['pair_index'];s,_=B.stitch(aa[0][i][1],aa[1][i][1]);data=B.bytes_from_dna(s);a=screen(s[:-8]);b=screen(s);duplicate=r['seed'] in seen;seen.add(r['seed']);payloads[r['seed']].add(r['payload_sha256']);counts['screen144_pass']+=a['pass'];counts['screen152_pass']+=b['pass'];counts['screen144_homopolymer_fail']+=a['homopolymer_fail'];counts['screen144_gc_fail']+=a['gc_fail'];counts['predicate_disagreement']+=a['pass']!=b['pass'];counts['duplicate_seed_pre_screen']+=duplicate;counts['unique_seed_screen144_pass']+=not duplicate and a['pass'];counts['unique_seed_screen144_fail']+=not duplicate and not a['pass'];rows.append({**r,'screen144':a,'screen152':b,'duplicate_seed_pre_screen':duplicate})
 files=['glass.pyx','droplet.pyx','utils.pyx'];j={'plan_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'prior_result_sha256':hashlib.sha256((ROOT/'results/dna_fountain_bounded_read_audit.json').read_bytes()).hexdigest(),'source_url':'https://github.com/TeamErlich/dna-fountain','source_commit':'8ee2777aa5e9e101e5d756f7e2449f6672b66f1f','source_file_sha256':{f:hashlib.sha256((Path(upstream)/f).read_bytes()).hexdigest() for f in files},'RS_exact_records':len(rows),'counts':dict(counts),'seed_payload_conflict_count':sum(len(s)>1 for s in payloads.values()),'records':rows,'decoding_performed':False,'limits':json.loads(p.read_bytes())['limits']};(ROOT/'results/dna_fountain_receiver_screen_audit.json').write_text(json.dumps(j,indent=2)+'\n');print({k:v for k,v in j.items() if k not in ['records','source_file_sha256']});return j
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--prefix',required=True);a.add_argument('--upstream',required=True);x=a.parse_args();compute(x.prefix,x.upstream)
