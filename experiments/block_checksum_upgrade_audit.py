"""Equal payload mutation library with low/full Fletcher block-check candidates.
Detection pilot only. Lengths differ, so not a matched-budget recovery result.
"""
import sys,json,itertools,hashlib
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from dnacell import encoder as E

def pack(payload,bi,full):
 checksum=E.checksum_trits(payload)
 return E.encode_trits_never_same(payload+(checksum if full else checksum[-6:]),E.BNS_SEEDS[bi%4])
def classify(seq,bi,original,full):
 try:trits=E.decode_never_same(seq,E.BNS_SEEDS[bi%4])
 except ValueError:return 'invalid_spacing'
 payload=trits[:32];checksum=E.checksum_trits(payload)
 if (checksum if full else checksum[-6:])!=trits[32:]:return 'checksum_reject'
 return 'accepted_wrong' if payload!=original else 'accepted_same'
def main():
 original=json.loads((ROOT/'results/block_checksum_dna_audit.json').read_text());msg=bytes.fromhex(original['message_hex']);dna=E.encode_message_v2(msg,1)[0]
 totals=Counter();rows=[];events=[];examples={}
 for bi in range(8):
  payload=E.decode_never_same(dna[bi*38:(bi+1)*38],E.BNS_SEEDS[bi%4])[:32];seqs={False:pack(payload,bi,False),True:pack(payload,bi,True)}
  assert seqs[False][:32]==seqs[True][:32]
  for order in [1,2]:
   count=Counter()
   for positions in itertools.combinations(range(32),order):
    opts=[[base for base in E.BASES if base!=seqs[False][pos]] for pos in positions]
    for repl in itertools.product(*opts):
     outcomes=[]
     for full in [False,True]:
      s=list(seqs[full])
      for pos,base in zip(positions,repl):s[pos]=base
      outcome=classify(''.join(s),bi,payload,full);outcomes.append(outcome)
      if full and outcome=='accepted_wrong' and str(order) not in examples:examples[str(order)]={'block':bi,'positions':list(positions),'replacements':''.join(repl),'original_dna':seqs[full],'changed_dna':''.join(s),'original_payload':payload,'changed_payload':E.decode_never_same(''.join(s),E.BNS_SEEDS[bi%4])[:32]}
     key='|'.join(outcomes);count[key]+=1;totals[f'{order}:{key}']+=1;events.append([bi,order,list(positions),''.join(repl),*outcomes])
   assert sum(count.values())==(96 if order==1 else 4464)
   rows.append({'block':bi,'order':order,'counts':dict(count)})
 raw=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in events);(ROOT/'results/block_checksum_upgrade_trials.jsonl').write_text(raw)
 out={'payload_source':'results/block_checksum_dna_audit.json','total_paired_events':len(events),'same_payload_prefix':True,'per_block_nt':{'low_byte':38,'full_fletcher':44},'length_increase_fraction':6/38,'asymptotic_three_copy_rate':{'low_byte':(8/6)*(32/38)/3,'full_fletcher':(8/6)*(32/44)/3},'counts':dict(totals),'rows':rows,'full_fletcher_misses':examples,'trial_sha256':hashlib.sha256(raw.encode()).hexdigest(),'limits':['Same exact payload position/letter mutation pairs for each candidate','Only payload positions mutate; check-region mutations and full-message effects not tested','Full Fletcher adds six check bases per block; not matched total budget or optimized density frontier','No consensus or physical read model; eight fixed block library','Candidate is audit-local, not a production encoder change']}
 (ROOT/'results/block_checksum_upgrade_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))
if __name__=='__main__':main()
