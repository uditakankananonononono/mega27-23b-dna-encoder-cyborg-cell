"""Exhaustive within-block strict one/two-DNA-substitution acceptance audit.
Finite enumerated library, not an independent biological substitution model.
"""
import hashlib,itertools,json,sys
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from dnacell import encoder as E

def classify(seq,bi,original):
 try:t=E.decode_never_same(seq,seed=E.BNS_SEEDS[bi%4])
 except ValueError:return 'invalid_spacing',None
 payload=t[:E.BNS_BLOCK]
 if E.checksum_trits(payload)[-E.BNS_CHK:]!=t[E.BNS_BLOCK:]:return 'checksum_reject',None
 return ('accepted_wrong_payload' if payload!=original else 'accepted_same_payload'),payload

def main():
 msg=bytes((17*i+31)%256 for i in range(128));dna=E.encode_message_v2(msg,copies=1)[0]
 aggregate=Counter();rows=[];examples={};records=[]
 for bi in range(8):
  block=dna[bi*38:(bi+1)*38];original=E.decode_never_same(block,seed=E.BNS_SEEDS[bi%4])[:32]
  for order in [1,2]:
   tally=Counter();distance=defaultdict(Counter);region=defaultdict(Counter)
   for positions in itertools.combinations(range(38),order):
    choices=[[base for base in E.BASES if base!=block[pos]] for pos in positions]
    for repl in itertools.product(*choices):
     changed=list(block)
     for pos,base in zip(positions,repl):changed[pos]=base
     outcome,payload=classify(''.join(changed),bi,original)
     tally[outcome]+=1;aggregate[f'{order}:{outcome}']+=1
     gap=positions[1]-positions[0] if order==2 else 0
     label='payload_only' if all(i<32 for i in positions) else ('checksum_only' if all(i>=32 for i in positions) else 'mixed')
     distance[gap][outcome]+=1;region[label][outcome]+=1
     records.append([bi,order,list(positions),''.join(repl),outcome])
     if outcome=='accepted_wrong_payload' and f'{order}:{label}' not in examples:
      examples[f'{order}:{label}']={'block':bi,'positions':list(positions),'replacement':''.join(repl),'original_dna':block,'changed_dna':''.join(changed),'original_trits':original,'changed_payload_trits':payload,'full_checksum_equal':E.checksum_trits(payload)==E.checksum_trits(original)}
   assert sum(tally.values())==(114 if order==1 else 6327)
   rows.append({'block':bi,'seed':E.BNS_SEEDS[bi%4],'order':order,'counts':dict(tally),'by_distance':{str(k):dict(v) for k,v in sorted(distance.items())},'by_region':{k:dict(v) for k,v in region.items()}})
 raw=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in records)
 (ROOT/'results/block_checksum_dna_trials.jsonl').write_text(raw)
 out={'message_hex':msg.hex(),'selected_blocks':list(range(8)),'per_block_nt':38,'total_trials':len(records),'trial_file_sha256':hashlib.sha256(raw.encode()).hexdigest(),'aggregate':dict(aggregate),'rows':rows,'examples':examples,'scope':['All strict single/two substitutions within each of eight fixed encoded payload blocks','Uniform enumeration over positions and three changed-letter choices, not a random physical read model','Only isolated block acceptance, not consensus or whole-file recovery; whole-message Fletcher and parity remain separate safeguards','The low-byte checksum is not improved or changed in this audit']}
 (ROOT/'results/block_checksum_dna_audit.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k] for k in ['total_trials','aggregate','trial_file_sha256','examples']},indent=2))
if __name__=='__main__':main()
