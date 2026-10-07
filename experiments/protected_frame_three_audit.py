"""Three-change adversarial prefix substitution sharpness audit."""
import hashlib,json,sys
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import protected_frame_audit as F
def prefix(n):return F.P.E.encode_trits_never_same([(n//3**p)%3 for p in range(3,-1,-1)],'A')*3
def compute():
 raw=(ROOT/'experiments/protected_frame_three_plan.json').read_bytes();pairs=[(a,b) for a in range(81) for b in range(81) if sum(x!=y for x,y in zip(prefix(a),prefix(b)))==3];rows=[]
 for payload_id,p in enumerate([[0]*32,[i%3 for i in range(32)]]):
  for a,b in pairs:
   original=F.frame(p,a);changed=prefix(b)+original[12:];positions=[i for i,(x,y) in enumerate(zip(original,changed)) if x!=y];assert len(positions)==3
   try:out=F.decode_frame(changed);outcome='accepted_unchanged' if out==p else 'accepted_wrong'
   except ValueError:out=None;outcome='rejected'
   rows.append({'payload_id':payload_id,'original_nonce':a,'changed_nonce':b,'positions':positions,'outcome':outcome,'original_frame_sha256':hashlib.sha256(original.encode()).hexdigest(),'changed_frame_sha256':hashlib.sha256(changed.encode()).hexdigest(),'decoded_payload':out})
 text=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in rows);return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'nonce_pair_count':len(pairs),'event_count':len(rows),'counts':dict(Counter(r['outcome'] for r in rows)),'ledger_sha256':hashlib.sha256(text.encode()).hexdigest(),'limits':json.loads(raw)['limits']},text
if __name__=='__main__':
 j,text=compute();(ROOT/'results/protected_frame_three_audit.json').write_text(json.dumps(j,indent=2)+'\n');(ROOT/'results/protected_frame_three.jsonl').write_text(text);print(j)
