"""Audit-local systematic [12,8,5] GF(81) detector, not production codec.
Four Vandermonde parity equations; each field symbol occupies four trits.
"""
import hashlib,itertools,json,sys
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from dnacell import encoder as E

def digits(x):return [(x//3**i)%3 for i in range(4)]
def number(a):return sum(x*3**i for i,x in enumerate(a))
def add(a,b):return number([(x+y)%3 for x,y in zip(digits(a),digits(b))])
def neg(a):return number([(-x)%3 for x in digits(a)])
# x^4+x+2. Irreducibility is verified independently below, not assumed.
MOD=[2,1,0,0,1]
def mul(a,b):
 c=[0]*7
 for i,x in enumerate(digits(a)):
  for j,y in enumerate(digits(b)):c[i+j]=(c[i+j]+x*y)%3
 for k in range(6,3,-1):
  q=c[k]
  for j in range(5):c[k-4+j]=(c[k-4+j]-q*MOD[j])%3
 return number(c[:4])
ADD=[[add(a,b) for b in range(81)] for a in range(81)]
MUL=[[mul(a,b) for b in range(81)] for a in range(81)]
INV={a:next((b for b in range(1,81) if MUL[a][b]==1),None) for a in range(1,81)}
def remainder(poly,divisor):
 p=poly[:]
 for k in range(len(p)-1,len(divisor)-2,-1):
  q=p[k]
  for j,c in enumerate(divisor):p[k-len(divisor)+1+j]=(p[k-len(divisor)+1+j]-q*c)%3
 return p[:len(divisor)-1]
def irreducible():
 return all(any(remainder(MOD,list(c)+[1])) for n in [1,2] for c in itertools.product(range(3),repeat=n))
def power(x,n):
 out=1
 for _ in range(n):out=MUL[out][x]
 return out
H=[[power(i,r) for i in range(12)] for r in range(4)]
def solve(matrix,rhs):
 a=[row[:]+[b] for row,b in zip(matrix,rhs)];n=len(a)
 for i in range(n):
  k=next(k for k in range(i,n) if a[k][i]);a[i],a[k]=a[k],a[i]
  inv=INV[a[i][i]];a[i]=[MUL[v][inv] for v in a[i]]
  for k in range(n):
   if k!=i:
    q=a[k][i];a[k]=[ADD[x][neg(MUL[q][y])] for x,y in zip(a[k],a[i])]
 return [row[-1] for row in a]
def syndrome(symbols):
 out=[]
 for row in H:
  s=0
  for h,x in zip(row,symbols):s=ADD[s][MUL[h][x]]
  out.append(s)
 return out
def symbols(trits):return [number(trits[i:i+4]) for i in range(0,len(trits),4)]
def encode(payload,bi):
 assert len(payload)==32
 u=symbols(payload);rhs=[neg(x) for x in syndrome(u+[0]*4)]
 parity=solve([row[8:] for row in H],rhs)
 return E.encode_trits_never_same(payload+sum([digits(x) for x in parity],[]),E.BNS_SEEDS[bi%4])
def classify(seq,bi,original):
 try:t=E.decode_never_same(seq,E.BNS_SEEDS[bi%4])
 except ValueError:return 'invalid_spacing'
 if len(t)!=48:raise ValueError('wrong block length')
 if any(syndrome(symbols(t))):return 'syndrome_reject'
 return 'accepted_same' if t[:32]==original else 'accepted_wrong'
def main():
 assert irreducible() and all(INV.values())
 # Every <=4-column dependence is impossible by Vandermonde determinant.
 # Compute a weight-five dependence as an explicit upper-bound witness.
 w=solve([row[:4] for row in H],[neg(row[4]) for row in H])+[1]+[0]*7
 assert sum(bool(x) for x in w)==5 and not any(syndrome(w))
 source=json.loads((ROOT/'results/block_checksum_dna_audit.json').read_text())
 dna=E.encode_message_v2(bytes.fromhex(source['message_hex']),1)[0]
 rows=[];totals=Counter();records=[];max_trits=0;max_symbols=0
 for bi in range(8):
  payload=E.decode_never_same(dna[bi*38:(bi+1)*38],E.BNS_SEEDS[bi%4])[:32]
  seq=encode(payload,bi);original=E.decode_never_same(seq,E.BNS_SEEDS[bi%4])
  assert classify(seq,bi,payload)=='accepted_same'
  for order in [1,2]:
   tally=Counter();region=Counter()
   for positions in itertools.combinations(range(48),order):
    opts=[[base for base in E.BASES if base!=seq[p]] for p in positions]
    for repl in itertools.product(*opts):
     s=list(seq)
     for p,b in zip(positions,repl):s[p]=b
     changed=''.join(s);outcome=classify(changed,bi,payload)
     tally[outcome]+=1;totals[f'{order}:{outcome}']+=1
     label='payload' if max(positions)<32 else ('parity' if min(positions)>=32 else 'mixed')
     region[f'{label}:{outcome}']+=1
     if outcome!='invalid_spacing':
      t=E.decode_never_same(changed,E.BNS_SEEDS[bi%4]);dt=sum(x!=y for x,y in zip(t,original));ds=sum(x!=y for x,y in zip(symbols(t),symbols(original)))
      assert 1<=ds<=dt<=2*order
      max_trits=max(max_trits,dt);max_symbols=max(max_symbols,ds)
     records.append([bi,order,list(positions),''.join(repl),outcome])
   assert sum(tally.values())==(144 if order==1 else 10152)
   rows.append({'block':bi,'order':order,'counts':dict(tally),'regions':dict(region)})
 raw=''.join(json.dumps(r,separators=(',',':'))+'\n' for r in records)
 (ROOT/'results/block_distance_code_trials.jsonl').write_text(raw)
 out={'field_modulus_low_to_high':MOD,'irreducibility_verified':True,'coordinates':list(range(12)),'field_code_n_k_distance':[12,8,5],'distance_upper_bound_witness':w,'parity_nt':16,'payload_nt':32,'block_nt':48,'three_copy_asymptotic_payload_bits_per_nt':(8/6)*(32/48)/3,'relative_length_vs_38_nt':48/38,'relative_length_vs_44_nt':48/44,'total_events':len(records),'counts':dict(totals),'rows':rows,'max_valid_trit_changes':max_trits,'max_valid_symbol_changes':max_symbols,'ledger_sha256':hashlib.sha256(raw.encode()).hexdigest(),'limits':['Eight fixed original blocks; all positions including parity region','Distance proof covers at most two strict substitutions per isolated seeded block, not indels','Detection only, no correction/consensus/outer parity/full-file benchmark','Neither matched-rate superiority nor physical validation; production codec unchanged']}
 (ROOT/'results/block_distance_code_audit.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))
if __name__=='__main__':main()
