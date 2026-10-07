import hashlib,json,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import bound_nonce_file_audit as F

def packed(msg,stripes):
 body=F.E.bytes_to_trits(msg);n=len(body);header=[(n//3**p)%3 for p in range(11,-1,-1)]
 t=F.E.scramble(header+body+F.E.checksum_trits(body))
 data=[(t[i:i+32]+[0]*32)[:32] for i in range(0,len(t),32)]
 parity=[[sum(data[b][k] for b in range(g,len(data),stripes))%3 for k in range(32)] for g in range(stripes)]
 return data+parity

def frames(blocks):return ''.join(F.B.frame(b,0) for b in blocks)

def test_saved_file_audit_reproduces():
 j=json.loads((ROOT/'results/bound_nonce_file_audit.json').read_text())
 assert F.compute()==j
 assert len(j['geometry'])==8 and len(j['trials'])==96 and not j['silent_wrong']
 assert all(c['outcome']==('exact_success' if c['case']=='single_data_erasure' else 'loud_failure') for c in j['controlled'])
 assert hashlib.sha256((ROOT/'experiments/bound_nonce_file_plan.json').read_bytes()).hexdigest()==j['plan_sha256']

@pytest.mark.parametrize('stripes',[1,4])
def test_clean_and_conflicting_valid_copies(stripes):
 s=frames(packed(b'abc',stripes));assert F.decode([s]*3,stripes)==b'abc'
 changed=F.B.frame([1]*32,0)+s[64:]
 with pytest.raises(ValueError,match='conflict'):F.decode([s,s,changed],stripes)

def test_parity_and_padding_are_validated():
 blocks=packed(b'abc',1);blocks[-1][0]=(blocks[-1][0]+1)%3
 with pytest.raises(ValueError,match='parity'):F.decode([frames(blocks)]*3,1)
 blocks=packed(b'abc',1);blocks[1][-1]=1
 blocks[-1]=[sum(b[k] for b in blocks[:-1])%3 for k in range(32)]
 with pytest.raises(ValueError,match='padding'):F.decode([frames(blocks)]*3,1)

def test_stripe_independent_erasure_recovery():
 s=frames(packed(bytes(range(20)),4))
 changed='A'*128+s[128:]
 assert F.decode([changed]*3,4)==bytes(range(20))

def test_missing_parity_with_intact_data():
 s=frames(packed(b'abc',1));assert F.decode([s[:-64]+'A'*64]*3,1)==b'abc'

@pytest.mark.parametrize('strands,stripes', [([],1),([''],1),(['A'*64],1),(['A'*128,'A'*127],1),(['N'*128],1),(['A'*128],True),(['A'*128],0),(['A'*128],1.0)])
def test_malformed_geometry(strands,stripes):
 with pytest.raises(ValueError):F.decode(strands,stripes)
