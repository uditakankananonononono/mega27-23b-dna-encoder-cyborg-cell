import sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import protected_nonce_audit as P
@pytest.mark.parametrize('strands,detector,stripes',[
 ([], 'distance5',1),([''],'distance5',1),(['A'*60],'distance5',4),(['A'*120],'distance5',0),(['A'*120],'distance5',-1),(['A'*120],'distance5',True),(['A'*120],'distance5',1.5),(['?'*120],'distance5',1),(['A'*120],'missing',1),(['A'*120,None],'distance5',1),(['A'*120,'A'*60],'distance5',1),(['A'*121],'distance5',1)])
def test_malformed_inputs_are_loud_valueerrors(strands,detector,stripes):
 with pytest.raises(ValueError):P.decode(strands,detector,stripes)
def test_clean_empty_payload_is_not_empty_strand():
 seq,ledger=P.encode(b'',1);assert seq is not None
 for d,s in seq.items():assert P.decode([s]*3,d,1)==b''
