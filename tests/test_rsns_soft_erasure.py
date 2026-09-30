from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'))
from rs_inner import descramble_soft,trits_to_syms,rsns_encode,rsns_decode
from dnacell.encoder import scramble

def test_invalid_marker_survives():
 raw=[3,0,0,0,0]
 assert descramble_soft(raw)[0]==3
 assert trits_to_syms(descramble_soft(raw))[1]==[True]

def test_valid_descramble_roundtrip():
 valid=[0,1,2]*80
 assert descramble_soft(scramble(valid))==valid

def test_rsns_clean_roundtrips():
 for n in [1,25,128,256]:
  msg=bytes(i%256 for i in range(n));assert rsns_decode(rsns_encode(msg))==msg
