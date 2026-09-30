import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from dnacell import encoder as E
def test_low_byte_checksum_can_miss_compensating_trit_change():
 a=[1]*32;b=a[:];b[0]=0;b[1]=2
 assert a!=b;assert E.checksum_trits(a)[-6:]==E.checksum_trits(b)[-6:]
 assert E.checksum_trits(a)!=E.checksum_trits(b)
def test_emitted_lengths_charge_byte_packing_and_copies():
 for L in [1,32,128,2048]:
  assert sum(map(len,E.encode_message(bytes(L))))==3*(6*L+24)
  assert sum(map(len,E.encode_message_v2(bytes(L))))==3*38*((6*L+24+31)//32+1)
