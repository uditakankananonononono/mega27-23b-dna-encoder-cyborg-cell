"""Exact implementation accounting and checksum image, not ideal collision theory."""
import sys,json,itertools,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from dnacell import encoder as E

def main():
 rows=[]
 for size in [1,8,32,128,256,2048,16384]:
  msg=bytes(size);a=E.encode_message(msg);b=E.encode_message_v2(msg);expected_v1=3*(6*size+24);expected_v2=3*38*(int(np.ceil((6*size+24)/32))+1)
  assert sum(map(len,a))==expected_v1;assert sum(map(len,b))==expected_v2
  rows.append({'payload_bytes':size,'v1_total_nt_3copies':expected_v1,'v2_total_nt_3copies':expected_v2,'v1_rate':8*size/expected_v1,'v2_rate':8*size/expected_v2})
 # Last six trits are the second checksum byte, i.e. s1, encoded as six trits.
 patterns={tuple(E.checksum_trits(list(t))[-6:]) for t in itertools.product(range(3),repeat=8)}
 a=[0]*32;b=a[:];b[0]=1;b[1]=2 # lifted s1 sum changes3, not collision.
 c=a[:];c[0]=1;c[1]=0 # find equal sum example on nonzero trits instead.
 a=[1]*32;b=a[:];b[0]=0;b[1]=2
 assert a!=b and E.checksum_trits(a)[-6:]==E.checksum_trits(b)[-6:]
 out={'rows':rows,'capacity_no_adjacent_repeat_bits_per_nt':float(np.log2(3)),'byte_mapping_bits_per_payload_trit':8/6,'v1_asymptotic_three_copy_payload_rate':8/18,'v2_asymptotic_three_copy_payload_rate':(8/6)*(32/38)/3,'v2_block_checksum':'last six trits of two-byte checksum encode only low byte s1, not full Fletcher16','checksum_possible_values_upper_bound':255,'length8_checksum_observed_distinct_values':len(patterns),'explicit_undetected_two_trit_error':{'before':a,'after':b,'checksum':E.checksum_trits(a)[-6:]},'proof':'checksum_trits writes (s2<<8)|s1 as two bytes, six trits per byte. Suffix six trits discard s2 and preserve only s1 modulo255. Two equal-and-opposite trit changes preserve s1. Collision probability depends on error distribution; no universal1/729 or worst-case exponential bound.','limits':['No noisy DNA probability law inferred from a trit counterexample','Three-copy rates exclude primers and index metadata','Historical false formulas preserved in repository history, withdrawn in current paper']}
 (ROOT/'results/codec_rate_checksum_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
