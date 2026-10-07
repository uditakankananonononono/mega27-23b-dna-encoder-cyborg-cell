import hashlib,json,sys,zlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from dna_fountain_bounded_read_audit import stitch,reverse_complement,bytes_from_dna,read_prefix
def test_stitch_and_byte_geometry():
 s='ACGT'*38
 joined,n=stitch(s[:151],reverse_complement(s[-151:]))
 assert joined==s and n==150
 assert bytes_from_dna('ACGT')==bytes([27])
 assert stitch('A'*151,'A'*151)==(None,0)
def test_snapshot_conserves_counts():
 j=json.loads((ROOT/'results/dna_fountain_bounded_read_audit.json').read_text());c=j['counts']
 assert j['pairs_checked']==10000==c['stitched']+c['no_exact_overlap']
 assert c['length152']==c['length152_RS_exact']+c['length152_RS_rejected']+c['length152_nonACGT']
 assert c['length152_RS_exact']==len(j['RS_exact_records'])==4873
 assert j['unique_RS_exact_seeds']==j['unique_RS_exact_sequences']==4677
 assert j['plan_sha256']==hashlib.sha256((ROOT/'research/dna_fountain_bounded_read_plan.json').read_bytes()).hexdigest()
 assert not j['decoding_performed']
def test_incomplete_gzip_prefix_keeps_complete_records(tmp_path):
 p=tmp_path/'data.gz';import gzip
 p.write_bytes(gzip.compress(b'@read/1\nACGT\n+\nIIII\n@tail\n')[:-4]);r,m=read_prefix(p)
 assert len(r)==1 and m['terminal_lines']==1 and not m['gzip_stream_complete']
