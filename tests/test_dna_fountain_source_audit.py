import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'))
from dna_fountain_source_audit import compute
def test_source_snapshot_replays():
 j=compute()
 assert j==json.loads((ROOT/'results/dna_fountain_source_audit.json').read_text())
 assert j['run_count']==16 and len({r['run_accession'] for r in j['runs']})==16
 assert j['compressed_fastq_bytes']==95237393879
 assert j['sum_ena_read_count']==479469649
 assert all(r['library_layout']=='PAIRED' for r in j['runs'])
 assert all(len(r['fastq_bytes'].split(';'))==2 for r in j['runs'])
 assert j['receiver']['padding_difference']==30208
 assert not j['decoding_performed'] and j['reads_downloaded']==0
