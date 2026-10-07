"""Frozen strict-overlap read-prefix audit; no PEAR, peeling or source reconstruction."""
import argparse,hashlib,json,zlib,sys
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def read_prefix(path):
 d=zlib.decompressobj(16+zlib.MAX_WBITS);raw=d.decompress(Path(path).read_bytes());lines=raw.decode('ascii').splitlines();n=len(lines)//4
 records=[]
 for i in range(n):
  name,seq,plus,qual=lines[4*i:4*i+4]
  if not name.startswith('@') or not plus.startswith('+') or len(seq)!=len(qual):raise ValueError('FASTQ geometry')
  records.append((name,seq,qual))
 return records,{'complete_records':n,'terminal_lines':len(lines)%4,'gzip_stream_complete':d.eof,'decompressed_bytes':len(raw)}
def reverse_complement(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def stitch(a,b):
 b=reverse_complement(b)
 for n in range(min(len(a),len(b)),19,-1):
  if a[-n:]==b[:n]:return a+b[n:],n
 return None,0
def bytes_from_dna(s):
 if len(s)%4 or set(s)-set('ACGT'):raise ValueError('non-ACGT or incomplete byte')
 tab={c:i for i,c in enumerate('ACGT')}
 return bytes(sum(tab[c]<<(6-2*k) for k,c in enumerate(s[i:i+4])) for i in range(0,len(s),4))
def compute(prefix):
 from reedsolo import RSCodec
 p=ROOT/'research/dna_fountain_bounded_read_plan.json';plan=json.loads(p.read_bytes());prefix=Path(prefix);ledger=json.loads((prefix/'ledger.json').read_text());rs=RSCodec(2)
 aa=[];meta=[]
 for mate,item in enumerate(ledger,1):
  f=prefix/f'mate{mate}.gz';b=f.read_bytes()
  if len(b)!=8388608 or hashlib.sha256(b).hexdigest()!=item['sha256']:raise ValueError('range digest')
  r,m=read_prefix(f);aa.append(r);meta.append(m)
 n=min(plan['record_budget'],*(len(a) for a in aa));hist=Counter();overlap=Counter();cnt=Counter();selected=[];lengths=[Counter(),Counter()]
 for i,(a,b) in enumerate(zip(aa[0][:n],aa[1][:n])):
  ida=a[0].split()[0].removesuffix('/1');idb=b[0].split()[0].removesuffix('/2')
  if ida!=idb:raise ValueError('pair-ID mismatch')
  for k,r in enumerate([a,b]):lengths[k][len(r[1])]+=1
  cnt['pairs_with_N']+=int('N' in a[1] or 'N' in b[1]);s,o=stitch(a[1],b[1])
  if s is None:cnt['no_exact_overlap']+=1;continue
  cnt['stitched']+=1;hist[len(s)]+=1;overlap[o]+=1
  if len(s)!=152:continue
  cnt['length152']+=1
  if set(s)-set('ACGT'):cnt['length152_nonACGT']+=1;continue
  data=bytes_from_dna(s)
  try:
   decoded=bytes(rs.decode(data)[0]);same=bytes(rs.encode(decoded))==data
  except Exception:same=False
  cnt['length152_RS_exact' if same else 'length152_RS_rejected']+=1
  if same:selected.append({'pair_index':i,'read_id':ida,'sequence_sha256':hashlib.sha256(s.encode()).hexdigest(),'seed':int.from_bytes(decoded[:4],'big'),'payload_sha256':hashlib.sha256(decoded[4:]).hexdigest()})
 return {'plan_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'source_run':'ERR1816980','prefix_ledger':ledger,'fastq_metadata':meta,'pairs_checked':n,'mate_length_histograms':[dict(sorted(c.items())) for c in lengths],'counts':dict(cnt),'stitched_length_histogram':dict(sorted(hist.items())),'overlap_length_histogram':dict(sorted(overlap.items())),'RS_exact_records':selected,'unique_RS_exact_seeds':len({r['seed'] for r in selected}),'unique_RS_exact_sequences':len({r['sequence_sha256'] for r in selected}),'reedsolo_version':'1.7.0','RS_settings':'RSCodec(2) defaults,zero tolerated corrected-byte differences; no source payload identity inferred','decoding_performed':False,'limits':plan['limits']}
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--prefix',required=True);x=a.parse_args();j=compute(x.prefix);(ROOT/'results/dna_fountain_bounded_read_audit.json').write_text(json.dumps(j,indent=2)+'\n');print({k:v for k,v in j.items() if k not in ['RS_exact_records','overlap_length_histogram']})
