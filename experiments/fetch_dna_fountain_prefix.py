"""Retrieve only preregistered16MiB ERR1816980 compressed prefix; no full reads."""
import argparse,urllib.request,json,hashlib
from pathlib import Path
a=argparse.ArgumentParser();a.add_argument('--out',required=True);x=a.parse_args();out=Path(x.out);out.mkdir(exist_ok=True,parents=True);ledger=[]
for mate in [1,2]:
 u=f'https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR181/000/ERR1816980/ERR1816980_{mate}.fastq.gz'
 r=urllib.request.urlopen(urllib.request.Request(u,headers={'Range':'bytes=0-8388607'}),timeout=45)
 if r.status!=206 or not r.headers.get('Content-Range','').startswith('bytes 0-8388607/'):raise ValueError('not requested range')
 b=r.read(8388609)
 if len(b)!=8388608:raise ValueError('byte budget')
 (out/f'mate{mate}.gz').write_bytes(b);ledger.append({'url':u,'status':r.status,'content_range':r.headers['Content-Range'],'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
(out/'ledger.json').write_text(json.dumps(ledger,indent=2)+'\n')
