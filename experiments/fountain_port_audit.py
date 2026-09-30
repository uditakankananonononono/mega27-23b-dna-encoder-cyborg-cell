"""Pinned external Python-3 port versus repaired local RSNS, strict channel.
Requires a separately installed unmodified port; one droplet read per encoded oligo.
"""
import json,sys,subprocess,argparse,hashlib,tempfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(Path(__file__).parent));sys.path.insert(0,str(ROOT/'src'))
from strict_substitution_benchmark import strict
from rs_inner import rsns_encode,rsns_decode

def run(port):
 port=Path(port).resolve();commit=subprocess.check_output(['git','-C',str(port),'rev-parse','HEAD']).decode().strip();assert commit=='f97c1b8a81f5c5b819209d5b5e26c9c8f4439495'
 rng=np.random.default_rng(511);rows=[];size=2048;messages=4
 with tempfile.TemporaryDirectory() as work:
  w=Path(work)
  for i in range(messages):
   msg=rng.bytes(size);inp=w/f'input_{i}.bin';inp.write_bytes(msg);dna=w/f'clean_{i}.fasta'
   cmd=[sys.executable,'encode.py','--file_in',str(inp),'--size','32','-m','3','--gc','.05','--rs','2','--delta','.001','--c_dist','.025','--out',str(dna),'--stop','180'];proc=subprocess.run(cmd,cwd=port,capture_output=True,text=True,timeout=20);assert proc.returncode==0,proc.stderr
   seqs=[s.strip() for s in dna.read_text().splitlines() if not s.startswith('>') and s.strip()];local=rsns_encode(msg);assert len(seqs)==180
   for rate in [0,.005,.01,.02,.03]:
    for seed in range(3):
     changed=0;noisy=[]
     for c,s in enumerate(seqs):z,n=strict(s,rate,91000+i*1000+seed*300+c);noisy.append(z);changed+=n
     nf=w/'noisy.fasta';nf.write_text(''.join(f'>droplet_{k}\n{s}\n' for k,s in enumerate(noisy)));out=w/'recover.bin'
     if out.exists():out.unlink()
     cmd=[sys.executable,'decode.py','--file_in',str(nf),'--fasta','--chunk_num','64','--size','32','--rs','2','--delta','.001','--c_dist','.025','--out',str(out),'--gc','.05','-m','3'];proc=subprocess.run(cmd,cwd=port,capture_output=True,text=True,timeout=20);ok=out.exists() and out.read_bytes()==msg
     row={'codec':'external_fountain_py3_port','message':i,'noise_seed':seed,'p':rate,'recovered':ok,'total_encoded_nt':sum(map(len,seqs)),'payload_bits_per_encoded_nt':8*size/sum(map(len,seqs)),'realized_substitution_fraction':changed/sum(map(len,seqs)),'process_exit':proc.returncode,'input_sha256':hashlib.sha256(msg).hexdigest()};rows.append(row)
     z,n=strict(local,rate,92000+i*1000+seed);err=None
     try:ok=rsns_decode(z)==msg
     except Exception as e:ok=False;err=type(e).__name__+': '+str(e)[:120]
     rows.append({'codec':'repaired_rsns','message':i,'noise_seed':seed,'p':rate,'recovered':ok,'total_encoded_nt':len(local),'payload_bits_per_encoded_nt':8*size/len(local),'realized_substitution_fraction':n/len(local),'error':err,'input_sha256':hashlib.sha256(msg).hexdigest()})
   print('message',i,'complete',flush=True)
 out={'source':'https://github.com/jdbrody/dna-fountain','upstream_original':'https://github.com/TeamErlich/dna-fountain','source_commit':commit,'payload_bytes':size,'n_messages':messages,'noise_seeds_per_message':3,'port_parameters':{'chunk_size':32,'rs_bytes':2,'delta':.001,'c_dist':.025,'oligos':180,'max_homopolymer':3,'gc_halfwidth':.05},'rows':rows,'limits':['Python-3 port of original, not experimental sequencing replication','One read per encoded oligo, no multiple-read consensus','Only four independent message identities; noise repeats not independent messages','Unequal encoded lengths consume noise RNG differently','Encoder budgets not optimized for either method','Density includes embedded seed/parity but excludes primers and molecules per copy; report this basis explicitly','Local RSNS different strand geometry and long single strand; physical oligo synthesis constraints not matched']}
 (ROOT/'results/fountain_port_audit.json').write_text(json.dumps(out,indent=2)+'\n')
 for codec in ['external_fountain_py3_port','repaired_rsns']:
  for p in [0,.005,.01,.02,.03]:
   rr=[r for r in rows if r['codec']==codec and r['p']==p];print(codec,p,sum(x['recovered'] for x in rr),'/',len(rr),rr[0]['payload_bits_per_encoded_nt'])
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--port',required=True);x=a.parse_args();run(x.port)
