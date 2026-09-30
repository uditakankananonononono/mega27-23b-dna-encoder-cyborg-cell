"""Clean-channel screen/peeling sensitivity. Original port unmodified.
GF(2) graph rank is a mathematical oracle, not the external decoder result.
"""
import argparse,hashlib,json,math,sys,subprocess
from collections import Counter
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];PIN='f97c1b8a81f5c5b819209d5b5e26c9c8f4439495'
def run(port,message_filter):
 port=Path(port).resolve();assert subprocess.check_output(['git','-C',str(port),'rev-parse','HEAD']).decode().strip()==PIN;sys.path.insert(0,str(port))
 from fountain import DNAFountain
 from glass import Glass
 from utils import prepare
 rng=np.random.default_rng(511);rows=[]
 for message in range(4):
  msg=rng.bytes(2048)
  if message!=message_filter:continue
  for rs in [2,8,16,24]:
   chunk=120-rs;k=math.ceil(len(msg)/chunk);padded=msg+b'\0'*(-len(msg)%chunk);chunks=[list(padded[i:i+chunk]) for i in range(0,len(padded),chunk)]
   for hp in [3,4,5]:
    prepare(hp)
    f=DNAFountain(file_in=chunks,file_size=len(padded),chunk_size=chunk,alpha=0,stop=63,rs=rs,c_dist=.025,delta=.001,max_homopolymer=hp,gc=.05)
    packed=[];degrees=[];singleton=[];masks=[];payload=[]
    while f.good<63:
     d=f.droplet()
     if f.screen(d):
      packed.append(d.to_human_readable_DNA());degrees.append(d.degree);singleton.extend(d.num_chunks if d.degree==1 else []);masks.append(sum(1<<j for j in d.num_chunks));payload.append(bytearray(d.data))
    g=Glass(k,out='unused',header_size=4,rs=rs,c_dist=.025,delta=.001,gc=.05,max_homopolymer=hp,chunk_size=chunk);rejected=0
    for dna in packed:
     if g.isDone():break
     seed,_=g.add_dna(dna);rejected+=int(seed==-1)
    peel_ok=bool(g.isDone()) and b''.join(bytes(x) for x in g.chunks)[:2048]==msg
    # Exact linear-system elimination with RHS payload bytes. No noisy oracle.
    basis={}
    for mask,rhs in zip(masks,payload):
     while mask:
      pivot=mask.bit_length()-1
      if pivot in basis:
       m,v=basis[pivot];mask^=m;rhs=bytearray(a^b for a,b in zip(rhs,v))
      else:basis[pivot]=(mask,rhs);break
    linear_ok=False
    if len(basis)==k:
     decoded={}
     for pivot in sorted(basis):
      mask,rhs=basis[pivot];rhs=bytearray(rhs)
      for lower in range(pivot):
       if mask>>lower&1:rhs=bytearray(a^b for a,b in zip(rhs,decoded[lower]))
      decoded[pivot]=rhs
     linear_ok=b''.join(bytes(decoded[i]) for i in range(k))[:2048]==msg
    rows.append({'message':message,'rs_bytes':rs,'chunk_bytes':chunk,'chunks':k,'max_homopolymer':hp,'gc_halfwidth':.05,'oligos':63,'oligo_nt':496,'total_encoded_nt':31248,'attempted_droplets':f.tries,'degree_histogram':dict(Counter(degrees)),'degree_one_oligos':degrees.count(1),'distinct_degree_one_chunks':len(set(singleton)),'gf2_rank':len(basis),'port_clean_recovered':peel_ok,'port_chunks_done':g.chunksDone(),'rejections':rejected,'clean_linear_oracle_recovered':linear_ok,'input_sha256':hashlib.sha256(msg).hexdigest(),'encoded_sha256':hashlib.sha256(''.join(packed).encode()).hexdigest()})
   print('message',message,'rs',rs,'done',flush=True)
 out={'source':'https://github.com/jdbrody/dna-fountain','source_commit':PIN,'rows':rows,'limits':['Four payloads, software-only','Default fixed LFSR stream; no independent encoding-seed study','Graph-rank and Gaussian solve are clean mathematical diagnostics, not the published peeling decoder','Relaxing max homopolymer changes the synthesis constraint, not a free recovery gain','Parameters retain c .025 and delta .001, no optimized degree distribution','No primers or copies; do not compare invalid clean configurations under noise']};(ROOT/'results'/f'fountain_clean_graph_message_{message_filter}.json').write_text(json.dumps(out,indent=2)+'\n')
 for rs in [2,8,16,24]:
  for hp in [3,4,5]:
   a=[r for r in rows if r['rs_bytes']==rs and r['max_homopolymer']==hp];print('rs',rs,'hp',hp,'peeling',sum(r['port_clean_recovered'] for r in a),'/',len(a),'linear',sum(r['clean_linear_oracle_recovered'] for r in a),'/',len(a))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--port',required=True);a.add_argument('--message',required=True,type=int);x=a.parse_args();run(x.port,x.message)
