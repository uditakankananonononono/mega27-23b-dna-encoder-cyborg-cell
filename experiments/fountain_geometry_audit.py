"""Fixed encoded-base and near-fixed oligo-length external parity sweep.
Directly calls compiled unmodified Glass decoder. Known source length allows
removal of zero padding, as required by the external chunk format.
"""
import argparse,hashlib,json,math,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from strict_substitution_benchmark import strict
PIN='f97c1b8a81f5c5b819209d5b5e26c9c8f4439495'
def run(port, message):
 port=Path(port).resolve();assert subprocess.check_output(['git','-C',str(port),'rev-parse','HEAD']).decode().strip()==PIN
 sys.path.insert(0,str(port));from glass import Glass
 from robust_solition import PRNG
 from utils import dna_to_int_array,int_to_four
 rng=np.random.default_rng(511);rows=[];size=2048
 with tempfile.TemporaryDirectory() as work:
  w=Path(work)
  for i in range(4):
   msg=rng.bytes(size)
   if i != message:continue
   inp=w/'input.bin';inp.write_bytes(msg)
   for rs in [2,8,16,24]:
    chunk=120-rs;nf=w/'clean.fasta';cmd=[sys.executable,'encode.py','--file_in',str(inp),'--size',str(chunk),'-m','3','--gc','.05','--rs',str(rs),'--delta','.001','--c_dist','.025','--out',str(nf),'--stop','63']
    proc=subprocess.run(cmd,cwd=port,capture_output=True,text=True,timeout=60);assert proc.returncode==0,proc.stderr
    seqs=[s.strip() for s in nf.read_text().splitlines() if s.strip() and not s.startswith('>')];assert len(seqs)==63 and all(len(s)==496 for s in seqs)
    nf.replace(ROOT/'results'/f'fountain_geometry_message_{i}_rs{rs}.fasta')
    graph=PRNG(K=math.ceil(size/chunk),delta=.001,c=.025);degrees=[];masks=[]
    for seq in seqs:
     data=dna_to_int_array(seq);seed=int.from_bytes(bytes(data[:4]),'big');graph.set_seed(seed);_,degree,indexes=graph.get_src_blocks_wrap();degrees.append(degree);masks.append(sum(1<<j for j in indexes))
    basis={}
    for mask in masks:
     while mask:
      pivot=mask.bit_length()-1
      if pivot in basis:mask^=basis[pivot]
      else:basis[pivot]=mask;break
    for p in [0,.005,.01,.02,.03]:
     for seed in range(3):
      g=Glass(math.ceil(size/chunk),out=str(w/'unused'),header_size=4,rs=rs,c_dist=.025,delta=.001,flag_correct=True,gc=.05,max_homopolymer=3,chunk_size=chunk);changes=0;reject=0;read=0;err=None
      for c,s in enumerate(seqs):
       z,n=strict(s,p,91000+i*1000+seed*300+c);changes+=n
       if g.isDone():continue
       try:k,_=g.add_dna(z);reject+=int(k==-1);read+=1
       except Exception as e:err=type(e).__name__+': '+str(e);break
      recovered=False
      if g.isDone():recovered=b''.join(bytes(x) for x in g.chunks)[:size]==msg
      rows.append({'degree_one_oligos':degrees.count(1),'graph_binary_rank':len(basis),'degree_histogram':{str(d):degrees.count(d) for d in sorted(set(degrees))},'message':i,'p':p,'noise_seed':seed,'rs_bytes':rs,'chunk_size':chunk,'chunks':math.ceil(size/chunk),'oligos':63,'oligo_nt':496,'total_encoded_nt':31248,'payload_bits_per_encoded_nt':size*8/31248,'input_sha256':hashlib.sha256(msg).hexdigest(),'encoded_sha256':hashlib.sha256(''.join(seqs).encode()).hexdigest(),'recovered':recovered,'chunks_done':g.chunksDone(),'reads_processed':read,'rejections':reject,'realized_substitution_fraction':changes/31248,'error':err})
    print('message',i,'rs',rs,'done',flush=True)
 out={'source':'https://github.com/jdbrody/dna-fountain','source_commit':PIN,'design':'63x496nt external oligos (31,248 nt) versus segmented RSNS 64x495nt (31,680 nt); 1.36% fewer external encoded bases','rows':rows,'limits':['Four message identities and three noise realizations, not 12 independent messages','Source payload length and codec parameters known to both decoders','External header indexes charged, local ideal oligo order uncharged','Near-matched strand length and encoded budget, not identical','GC screening and no-homopolymer constraints differ','No primers, copies, consensus, synthesis or physical errors','External RS setting selected within evaluated grid; not optimized universal frontier']}
 (ROOT/'results'/f'fountain_geometry_message_{message}.json').write_text(json.dumps(out,indent=2)+'\n')
 for rs in [2,8,16,24]:
  print('rs',rs,[(p,sum(r['recovered'] for r in rows if r['rs_bytes']==rs and r['p']==p)) for p in [0,.005,.01,.02,.03]])
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--port',required=True);a.add_argument('--message',required=True,type=int);x=a.parse_args();run(x.port,x.message)
