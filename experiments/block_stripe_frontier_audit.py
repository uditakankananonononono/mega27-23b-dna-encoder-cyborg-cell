"""Test added outer parity stripes vs copy allocation with all nucleotide costs charged."""
import hashlib,json,sys,argparse
from pathlib import Path
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import block_matched_geometry_audit as M
F=M.F;E=M.E

def encode(msg,copies,detector,stripes):
 body=E.bytes_to_trits(msg);n=len(body);header=[(n//3**p)%3 for p in range(11,-1,-1)]
 t=E.scramble(header+body+E.checksum_trits(body));blocks=[(t[i:i+32]+[0]*32)[:32] for i in range(0,len(t),32)]
 assert stripes<=len(blocks)
 parity=[[sum(blocks[b][k] for b in range(g,len(blocks),stripes))%3 for k in range(32)] for g in range(stripes)]
 seq=''.join(M.pack(b,i,detector) for i,b in enumerate(blocks+parity));return [seq]*copies

def decode(strands,detector,stripes):
 if not strands or len(set(map(len,strands)))!=1 or len(strands[0])%48:raise ValueError('geometry')
 n=len(strands[0])//48;data_count=n-stripes
 if data_count<stripes:raise ValueError('stripe geometry')
 voted=''.join(max(E.BASES,key=lambda b:sum(s[i]==b for s in strands)) for i in range(n*48));blocks=[]
 for bi in range(n):
  p=M.block_decode(voted[bi*48:(bi+1)*48],bi,detector)
  if p is None:
   for seq in strands:
    p=M.block_decode(seq[bi*48:(bi+1)*48],bi,detector)
    if p is not None:break
  blocks.append(p)
 for g in range(stripes):
  ids=list(range(g,data_count,stripes));bad=[i for i in ids if blocks[i] is None];parity=blocks[data_count+g]
  if len(bad)==1 and parity is not None:
   blocks[bad[0]]=[(parity[k]-sum(blocks[i][k] for i in ids if blocks[i] is not None))%3 for k in range(32)]
  elif bad:raise ValueError('unrecoverable stripe')
 t=E.descramble(sum(blocks[:data_count],[]));length=sum(x*3**(11-i) for i,x in enumerate(t[:12]))
 if length%6 or length>len(t)-24:raise ValueError('header')
 body=t[12:12+length]
 if t[12+length:24+length]!=E.checksum_trits(body):raise ValueError('global checksum')
 if any(sum(body[i+j]*3**(5-j) for j in range(6))>255 for i in range(0,length,6)):raise ValueError('byte')
 return E.trits_to_bytes(body)

def gc_range(s,width):
 a=np.array([b in 'GC' for b in s],dtype=np.int32);c=np.r_[0,np.cumsum(a)];vals=(c[width:]-c[:-width])/width
 return [float(vals.min()),float(vals.max())]

def main(size_filter=None,stripe_filter=None):
 raw=(ROOT/'experiments/block_stripe_frontier_plan.json').read_bytes();p=json.loads(raw);rng=np.random.default_rng(p['seed']);messages=[rng.bytes(n) for n in p['sizes'] for _ in range(p['identities_per_size'])]
 geometry=[];packed={};excluded=[]
 for i,msg in enumerate(messages):
  for stripes in p['parity_stripes']:
   cand={d:encode(msg,1,d,stripes) for d in p['detectors']};met={}
   for d,ss in cand.items():
    s=ss[0];block_gc=[E.gc_content(s[k:k+48]) for k in range(0,len(s),48)];window=gc_range(s,64)
    met[d]={'nt':len(s),'hp':E.max_homopolymer(s),'gc':E.gc_content(s),'block48_gc_minmax':[min(block_gc),max(block_gc)],'sliding64_gc_minmax':window,'block_gc_gate':min(block_gc)>=.35 and max(block_gc)<=.65,'sliding_gc_gate':window[0]>=.35 and window[1]<=.65}
   good=all(v['hp']<=2 and .45<=v['gc']<=.55 for v in met.values());geometry.append({'message':i,'bytes':len(msg),'stripes':stripes,'included':good,'source_sha256':hashlib.sha256(msg).hexdigest(),'detectors':met})
   if not good:excluded.append([i,stripes]);continue
   for d,ss in cand.items():
    assert decode(ss,d,stripes)==msg
    for copies in p['copies']:packed[i,stripes,d,copies]=ss*copies
 rows=[];wrong=[]
 for size in p['sizes']:
  if size_filter is not None and size!=size_filter:continue
  for stripes in p['parity_stripes']:
   if stripe_filter is not None and stripes!=stripe_filter:continue
   ids=[i for i,m in enumerate(messages) if len(m)==size and (i,stripes,p['detectors'][0],1) in packed]
   for copies in p['copies']:
    for rate in p['strict_p']:
     for d in p['detectors']:
      counts=Counter();changed=total=0
      for i in ids:
       clean=packed[i,stripes,d,copies]
       for rep in range(p['noise_repeats']):
        noisy=[]
        for c,s in enumerate(clean):
         z,k=F.strict(s,rate,800000+i*1000+rep*20+c);noisy.append(z);changed+=k;total+=len(s)
        try:out=decode(noisy,d,stripes);outcome='exact_success' if out==messages[i] else 'silent_wrong'
        except ValueError:outcome='loud_failure'
        counts[outcome]+=1
        if outcome=='silent_wrong':wrong.append({'message':i,'stripes':stripes,'detector':d,'copies':copies,'p':rate,'repeat':rep,'expected_hex':messages[i].hex(),'actual_hex':out.hex(),'noisy_strands':noisy})
      nt=len(packed[ids[0],stripes,d,copies][0])*copies if ids else None
      row={'bytes':size,'stripes':stripes,'copies':copies,'detector':d,'p':rate,'n_identities':len(ids),'n_trials':len(ids)*p['noise_repeats'],'counts':dict(counts),'total_nt_per_file':nt,'payload_bits_per_nt':size*8/nt if nt else None,'mutated_fraction':changed/total if total else None};rows.append(row)
      print(size,stripes,copies,rate,d,dict(counts),flush=True)
 frontier=[]
 for size in p['sizes']:
  for rate in p['strict_p']:
   choices=[r for r in rows if r['bytes']==size and r['p']==rate and r['n_trials']]
   for r in choices:
    if not any(z['total_nt_per_file']<=r['total_nt_per_file'] and z['counts'].get('exact_success',0)>=r['counts'].get('exact_success',0) and z['counts'].get('silent_wrong',0)<=r['counts'].get('silent_wrong',0) and (z['total_nt_per_file']<r['total_nt_per_file'] or z['counts'].get('exact_success',0)>r['counts'].get('exact_success',0) or z['counts'].get('silent_wrong',0)<r['counts'].get('silent_wrong',0)) for z in choices):frontier.append(r)
 out={'plan_sha256':hashlib.sha256(raw).hexdigest(),'geometry':geometry,'excluded_geometry':excluded,'rows':rows,'empirical_nondominated_tested_choices':frontier,'limits':p['limits'],'frontier_scope':p['frontier'],'constraint_scope':p['additional_constraint_audit']}
 suffix=f'_{size_filter}_{stripe_filter}' if size_filter is not None else ''
 (ROOT/f'results/block_stripe_frontier_audit{suffix}.json').write_text(json.dumps(out,indent=2)+'\n');(ROOT/f'results/block_stripe_frontier_wrong{suffix}.json').write_text(json.dumps(wrong,indent=2)+'\n')
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--size',type=int);parser.add_argument('--stripes',type=int);a=parser.parse_args();main(a.size,a.stripes)
