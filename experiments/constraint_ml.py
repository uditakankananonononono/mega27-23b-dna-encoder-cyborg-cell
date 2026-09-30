"""Optional verdict #15: learn sequence-origin constraints, NOT synthesis difficulty.
Matched 300-nt windows; grouped holdouts prevent splitting one natural accession.
Ablation asks whether easy homopolymer leakage explains discrimination.
"""
from pathlib import Path
import sys,json,hashlib,itertools
import numpy as np
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score,balanced_accuracy_score
sys.path.insert(0,str(Path(__file__).parent))
from rs_inner import rsns_encode
from natural_genome_realism import stats
ROOT=Path(__file__).resolve().parents[1]
kmers=[''.join(k) for k in itertools.product('ACGT',repeat=3)]
def feat(s):
 d=stats(s);return [d['gc'],d['entropy'],d['max_homopolymer']]+[d['kspec'].get(k,0) for k in kmers]
def main():
 rng=np.random.default_rng(20260930);X=[];y=[];groups=[];manifest=[]
 for f in sorted((ROOT/'data/payloads').glob('*.fasta')):
  seq=''.join(x.strip().upper() for x in f.read_text().splitlines() if not x.startswith('>'))
  if len(seq)<300 or set(seq)-set('ACGT'):continue
  # One window per accession removes dependence between windows within a gene.
  start=int(rng.integers(len(seq)-299));s=seq[start:start+300]
  X.append(feat(s));y.append(1);groups.append(f.name)
  manifest.append({'header':f.read_text().splitlines()[0],'file':f.name,'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'window_start':start})
 n=len(y)
 for i in range(n):
  s=rsns_encode(rng.bytes(128))[:300];assert len(s)==300
  X.append(feat(s));y.append(0);groups.append('synthetic-'+str(i))
 X=np.array(X);y=np.array(y);groups=np.array(groups);out={}
 views={'full':list(range(X.shape[1])),'no_homopolymer':[0,1]+list(range(3,X.shape[1])),'composition_only':[0,1],'kmer_only':list(range(3,X.shape[1]))}
 splits=list(StratifiedGroupKFold(5,shuffle=True,random_state=73).split(X,y,groups))
 for name,cols in views.items():
  rows=[];pred=np.zeros(len(y))
  for fold,(tr,te) in enumerate(splits):
   m=RandomForestClassifier(n_estimators=200,max_depth=4,min_samples_leaf=5,random_state=73+fold,n_jobs=1)
   m.fit(X[tr][:,cols],y[tr]);p=m.predict_proba(X[te][:,cols])[:,1];pred[te]=p
   rows.append({'auc':float(roc_auc_score(y[te],p)),'balanced_accuracy':float(balanced_accuracy_score(y[te],p>=.5)),'test_n':len(te)})
  boot=[]
  for _ in range(1000):
   a=np.r_[rng.choice(np.where(y==0)[0],n,replace=True),rng.choice(np.where(y==1)[0],n,replace=True)]
   boot.append(roc_auc_score(y[a],pred[a]))
  out[name]={'folds':rows,'pooled_oof_auc':float(roc_auc_score(y,pred)),'bootstrap95_oof_auc':np.quantile(boot,[.025,.975]).tolist()}
 result={'task':'natural-accession vs RSNS storage-sequence origin discrimination','n_per_class':n,'window_nt':300,'seed':20260930,'views':out,'natural_manifest':manifest,'natural_source_caveat':'Accessions share one WGS project/strain, so holdout genes are not independent genomes. Sequence-origin recognition can exploit codec construction and assembly patterns.','limits':['NOT a synthesis-difficulty predictor: no experimental synthesis labels','E. coli K-12 strain C3 WGS accession windows versus one codec only, not whole-genome or organism transfer','OOF bootstrap conditions on fitted fold models, not full refit uncertainty','Homopolymer ablation still includes correlated k-mers; high accuracy does not prove a new biological mechanism']}
 (ROOT/'results/constraint_ml.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:{'auc':v['pooled_oof_auc'],'ci':v['bootstrap95_oof_auc']} for k,v in out.items()},indent=2))
if __name__=='__main__':main()
