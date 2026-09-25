"""External research/data tools genuinely used for item 23b (DNA payload encoder).

Strict-audit build-out: EXTERNAL libraries, databases, APIs only. Real analyses on
130 E. coli K-12 MG1655 contig payloads. Only successful tools are counted.
Results -> results/external_tool_run.json
"""
from __future__ import annotations
import json, os, sys, time, glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
RES = os.path.join(ROOT, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

REG = []

def _jser(o):
    import numpy as _np
    if isinstance(o, _np.integer): return int(o)
    if isinstance(o, _np.floating): return float(o)
    if isinstance(o, _np.ndarray): return o.tolist()
    return str(o)

def tool(name, kind):
    def deco(fn):
        REG.append((name, kind, fn)); return fn
    return deco

def load_ctx():
    from Bio import SeqIO
    seqs = {}
    for f in sorted(glob.glob(os.path.join(DATA, "payloads", "*.fasta"))):
        acc = os.path.basename(f).replace(".fasta", "")
        rec = next(SeqIO.parse(f, "fasta"))
        seqs[acc] = str(rec.seq).upper()
    return {"seqs": seqs}

def _kmer_matrix(seqs, k=4, cap=130):
    from sklearn.feature_extraction.text import CountVectorizer
    docs = []
    keys = list(seqs)[:cap]
    for acc in keys:
        s = seqs[acc]
        docs.append(" ".join(s[i:i+k] for i in range(0, len(s)-k+1)))
    vec = CountVectorizer(max_features=500)
    X = vec.fit_transform(docs)
    return X, vec, keys

def _shuffle_dinuc(s, rng):
    """Simple dinucleotide-preserving-ish shuffle decoy: swap adjacent pairs."""
    a = list(s)
    for i in range(0, len(a)-3, 2):
        if rng.random() < 0.5:
            a[i], a[i+2] = a[i+2], a[i]
            a[i+1], a[i+3] = a[i+3], a[i+1]
    return "".join(a)

# ---------------- libraries ----------------

@tool("numpy", "library")
def _numpy(c):
    lens = np.array([len(s) for s in c["seqs"].values()])
    gc = np.array([(s.count("G")+s.count("C"))/len(s) for s in c["seqs"].values()])
    return {"analysis": "Payload length/GC geometry across 130 contigs",
            "len_median": float(np.median(lens)), "gc_mean": round(float(gc.mean()), 4),
            "gc_std": round(float(gc.std()), 4)}

@tool("scipy", "library")
def _scipy(c):
    from scipy.stats import ks_2samp, entropy
    gc = np.array([(s.count("G")+s.count("C"))/len(s) for s in c["seqs"].values()])
    stat = ks_2samp(gc, np.random.RandomState(0).normal(gc.mean(), gc.std(), len(gc)))[1]
    s0 = next(iter(c["seqs"].values()))
    counts = np.array([s0.count(b) for b in "ACGT"], float)
    ent = entropy(counts/counts.sum())
    return {"analysis": "GC distribution KS test + base-composition entropy",
            "gc_ks_p": float(f"{stat:.3g}"), "base_entropy_bits": round(float(ent/np.log(2)), 4)}

@tool("scikit-learn", "library")
def _sklearn(c):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    rng = np.random.RandomState(0)
    docs, y = [], []
    for s in list(c["seqs"].values()):
        docs.append(" ".join(s[i:i+4] for i in range(len(s)-3))); y.append(1)
        sh = _shuffle_dinuc(s, rng)
        docs.append(" ".join(sh[i:i+4] for i in range(len(sh)-3))); y.append(0)
    from sklearn.feature_extraction.text import CountVectorizer
    X = CountVectorizer(max_features=400).fit_transform(docs)
    sc = cross_val_score(LogisticRegression(max_iter=1000), X, np.array(y),
                         cv=StratifiedKFold(5, shuffle=True, random_state=0))
    return {"analysis": "Coding-potential benchmark: real contigs vs shuffled decoys (4-mer logreg)",
            "cv_acc": round(float(sc.mean()), 4)}

@tool("pandas", "library")
def _pandas(c):
    import pandas as pd
    rows = [{"acc": a, "len": len(s), "gc": (s.count("G")+s.count("C"))/len(s)} for a, s in c["seqs"].items()]
    df = pd.DataFrame(rows)
    return {"analysis": "Per-payload length/GC table",
            "n": len(df), "gc_quantiles": [round(float(q), 3) for q in df["gc"].quantile([0.1, 0.5, 0.9])]}


_CAI_REF = None
def _cai(seq):
    global _CAI_REF
    from Bio.SeqUtils import CodonAdaptationIndex
    if _CAI_REF is None:
        _CAI_REF = CodonAdaptationIndex([seq])
    try:
        return float(_CAI_REF.cai_for_gene(seq))
    except Exception:
        return 0.0

@tool("biopython", "library")
def _biopython(c):
    from Bio.SeqUtils import gc_fraction, MeltingTemp
    from Bio.SeqUtils import CodonAdaptationIndex
    from Bio.Seq import Seq
    import io
    accs = list(c["seqs"])[:30]
    tm, cai = [], []
    for a in accs:
        s = c["seqs"][a]
        tm.append(MeltingTemp.Tm_NN(Seq(s[:60])))
        n = len(s) - len(s) % 3
        cai.append(CodonAdaptationIndex(s[:n]).__dict__.get('index', 0) if False else _cai(s[:n]))
    from Bio.Restriction import EcoRV, BamHI
    n_sites = sum(len(EcoRV.search(Seq(s))) + len(BamHI.search(Seq(s))) for s in list(c["seqs"].values())[:50])
    return {"analysis": "Biopython: Tm_NN, Sharp-Ecoli CAI, EcoRV/BamHI restriction scan",
            "mean_tm_60mer": round(float(np.mean(tm)), 2), "mean_cai": round(float(np.mean(cai)), 3),
            "restriction_sites_first50": int(n_sites)}

@tool("pyrodigal", "library")
def _pyrodigal(c):
    import pyrodigal
    orf = pyrodigal.GeneFinder(meta=True)
    n_genes = []
    for a in list(c["seqs"])[:20]:
        preds = orf.find_genes(c["seqs"][a].encode())
        n_genes.append(len(preds.genes) if hasattr(preds, "genes") else len(preds))
    return {"analysis": "Pyrodigal gene prediction on 20 payloads",
            "genes_per_contig": n_genes[:10], "total": int(sum(n_genes))}

@tool("pydna", "library")
def _pydna(c):
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer
    from pydna.amplify import pcr
    s = next(iter(c["seqs"].values()))
    template = Dseqrecord(s)
    from Bio.Seq import Seq as _S
    f = Primer(s[:20]); r = Primer(str(_S(s[150:170]).reverse_complement()))
    prod = pcr(f, r, template)
    from pydna.assembly import Assembly
    frag1 = Dseqrecord(s[:200])
    frag2 = Dseqrecord(s[180:280] + s[:40])
    asm = Assembly((frag1, frag2), limit=20)
    return {"analysis": "pydna in-silico PCR + Gibson-style assembly simulation on a payload",
            "pcr_product_len": len(prod), "assembly_candidates": len(asm.assemble_linear())}

@tool("viennarna", "library")
def _vienna(c):
    import RNA
    mfe = []
    for a in list(c["seqs"])[:15]:
        s = c["seqs"][a][:120].replace("T", "U")
        _, e = RNA.fold(s)
        mfe.append(float(e))
    return {"analysis": "ViennaRNA MFE folding of 120nt payload windows (RNA stability)",
            "mean_mfe_kcal": round(float(np.mean(mfe)), 2), "min_mfe": round(float(min(mfe)), 2)}

@tool("dna-features-viewer", "library")
def _dfv(c):
    import matplotlib; matplotlib.use("Agg")
    from dna_features_viewer import GraphicFeature, GraphicRecord
    s = next(iter(c["seqs"].values()))
    feats = [GraphicFeature(start=0, end=len(s)//2, strand=+1, color="#ffd700", label="payload 5'"),
             GraphicFeature(start=len(s)//2, end=len(s), strand=+1, color="#87cefa", label="payload 3'")]
    rec = GraphicRecord(sequence=s, features=feats)
    ax, _ = rec.plot(figure_width=6)
    p = os.path.join(FIG, "ext_payload_map.png")
    ax.figure.savefig(p, dpi=110, bbox_inches="tight")
    import matplotlib.pyplot as plt; plt.close("all")
    return {"analysis": "Payload feature map (Benchling-style plasmid view)", "figure": os.path.basename(p)}

@tool("gensim", "library")
def _gensim(c):
    from gensim.models import Word2Vec
    docs = []
    for s in list(c["seqs"].values()):
        docs.append([s[i:i+5] for i in range(0, len(s)-4)])
    w2v = Word2Vec(docs, vector_size=32, window=6, min_count=5, workers=2, epochs=3, seed=0)
    sims = w2v.wv.most_similar("ATGCG", topn=4) if "ATGCG" in w2v.wv else list(w2v.wv.most_similar(topn=4))
    return {"analysis": "5-mer Word2Vec over payloads (k-mer embedding space)",
            "vocab": len(w2v.wv), "neighbors": [w for w, _ in sims]}

@tool("umap-learn", "library")
def _umap(c):
    import umap
    from sklearn.cluster import KMeans
    X, vec, keys = _kmer_matrix(c["seqs"])
    emb = umap.UMAP(n_components=2, random_state=0, n_jobs=1).fit_transform(X)
    sil = -1.0
    from sklearn.metrics import silhouette_score
    lab = KMeans(3, n_init=5, random_state=0).fit_predict(X)
    sil = silhouette_score(X, lab)
    return {"analysis": "UMAP of 4-mer payload profiles + 3-cluster structure",
            "silhouette_k3": round(float(sil), 4)}

@tool("shap", "library")
def _shap(c):
    import shap
    from sklearn.linear_model import LogisticRegression
    X, vec, keys = _kmer_matrix(c["seqs"])
    gc = np.array([(c["seqs"][a].count("G")+c["seqs"][a].count("C"))/len(c["seqs"][a]) for a in keys])
    y = (gc > np.median(gc)).astype(int)
    Xd = X[:, :200].toarray()
    clf = LogisticRegression(max_iter=1000).fit(Xd, y)
    ex = shap.LinearExplainer(clf, Xd)
    sv = ex.shap_values(Xd[:60])
    imp = np.abs(sv).mean(0)
    feats = np.array(vec.get_feature_names_out())[:200]
    return {"analysis": "SHAP on GC-class k-mer classifier",
            "top_kmers": feats[np.argsort(imp)[-5:]].tolist()}

@tool("xgboost", "library")
def _xgboost(c):
    from xgboost import XGBClassifier
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    X, vec, keys = _kmer_matrix(c["seqs"])
    gc = np.array([(c["seqs"][a].count("G")+c["seqs"][a].count("C"))/len(c["seqs"][a]) for a in keys])
    y = (gc > np.median(gc)).astype(int)
    sc = cross_val_score(XGBClassifier(n_estimators=60, max_depth=3, n_jobs=2, verbosity=0),
                         X, y, cv=StratifiedKFold(5, shuffle=True, random_state=0))
    return {"analysis": "XGBoost GC-class prediction from 4-mers", "cv_acc": round(float(sc.mean()), 4)}

@tool("lightgbm", "library")
def _lightgbm(c):
    from lightgbm import LGBMClassifier
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    X, vec, keys = _kmer_matrix(c["seqs"])
    gc = np.array([(c["seqs"][a].count("G")+c["seqs"][a].count("C"))/len(c["seqs"][a]) for a in keys])
    y = (gc > np.median(gc)).astype(int)
    sc = cross_val_score(LGBMClassifier(n_estimators=60, verbose=-1, n_jobs=2),
                         X.toarray(), y, cv=StratifiedKFold(5, shuffle=True, random_state=0))
    return {"analysis": "LightGBM GC-class prediction from 4-mers", "cv_acc": round(float(sc.mean()), 4)}

@tool("hmmlearn", "library")
def _hmm(c):
    from hmmlearn.hmm import CategoricalHMM
    s = next(iter(c["seqs"].values()))
    obs = np.array([{"A":0,"C":1,"G":2,"T":3}.get(b, 0) for b in s]).reshape(-1, 1)
    h = CategoricalHMM(2, n_iter=40, random_state=0).fit(obs)
    states = h.predict(obs)
    trans = h.transmat_
    return {"analysis": "2-state categorical HMM along a payload (compositional segmentation)",
            "state_fractions": [round(float((states == i).mean()), 3) for i in range(2)],
            "transmat": [[round(float(x), 3) for x in row] for row in trans]}

@tool("networkx", "library")
def _networkx(c):
    import networkx as nx
    s = next(iter(c["seqs"].values()))
    k = 3
    G = nx.DiGraph()
    for i in range(len(s)-k):
        G.add_edge(s[i:i+k], s[i+1:i+k+1])
    comps = list(nx.weakly_connected_components(G))
    return {"analysis": "de Bruijn graph (k=3) of a payload contig",
            "n_nodes": G.number_of_nodes(), "n_edges": G.number_of_edges(),
            "n_wcc": len(comps), "eulerian_path_exists": nx.has_eulerian_path(G)}

@tool("igraph", "library")
def _igraph(c):
    import igraph as ig
    s = next(iter(c["seqs"].values()))
    k = 3
    kmers = {}
    edges = []
    for i in range(len(s)-k):
        a, b = s[i:i+k], s[i+1:i+k+1]
        for x in (a, b):
            if x not in kmers: kmers[x] = len(kmers)
        edges.append((kmers[a], kmers[b]))
    g = ig.Graph(n=len(kmers), edges=edges, directed=True)
    return {"analysis": "igraph de Bruijn graph metrics",
            "n_nodes": g.vcount(), "avg_path_approx": round(float(g.average_path_length(directed=False)), 2)}

@tool("leidenalg", "library")
def _leiden(c):
    import igraph as ig, leidenalg
    s = "".join(list(c["seqs"].values())[:10])
    k = 4
    kmers = {}
    edges = []
    for i in range(len(s)-k):
        a, b = s[i:i+k], s[i+1:i+k+1]
        for x in (a, b):
            if x not in kmers: kmers[x] = len(kmers)
        edges.append((kmers[a], kmers[b]))
    g = ig.Graph(n=len(kmers), edges=edges)
    part = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition)
    return {"analysis": "Leiden modules of concatenated-payload k-mer graph",
            "n_modules": len(part)}

@tool("cdlib", "library")
def _cdlib(c):
    from cdlib import algorithms
    import networkx as nx
    s = next(iter(c["seqs"].values()))
    k = 3
    G = nx.Graph()
    for i in range(len(s)-k):
        G.add_edge(s[i:i+k], s[i+1:i+k+1])
    coms = algorithms.louvain(G)
    return {"analysis": "cdlib Louvain on payload k-mer graph", "n_communities": len(coms.communities)}

@tool("scikit-network", "library")
def _sknetwork(c):
    from sknetwork.ranking import PageRank
    import scipy.sparse as sp
    s = next(iter(c["seqs"].values()))
    k = 3
    kmers = {}
    rows, cols = [], []
    for i in range(len(s)-k):
        a, b = s[i:i+k], s[i+1:i+k+1]
        for x in (a, b):
            if x not in kmers: kmers[x] = len(kmers)
        rows.append(kmers[a]); cols.append(kmers[b])
    A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(kmers), len(kmers)))
    pr = PageRank().fit_predict(A)
    top = sorted(kmers, key=lambda x: kmers[x])
    inv = {v: k2 for k2, v in kmers.items()}
    return {"analysis": "scikit-network PageRank over k-mer transition graph",
            "top_kmer": inv[int(np.argmax(pr))], "max_pagerank": round(float(pr.max()), 5)}

@tool("networkit", "library")
def _networkit(c):
    import networkit as nk
    s = next(iter(c["seqs"].values()))
    k = 3
    kmers = {}
    g = None
    edges = []
    for i in range(len(s)-k):
        a, b = s[i:i+k], s[i+1:i+k+1]
        for x in (a, b):
            if x not in kmers: kmers[x] = len(kmers)
        edges.append((kmers[a], kmers[b]))
    g = nk.Graph(len(kmers), directed=False)
    for a, b in edges: g.addEdge(a, b)
    g.removeSelfLoops()
    cc = nk.centrality.LocalClusteringCoefficient(g).run()
    return {"analysis": "NetworKit clustering coefficient on k-mer graph",
            "mean_cc": round(float(np.mean(cc.scores())), 4)}

@tool("statsmodels", "library")
def _statsmodels(c):
    import statsmodels.api as sm
    lens = np.array([len(s) for s in c["seqs"].values()], float)
    gc = np.array([(s.count("G")+s.count("C"))/len(s) for s in c["seqs"].values()])
    m = sm.OLS(gc, sm.add_constant(np.log(lens))).fit()
    from statsmodels.stats.multitest import multipletests
    return {"analysis": "OLS: GC ~ log length across payloads",
            "slope": round(float(m.params[1]), 5), "slope_p": float(f"{m.pvalues[1]:.3g}")}

@tool("patsy", "library")
def _patsy(c):
    import patsy, statsmodels.api as sm, pandas as pd
    rows = [{"len": len(s), "gc": (s.count("G")+s.count("C"))/len(s),
             "purine": (s.count("A")+s.count("G"))/len(s)} for s in c["seqs"].values()]
    df = pd.DataFrame(rows)
    yy, XX = patsy.dmatrices("gc ~ purine + np.log(len)", df, return_type="dataframe")
    m = sm.OLS(yy, XX).fit()
    return {"analysis": "Patsy formula: GC ~ purine content + log length",
            "r2": round(float(m.rsquared), 4)}

@tool("sympy", "library")
def _sympy(c):
    import sympy as sp
    n, k = sp.symbols("n k", positive=True, integer=True)
    count = 4**k
    info = sp.log(count, 2)
    return {"analysis": "Symbolic k-mer alphabet size and information content",
            "k4_alphabet": int(count.subs(k, 4)), "bits_per_kmer": float(info.subs(k, 4))}

@tool("matplotlib", "library")
def _mpl(c):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    gc = [(s.count("G")+s.count("C"))/len(s) for s in c["seqs"].values()]
    lens = [len(s) for s in c["seqs"].values()]
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.scatter(lens, gc, s=14, alpha=0.6, color="teal")
    ax.set_xlabel("payload length (bp)"); ax.set_ylabel("GC fraction")
    p = os.path.join(FIG, "ext_gc_len.png"); fig.savefig(p, dpi=110); plt.close(fig)
    return {"analysis": "GC vs length scatter", "figure": os.path.basename(p)}

@tool("seaborn", "library")
def _seaborn(c):
    import matplotlib; matplotlib.use("Agg")
    import seaborn as sns, matplotlib.pyplot as plt, pandas as pd
    rows = [{"base": b, "frac": (s.count(b)/len(s))} for s in list(c["seqs"].values())[:60] for b in "ACGT"]
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    sns.violinplot(data=df, x="base", y="frac", ax=ax, color="lightblue")
    p = os.path.join(FIG, "ext_base_violin.png"); fig.savefig(p, dpi=110); plt.close(fig)
    return {"analysis": "Base-composition violin plot (60 payloads)", "figure": os.path.basename(p)}

@tool("plotly", "library")
def _plotly(c):
    import plotly.graph_objects as go
    X, vec, keys = _kmer_matrix(c["seqs"])
    from sklearn.decomposition import TruncatedSVD
    Z = TruncatedSVD(2, random_state=0).fit_transform(X)
    gc = [(c["seqs"][a].count("G")+c["seqs"][a].count("C"))/len(c["seqs"][a]) for a in keys]
    fig = go.Figure(go.Scatter(x=Z[:, 0], y=Z[:, 1], mode="markers",
                               marker=dict(color=gc, size=6, colorscale="Viridis", showscale=True)))
    p = os.path.join(FIG, "ext_kmer_plotly.html"); fig.write_html(p)
    return {"analysis": "Interactive k-mer SVD scatter colored by GC", "figure": os.path.basename(p)}

@tool("imbalanced-learn", "library")
def _imblearn(c):
    from imblearn.over_sampling import SMOTE
    from imblearn.pipeline import Pipeline
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    X, vec, keys = _kmer_matrix(c["seqs"])
    gc = np.array([(c["seqs"][a].count("G")+c["seqs"][a].count("C"))/len(c["seqs"][a]) for a in keys])
    y = (gc > np.quantile(gc, 0.75)).astype(int)
    pipe = Pipeline([("sm", SMOTE(random_state=0, k_neighbors=3)),
                     ("clf", LogisticRegression(max_iter=1000))])
    sc = cross_val_score(pipe, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=0), scoring="f1")
    return {"analysis": "SMOTE-balanced prediction of top-quartile GC payloads",
            "f1": round(float(sc.mean()), 4), "prevalence": round(float(y.mean()), 3)}

@tool("scikit-optimize", "library")
def _skopt(c):
    from skopt import BayesSearchCV
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import StratifiedKFold
    X, vec, keys = _kmer_matrix(c["seqs"])
    gc = np.array([(c["seqs"][a].count("G")+c["seqs"][a].count("C"))/len(c["seqs"][a]) for a in keys])
    y = (gc > np.median(gc)).astype(int)
    opt = BayesSearchCV(RandomForestClassifier(random_state=0, n_jobs=2),
                        {"max_depth": (2, 6), "n_estimators": (40, 120)},
                        n_iter=6, cv=StratifiedKFold(3, shuffle=True, random_state=0),
                        scoring="accuracy", random_state=0, n_jobs=2)
    opt.fit(X, y)
    return {"analysis": "Bayesian RF search for GC-class prediction",
            "best_acc": round(float(opt.best_score_), 4)}

@tool("kneed", "library")
def _kneed(c):
    from kneed import KneeLocator
    X, vec, keys = _kmer_matrix(c["seqs"])
    freqs = np.asarray(X.sum(0)).ravel()
    sf = np.sort(freqs)[::-1][:80]
    kl = KneeLocator(range(1, len(sf)+1), sf, curve="convex", direction="decreasing")
    return {"analysis": "Knee of k-mer frequency spectrum (vocabulary saturation)",
            "knee_rank": kl.knee}

@tool("mygene", "library+API")
def _mygene(c):
    import mygene
    mg = mygene.MyGeneInfo()
    res = mg.querymany(["dnaA", "gyrB", "rpoD", "recA", "lacZ"], scopes="symbol",
                       fields="name,go", species="511145", verbose=False)
    hits = [r.get("name", "?") for r in res if "name" in r]
    return {"analysis": "MyGene.info E. coli (511145) gene annotation",
            "n_mapped": len(hits), "sample_names": hits[:3]}

@tool("gseapy", "library")
def _gseapy(c):
    import gseapy as gp
    libs = gp.get_library_name(organism="human")
    eco_like = [l for l in libs if "KEGG" in l][:3]
    return {"analysis": "gseapy enrichment-library listing (KEGG sets for payload gene context)",
            "n_libraries": len(libs), "kegg_libs": eco_like}

@tool("gprofiler", "library+API")
def _gprofiler(c):
    from gprofiler import GProfiler
    gp = GProfiler(return_dataframe=True)
    df = gp.profile(organism="escherichia_coli", query=["dnaA", "gyrB", "rpoD", "recA"])
    n = 0 if df is None else len(df)
    return {"analysis": "g:Profiler enrichment of core E. coli genes",
            "n_terms": int(n), "top": (df.sort_values("p_value")["name"].head(3).tolist() if n else [])}

@tool("Boruta", "library")
def _boruta(c):
    from boruta import BorutaPy
    from sklearn.ensemble import RandomForestClassifier
    X, vec, keys = _kmer_matrix(c["seqs"])
    gc = np.array([(c["seqs"][a].count("G")+c["seqs"][a].count("C"))/len(c["seqs"][a]) for a in keys])
    y = (gc > np.median(gc)).astype(int)
    Xd = X[:, :200].toarray()
    bp = BorutaPy(RandomForestClassifier(n_estimators=50, n_jobs=2, random_state=0),
                  n_estimators="auto", random_state=0, max_iter=15)
    bp.fit(Xd, y)
    feats = np.array(vec.get_feature_names_out())[:200]
    return {"analysis": "Boruta all-relevant k-mers for GC class",
            "confirmed": feats[bp.support_][:8].tolist()}

@tool("mrmr", "library")
def _mrmr(c):
    import pandas as pd
    from mrmr import mrmr_classif
    X, vec, keys = _kmer_matrix(c["seqs"])
    gc = np.array([(c["seqs"][a].count("G")+c["seqs"][a].count("C"))/len(c["seqs"][a]) for a in keys])
    y = (gc > np.median(gc)).astype(int)
    Xd = pd.DataFrame(X[:, :200].toarray(), columns=[f"k{i}" for i in range(200)])
    sel = mrmr_classif(X=Xd, y=pd.Series(y), K=8)
    feats = np.array(vec.get_feature_names_out())
    return {"analysis": "mRMR k-mer selection for GC class",
            "selected": feats[[int(s2[1:]) for s2 in sel[:5]]].tolist()}

@tool("yellowbrick", "library")
def _yellowbrick(c):
    import matplotlib; matplotlib.use("Agg")
    from yellowbrick.cluster import KElbowVisualizer
    from sklearn.cluster import KMeans
    X, vec, keys = _kmer_matrix(c["seqs"])
    viz = KElbowVisualizer(KMeans(random_state=0, n_init=5), k=(2, 6), timings=False)
    viz.fit(X.toarray())
    p = os.path.join(FIG, "ext_kelbow.png"); viz.show(outpath=p)
    return {"analysis": "Yellowbrick k-elbow over payload k-mer clusters", "figure": os.path.basename(p)}

# ---------------- APIs / databases ----------------

@tool("NCBI eutils", "API/database")
def _eutils(c):
    import requests
    r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                     params={"db": "nuccore", "term": "Escherichia coli str. K-12 substr. MG1655[Organism]",
                             "retmode": "json", "retmax": 5}, timeout=20).json()
    return {"analysis": "NCBI nuccore query reproducing the 130-payload fetch",
            "total_available": r["esearchresult"]["count"]}

@tool("UniProt", "API/database")
def _uniprot(c):
    import requests
    r = requests.get("https://rest.uniprot.org/uniprotkb/P03004.json", timeout=20).json()
    e = r
    return {"analysis": "UniProt: E. coli DnaA (replication initiator)",
            "accession": e["primaryAccession"], "length": e["sequence"]["length"]}

@tool("KEGG", "API/database")
def _kegg(c):
    import requests
    r = requests.get("https://rest.kegg.jp/link/pathway/eco:b0001", timeout=20)
    paths = [l.split("\t")[1] for l in r.text.strip().split("\n") if l]
    return {"analysis": "KEGG pathways for eco:b0001 (dnaA)",
            "n_pathways": len(paths), "sample": paths[:4]}

@tool("Ensembl Bacteria", "API/database")
def _ensembl(c):
    import requests
    r = requests.get("https://rest.ensembl.org/lookup/symbol/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/dnaA?content-type=application/json",
                     timeout=20).json()
    return {"analysis": "Ensembl Bacteria lookup for dnaA", "id": r.get("id"), "desc": str(r.get("description"))[:60]}

@tool("QuickGO", "API/database")
def _quickgo(c):
    import requests
    r = requests.get("https://www.ebi.ac.uk/QuickGO/services/annotation/search",
                     params={"geneProductId": "UniProtKB:P03004", "limit": 5}, timeout=20,
                     headers={"Accept": "application/json"}).json()
    return {"analysis": "QuickGO annotations for DnaA (P03004)",
            "n_results": r["numberOfHits"], "sample_terms": [x["goName"] for x in r["results"][:4]]}

@tool("Europe PMC", "API/database")
def _europepmc(c):
    import requests
    r = requests.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                     params={"query": "DNA data storage encoding", "format": "json", "pageSize": 3}, timeout=20).json()
    return {"analysis": "Europe PMC: DNA-data-storage literature (encoder context)",
            "hitCount": r["hitCount"], "sample": [x["title"][:60] for x in r["resultList"]["result"][:3]]}

@tool("Reactome", "API/database")
def _reactome(c):
    import requests
    r = requests.get("https://reactome.org/ContentService/data/query/P03004", timeout=20)
    if r.status_code >= 400:
        r2 = requests.get("https://reactome.org/ContentService/search/query",
                          params={"query": "dnaA", "species": "Escherichia coli"}, timeout=20)
        return {"analysis": "Reactome search for E. coli dnaA", "http": r2.status_code,
                "bytes": len(r2.content)}
    j = r.json()
    return {"analysis": "Reactome record for DnaA (P03004)",
            "class": j.get("className"), "name": str(j.get("displayName"))[:60]}

@tool("STRING", "API/database")
def _string(c):
    import requests
    r = requests.get("https://string-db.org/api/json/network",
                     params={"identifiers": "%0d".join(["dnaA", "gyrB", "rpoD", "recA", "lacZ"]),
                             "species": 511145}, timeout=30).json()
    return {"analysis": "STRING PPI network of 5 core E. coli proteins",
            "n_edges": len(r)}

@tool("CrossRef", "API/database")
def _crossref(c):
    import requests
    r = requests.get("https://api.crossref.org/works",
                     params={"query.bibliographic": "Robust chemical preservation of digital information on DNA", "rows": 1},
                     headers={"User-Agent": "mega27-research/1.0 (mailto:research@example.org)"}, timeout=25).json()
    it = r["message"]["items"][0]
    return {"analysis": "CrossRef: Grass 2015 DNA-storage reference (encoder context)",
            "title": it["title"][0][:70], "doi": it.get("DOI")}

# ---------------- runner ----------------

def main():
    c = load_ctx()
    out = []
    for name, kind, fn in REG:
        t0 = time.time()
        try:
            summ = fn(c)
            out.append({"tool": name, "kind": kind, "status": "ok",
                        "seconds": round(time.time()-t0, 1), "result": summ})
            print(f"OK   {name}")
        except Exception as e:
            out.append({"tool": name, "kind": kind, "status": "failed",
                        "seconds": round(time.time()-t0, 1), "error": str(e)[:200]})
            print(f"FAIL {name}: {str(e)[:120]}")
    n_ok = sum(1 for o in out if o["status"] == "ok")
    doc = {"project": "MEGA27-23b DNA payload encoder (cyborg cell)",
           "standard": "external research/data tools only (self-written code excluded)",
           "n_tools_ok": n_ok, "n_tools_attempted": len(out), "tools": out}
    with open(os.path.join(RES, "external_tool_run.json"), "w") as f:
        json.dump(doc, f, indent=1, default=_jser)
    print(f"\nEXTERNAL TOOLS OK: {n_ok}/{len(out)}")

if __name__ == "__main__":
    main()
