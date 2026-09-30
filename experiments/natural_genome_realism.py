#!/usr/bin/env python3
"""Verdict #9/#17: natural-genome realism. Compare our encoded DNA payloads
against E. coli K-12 strain C3 WGS accession sequences (data/payloads/) on the constraint
statistics biology 'chooses': GC content, max homopolymer, Shannon entropy,
k-mer (k=3) spectrum distance, dinucleotide odds. Question: do our synthetic
storage sequences sit inside or outside the natural distribution, and which
constraints separate synthetic from natural DNA?
"""
import json, pathlib, sys, os
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from rs_inner import rsns_encode

BASES = "ACGT"

def stats(seq):
    n = len(seq)
    c = Counter(seq)
    gc = (c["G"] + c["C"]) / n
    hp = 1; cur = 1
    for i in range(1, n):
        cur = cur + 1 if seq[i] == seq[i-1] else 1
        hp = max(hp, cur)
    p = np.array([c[b] / n for b in BASES]); p = p[p > 0]
    ent = float(-(p * np.log2(p)).sum())
    km = Counter(seq[i:i+3] for i in range(n - 2))
    tot = sum(km.values())
    kspec = {k: v / tot for k, v in km.items()}
    din = Counter(seq[i:i+2] for i in range(n - 1))
    dtot = sum(din.values())
    odds = {}
    for d_, v in din.items():
        exp = (c[d_[0]] / n) * (c[d_[1]] / n) * dtot
        odds[d_] = v / exp if exp else 0.0
    return {"n": n, "gc": gc, "max_homopolymer": hp, "entropy": ent,
            "kspec": kspec, "dinuc_odds": odds}

def main():
    rng = np.random.default_rng(21)
    nat = []
    for f in sorted(pathlib.Path("data/payloads").glob("*.fasta")):
        seq = "".join(l.strip() for l in open(f) if not l.startswith(">"))
        if len(seq) >= 300:
            nat.append((f.name, stats(seq)))
    syn = []
    for t in range(131):
        msg = rng.bytes(256)
        syn.append(stats(rsns_encode(msg)))
    def agg(rows, key):
        v = np.array([r[1][key] if isinstance(r, tuple) else r[key] for r in rows], dtype=float)
        return {"mean": float(v.mean()), "sd": float(v.std(ddof=1)),
                "min": float(v.min()), "max": float(v.max())}
    out = {"natural_n": len(nat), "synthetic_n": len(syn),
           "natural": {k: agg(nat, k) for k in ("gc", "max_homopolymer", "entropy")},
           "synthetic": {k: agg(syn, k) for k in ("gc", "max_homopolymer", "entropy")}}
    # k-mer spectrum distance (L1) between mean spectra
    keys = sorted(set().union(*[set(r[1]["kspec"]) for r in nat]))
    nat_mean = np.array([np.mean([r[1]["kspec"].get(k, 0) for r in nat]) for k in keys])
    syn_mean = np.array([np.mean([r["kspec"].get(k, 0) for r in syn]) for k in keys])
    out["kmer3_L1_natural_vs_synthetic"] = float(np.abs(nat_mean - syn_mean).sum())
    # GC z-score of synthetic vs natural distribution
    gcn = np.array([r[1]["gc"] for r in nat])
    gcs = np.array([r["gc"] for r in syn])
    out["gc_z_synthetic_vs_natural"] = float((gcs.mean() - gcn.mean()) / gcn.std(ddof=1))
    # max homopolymer: fraction of synthetic exceeding natural max
    natmax = max(r[1]["max_homopolymer"] for r in nat)
    out["natural_max_homopolymer_observed"] = int(natmax)
    out["synthetic_frac_exceeding_natural_max_hp"] = float(np.mean([r["max_homopolymer"] > natmax for r in syn]))
    # per-sequence z of entropy
    entn = np.array([r[1]["entropy"] for r in nat])
    ents = np.array([r["entropy"] for r in syn])
    out["entropy_z_synthetic_vs_natural"] = float((ents.mean() - entn.mean()) / entn.std(ddof=1))
    out["reading"] = ("synthetic sequences near the natural GC band but with distinct "
                      "k-mer/homopolymer structure identify which biological constraints "
                      "(not base composition) separate storage DNA from genome DNA")
    pathlib.Path("results/natural_genome_realism.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1)[:1200])

if __name__ == "__main__":
    main()
