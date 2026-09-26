#!/usr/bin/env python3
"""RS-ns capacity scaling: does the correcting inner code hold as message size grows?

Known: 256B at density 0.813 bits/base, recovery 1.0/1.0/0.583/0.083 at 0/1/2/3%
substitution. Here: sizes 256/512/1024/2048 bytes, same RS(45,30) body blocks +
RS(9,3) header, 16 trials per (size, rate), exact-recovery fraction and density.
Question: does per-message recovery degrade with block count (header stays fixed,
body blocks scale), and does density approach the asymptotic rate?
"""
import json, pathlib, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from rs_inner import rsns_encode, rsns_decode

BASES = "ACGT"
def substitute(dna, p, rng):
    a = list(dna)
    for i, ch in enumerate(a):
        if rng.random() < p:
            a[i] = rng.choice([b for b in BASES if b != ch])
    return "".join(a)

def main():
    rng = np.random.default_rng(5)
    out = {"trials_per_cell": 16, "cells": []}
    for size in (256, 512, 1024, 2048):
        for rate in (0.0, 0.01, 0.02, 0.03):
            ok, dens = 0, []
            for t in range(16):
                msg = rng.bytes(size)
                dna = rsns_encode(msg)
                dens.append(8.0 * size / len(dna))
                try:
                    dec = rsns_decode(substitute(dna, rate, rng))
                    if dec == msg:
                        ok += 1
                except Exception:
                    pass
            cell = {"bytes": size, "subst_rate": rate,
                    "recovery": round(ok / 16, 4),
                    "density_bits_per_base": round(float(np.mean(dens)), 4)}
            out["cells"].append(cell)
            print(cell, flush=True)
    pathlib.Path("results/rsns_capacity_sweep.json").write_text(json.dumps(out, indent=1))
    print("saved results/rsns_capacity_sweep.json")

if __name__ == "__main__":
    main()
