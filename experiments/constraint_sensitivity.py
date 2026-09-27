#!/usr/bin/env python3
"""Verdict #14: robustness under varied constraints. For the v2 codec:
- GC-window deviation distribution (windows 20/30/50 nt) over 64 messages
- homopolymer guarantee check (never-same => max run 1 by construction)
- parity-size sensitivity: BNS_CHK swept, density + recovery at 1/2/3% subst
- copies swept (replication redundancy cost)
"""
import json, sys
sys.path.insert(0, "src")
import numpy as np
import dnacell.encoder as enc
from dnacell.encoder import encode_message_v2, decode_message_v2, introduce_errors

def gc_window_dev(dna, w):
    devs = []
    for i in range(0, len(dna) - w + 1, w // 2):
        seg = dna[i:i+w]
        gc = (seg.count("G") + seg.count("C")) / w
        devs.append(abs(gc - 0.5))
    return max(devs), float(np.mean(devs))

def max_hp(dna):
    hp = cur = 1
    for i in range(1, len(dna)):
        cur = cur + 1 if dna[i] == dna[i-1] else 1
        hp = max(hp, cur)
    return hp

def recovery(chk, copies, rate, trials=16, msg_len=64):
    old = enc.BNS_CHK
    enc.BNS_CHK = chk
    try:
        rng = np.random.default_rng(11); ok = 0; dens = []
        for t in range(trials):
            msg = bytes(rng.integers(256, size=msg_len, dtype=np.uint8))
            strands = encode_message_v2(msg, copies=copies)
            dens.append(8.0 * msg_len / len(strands[0]))
            noisy = [introduce_errors(s, sub_rate=rate, seed=500+31*t+c) for c, s in enumerate(strands)]
            try: ok += decode_message_v2(noisy) == msg
            except Exception: pass
        return ok / trials, float(np.mean(dens))
    finally:
        enc.BNS_CHK = old

def main():
    rng = np.random.default_rng(7)
    out = {"gc_windows": {}, "homopolymer": {}, "parity_sweep": {}, "copies_sweep": {}}
    for w in (20, 30, 50):
        mx, mn = [], []
        for t in range(32):
            msg = bytes(rng.integers(256, size=128, dtype=np.uint8))
            dna = encode_message_v2(msg)[0]
            a, b = gc_window_dev(dna, w)
            mx.append(a); mn.append(b)
        out["gc_windows"][str(w)] = {"max_dev_worst": float(np.max(mx)), "max_dev_mean": float(np.mean(mx)), "mean_dev": float(np.mean(mn))}
        print("gc", w, out["gc_windows"][str(w)], flush=True)
    hps = [max_hp(encode_message_v2(bytes(rng.integers(256, size=128, dtype=np.uint8)))[0]) for _ in range(32)]
    out["homopolymer"] = {"max_run_over_32_messages": int(max(hps)), "guarantee": "never-same encoding: run length 1 by construction"}
    print("hp", out["homopolymer"], flush=True)
    for chk in (4, 8, 16):
        row = {}
        for rate in (0.01, 0.02, 0.03):
            r, d = recovery(chk, 3, rate)
            row[str(rate)] = {"recovery": round(r, 4), "density": round(d, 4)}
        out["parity_sweep"][f"chk_{chk}"] = row
        print("chk", chk, row, flush=True)
    for c in (1, 3, 5):
        row = {}
        for rate in (0.01, 0.02):
            r, d = recovery(8, c, rate)
            row[str(rate)] = {"recovery": round(r, 4), "effective_density": round(d / c, 4)}
        out["copies_sweep"][str(c)] = row
        print("copies", c, row, flush=True)
    out["reading"] = "homopolymer constraint is structural (free); GC balance via scrambling has measurable worst-case deviation; parity size and copy count buy recovery at exact density cost - the cellular-resource ledger in bits/base"
    json.dump(out, open("results/constraint_sensitivity.json", "w"), indent=1)
    print("saved results/constraint_sensitivity.json")

if __name__ == "__main__":
    main()
