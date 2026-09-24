"""Sensitivity sweeps: block size B, replicate count, k_sen, P/O ratio."""
import json, os, sys
import numpy as np
sys.path.insert(0, "src")
import dnacell.encoder as enc
from dnacell.encoder import (encode_message_v2, decode_message_v2, introduce_errors)
from dnacell.metabolism import solve_fba, RXN
import dnacell.metabolism as mb

def codec_recovery(B, copies, rate, trials=16, msg_len=64):
    enc.BNS_BLOCK = B
    rng = np.random.default_rng(11)
    ok = 0
    for t in range(trials):
        msg = bytes(rng.integers(256, size=msg_len, dtype=np.uint8))
        strands = encode_message_v2(msg, copies=copies)
        noisy = [introduce_errors(s, sub_rate=rate, seed=500 + 31 * t + c)
                 for c, s in enumerate(strands)]
        try:
            ok += decode_message_v2(noisy) == msg
        except Exception:
            pass
    enc.BNS_BLOCK = 32
    return ok / trials

def main():
    out = {}
    out["block_size"] = {str(B): {str(r): codec_recovery(B, 3, r) for r in (0.01, 0.02, 0.03)}
                         for B in (16, 32, 64)}
    out["copies"] = {str(c): {str(r): codec_recovery(32, c, r) for r in (0.01, 0.02, 0.03)}
                     for c in (3, 5)}
    # metabolism sweeps
    ks, out["k_sen"] = mb.K_SEN, {}
    for k in (1.0, 2.0, 4.0):
        mb.K_SEN = k
        v0 = solve_fba(0.0)[0][RXN.index("sen")]
        v1 = solve_fba(1.0)[0][RXN.index("sen")]
        out["k_sen"][str(k)] = dict(sen_g0=v0, sen_g1=v1, delta=v0 - v1)
    mb.K_SEN = ks
    json.dump(out, open("results/sensitivity.json", "w"), indent=1)
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
