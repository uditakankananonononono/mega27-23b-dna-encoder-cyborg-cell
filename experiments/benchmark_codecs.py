"""Head-to-head: our codec vs Goldman et al. 2013 and DNA Fountain (Erlich 2017).

All three encode the same messages; structural stats (homopolymer, GC window)
and recovery under substitution noise are measured on the full strand sets.
Density = message bits / total synthesized bases (all copies counted).
Goldman: base-3 fixed 6 trits/byte, never-same rotating code, no scrambler,
same 3-fold replication (their overlap redundancy replaced so encodings are
compared, not redundancy schemes). Fountain: Luby droplets, 2 bits/base direct map, GC/homopolymer screening, CRC8 inner
error-detection per droplet. FAIRNESS FIX (2026-09-26): the original harness had NO inner
error detection and claimed 'no inner ECC, faithful to the paper' - that claim was wrong
(Erlich & Zielinski 2017 use a Reed-Solomon inner code to detect/correct bad droplets),
so every substitution silently corrupted chunks and cascaded through XOR. Prior fountain
recovery numbers (0% at >=1% substitution) were a harness artifact and are preserved in
git history as a documented mistake; they must not be cited as a fountain weakness.
"""
import json, os, sys
import numpy as np
sys.path.insert(0, "src")
from dnacell.encoder import (encode_message, decode_message, encode_message_v2, decode_message_v2, introduce_errors,
                             max_homopolymer, gc_content, BASES, B2I,
                             bytes_to_trits, trits_to_bytes, checksum_trits,
                             HEADER_TRITS, CHECKSUM_TRITS,
                             encode_trits_never_same, decode_never_same)

# ---------- Goldman-style baseline ----------
def goldman_encode(data: bytes, copies: int = 3):
    body = bytes_to_trits(data)
    L = len(body)
    header = [(L // 3**p) % 3 for p in range(HEADER_TRITS - 1, -1, -1)]
    dna = encode_trits_never_same(header + body + checksum_trits(body))
    return [dna] * copies

def goldman_decode(strands):
    n = len(strands[0])
    voted = "".join(max({s[i] for s in strands},
                        key=[s[i] for s in strands].count) for i in range(n))
    trits = decode_never_same(voted)
    L = 0
    for t in trits[:HEADER_TRITS]:
        L = L * 3 + t
    body = trits[HEADER_TRITS:HEADER_TRITS + L]
    chk = trits[HEADER_TRITS + L:HEADER_TRITS + L + CHECKSUM_TRITS]
    if checksum_trits(body) != chk:
        raise ValueError("checksum mismatch")
    return trits_to_bytes(body)

# ---------- DNA Fountain baseline ----------
B4 = ["00", "01", "10", "11"]
I2B = {"00": "A", "01": "C", "10": "G", "11": "T"}
B2I4 = {v: k for k, v in I2B.items()}
CHUNK = 16  # bytes per chunk

def _dbytes_to_dna(bs: bytes) -> str:
    return "".join(I2B[f"{b:08b}"[i:i+2]] for b in bs for i in range(0, 8, 2))

def _dna_to_bytes(dna: str, nbytes: int) -> bytes:
    bits = "".join(B2I4[c] for c in dna)
    return bytes(int(bits[i*8:(i+1)*8], 2) for i in range(nbytes))

def _screen_ok(dna: str) -> bool:
    return 0.44 <= gc_content(dna) <= 0.56 and max_homopolymer(dna) <= 3

def _crc8(data: bytes) -> int:
    crc = 0
    for b in data:
        crc ^= b
        for _ in range(8):
            crc = ((crc << 1) ^ 0x07) & 0xFF if crc & 0x80 else (crc << 1) & 0xFF
    return crc

def fountain_encode(data: bytes, seed0: int = 0, over: float = 2.5):
    chunks = [data[i:i+CHUNK].ljust(CHUNK, b"\0") for i in range(0, len(data), CHUNK)]
    k = len(chunks)
    n = int(np.ceil(k * over)) + 2
    strands, meta = [], []
    seed = seed0
    while len(strands) < n:
        rng = np.random.default_rng(seed)
        deg = int(rng.integers(1, min(4, k) + 1))
        idx = sorted(rng.choice(k, size=deg, replace=False).tolist())
        d = bytearray(CHUNK)
        for i in idx:
            d = bytearray(a ^ b for a, b in zip(d, chunks[i]))
        dna = _dbytes_to_dna(bytes(d) + bytes([_crc8(bytes(d))]))  # CRC8 inner check, standing in for the paper's RS inner code
        if _screen_ok(dna):
            strands.append(dna)
            meta.append((seed, idx))
        seed += 1
    return strands, (k, meta)

def fountain_decode(strands, packing):
    k, meta = packing
    known = {}
    droplets = []
    for dna, (seed, idx) in zip(strands, meta):
        try:
            raw = _dna_to_bytes(dna, CHUNK + 1)
        except KeyError:
            continue
        d, chk = raw[:CHUNK], raw[CHUNK]
        if _crc8(d) != chk:
            continue  # CRC8 detects corrupted droplet; fountain redundancy absorbs the loss
        droplets.append((set(idx), d))
    changed = True
    while changed and len(known) < k:
        changed = False
        nxt = []
        for idx, d in droplets:
            idx = set(idx)
            for i in list(idx):
                if i in known:
                    d = bytes(a ^ b for a, b in zip(d, known[i]))
                    idx.discard(i)
            if len(idx) == 1:
                i = idx.pop()
                if i not in known:
                    known[i] = d
                    changed = True
            else:
                nxt.append((idx, d))
        droplets = nxt
    if len(known) < k:
        raise ValueError(f"fountain decode failed: {len(known)}/{k} chunks")
    return b"".join(known[i] for i in range(k))

# ---------- never-same fountain (v3 candidate, 2026-09-26) ----------
# LT droplet structure as in fountain_encode, but each droplet's DNA is the
# never-same rotating trit code (homopolymer 1 by construction, no screening)
# with per-droplet CRC8. Substitutions desynchronize the trit stream, so the
# CRC turns them into droplet ERASURES absorbed by the oversampling.
NSF_TRITS = None  # trit count comes from bytes_to_trits directly

def nsfountain_encode(data: bytes, seed0: int = 0, over: float = 2.5):
    chunks = [data[i:i+CHUNK].ljust(CHUNK, b"\0") for i in range(0, len(data), CHUNK)]
    k = len(chunks)
    n = int(np.ceil(k * over)) + 2
    strands, meta = [], []
    seed = seed0
    while len(strands) < n:
        rng = np.random.default_rng(seed)
        deg = int(rng.integers(1, min(4, k) + 1))
        idx = sorted(rng.choice(k, size=deg, replace=False).tolist())
        d = bytearray(CHUNK)
        for i in idx:
            d = bytearray(a ^ b for a, b in zip(d, chunks[i]))
        payload = bytes(d) + bytes([_crc8(bytes(d))])
        strands.append(encode_trits_never_same(bytes_to_trits(payload)))
        meta.append((seed, idx))
        seed += 1
    return strands, (k, meta)

def nsfountain_decode(strands, packing):
    k, meta = packing
    known, droplets = {}, []
    for dna, (seed, idx) in zip(strands, meta):
        try:
            raw = trits_to_bytes(decode_never_same(dna))
        except (ValueError, AssertionError):
            continue
        d, chk = raw[:CHUNK], raw[CHUNK]
        if _crc8(d) != chk:
            continue  # corrupted droplet -> erasure
        droplets.append((set(idx), d))
    changed = True
    while changed and len(known) < k:
        changed = False
        nxt = []
        for idx, d in droplets:
            idx = set(idx)
            for i in list(idx):
                if i in known:
                    d = bytes(a ^ b for a, b in zip(d, known[i]))
                    idx.discard(i)
            if len(idx) == 1:
                i = idx.pop()
                if i not in known:
                    known[i] = d
                    changed = True
            else:
                nxt.append((idx, d))
        droplets = nxt
    if len(known) < k:
        raise ValueError(f"nsfountain decode failed: {len(known)}/{k} chunks")
    return b"".join(known[i] for i in range(k))


# ---------- evaluation ----------
def windowed_gc_dev(seq: str, w: int = 50) -> float:
    devs = [abs(gc_content(seq[i:i+w]) - 0.5) for i in range(0, len(seq) - w + 1, w)]
    return float(max(devs)) if devs else 0.0

def run():
    rng = np.random.default_rng(7)
    results = {}
    trials_log = []
    for msg_len in (64, 256):
        msgs = [bytes(rng.integers(256, size=msg_len, dtype=np.uint8)) for _ in range(24)]
        for name, enc in (("ours", lambda m: (encode_message(m), None)),
                          ("ours_v2", lambda m: (encode_message_v2(m), None)),
                          ("ours_2x", lambda m: (encode_message(m, copies=2), None)),
                          ("ours_v2_2x", lambda m: (encode_message_v2(m, copies=2), None)),
                          ("goldman", lambda m: (goldman_encode(m), None)),
                          ("fountain", fountain_encode),
                          ("nsfountain", nsfountain_encode)):
            dens, hp, gcdev, rec = [], [], [], {0.0: 0, 0.01: 0, 0.02: 0, 0.03: 0}
            trials = 0
            for mi, m in enumerate(msgs):
                strands, packing = enc(m)
                total_bases = sum(len(s) for s in strands)
                dens.append(8 * len(m) / total_bases)
                hp.append(max(max_homopolymer(s) for s in strands))
                gcdev.append(max(windowed_gc_dev(s) for s in strands))
                for rate in (0.0, 0.01, 0.02, 0.03):
                    noisy = [introduce_errors(s, sub_rate=rate, seed=1000 + 97 * mi + ci) for ci, s in enumerate(strands)]
                    try:
                        if name == "fountain":
                            pass
                        if name == "nsfountain":
                            dec = nsfountain_decode(noisy, packing)[:len(m)]
                        elif name == "fountain":
                            dec = fountain_decode(noisy, packing)[:len(m)]
                        elif name == "goldman":
                            dec = goldman_decode(noisy)
                        elif name in ("ours_v2", "ours_v2_2x"):
                            dec = decode_message_v2(noisy)
                        else:
                            dec = decode_message(noisy)
                        ok = dec == m
                    except Exception:
                        ok = False
                    rec[rate] += ok
                    trials += 1 if rate == 0 else 0
                    trials_log.append(dict(codec=name, size=msg_len, msg=mi, rate=rate, recovered=bool(ok)))
            n_msg = len(msgs)
            results[f"{name}_{msg_len}B"] = dict(
                density_bits_per_base=float(np.mean(dens)),
                max_homopolymer=int(max(hp)),
                max_window_gc_dev=float(max(gcdev)),
                recovery={str(r): rec[r] / n_msg for r in rec})
            print(name, msg_len, results[f"{name}_{msg_len}B"])
    os.makedirs("results", exist_ok=True)
    json.dump(results, open("results/codec_benchmark.json", "w"), indent=1)
    json.dump(trials_log, open("results/codec_trials.json", "w"))
    print("saved results/codec_benchmark.json")

if __name__ == "__main__":
    run()
