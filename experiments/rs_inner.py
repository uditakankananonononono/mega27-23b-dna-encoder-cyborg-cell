"""RS-over-GF(3^5) correcting inner code for the never-same DNA codec.

Frontier conclusion (results/CODEC_FRONTIER.md, 2026-09-26): density closure
needs a CORRECTING inner code over trits (RS-style), not detection. This module
implements it: message trits -> RS(n,k) over GF(243), 5 trits per symbol ->
scramble -> never-same DNA. Homopolymer 1 by construction; substitutions
corrupt 2 adjacent trits and often surface as out-of-range trit 3 -> symbol
ERASURE, which RS decodes at half the redundancy cost of an error.

No external deps beyond numpy (already required by the repo).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from dnacell.encoder import (bytes_to_trits, trits_to_bytes, scramble, descramble, keystream,
                             encode_trits_never_same, TRITS_PER_BYTE)

# ---------------- GF(243) = GF(3^5) ----------------
# Element = int 0..242, base-3 digits = polynomial coefficients (little-endian).

def _add(a, b):
    r, m, carry = 0, 1, 0
    while m <= 243 * 3:
        da, db = (a // m) % 3, (b // m) % 3
        r += ((da + db) % 3) * m
        m *= 3
        if m > 3**5:
            break
    return r

def _mul_raw(a, b):
    """polynomial multiply in GF(3)[x], no reduction"""
    da = [(a // 3**i) % 3 for i in range(5)]
    db = [(b // 3**i) % 3 for i in range(5)]
    c = [0] * 9
    for i in range(5):
        for j in range(5):
            c[i + j] = (c[i + j] + da[i] * db[j]) % 3
    return sum(d * 3**i for i, d in enumerate(c))

def _mod_p(a, p):
    """reduce polynomial a (up to deg 8) modulo p (deg 5), coefficients base 3"""
    digs = [(a // 3**i) % 3 for i in range(9)]
    pd = [(p // 3**i) % 3 for i in range(6)]
    assert pd[5] == 1
    for i in range(8, 4, -1):
        if digs[i]:
            f = digs[i]
            for j in range(6):
                digs[i - 5 + j] = (digs[i - 5 + j] - f * pd[j]) % 3
    return sum(d * 3**i for i, d in enumerate(digs[:5]))

def _pow_p(x, e, p):
    r = 1
    while e:
        if e & 1:
            r = _mod_p(_mul_raw(r, x), p)
        x = _mod_p(_mul_raw(x, x), p)
        e >>= 1
    return r

def _find_primitive_poly():
    """monic deg-5 poly over GF(3): irreducible (prime degree: no root in GF(3)
    and x^(3^5)==x mod p) and x has order 242 (=2*11^2)."""
    candidates = []
    for c in range(1, 3**5):
        digs = [(c // 3**i) % 3 for i in range(5)]
        if digs[0] == 0:
            continue
        p = c + 3**5
        if any(p % 3**1 == 0 for _ in ()):
            pass
        # no root in GF(3)
        def ev(v):
            return sum(((p // 3**i) % 3) * pow(v, i, 3) for i in range(6)) % 3
        if ev(0) == 0 or ev(1) == 0 or ev(2) == 0:
            continue
        if _pow_p(3, 3**5, p) != 3:  # x^(3^5) == x ?
            continue
        if _pow_p(3, 242, p) != 1:
            continue
        if _pow_p(3, 121, p) == 1 or _pow_p(3, 22, p) == 1:
            continue
        candidates.append(p)
    assert candidates, "no primitive poly found"
    return candidates[0]

PRIM = _find_primitive_poly()
EXP = [0] * 243
LOG = [0] * 243
_v = 1
for _i in range(242):
    EXP[_i] = _v
    LOG[_v] = _i
    _v = _mod_p(_mul_raw(_v, 3), PRIM)  # multiply by x (poly '3' = x)
assert _v == 1, "generator order != 242"
EXP[242] = EXP[0]

def gadd(a, b):
    return _add(a, b) if a and b else (a or b)

def gmul(a, b):
    if a == 0 or b == 0:
        return 0
    return EXP[(LOG[a] + LOG[b]) % 242]

def ginv(a):
    return EXP[(242 - LOG[a]) % 242]

def gpow(a, e):
    if a == 0:
        return 0
    return EXP[(LOG[a] * e) % 242]

# polynomial ops over GF(243), lists little-endian
def pmul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] = gadd(r[i + j], gmul(a, b))
    return r

def peval(p, x):
    r = 0
    for c in reversed(p):
        r = gadd(gmul(r, x), c)
    return r

# ---------------- Reed-Solomon over GF(243) ----------------
def rs_generator(nsym):
    # char 3: (x - a^i) = (x + 2*a^i)
    g = [1]
    for i in range(1, nsym + 1):
        g = pmul(g, [gmul(2, EXP[i]), 1])
    return g

def rs_encode(msg, nsym):
    """systematic: msg (k symbols) -> msg + nsym parity"""
    g = rs_generator(nsym)
    mr = msg[:] + [0] * nsym             # big-endian: msg * x^nsym
    gr = g[::-1]
    for i in range(len(msg)):
        coef = mr[i]
        if coef:
            for j in range(1, len(gr)):
                mr[i + j] = gadd(mr[i + j], gmul(gmul(2, coef), gr[j]))  # subtract
    # char 3: codeword tail is -remainder; big-endian order, no reversal
    return msg + [gmul(2, c) if c else 0 for c in mr[len(msg):len(msg) + nsym]]

def _syndromes(r, nsym):
    return [peval(r[::-1], EXP[j]) for j in range(1, nsym + 1)]

def _berlekamp_massey(s):
    C = [1]; B = [1]; L = 0; m = 1; b = 1
    for n in range(len(s)):
        d = s[n]
        for i in range(1, L + 1):
            d = gadd(d, gmul(C[i], s[n - i]))
        if d == 0:
            m += 1
        elif 2 * L <= n:
            T = C[:]
            coef = gmul(d, ginv(b))
            while len(C) < len(B) + m:
                C.append(0)
            for i in range(len(B)):
                C[i + m] = gadd(C[i + m], gmul(gmul(2, coef), B[i]))  # subtract
            L = n + 1 - L; B = T; b = d; m = 1
        else:
            coef = gmul(d, ginv(b))
            while len(C) < len(B) + m:
                C.append(0)
            for i in range(len(B)):
                C[i + m] = gadd(C[i + m], gmul(gmul(2, coef), B[i]))
            m += 1
    return C, L

def rs_decode(r, nsym, erasures=()):
    """r: received n symbols; erasures: known-bad positions. Returns corrected
    list or None. 2*errors + len(erasures) <= nsym required."""
    n = len(r)
    S = _syndromes(r, nsym)
    if max(S) == 0 and not erasures:
        return list(r)
    # errata locator: erasure part first, BM on Forney syndrome for errors
    eloc = [1]
    for pos in erasures:
        X = EXP[n - 1 - pos]
        eloc = pmul(eloc, [1, gmul(2, X)])  # little-endian (1 - X x): const 1, coef -X
    # Forney syndrome: S(x) * eloc(x), take low nsym terms (S_j at index j-1)
    Spoly = S[:]  # S_j at index j-1, little-endian in j
    fs_full = pmul(Spoly, eloc)
    fsyn = fs_full[:nsym]
    # strip leading zero terms? BM on fsyn starting after len(erasures) terms
    C, L = _berlekamp_massey(fsyn[len(erasures):]) if erasures else _berlekamp_massey(S)
    # combined errata locator
    Lambda = pmul(C, eloc)
    # Chien: find positions where Lambda( X^{-1} ) == 0 with X = a^{n-1-pos}
    errata = []
    for pos in range(n):
        Xinv = EXP[(242 - (n - 1 - pos)) % 242]
        if peval(Lambda, Xinv) == 0:
            errata.append(pos)
    if len(errata) != len(Lambda) - 1 or len(errata) > nsym:
        return None
    # solve syndromes S_j = sum_m Y_m X_m^j, j=1..len(errata), for magnitudes
    v = len(errata)
    Xs = [EXP[n - 1 - p] for p in errata]
    A = [[gpow(Xs[m], j + 1) for m in range(v)] for j in range(v)]
    bb = S[:v]
    # Gaussian elimination over GF(243)
    for col in range(v):
        piv = next((rr for rr in range(col, v) if A[rr][col]), None)
        if piv is None:
            return None
        A[col], A[piv] = A[piv], A[col]
        bb[col], bb[piv] = bb[piv], bb[col]
        inv = ginv(A[col][col])
        A[col] = [gmul(x, inv) for x in A[col]]
        bb[col] = gmul(bb[col], inv)
        for rr in range(v):
            if rr != col and A[rr][col]:
                f = A[rr][col]
                A[rr] = [gadd(x, gmul(gmul(2, f), y)) for x, y in zip(A[rr], A[col])]
                bb[rr] = gadd(bb[rr], gmul(gmul(2, f), bb[col]))
    out = list(r)
    for m, pos in enumerate(errata):
        out[pos] = gadd(out[pos], gmul(2, bb[m]))  # subtract error magnitude
    if max(_syndromes(out, nsym)) != 0:
        return None
    return out

# ---------------- symbol <-> trits ----------------
def syms_to_trits(syms):
    out = []
    for s in syms:
        for i in range(5):
            out.append((s // 3**i) % 3)
    return out

def trits_to_syms(trits):
    """returns (symbols, erasure_flags); any trit==3 marks the symbol erased"""
    assert len(trits) % 5 == 0
    syms, er = [], []
    for i in range(0, len(trits), 5):
        seg = trits[i:i + 5]
        if any(t == 3 for t in seg):
            syms.append(0)
            er.append(True)
        else:
            v = sum(t * 3**j for j, t in enumerate(seg))
            syms.append(v)
            er.append(False)
    return syms, er

# ---------------- codec ----------------
HDR_K, HDR_NSYM = 3, 6      # RS(9,3): 12+3 pad trits of length header
BLK_K, BLK_NSYM = 30, 15    # RS(45,30) body blocks, rate 2/3

def decode_never_same_soft(seq, seed="A"):
    """like decode_never_same but keeps out-of-range trit 3 as an erasure flag"""
    from dnacell.encoder import B2I
    trits = []
    prev = seed
    for ch in seq:
        trits.append((B2I[ch] - B2I[prev] - 1) % 4)
        prev = ch
    return trits

def rsns_encode(data: bytes):
    trits = bytes_to_trits(data)
    L = len(trits)
    hdr = [(L // 3**p) % 3 for p in range(14, -1, -1)]  # 15 trits = 3 symbols
    hsyms, _ = trits_to_syms(hdr)
    hcw = rs_encode(hsyms, HDR_NSYM)
    # body blocks
    body = trits[:]
    while len(body) % (BLK_K * 5):
        body.append(0)
    blocks = [body[i:i + BLK_K * 5] for i in range(0, len(body), BLK_K * 5)]
    bsyms = []
    for b in blocks:
        s, _ = trits_to_syms(b)
        bsyms.extend(rs_encode(s, BLK_NSYM))
    allsyms = hcw + bsyms
    return encode_trits_never_same(scramble(syms_to_trits(allsyms)))

def descramble_soft(raw):
    """Preserve invalid trit 3 as an erasure instead of reducing modulo three."""
    return [3 if t == 3 else (t-k) % 3 for t,k in zip(raw,keystream(len(raw)))]

def rsns_decode(dna: str, nsyms_hint=None):
    trits = descramble_soft(decode_never_same_soft(dna))
    n_hdr_trits = (HDR_K + HDR_NSYM) * 5
    hsyms, her = trits_to_syms(trits[:n_hdr_trits])
    hdec = rs_decode(hsyms, HDR_NSYM, erasures=[i for i, e in enumerate(her) if e])
    if hdec is None:
        raise ValueError("header RS decode failed")
    L = 0
    htr = syms_to_trits(hdec[:HDR_K])
    for t in htr:
        L = L * 3 + t
    n_body_trits = ((L + BLK_K * 5 - 1) // (BLK_K * 5)) * (BLK_K + BLK_NSYM) * 5
    btrits = trits[n_hdr_trits:n_hdr_trits + n_body_trits]
    if len(btrits) < n_body_trits:
        raise ValueError("strand truncated")
    bsyms, ber = trits_to_syms(btrits)
    out_syms = []
    nblk = n_body_trits // ((BLK_K + BLK_NSYM) * 5)
    for i in range(nblk):
        seg = bsyms[i * (BLK_K + BLK_NSYM):(i + 1) * (BLK_K + BLK_NSYM)]
        er = [j for j, e in enumerate(ber[i * (BLK_K + BLK_NSYM):(i + 1) * (BLK_K + BLK_NSYM)]) if e]
        dec = rs_decode(seg, BLK_NSYM, erasures=er)
        if dec is None:
            raise ValueError(f"body block {i} RS decode failed")
        out_syms.extend(dec[:BLK_K])
    body_trits = syms_to_trits(out_syms)[:L]
    return trits_to_bytes(body_trits)
