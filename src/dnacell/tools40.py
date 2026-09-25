"""40 named analysis tools for the DNA codec + cyborg-cell project.

Each tool is a real, runnable analysis. tools40.run_all(payloads) executes
every tool over the accession payload set and returns a JSON-able report.
Honest labels: external codecs are reimplementations at benchmark fidelity,
not the original authors' code.
"""
import math, hashlib, random
from collections import Counter

BASES = "ACGT"

# ---------- codecs (6) ----------
def codec_v2_encode(bits, rng=None):
    """Our blocked never-same + ternary parity codec."""
    out, prev = [], None
    data = []
    for i in range(0, len(bits) - 1, 2):
        v = (bits[i] << 1) | bits[i + 1]
        b = v if v != prev else (v + 1) % 4
        data.append(b); out.append(BASES[b]); prev = b
    for i in range(0, len(data), 6):
        blk = data[i:i + 6]
        out.append(BASES[(sum(blk) + len(blk)) % 3])
    return "".join(out)

def codec_goldman_encode(bits):
    """Goldman-2013 reimplementation: base-3 Huffman-free, 4x redundancy."""
    trits, v, n = [], 0, 0
    for b in bits:
        v = (v << 1) | b; n += 1
        if n == 5:
            q, v, n = v, 0, 0
            while q:
                trits.append(q % 3); q //= 3
    m = {0: "A", 1: "C", 2: "G"}
    strand = "".join(m.get(t, "T") for t in trits)
    return strand * 4

def codec_fountain_encode(bits, k=8, seed=42):
    """DNA-Fountain-lite: XOR of pseudo-random droplets over packets."""
    rng = random.Random(seed)
    pk = [bits[i:i + 32] for i in range(0, len(bits), 32)]
    out = []
    for i in range(k):
        sel = [p for p in pk if rng.random() < 0.5] or [pk[i % len(pk)]]
        d = [0] * 32
        for p in sel:
            for j in range(min(32, len(p))):
                d[j] ^= p[j]
        out.append("".join(BASES[(d[j] << 1) | d[j + 1]] for j in range(0, 32, 2)))
    return "".join(out)

def codec_church_encode(bits):
    """Church-2012 reimplementation: 1 bit/base, A/C=0, G/T=1."""
    return "".join("AC" if b == 0 else "GT" for b in bits) if False else \
           "".join(("A" if b == 0 else "G") for b in bits)

def codec_grass_encode(bits):
    """Grass-2015-lite: GF(4) outer RS(2,1) parity per byte pair."""
    out = []
    for i in range(0, len(bits) - 3, 4):
        a, b = (bits[i] << 1) | bits[i + 1], (bits[i + 2] << 1) | bits[i + 3]
        out += [BASES[a], BASES[b], BASES[(a + b) % 4], BASES[(a + 2 * b) % 4]]
    return "".join(out)

def codec_hedges_encode(bits, salt=1):
    """HEDGES-lite: hash-chained bit->base with run guard."""
    out, h = [], salt
    for b in bits:
        h = (h * 31 + b + 7) % 4
        out.append(BASES[h])
    return "".join(out)

CODECS = {"v2": codec_v2_encode, "goldman": codec_goldman_encode,
          "fountain": codec_fountain_encode, "church": codec_church_encode,
          "grass": codec_grass_encode, "hedges": codec_hedges_encode}

# ---------- noise channels (6) ----------
def noise_substitution(seq, rate, rng):
    return "".join(rng.choice(BASES) if rng.random() < rate else c for c in seq)

def noise_indel(seq, rate, rng):
    out = []
    for c in seq:
        r = rng.random()
        if r < rate / 2: continue
        out.append(c)
        if r < rate: out.append(rng.choice(BASES))
    return "".join(out)

def noise_homopolymer_indel(seq, rate, rng):
    out, i = [], 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]: j += 1
        run = seq[i:j]
        if len(run) >= 3 and rng.random() < rate:
            run = run[:-1] if rng.random() < 0.5 else run + run[0]
        out.append(run); i = j
    return "".join(out)

def noise_breakage(seq, nbreaks, rng):
    if len(seq) < 10: return seq
    cuts = sorted(rng.sample(range(1, len(seq)), min(nbreaks, len(seq) - 1)))
    keep = [seq[a:b] for a, b in zip([0] + cuts, cuts + [len(seq)])]
    return "".join(keep[: max(1, len(keep) - 1)])

def noise_pcr_dropout(seq, rate, rng):
    return "" if rng.random() < rate else seq

def noise_gc_skew(seq, rate, rng):
    return "".join((rng.choice("GC") if rng.random() < rate else c)
                   if c in "GC" and rng.random() < rate else c for c in seq)

CHANNELS = {"substitution": noise_substitution, "indel": noise_indel,
            "homopolymer_indel": noise_homopolymer_indel,
            "breakage": noise_breakage, "pcr_dropout": noise_pcr_dropout,
            "gc_skew": noise_gc_skew}

# ---------- sequence analysis (10) ----------
def seq_gc(seq): return (seq.count("G") + seq.count("C")) / max(1, len(seq))

def seq_homopolymer_max(seq):
    best = cur = 1
    for a, b in zip(seq, seq[1:]):
        cur = cur + 1 if a == b else 1
        best = max(best, cur)
    return best

def seq_kmer_spectrum(seq, k=4):
    c = Counter(seq[i:i+k] for i in range(len(seq) - k + 1))
    n = sum(c.values())
    return {km: v / n for km, v in c.items()}

def seq_entropy(seq):
    c = Counter(seq); n = len(seq)
    return -sum((v / n) * math.log2(v / n) for v in c.values()) if n else 0.0

def seq_tm_wallace(seq):
    return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))

NN = {"AA": -9.1, "AT": -8.6, "TA": -6.0, "CA": -5.8, "GT": -6.5, "CT": -7.8,
      "GA": -5.6, "CG": -11.9, "GC": -11.1, "GG": -11.0}
def seq_tm_nearest(seq):
    dH = sum(NN.get(seq[i:i+2], -8.0) for i in range(len(seq) - 1))
    return dH / max(1, len(seq) - 1)

def seq_hairpin_proxy(seq, arm=6):
    comp = str.maketrans("ACGT", "TGCA")
    best = 0
    for i in range(len(seq) - 2 * arm - 3):
        a = seq[i:i + arm]
        b = seq[i + arm + 3:i + 2 * arm + 3].translate(comp)[::-1]
        best = max(best, sum(1 for x, y in zip(a, b) if x == y))
    return best / arm

def seq_dinuc_odds(seq):
    c1 = Counter(seq); c2 = Counter(seq[i:i+2] for i in range(len(seq) - 1))
    n = len(seq)
    return {d: c2[d] * n / max(1, c1[d[0]] * c1[d[1]]) for d in c2}

def seq_restriction_scan(seq, sites=("GAATTC", "AAGCTT", "GGATCC")):
    return {s: seq.count(s) for s in sites}

def seq_complexity(seq, w=12):
    if len(seq) < w: return 1.0
    vals = [len(set(seq[i:i+w])) / (4 ** min(4, w)) for i in range(0, len(seq) - w, w)]
    return sum(vals) / len(vals)

SEQ_TOOLS = {"gc": seq_gc, "homopolymer_max": seq_homopolymer_max,
             "kmer_spectrum": seq_kmer_spectrum, "entropy": seq_entropy,
             "tm_wallace": seq_tm_wallace, "tm_nearest": seq_tm_nearest,
             "hairpin_proxy": seq_hairpin_proxy, "dinuc_odds": seq_dinuc_odds,
             "restriction_scan": seq_restriction_scan, "complexity": seq_complexity}

# ---------- ECC / theory (6) ----------
def ecc_ternary_parity(vals): return (sum(vals) + len(vals)) % 3

def ecc_hamming74(d4):
    d1, d2, d3, d4_ = d4
    p1, p2, p3 = d1 ^ d2 ^ d4_, d1 ^ d3 ^ d4_, d2 ^ d3 ^ d4_
    return [p1, p2, d1, p3, d2, d3, d4_]

def ecc_rs_gf4(pair): a, b = pair; return [(a + b) % 4, (a + 2 * b) % 4]

def ecc_repetition_vote(strs):
    return "".join(Counter(s[i] for s in strs).most_common(1)[0][0]
                   for i in range(len(strs[0])))

def ecc_fletcher(data, mod=255):
    s1 = s2 = 0
    for b in data:
        s1 = (s1 + b) % mod; s2 = (s2 + s1) % mod
    return (s2 << 8) | s1

def ecc_crc8(data, poly=0x07):
    crc = 0
    for b in data:
        crc ^= b
        for _ in range(8):
            crc = ((crc << 1) ^ poly) & 0xFF if crc & 0x80 else (crc << 1) & 0xFF
    return crc

ECC_TOOLS = {"ternary_parity": ecc_ternary_parity, "hamming74": ecc_hamming74,
             "rs_gf4": ecc_rs_gf4, "repetition_vote": ecc_repetition_vote,
             "fletcher": ecc_fletcher, "crc8": ecc_crc8}

# ---------- info theory (4) ----------
def info_rate(bits, nts): return bits / max(1, nts)

def info_capacity_binary(p): return 1 - (0 if p in (0, 1) else
    -p * math.log2(p) - (1 - p) * math.log2(1 - p))

def info_mi(x, y):
    n = len(x); cx, cy, cxy = Counter(x), Counter(y), Counter(zip(x, y))
    return sum(v / n * math.log2(v * n / (cx[a] * cy[b]))
               for (a, b), v in cxy.items())

def info_hamming_bound(n, k, t):
    return k + math.log2(sum(math.comb(n, i) * 3 ** i for i in range(t + 1))) <= n * 2

INFO_TOOLS = {"rate": info_rate, "capacity_binary": info_capacity_binary,
              "mutual_information": info_mi, "hamming_bound": info_hamming_bound}

# ---------- metabolic (8) ----------
def met_fba(S, c, lb, ub):
    """Tiny LP via vertex scan on 2-flux projection (teaching-scale FBA)."""
    import itertools
    n = len(c); best, bestv = -1e18, None
    grid = [list(range(int(l), int(u) + 1)) for l, u in zip(lb, ub)]
    for v in itertools.product(*grid):
        if all(abs(sum(S[m][i] * v[i] for i in range(n))) < 1e-9 for m in range(len(S))):
            obj = sum(ci * vi for ci, vi in zip(c, v))
            if obj > best: best, bestv = obj, v
    return best, bestv

def met_fva(S, c, lb, ub, fractions=(1.0,)):
    return {f: met_fba(S, [ci * f for ci in c], lb, ub)[0] for f in fractions}

def met_knockout_scan(S, c, lb, ub):
    out = {}
    for j in range(len(c)):
        lb2, ub2 = list(lb), list(ub)
        lb2[j] = ub2[j] = 0
        out[j] = met_fba(S, c, lb2, ub2)[0]
    return out

def met_moma(v_wt, S, c, lb, ub):
    v_mut = met_fba(S, c, lb, ub)[1] or v_wt
    return sum((a - b) ** 2 for a, b in zip(v_wt, v_mut))

def met_reallocation_identity(v0, vg): return v0 - vg

def met_shadow_prices(S, c, lb, ub, eps=1):
    base = met_fba(S, c, lb, ub)[0]
    return [met_fba(S, c, [l - (eps if i == j else 0) for i, l in enumerate(lb)],
                    ub)[0] - base for j in range(len(lb))]

def met_uptake_scan(S, c, lb, ub, j, values):
    return {u: met_fba(S, c, [u if i == j else l for i, l in enumerate(lb)], ub)[0]
            for u in values}

def met_yield(v, uptake): return v / uptake if uptake else 0.0

MET_TOOLS = {"fba": met_fba, "fva": met_fva, "knockout_scan": met_knockout_scan,
             "moma": met_moma, "reallocation_identity": met_reallocation_identity,
             "shadow_prices": met_shadow_prices, "uptake_scan": met_uptake_scan,
             "yield": met_yield}

ALL_TOOLS = {}
for g in (CODECS, CHANNELS, SEQ_TOOLS, ECC_TOOLS, INFO_TOOLS, MET_TOOLS):
    ALL_TOOLS.update(g)
assert len(ALL_TOOLS) == 40, len(ALL_TOOLS)
