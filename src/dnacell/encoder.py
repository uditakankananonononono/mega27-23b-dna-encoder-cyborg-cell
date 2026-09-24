"""DNA data storage codec with biological constraints.

Pipeline: bytes -> base-3 trits -> deterministic scrambler -> never-same
rotating base code -> DNA -> (synthesis errors) -> majority-vote decode.

Constraints (cf. Goldman et al. 2013, Grass et al. 2015): the never-same
mapping makes homopolymer runs > 1 impossible by construction; the
keystream scrambler decorrelates payload statistics so windowed GC content
stays near 0.5. Integrity: Fletcher-16 checksum; errors: 3-fold replication
with per-position majority vote.
"""
from __future__ import annotations
import hashlib
import numpy as np

BASES = "ACGT"
B2I = {b: i for i, b in enumerate(BASES)}
TRITS_PER_BYTE = 6          # 3^6 = 729 >= 256
HEADER_TRITS = 12           # message length in trits
CHECKSUM_TRITS = 12         # Fletcher-16 as 12 trits


def bytes_to_trits(data: bytes) -> list[int]:
    out = []
    for byte in data:
        for p in range(TRITS_PER_BYTE - 1, -1, -1):
            out.append((byte // 3**p) % 3)
    return out


def trits_to_bytes(trits: list[int]) -> bytes:
    assert len(trits) % TRITS_PER_BYTE == 0
    out = bytearray()
    for i in range(0, len(trits), TRITS_PER_BYTE):
        v = 0
        for t in trits[i:i + TRITS_PER_BYTE]:
            v = v * 3 + t
        out.append(v & 0xFF)
    return bytes(out)


def keystream(n: int) -> list[int]:
    """Deterministic trit keystream from SHA-256 counter mode."""
    out, counter = [], 0
    while len(out) < n:
        block = hashlib.sha256(f"mega27-dna-{counter}".encode()).digest()
        out.extend(b % 3 for b in block)
        counter += 1
    return out[:n]


def scramble(trits: list[int]) -> list[int]:
    ks = keystream(len(trits))
    return [(t + k) % 3 for t, k in zip(trits, ks)]


def descramble(trits: list[int]) -> list[int]:
    ks = keystream(len(trits))
    return [(t - k) % 3 for t, k in zip(trits, ks)]


def encode_trits_never_same(trits: list[int], seed: str = "A") -> str:
    """trit t -> the (t+1)-th next base cyclically; never equals previous."""
    prev = seed
    out = []
    for t in trits:
        assert t in (0, 1, 2)
        nxt = BASES[(B2I[prev] + 1 + t) % 4]
        out.append(nxt)
        prev = nxt
    return "".join(out)


def decode_never_same(seq: str, seed: str = "A") -> list[int]:
    trits = []
    prev = seed
    for ch in seq:
        trits.append((B2I[ch] - B2I[prev] - 1) % 4)
        prev = ch
    if any(t not in (0, 1, 2) for t in trits):
        raise ValueError("corrupt base spacing - unrepairable strand")
    return trits


def max_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best


def gc_content(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / max(len(seq), 1)


def checksum_trits(trits: list[int]) -> list[int]:
    s1 = s2 = 0
    for t in trits:
        s1 = (s1 + t + 1) % 255
        s2 = (s2 + s1) % 255
    val = (s2 << 8) | s1
    return bytes_to_trits(val.to_bytes(2, "big"))


def encode_message(data: bytes, copies: int = 3) -> list[str]:
    body_trits = bytes_to_trits(data)
    L = len(body_trits)
    header = [(L // 3**p) % 3 for p in range(HEADER_TRITS - 1, -1, -1)]
    scrambled = scramble(header + body_trits + checksum_trits(body_trits))
    dna = encode_trits_never_same(scrambled)
    return [dna for _ in range(copies)]


def introduce_errors(seq: str, sub_rate: float = 0.01, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    return "".join(BASES[rng.integers(4)] if rng.random() < sub_rate else ch
                   for ch in seq)


def decode_message(strands: list[str]) -> bytes:
    assert strands, "no strands"
    n = len(strands[0])
    assert all(len(s) == n for s in strands), "indel handling out of scope"
    voted = "".join(max({s[i] for s in strands},
                        key=[s[i] for s in strands].count) for i in range(n))
    trits = descramble(decode_never_same(voted))
    L = 0
    for t in trits[:HEADER_TRITS]:
        L = L * 3 + t
    body = trits[HEADER_TRITS:HEADER_TRITS + L]
    chk = trits[HEADER_TRITS + L:HEADER_TRITS + L + CHECKSUM_TRITS]
    if checksum_trits(body) != chk:
        raise ValueError("checksum mismatch - corruption beyond repair")
    return trits_to_bytes(body)
