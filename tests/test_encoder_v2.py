import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from dnacell.encoder import (encode_message_v2, decode_message_v2,
                             introduce_errors, max_homopolymer, gc_content)


def test_roundtrip_no_noise():
    msg = b"blocked codec with parity repair" * 3
    assert decode_message_v2(encode_message_v2(msg)) == msg


def test_odd_length_roundtrip():
    for n in (1, 5, 31, 32, 33, 100):
        msg = bytes(range(n))
        assert decode_message_v2(encode_message_v2(msg)) == msg


def test_constraints():
    dna = encode_message_v2(os.urandom(120))[0]
    assert max_homopolymer(dna) <= 2
    assert abs(gc_content(dna) - 0.5) < 0.15


def test_low_noise_recovers():
    msg = b"noise resilience check" * 4
    strands = encode_message_v2(msg)
    noisy = [introduce_errors(s, sub_rate=0.005, seed=i) for i, s in enumerate(strands)]
    assert decode_message_v2(noisy) == msg


def test_single_block_erasure_repaired_by_parity():
    msg = b"P" * 200
    strands = encode_message_v2(msg)
    # corrupt one block in ALL copies identically so voting cannot fix it;
    # tier-3 parity must repair it
    blen = 38
    bad = [s[:blen] + "A" * blen + s[2 * blen:] for s in strands]
    # force all-three-agree corruption: same base at same positions
    assert decode_message_v2(bad) == msg
