from dnacell.encoder import (bytes_to_trits, trits_to_bytes, encode_message,
                             decode_message, introduce_errors, max_homopolymer,
                             gc_content, encode_trits_never_same,
                             decode_never_same)


def test_trits_roundtrip():
    data = b"Hello, cyborg cell! \x00\xff"
    assert trits_to_bytes(bytes_to_trits(data)) == data


def test_never_same_roundtrip_and_no_homopolymers():
    trits = [1] * 400  # worst case: constant trit stream
    seq = encode_trits_never_same(trits)
    assert decode_never_same(seq) == trits
    assert max_homopolymer(seq) == 1


def test_end_to_end_clean():
    msg = b"MEGA27 cyborg cell logic loop payload"
    strands = encode_message(msg)
    assert decode_message(strands) == msg


def test_end_to_end_with_errors():
    msg = b"error tolerance test payload 0123456789"
    strands = [introduce_errors(s, sub_rate=0.02, seed=i)
               for i, s in enumerate(encode_message(msg))]
    assert decode_message(strands) == msg


def test_constraints_gc_and_homopolymer():
    msg = bytes(range(256))
    dna = encode_message(msg)[0]
    assert max_homopolymer(dna) <= 3
    assert 0.35 <= gc_content(dna) <= 0.65


def test_checksum_detects_heavy_corruption():
    msg = b"checksum detection payload"
    strands = [introduce_errors(s, sub_rate=0.35, seed=i)
               for i, s in enumerate(encode_message(msg, copies=3))]
    try:
        out = decode_message(strands)
        assert out == msg  # either repaired...
    except ValueError:
        pass  # ...or loudly rejected; never silently wrong
