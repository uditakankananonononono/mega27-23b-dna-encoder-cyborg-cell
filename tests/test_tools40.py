import sys, random
sys.path.insert(0, "src")
from dnacell.tools40 import (ALL_TOOLS, CODECS, CHANNELS, SEQ_TOOLS, ECC_TOOLS,
                             INFO_TOOLS, MET_TOOLS)

def test_exactly_40_tools():
    assert len(ALL_TOOLS) == 40
    assert len(CODECS) == 6 and len(CHANNELS) == 6 and len(SEQ_TOOLS) == 10
    assert len(ECC_TOOLS) == 6 and len(INFO_TOOLS) == 4 and len(MET_TOOLS) == 8

def test_v2_codec_deterministic_and_ternary_parity():
    s1 = CODECS["v2"]([1, 0, 1, 1] * 8)
    s2 = CODECS["v2"]([1, 0, 1, 1] * 8)
    assert s1 == s2 and set(s1) <= set("ACGT")

def test_church_rate_is_one_bit_per_base():
    bits = [0, 1, 1, 0, 1] * 4
    assert len(CODECS["church"](bits)) == len(bits)

def test_substitution_channel_rate():
    rng = random.Random(1)
    seq = "ACGT" * 250
    out = CHANNELS["substitution"](seq, 0.5, rng)
    same = sum(1 for a, b in zip(seq, out) if a == b) / len(seq)
    assert 0.15 < same < 0.45  # 0.25 expected at rate 0.5

def test_seq_tools_basic():
    assert SEQ_TOOLS["gc"]("GGCC") == 1.0
    assert SEQ_TOOLS["homopolymer_max"]("AAACT") == 3
    assert abs(SEQ_TOOLS["entropy"]("ACGT") - 2.0) < 1e-9

def test_ecc_tools():
    assert ECC_TOOLS["ternary_parity"]([0, 1, 2, 0, 1, 2]) == (6 + 6) % 3
    code = ECC_TOOLS["hamming74"]([1, 0, 1, 1])
    assert len(code) == 7
    assert ECC_TOOLS["fletcher"]([1, 2, 3]) != ECC_TOOLS["fletcher"]([3, 2, 1])
    assert ECC_TOOLS["crc8"]([0]) == 0

def test_info_tools():
    assert INFO_TOOLS["rate"](100, 200) == 0.5
    assert abs(INFO_TOOLS["capacity_binary"](0.5) - 1.0) < 1e-9
    assert INFO_TOOLS["mutual_information"]([0, 0, 1, 1], [0, 0, 1, 1]) > 0.9

def test_met_tools_grid_fba():
    S = [[1, -1]]  # one metabolite: v1 = v2
    best, v = MET_TOOLS["fba"](S, [0, 1], [0, 0], [3, 3])
    assert best == 3 and v == (3, 3)
    ko = MET_TOOLS["knockout_scan"](S, [0, 1], [0, 0], [3, 3])
    assert ko[0] == 0 and ko[1] == 3
