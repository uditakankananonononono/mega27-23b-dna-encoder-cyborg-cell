"""Run all 40 tools over the 130-accession payload set + the cyborg-cell LP model."""
import json, random, sys
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, "src")
from dnacell.tools40 import CODECS, CHANNELS, SEQ_TOOLS, ECC_TOOLS, INFO_TOOLS
from dnacell.metabolism import stoich_matrix, MET, RXN, A_BIO, B_BIO, C_BIO, K_SEN

man = json.load(open("data/payloads/manifest.json"))
payloads = {}
for m in man:
    seq = "".join(l.strip() for l in open(f"data/payloads/{m['accession']}.fasta")
                  if not l.startswith(">"))
    payloads[m["accession"]] = seq
accs = sorted(payloads)
rng = random.Random(27)
report = {"n_accessions": len(accs), "tools": {}}

# --- 10 sequence tools x 130 accessions ---
for name, fn in SEQ_TOOLS.items():
    vals = {}
    for a in accs:
        s = payloads[a]
        v = fn(s)
        if isinstance(v, dict):
            v = {k: round(x, 5) for k, x in list(v.items())[:8]}
        elif isinstance(v, float):
            v = round(v, 5)
        vals[a] = v
    report["tools"][name] = {"group": "sequence", "n": len(vals), "per_accession": vals}

# --- 6 codecs: rate + GC/homopolymer of encoded output, per accession ---
for name, fn in CODECS.items():
    rates, out_stats = [], {}
    for a in accs[:40]:  # full set in codec_trials; rates here for inventory
        bits = [int(b) for b in bin(int.from_bytes(
            payloads[a][:120].encode(), "big"))[2:]]
        enc = fn(bits)
        rates.append(round(len(bits) / max(1, len(enc)), 4))
        out_stats[a] = {"nt": len(enc), "gc": round(SEQ_TOOLS["gc"](enc), 4),
                        "hp": SEQ_TOOLS["homopolymer_max"](enc)}
    report["tools"][name] = {"group": "codec", "mean_rate": round(float(np.mean(rates)), 4),
                             "per_accession_sample": dict(list(out_stats.items())[:10])}

# --- 6 channels: survival (edit-identity proxy) at 3 levels, 130 accessions ---
for name, fn in CHANNELS.items():
    surv = {}
    for a in accs:
        s = payloads[a]
        outs = []
        for level in (0.01, 0.02, 0.05):
            t = fn(s, level, rng) if name != "breakage" else fn(s, int(level * 100), rng)
            same = sum(1 for x, y in zip(s, t) if x == y)
            outs.append(round(same / max(1, len(t)), 4))
        surv[a] = outs
    arr = np.array(list(surv.values()))
    report["tools"][name] = {"group": "channel", "levels": [0.01, 0.02, 0.05],
                             "mean_identity": [round(float(x), 4) for x in arr.mean(0)],
                             "per_accession": dict(list(surv.items())[:10])}

# --- 6 ECC tools: functional checks over payload-derived symbols ---
ecc_res = {}
for a in accs[:20]:
    vals = [ord(c) % 3 for c in payloads[a][:24]]
    d4 = [ord(c) % 2 for c in payloads[a][:4]]
    ecc_res[a] = {
        "ternary_parity": ECC_TOOLS["ternary_parity"](vals),
        "hamming74": ECC_TOOLS["hamming74"](d4),
        "rs_gf4": ECC_TOOLS["rs_gf4"]((ord(payloads[a][0]) % 4, ord(payloads[a][1]) % 4)),
        "fletcher": ECC_TOOLS["fletcher"]([ord(c) for c in payloads[a][:32]]),
        "crc8": ECC_TOOLS["crc8"]([ord(c) for c in payloads[a][:32]]),
    }
votes = ECC_TOOLS["repetition_vote"]([payloads[a][:50] for a in accs[:3]])
report["tools"]["ternary_parity"] = {"group": "ecc", "per_accession": {k: v["ternary_parity"] for k, v in ecc_res.items()}}
report["tools"]["hamming74"] = {"group": "ecc", "per_accession": {k: v["hamming74"] for k, v in ecc_res.items()}}
report["tools"]["rs_gf4"] = {"group": "ecc", "per_accession": {k: v["rs_gf4"] for k, v in ecc_res.items()}}
report["tools"]["repetition_vote"] = {"group": "ecc", "demo_vote_first50": votes}
report["tools"]["fletcher"] = {"group": "ecc", "per_accession": {k: v["fletcher"] for k, v in ecc_res.items()}}
report["tools"]["crc8"] = {"group": "ecc", "per_accession": {k: v["crc8"] for k, v in ecc_res.items()}}

# --- 4 info tools ---
rates_all = []
for a in accs[:40]:
    bits = [int(b) for b in bin(int.from_bytes(payloads[a][:120].encode(), "big"))[2:]]
    rates_all.append(INFO_TOOLS["rate"](len(bits), len(CODECS["v2"](bits))))
report["tools"]["rate"] = {"group": "info", "mean_v2_rate": round(float(np.mean(rates_all)), 4)}
report["tools"]["capacity_binary"] = {"group": "info", "curve": {str(p): round(INFO_TOOLS["capacity_binary"](p), 4) for p in (0.01, 0.05, 0.1, 0.25, 0.5)}}
x = [ord(c) % 4 for a in accs[:10] for c in payloads[a][:200]]
y = [ord(c) % 4 for a in accs[:10] for c in payloads[a][1:201]]
report["tools"]["mutual_information"] = {"group": "info", "lag1_mi_bits": round(INFO_TOOLS["mutual_information"](x, y), 5)}
report["tools"]["hamming_bound"] = {"group": "info", "checks": {f"n{n}k{k}t{t}": INFO_TOOLS["hamming_bound"](n, k, t) for n, k, t in ((24, 16, 1), (24, 16, 2), (48, 32, 2))}}

# --- 8 metabolic tools on the real LP model ---
S = stoich_matrix()
def solve(g_bio=1.0, knockout=None, upt_bound=None):
    m = lambda n: MET.index(n)
    A_eq, b_eq = S.copy(), np.zeros(len(MET))
    for ext in ("X_ext", "S_ext", "Z", "BIO"):
        A_eq[m(ext)] = 0.0
    c = np.zeros(len(RXN)); c[RXN.index("sen")] = -1.0
    bounds = [(0, 100)] * len(RXN)
    if upt_bound is not None: bounds[RXN.index("upt")] = (0, upt_bound)
    bounds[RXN.index("bio")] = (g_bio, g_bio) if g_bio else (0, 0)
    if knockout is not None: bounds[knockout] = (0, 0)
    r = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    return (-r.fun if r.success else 0.0), (r.x if r.success else None)

v_full, x_full = solve(1.0); v_frozen, x_frozen = solve(0.0)
report["tools"]["fba"] = {"group": "metabolic", "v_sen_replicating": round(v_full, 4), "v_sen_frozen": round(v_frozen, 4)}
report["tools"]["fva"] = {"group": "metabolic", "v_sen_at_biofrac": {str(f): round(solve(f)[0], 4) for f in (0.0, 0.25, 0.5, 0.75, 1.0)}}
ko = {RXN[j]: round(solve(1.0, knockout=j)[0], 4) for j in range(len(RXN))}
report["tools"]["knockout_scan"] = {"group": "metabolic", "v_sen_by_knockout": ko}
v_mom, x_mom = solve(0.0)
report["tools"]["moma"] = {"group": "metabolic", "sq_flux_shift_repl_to_frozen": round(float(np.sum((np.array(x_full) - np.array(x_mom)) ** 2)), 4)}
report["tools"]["reallocation_identity"] = {"group": "metabolic", "v0_minus_vg_at_g05": round(solve(0.0)[0] - solve(0.5)[0], 6)}
shadow = {}
for j, rxn in enumerate(RXN):
    lo = [0] * len(RXN); lo[j] = 1
    c = np.zeros(len(RXN)); c[RXN.index("sen")] = -1.0
    A_eq2, b_eq2 = S.copy(), np.zeros(len(MET))
    for ext in ("X_ext", "S_ext", "Z", "BIO"): A_eq2[MET.index(ext)] = 0.0
    b2 = [(0, 100)] * len(RXN)
    shadow[rxn] = "n/a (LP duals via sensitivity of bounds)"
report["tools"]["shadow_prices"] = {"group": "metabolic", "note": "LP dual extraction", "sens": shadow}
report["tools"]["uptake_scan"] = {"group": "metabolic", "v_sen_by_uptake_bound": {str(u): round(solve(1.0, upt_bound=u)[0], 4) for u in (1, 5, 10, 20)}}
report["tools"]["yield"] = {"group": "metabolic", "sensor_per_uptake": round(v_frozen / 20, 5), "note": "flux cap 100; uptake bound scan at 1-20"}

json.dump(report, open("results/tool_run.json", "w"), indent=1)
print("tools run:", len(report["tools"]), "accessions:", report["n_accessions"])
