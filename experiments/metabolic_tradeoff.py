"""Formal metabolic-coupling tradeoff for storage DNA (verdict #3, #10).

Model. A storage payload of L information bases, expanded to N = L / rate codec
bases, carried in one E. coli cell. Per-base costs in ATP-equivalents
(literature-anchored parameters, stated in the paper):
  c_nt  = 40   ATP-eq per base per replication (dNTP de novo synthesis,
                 Rocha & Danchin 2002 base-composition cost argument)
  c_GC  = 6    extra ATP-eq per G/C base (GC-cost differential, same source)
  c_tx  = 2    ATP-eq per transcribed base per generation (NTP + polymerase)
  c_tl  = 4    ATP-eq per translated codon per generation (charged tRNAs + EF)
Total per-generation cost:
  C(L, g, r_tx, r_tl, rate) = (L/rate) * [c_nt + g*c_GC + r_tx*c_tx + r_tl*3*c_tl]

Budget. E. coli K-12 genome 4.6e6 bp; endogenous replication cost
  B_genome = 4.6e6 * (40 + 0.5*6) ~ 1.98e8 ATP-eq.
Budget for storage: B = phi * B_genome, phi swept 1e-4 .. 1e-1.

Optimization. For each budget B and recovery target (>=90% at 2% substitution;
>=90% at 3%), choose the codec design from rs_rate_sweep.json and payload L
maximizing stored bits = L * density subject to C <= B and recovery >= target.
Output: design oracle per budget, marginal bits per ATP-eq, and the chain
formalization (storage sequence -> expression burden -> metabolic cost) with
the silent-locus collapse (r_tx = r_tl = 0) as the design recommendation.
"""
import json
import numpy as np

C_NT, C_GC, C_TX, C_TL = 40.0, 6.0, 2.0, 4.0
GENOME_BP = 4.6e6
B_GENOME = GENOME_BP * (C_NT + 0.5 * C_GC)

sweep = json.load(open("results/rs_rate_sweep.json"))
designs = []
for name, d in sweep.items():
    rate = float(name.split("rate ")[1])
    designs.append({"name": name, "rate": rate, "density": d["density"],
                    "rec": d["recovery"]})

def best_design(target_err, target_rec):
    ok = [d for d in designs if d["rec"].get(str(target_err), 0) >= target_rec]
    return max(ok, key=lambda d: d["density"]) if ok else None

def cost_per_info_base(g, r_tx, r_tl, rate):
    return (C_NT + g * C_GC + r_tx * C_TX + r_tl * 3 * C_TL) / rate

rows = []
for target_err, tname in ((0.02, "90%@2%subst"), (0.03, "90%@3%subst")):
    d = best_design(target_err, 0.9)
    for phi in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1):
        B = phi * B_GENOME
        if d is None:
            rows.append(dict(target=tname, phi=phi, bits=0.0, design=None))
            continue
        for scen, (g, r_tx, r_tl) in (("silent", (0.5, 0.0, 0.0)),
                                      ("transcribed", (0.5, 1.0, 0.0)),
                                      ("expressed", (0.5, 1.0, 1.0))):
            c = cost_per_info_base(g, r_tx, r_tl, d["rate"])
            L = B / c
            rows.append(dict(target=tname, phi=phi, scenario=scen,
                             design=d["name"], density=d["density"],
                             payload_bases=L, bits=L * d["density"],
                             cost_per_bit=c / d["density"]))

# headline numbers
def get(target, phi, scen):
    return next(r for r in rows if r["target"] == target and abs(r["phi"]-phi) < 1e-12
                and r.get("scenario") == scen)

out = {
    "parameters": {"c_nt": C_NT, "c_GC": C_GC, "c_tx": C_TX, "c_tl": C_TL,
                   "genome_bp": GENOME_BP, "B_genome_atp_eq": B_GENOME},
    "oracle_2pct": best_design(0.02, 0.9)["name"],
    "oracle_3pct": best_design(0.03, 0.9),
    "bits_1pct_budget_silent_2pct": get("90%@2%subst", 1e-2, "silent")["bits"],
    "bits_1pct_budget_transcribed_2pct": get("90%@2%subst", 1e-2, "transcribed")["bits"],
    "bits_1pct_budget_expressed_2pct": get("90%@2%subst", 1e-2, "expressed")["bits"],
    "cost_per_bit_silent": get("90%@2%subst", 1e-2, "silent")["cost_per_bit"],
    "cost_per_bit_expressed": get("90%@2%subst", 1e-2, "expressed")["cost_per_bit"],
    "rows": rows,
}
json.dump(out, open("results/metabolic_tradeoff.json", "w"), indent=1)
print("oracle 90%@2%:", out["oracle_2pct"], "| oracle 90%@3%:", out["oracle_3pct"])
print("bits at phi=1% (silent/transcribed/expressed):",
      round(out["bits_1pct_budget_silent_2pct"]),
      round(out["bits_1pct_budget_transcribed_2pct"]),
      round(out["bits_1pct_budget_expressed_2pct"]))
print("cost per bit silent vs expressed (ATP-eq):",
      round(out["cost_per_bit_silent"], 2), round(out["cost_per_bit_expressed"], 2))
