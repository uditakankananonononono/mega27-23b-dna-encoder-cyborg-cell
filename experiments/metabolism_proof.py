"""Numerical verification of the Reallocation Identity (cyborg-cell theorem).

THEOREM (Reallocation Identity). In the steady-state LP of
dnacell.metabolism with objective max v_sen, let v*(g) be the optimum with
biomass floor g. If the optimal basis is constant on [0, g], then

    v*(0) - v*(g) = g * (a_bio + b_bio*e_P + c_bio*e_H) / k_sen

where e_P = 1 + R*(P/O) is the ATP yield of one pyruvate (TCA + ETC,
R = 4 NADH/TCA turn, P/O = 2.5) and e_H = P/O is the ATP yield of NADH.

Proof sketch (full version in the paper): the biomass reaction consumes
a_bio ATP, b_bio P, c_bio NADH per unit. Fixing its flux to zero frees those
resources; under a constant optimal basis, each freed P is oxidized for
e_P ATP and each freed NADH for e_H ATP; all freed ATP flows to v_sen at
cost k_sen per unit (LP strong duality / shadow-price constancy on the
basis interval). QED.

COROLLARY (100% reallocation). With X_ext unbounded and maintenance fixed,
100% of the freed ATP-equivalents are reallocated to complex-input
processing. Counterexample when it fails: X_ext bounded - excess freed
energy has no sink and the reallocated fraction is x_max*k_sen / (freed ATP).
"""
import json, os, sys
import numpy as np
sys.path.insert(0, "src")
from dnacell.metabolism import solve_fba, A_BIO, B_BIO, C_BIO, K_SEN, RXN

R, PO = 4.0, 2.5
E_P = 1 + R * PO      # 11 ATP per pyruvate
E_H = PO              # 2.5 ATP per NADH

def identity_rhs(g):
    return g * (A_BIO + B_BIO * E_P + C_BIO * E_H) / K_SEN

def main():
    rows = []
    for g in (0.25, 0.5, 1.0, 1.5, 2.0):
        vg, _ = solve_fba(g)
        v0, _ = solve_fba(0.0)
        lhs = v0[RXN.index("sen")] - vg[RXN.index("sen")]
        rhs = identity_rhs(g)
        rows.append(dict(g=g, delta_sen=float(lhs), identity=float(rhs),
                         abs_err=float(abs(lhs - rhs))))
        print(f"g={g}: delta_sen={lhs:.6f} identity={rhs:.6f} |err|={abs(lhs-rhs):.2e}")

    # 100% reallocation check: freed ATP equivalents vs extra sensor ATP spend
    v0, _ = solve_fba(0.0)
    vg, _ = solve_fba(1.0)
    freed = A_BIO + B_BIO * E_P + C_BIO * E_H
    spent = (v0[RXN.index("sen")] - vg[RXN.index("sen")]) * K_SEN
    pct = 100 * spent / freed
    print(f"freed ATP-eq={freed:.3f} extra sensor ATP={spent:.3f} reallocated={pct:.4f}%")

    # counterexample: bounded complex input
    vcap, _ = solve_fba(0.0, x_max=v0[RXN.index("sen")] - 10)
    capped_gain = (vcap[RXN.index("sen")] - vg[RXN.index("sen")]) * K_SEN
    print(f"counterexample: x_max caps reallocation at {100*capped_gain/freed:.2f}% (<100%)")

    os.makedirs("results", exist_ok=True)
    json.dump(dict(rows=rows, freed_atp_eq=freed, reallocated_pct=pct),
              open("results/metabolism_identity.json", "w"), indent=1)
    print("saved results/metabolism_identity.json")

if __name__ == "__main__":
    main()
