"""FBA-style energy reallocation model for the cyborg-cell logic loop.

Minimal but stoichiometrically real core network (lumped reactions in the
spirit of Orth et al. 2010 core E. coli):

  v_upt : S_ext -> S                     (substrate uptake, bounded)
  v_gly : S + 2 ADP + 2 NAD+ -> 2 P + 2 ATP + 2 NADH       (glycolysis)
  v_tca : P + 4 NAD+ + ADP -> 4 NADH + ATP                 (TCA lump, CO2 omitted)
  v_etc : NADH + 2.5 ADP -> NAD+ + 2.5 ATP                 (ETC, P/O = 2.5)
  v_atm : ATP -> ADP                     (maintenance, fixed >= atm)
  v_bio : a_bio ATP + b_bio P + c_bio NADH -> biomass      (replication drain)
  v_sen : X_ext + k_sen ATP -> Z         (complex-input processing = sensor)

Per glucose at full oxidation: 2 (gly) + 2x(1 + 4x2.5) = 24 ATP - within the
canonical 26-30 range for lumped core models.

Steady state S v = 0. Objective: maximize v_sen.
The "cyborg" condition is v_bio = 0 (replication frozen).
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import linprog

# metabolites: S, P, ATP, ADP, NADH, NAD+, X_ext, S_ext, Z, biomass
MET = ["S", "P", "ATP", "ADP", "NADH", "NAD", "X_ext", "S_ext", "Z", "BIO"]
RXN = ["upt", "gly", "tca", "etc", "atm", "bio", "sen"]
A_BIO, B_BIO, C_BIO = 30.0, 4.0, 6.0   # ATP / pyruvate / NADH cost per unit biomass
K_SEN = 2.0                            # ATP cost per unit complex input processed


def stoich_matrix() -> np.ndarray:
    S = np.zeros((len(MET), len(RXN)))
    def m(name): return MET.index(name)
    S[m("S_ext"), 0], S[m("S"), 0] = -1, 1                       # upt
    S[m("S"), 1], S[m("P"), 1] = -1, 2                           # gly
    S[m("ADP"), 1], S[m("ATP"), 1] = -2, 2
    S[m("NAD"), 1], S[m("NADH"), 1] = -2, 2
    S[m("P"), 2], S[m("NAD"), 2] = -1, -4                        # tca
    S[m("ADP"), 2], S[m("ATP"), 2] = -1, 1
    S[m("NADH"), 2] = 4
    S[m("NADH"), 3], S[m("NAD"), 3] = -1, 1                      # etc
    S[m("ADP"), 3], S[m("ATP"), 3] = -2.5, 2.5
    S[m("ATP"), 4], S[m("ADP"), 4] = -1, 1                       # atm
    S[m("ATP"), 5], S[m("P"), 5] = -A_BIO, -B_BIO                # bio
    S[m("NADH"), 5], S[m("BIO"), 5] = -C_BIO, 1
    S[m("NAD"), 5] = C_BIO  # reductant is returned oxidized
    S[m("ADP"), 5] = A_BIO  # ATP hydrolysis returns ADP
    S[m("X_ext"), 6], S[m("Z"), 6] = -1, 1                       # sen
    S[m("ATP"), 6] = -K_SEN
    S[m("ADP"), 6] = K_SEN
    return S


def solve_fba(g_min: float, uptake_max: float = 10.0, atm: float = 1.0,
              x_max: float = 1e6):
    """Maximize sensor flux s.t. S v = 0, biomass >= g_min. Returns (v, result)."""
    S = stoich_matrix()
    c = np.zeros(len(RXN)); c[RXN.index("sen")] = -1.0
    A_eq, b_eq = S.copy(), np.zeros(len(MET))
    # external metabolites are not balanced
    for ext in ("S_ext", "X_ext", "Z", "BIO"):
        A_eq[MET.index(ext)] = 0.0
    bounds = [(0, uptake_max), (0, None), (0, None), (0, None),
              (atm, None), (g_min, None), (0, x_max)]
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    assert res.status == 0, res.message
    return res.x, res


def replication_tax(g: float) -> float:
    """ATP equivalents drained per unit time by replication at rate g."""
    return A_BIO * g
