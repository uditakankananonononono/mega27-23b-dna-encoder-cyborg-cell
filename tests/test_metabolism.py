import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from dnacell.metabolism import (solve_fba, stoich_matrix, RXN, MET,
                                A_BIO, B_BIO, C_BIO, K_SEN)

E_P, E_H = 11.0, 2.5

def test_steady_state_residual():
    S = stoich_matrix()
    v, _ = solve_fba(0.5)
    internal = [i for i, m in enumerate(MET) if m not in ("S_ext", "X_ext", "Z", "BIO")]
    assert np.allclose((S @ v)[internal], 0, atol=1e-8)

def test_freezing_replication_increases_sensor_flux():
    assert solve_fba(0.0)[0][RXN.index("sen")] > solve_fba(1.0)[0][RXN.index("sen")]

def test_reallocation_identity_exact():
    v0 = solve_fba(0.0)[0][RXN.index("sen")]
    for g in (0.5, 1.0, 2.0):
        vg = solve_fba(g)[0][RXN.index("sen")]
        assert abs((v0 - vg) - g * (A_BIO + B_BIO * E_P + C_BIO * E_H) / K_SEN) < 1e-9

def test_bounded_input_breaks_full_reallocation():
    v0 = solve_fba(0.0)[0][RXN.index("sen")]
    vc = solve_fba(0.0, x_max=v0 - 10)[0][RXN.index("sen")]
    vg = solve_fba(1.0)[0][RXN.index("sen")]
    assert (vc - vg) * K_SEN < A_BIO + B_BIO * E_P + C_BIO * E_H - 1e-6
