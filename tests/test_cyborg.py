import numpy as np
from dnacell.cyborg import steady_state, sensing_gain_bound, simulate


def test_steady_state_feasible_region():
    s = steady_state(mu_max=0.8, Ks=0.2, Y=0.5, D=0.3, S0=1.0, r=0.5)
    assert s["feasible"]
    # washout balance: (1-r) mu(S*) = D
    mu = 0.8 * s["S_star"] / (0.2 + s["S_star"]) * 0.5
    assert abs(mu - 0.3) < 1e-9


def test_steady_state_infeasible_when_growth_below_dilution():
    s = steady_state(mu_max=0.4, Ks=0.2, Y=0.5, D=0.5, S0=1.0, r=0.0)
    assert not s["feasible"]


def test_simulation_converges_to_analytic_steady_state():
    T, X, S = simulate(r=0.5, t_end=400)
    s = steady_state(mu_max=0.8, Ks=0.2, Y=0.5, D=0.3, S0=1.0, r=0.5)
    assert abs(X[-1] - s["x_star"]) < 0.02
    assert abs(S[-1] - s["S_star"]) < 0.02


def test_sensing_gain_bound_holds_for_high_alpha():
    realloc, base, holds = sensing_gain_bound(0.8, 0.2, 0.5, 0.3, 1.0,
                                              r=0.5, alpha=3.0)
    assert holds and realloc > base
    # and fails for negligible sensor conversion efficiency
    _, _, holds_low = sensing_gain_bound(0.8, 0.2, 0.5, 0.3, 1.0,
                                         r=0.5, alpha=0.01)
    assert not holds_low


def test_sensing_gain_bound_monotone_in_alpha():
    alphas = np.linspace(0.01, 4, 12)
    holds = [sensing_gain_bound(0.8, 0.2, 0.5, 0.3, 1.0, 0.5, a)[2] for a in alphas]
    assert holds == sorted(holds)  # once true, stays true
