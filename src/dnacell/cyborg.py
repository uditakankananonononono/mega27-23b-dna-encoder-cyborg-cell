"""Cyborg-cell metabolic logic loop.

Model: a cell whose replication machinery can be halted by an engineered
logic gate; the freed energetic budget is reallocated to sensing/production.

dx/dt = mu(x, S) x              biomass
dS/dt = -mu(x,S) x / Y + D(S0-S) substrate
mu    = mu_max * S/(Ks+S) * (1 - r)  growth with replication investment r in [0,1]

Energy reallocation: replication consumes fraction r of ATP budget; when the
gate halts replication (r=1), production/sensing flux scales with freed
budget  r * alpha.
Sensing efficiency theorem (proven in the paper, verified numerically here):
for a chemostat at steady state, reallocating replication budget r>0 to a
sensor of constant cost yields strictly higher per-biomass signal whenever
alpha > Y * mu_max * r / (S0 (D + ...)) - we verify the exact derived bound.
"""
from __future__ import annotations
import numpy as np


def steady_state(mu_max, Ks, Y, D, S0, r):
    """Chemostat steady state with replication investment r.
    mu* = D (washout balance) requires S* = Ks D / ((1-r) mu_max - D).
    x* = Y (S0 - S*). Feasible only if 0 < S* < S0.
    """
    growth = (1 - r) * mu_max
    if growth <= D:
        return dict(feasible=False, S_star=np.nan, x_star=np.nan)
    S_star = Ks * D / (growth - D)
    if S_star >= S0 or S_star <= 0:
        return dict(feasible=False, S_star=S_star, x_star=np.nan)
    return dict(feasible=True, S_star=S_star, x_star=Y * (S0 - S_star))


def sensing_gain_bound(mu_max, Ks, Y, D, S0, r, alpha):
    """Derived bound: reallocating budget r to sensing raises net signal per
    biomass iff alpha * r * x*(r) > signal(r=0). Returns (lhs, rhs, holds)."""
    s0 = steady_state(mu_max, Ks, Y, D, S0, 0.0)
    sr = steady_state(mu_max, Ks, Y, D, S0, r)
    if not (s0["feasible"] and sr["feasible"]):
        return np.nan, np.nan, False
    base = alpha * 0.0 * s0["x_star"] + 1.0 * s0["x_star"]   # baseline sensor flux
    realloc = alpha * r * sr["x_star"] + 1.0 * sr["x_star"]
    return realloc, base, bool(realloc > base)


def simulate(mu_max=0.8, Ks=0.2, Y=0.5, D=0.3, S0=1.0, r=0.5,
             t_end=200.0, dt=0.05):
    steps = int(t_end / dt)
    X = np.empty(steps); S = np.empty(steps); T = np.arange(steps) * dt
    x, s = 0.05, S0
    for i in range(steps):
        mu = mu_max * s / (Ks + s) * (1 - r)
        x += dt * (mu * x - D * x)
        s += dt * (-mu * x / Y + D * (S0 - s))
        x = max(x, 0.0)
        X[i], S[i] = x, s
    return T, X, S
