"""The full model (surface forcing + firn), its net balance and the equilibrium-line altitude."""

import numpy as np
from scipy.optimize import brentq

from .melt import DAYS_IN_MONTH, surface_forcing
from .params import resolve
from .refreeze import firn

T_RANGE = (-100.0, 40.0)            # C, annual mean temperatures searched for the ELA
_Z_10M = np.array([0.0, 10.0])      # enough depth grid for bdot


def run(T_m, P_m, SW_m=None, p=None, z=None, regime=None, **params):
    """
    SnowDRIP for one climatological year: monthly mean temperature T_m [C], precipitation P_m [m w.e./yr] and
    shortwave radiation SW_m [W m-2] at the elevation of interest (12 values each). Returns the outputs of
    surface_forcing() and firn() in one dict. Parameters: p (default Params()), with keyword overrides,
    e.g. run(T_m, P_m, SW_m, p_max=0.55).
    """
    p = resolve(p, params)
    s = surface_forcing(T_m, P_m, SW_m, p)
    return s | firn(s["S"], s["M"], s["T_sfc"], s["f_c"], s["R"], p, z, regime)


def _full_balance(s, p):
    if s["S"] <= 1e-9:
        return -s["M"]
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        return firn(s["S"], s["M"], s["T_sfc"], s["f_c"], s["R"], p, z=_Z_10M)["bdot"]


def _cold_balance(s, p):
    return s["S"] - s["M"] + min(s["M"] + s["R"], p.p_max * s["S"])


def _ela_root(forcing_at, lo, hi, step, p, xtol):
    """Root of bdot in x for forcing that warms with x. The cold branch retains at least as much as the full
    model, so its (cheap) root bounds the full-model root from above. NaN if there is no root in [lo, hi]."""
    cold = lambda x: _cold_balance(forcing_at(x), p)
    full = lambda x: _full_balance(forcing_at(x), p)
    if cold(lo) <= 0.0 or cold(hi) > 0.0:
        return np.nan
    hi = brentq(cold, lo, hi, xtol=xtol)
    if full(hi) >= 0.0:
        return hi
    a = max(hi - step, lo)
    while full(a) <= 0.0:
        if a <= lo:
            return np.nan
        a = max(a - step, lo)
    return brentq(full, a, hi, xtol=xtol)


def net_balance(T_m, P_m, SW_m=None, p=None, **params):
    """Net balance bdot [m w.e./yr] for a climatological year (arguments as run())."""
    p = resolve(p, params)
    return _full_balance(surface_forcing(T_m, P_m, SW_m, p), p)


def ela_offset(T_m, P_m, SW_m=None, p=None, xtol=0.01, **params):
    """Uniform temperature shift dT [K] of the climatological year that gives bdot = 0."""
    p = resolve(p, params)
    T_m = np.asarray(T_m, float)
    if not np.sum(P_m) > 0.0:
        return np.nan
    T_mu = np.sum(DAYS_IN_MONTH * T_m) / DAYS_IN_MONTH.sum()
    at = lambda dT: surface_forcing(T_m + dT, P_m, SW_m, p)
    return _ela_root(at, T_RANGE[0] - T_mu, T_RANGE[1] - T_mu, 5.0, p, xtol)


def ela_elevation(T_m, z_ref, lapse_rate, P_m, SW_m=None, p=None, xtol=1.0, **params):
    """
    ELA [m] from monthly temperatures T_m at z_ref [m], moved with lapse_rate [K/km, positive; one or 12
    values]. P_m is constant with elevation; SW_m is 12 values or a function of elevation returning them.
    """
    p = resolve(p, params)
    T_m = np.asarray(T_m, float)
    gam = np.broadcast_to(np.asarray(lapse_rate, float), T_m.shape) / 1000.0
    if not np.sum(P_m) > 0.0 or np.any(gam <= 0.0):
        return np.nan
    w = DAYS_IN_MONTH / DAYS_IN_MONTH.sum()
    T_mu, g = np.sum(w * T_m), np.sum(w * gam)
    sw = SW_m if callable(SW_m) else (lambda z: SW_m)
    at = lambda x: surface_forcing(T_m + gam * (x + z_ref), P_m, sw(-x), p)   # x = -z
    x_lo, x_hi = -(z_ref + (T_mu - T_RANGE[0]) / g), -(z_ref + (T_mu - T_RANGE[1]) / g)
    return -_ela_root(at, x_lo, x_hi, 5.0 / g, p, xtol)
