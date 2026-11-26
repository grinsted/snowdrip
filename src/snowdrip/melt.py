"""
Temperature-index melt over a climatological year.

Daily melt on surface s (snow, ice) is m_s = max(0, beta_T_s T + beta_SW_s SW + beta_0_s): degree days above
the threshold T*_s = -(beta_SW_s SW + beta_0_s) / beta_T_s. Over Gaussian daily weather (p.sigma_w) about each
monthly mean the melt potential is M_s = beta_T_s PDD(T*_s). The surface melts as snow until melt and rain fill
the pores of the year's snow, A = max(0, p_max S - R), and as ice after that:

    M = min(M_snow, A) + max(0, 1 - A / M_snow) M_ice.
"""

import numpy as np
from scipy.special import ndtr

from .params import resolve

DAYS_IN_MONTH = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
_SQRT_2PI = np.sqrt(2 * np.pi)


def _pdf(x):
    return np.exp(-x ** 2 / 2.0) / _SQRT_2PI


def _scalar(x):
    return float(x) if np.ndim(x) == 0 else x


def _shortwave_or_zero(SW, p):
    if SW is None:
        if p.beta_SW_snow != 0.0 or p.beta_SW_ice != 0.0:
            raise ValueError("shortwave radiation is required when the shortwave coefficients are non-zero "
                             "(or use degree_day_params() for the degree-day form)")
        return 0.0
    return np.asarray(SW, float)


def degree_days(T_m, sigma_w, T_star=0.0):
    """Annual degree days [K day/yr] above T_star [C] for Gaussian daily temperatures about the monthly
    means T_m (..., 12)."""
    d = np.asarray(T_m, float) - T_star
    y = d / sigma_w
    return _scalar(np.sum(DAYS_IN_MONTH * (sigma_w * _pdf(y) + d * ndtr(y)), axis=-1))


def _combine(M_snow, M_ice, S, R, p):
    A = np.maximum(0.0, p.p_max * np.asarray(S, float) - R)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(M_snow <= A, M_snow, A + (1.0 - A / M_snow) * M_ice)


def _melt_monthly(T_m, S, R, SW_m, p):
    """Annual melt [m w.e./yr] from monthly means T_m [C] and SW_m [W m-2], (..., 12), snowfall S and rain R."""
    SW_m = _shortwave_or_zero(SW_m, p)
    M_snow = p.beta_T_snow * degree_days(T_m, p.sigma_w, -(p.beta_SW_snow * SW_m + p.beta_0_snow) / p.beta_T_snow)
    M_ice = p.beta_T_ice * degree_days(T_m, p.sigma_w, -p.beta_SW_ice * SW_m / p.beta_T_ice)
    return _scalar(_combine(M_snow, M_ice, S, R, p))


def melt_daily(T_d, S, R=0.0, SW_d=None, p=None, n_years=1.0, **params):
    """Annual melt [m w.e./yr] from daily temperatures T_d [C] and shortwave SW_d [W m-2] covering n_years, with
    annual snowfall S and rain R [m w.e./yr]."""
    p = resolve(p, params)
    T_d, SW_d = np.asarray(T_d, float), _shortwave_or_zero(SW_d, p)
    M_snow = np.sum(np.maximum(p.beta_T_snow * T_d + p.beta_SW_snow * SW_d + p.beta_0_snow, 0.0)) / n_years
    M_ice = np.sum(np.maximum(p.beta_T_ice * T_d + p.beta_SW_ice * SW_d, 0.0)) / n_years
    return float(_combine(M_snow, M_ice, S, R, p))


def surface_forcing(T_m, P_m, SW_m=None, p=None, **params):
    """
    Inputs of firn() from a climatological year: monthly mean temperature T_m [C], precipitation P_m
    [m w.e./yr] and shortwave radiation SW_m [W m-2] at the elevation of interest, (..., 12), with daily weather
    of spread sigma_w about each monthly mean. Returns snowfall S, rain R and melt M [m w.e./yr],
    T_sfc = mean(min(T, 0)) [C], the cold-season fraction f_c and the annual mean T_mu [C].
    Parameters: p (default Params()), with keyword overrides, e.g. surface_forcing(T_m, P_m, SW_m, sigma_w=2.5).
    """
    p = resolve(p, params)
    T_m, P_m = np.asarray(T_m, float), np.asarray(P_m, float)
    w = DAYS_IN_MONTH / DAYS_IN_MONTH.sum()
    sigma_w = p.sigma_w
    x = T_m / sigma_w
    T_sfc = np.sum(w * (T_m * ndtr(-x) - sigma_w * _pdf(x)), axis=-1)
    f_c = np.sum(w * ndtr(-x), axis=-1)
    R = np.sum(w * P_m * ndtr(-((p.rain_threshold - T_m) / sigma_w)), axis=-1)
    S = np.sum(w * P_m, axis=-1) - R
    M = _melt_monthly(T_m, S, R, SW_m, p)
    return dict(S=S, R=R, M=M, T_sfc=T_sfc, f_c=f_c, T_mu=np.sum(w * T_m, axis=-1))
