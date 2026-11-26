"""
Steady-state refreezing, temperature and density of the near-surface firn (see the SnowDRIP v1.0 paper).

Mass
    Retained water W = min(M + R, p_max S): limited by the water or by the pore space of the year's snow.
    Net accumulation bdot = S - M + W; runoff = M + R - W.

Cold firn
    Water refreezes in two reservoirs: shallow, at the weather-scale skin depth ell_s (refrozen by the next
    cold spell), and deep, at ell_d = H_snow / 2, beyond the reach of the winter cold wave (skin depth ell_c).
    The deep share is s(h) = 1 - (1 - exp(-h))/h, h = H_snow / ell_c, scaled towards the hard cut-off
    max(0, 1 - 1/h) by the open porosity phi of the firn below. Each source warms the firn as

        T(z) - T_sfc = dT (1 - exp(-z/ell)),  dT = W_i (L/c + |T_sfc|) / (bdot + G(ell)),

    G(ell) = k/ell yr / (rho_w c) the conductive cold supply as an equivalent burial rate.

Temperate firn
    0 C below ell_d, T_sfc exp(-z/ell_d) above. Retention is what offsets the cold supply,
    c|T_sfc|/L (bdot + G(ell_d)), plus the capillary water carried down; the rest drains.

Between the two
    The water excess of each state, E_C (cold) and E_T (temperate), decides the regime; where they disagree
    the branches are blended with weight max(E_C, E_T) / (|E_C| + |E_T|).

Density
    rho_deep = rho_snow + (rho_ice - rho_snow) W / bdot, approached exponentially with e-folding depth ell_d.
    Dry densification is left to the user.
"""

import warnings

import numpy as np

from .params import CP, LF, RHO_CO, RHO_I, RHO_SNOW, RHO_WATER, SEC_PER_YEAR, resolve

W_TOL = 1e-6                # m w.e./yr, convergence of the temperate retention


def equivalent_burial(k, ell):
    """Burial rate [m w.e./yr] that brings the same cold as conduction across depth ell."""
    return k / ell * SEC_PER_YEAR / (RHO_WATER * CP)


def skin_depth(k, f_c=1.0):
    """Depth [m] winter cold reaches into snow during a cold season of f_c years."""
    return np.sqrt(k / (RHO_SNOW * CP) * f_c * SEC_PER_YEAR / np.pi)


def shallow_depth(tau_w_days, k_snow):
    """Depth [m] a cold spell of tau_w days reaches into snow."""
    return np.sqrt(k_snow / (RHO_SNOW * CP) * tau_w_days * 86400.0 / np.pi)


def snow_conductivity(rho):
    """Conductivity of snow [W/m/K], Calonne et al. (2011)."""
    return 0.024 - 1.23e-4 * rho + 2.5e-6 * rho**2


def capillary_export(rho, bdot):
    """Irreducible water carried down with the firn [m w.e./yr], Coleou and Lesaffre (1998)."""
    wmi = 0.057 * (RHO_I - rho) / rho + 0.017
    return wmi / (1.0 - wmi) * bdot


def refrozen_share(W, bdot):
    """W / bdot, capped at 1; 1 where nothing survives."""
    return min(W / bdot, 1.0) if bdot > 0.0 else 1.0


def deep_share(h, phi, f_pf=0.0):
    """Share of the water held beyond the reach of the winter cold, for snowpack thickness h (in skin depths)
    and open porosity phi of the firn below; f_pf lets a further fraction bypass the cold as preferential flow."""
    s_closed = max(0.0, 1.0 - 1.0 / h)
    s_open = 1.0 - (1.0 - np.exp(-h)) / h
    s = (1.0 - phi) * s_closed + phi * s_open
    return s + f_pf * (1.0 - s)


def temperate_weight(E_C, E_T):
    """Weight of the temperate branch from the water excess of the cold and the temperate state."""
    if E_C >= 0 and E_T >= 0:
        return 1.0
    if E_C < 0 and E_T < 0:
        return 0.0
    return max(E_C, E_T) / (abs(E_C) + abs(E_T))


def ablation_column_density(z):
    """Density [kg/m3] below the ELA: solid ice."""
    return np.full_like(np.asarray(z, dtype=float), RHO_I)


def firn(S, M, T_sfc, f_c, rain=0.0, p=None, z=None, regime=None, **params):
    """
    The firn model for snowfall S, melt M and rain [m w.e./yr], surface temperature T_sfc = mean(min(T, 0)) [C]
    and cold-season fraction f_c (surface_forcing() gives them). z: depth grid [m] (default 0-20 m); regime:
    None (blend), "cold" or "temperate". Parameters: p (default Params()), with keyword overrides.

    Returns a dict: retained water W, runoff, net accumulation bdot, ablation, T_deep, T10m, rho_deep, the
    profiles z, T(z), rho(z), the depth scales ell_c, ell_s, ell_d, H_snow, xi = c|T_sfc|/L, and the regime
    diagnostics w_T, temperate, E_C, E_T, phi, deep_share. Below the ELA (firn False) retention follows the
    seasonal-snow rule, density is ice and temperature and the regime diagnostics are NaN.
    """
    p = resolve(p, params)
    z = np.arange(0.0, 20.0 + 1e-9, 0.1) if z is None else z
    H_snow = S * RHO_WATER / RHO_SNOW
    k_snow = snow_conductivity(RHO_SNOW)
    ell_c = skin_depth(k_snow, f_c)
    h = H_snow / ell_c
    cold_per_water = CP * -T_sfc / LF
    ell_s = shallow_depth(p.tau_w, k_snow)
    ell_d = max(0.5 * H_snow, ell_s)

    water = M + rain
    W_c = min(water, p.p_max * S)
    bdot_c = S - M + W_c
    rho_deep_c = RHO_SNOW + (RHO_I - RHO_SNOW) * refrozen_share(W_c, bdot_c)
    phi = float(np.clip((RHO_CO - rho_deep_c) / (RHO_CO - RHO_SNOW), 0.0, 1.0))

    def density(W, bdot):
        rho_deep = RHO_SNOW + (RHO_I - RHO_SNOW) * refrozen_share(W, bdot)
        return rho_deep, rho_deep - (rho_deep - RHO_SNOW) * np.exp(-np.asarray(z, float) / ell_d)

    # cold branch; retained rain refreezes deep when the snowpack exceeds the cold wave
    Lam = LF / CP - T_sfc
    share = deep_share(h, phi, p.f_pf)
    W_rain = W_c * rain / water if water > 0 else 0.0
    W_deep = share * (W_c - W_rain) + (W_rain if h > 1 else share * W_rain)
    dT_s = (W_c - W_deep) * Lam / (bdot_c + equivalent_burial(k_snow, ell_s))
    dT_d = W_deep * Lam / (bdot_c + equivalent_burial(p.k, ell_d))
    T_deep_c = min(T_sfc + dT_s + dT_d, 0.0)
    T_c = np.minimum(T_deep_c - dT_s * np.exp(-z / ell_s) - dT_d * np.exp(-z / ell_d), 0.0)
    cold = dict(W=W_c, bdot=bdot_c, T_deep=T_deep_c, T=T_c)

    # temperate branch
    G_T = equivalent_burial(p.k, ell_d)
    W, bdot = W_c, bdot_c
    for _ in range(100):
        rho_deep = RHO_SNOW + (RHO_I - RHO_SNOW) * refrozen_share(W, bdot)
        burial = max(bdot, 0.0)
        W_new = min(water, p.p_max * S, cold_per_water * (burial + G_T) + capillary_export(rho_deep, burial))
        bdot = S - M + W_new
        if abs(W_new - W) < W_TOL:
            break
        W = W_new
    else:
        warnings.warn(f"temperate retention did not converge (S={S}, M={M}, T_sfc={T_sfc})")
    temp = dict(W=W_new, bdot=bdot, T_deep=0.0, T=T_sfc * np.exp(-z / ell_d))

    # water excess of each state, on the water that reaches the firn within ell_c
    f_new = min(1.0, h)
    f_a = (f_new + (1.0 - f_new) * phi) * (W_c / water if water > 0 else 1.0)
    E_C = f_a * W_c / (cold_per_water * (max(bdot_c, 0.0) + equivalent_burial(p.k, ell_c))) - 1.0
    E_T = f_a * W_c / (cold_per_water * (max(temp["bdot"], 0.0) + G_T)) - 1.0
    w_T = {"cold": 0.0, "temperate": 1.0}.get(regime, temperate_weight(E_C, E_T))

    def mix(a, b):
        return (1.0 - w_T) * a + w_T * b

    rho_deep_c, rho_c = density(cold["W"], cold["bdot"])
    rho_deep_t, rho_t = density(temp["W"], temp["bdot"])
    W, T, bdot, T_deep = (mix(cold[k], temp[k]) for k in ("W", "T", "bdot", "T_deep"))
    rho, rho_deep = mix(rho_c, rho_t), mix(rho_deep_c, rho_deep_t)

    firn = bdot > 0.0
    if not firn:
        # seasonal snow on ice: the temperate retention at bdot = 0 (continuous at the ELA)
        W = min(water, p.p_max * S, cold_per_water * G_T)
        bdot = S - M + W
        rho, rho_deep = ablation_column_density(z), RHO_I
        T, T_deep = np.full_like(z, np.nan, dtype=float), np.nan
        w_T = E_C = E_T = phi = share = np.nan
    return dict(firn=firn, temperate=bool(firn and w_T >= 0.5), w_T=w_T, E_C=E_C, E_T=E_T, phi=phi,
                deep_share=share, W=W, runoff=water - W, bdot=bdot, ablation=max(-bdot, 0.0), H_snow=H_snow,
                xi=cold_per_water, ell_c=ell_c, ell_s=ell_s, ell_d=ell_d, T_deep=T_deep,
                T10m=np.interp(10.0, z, T), rho_deep=rho_deep, z=z, T=T, rho=rho)
