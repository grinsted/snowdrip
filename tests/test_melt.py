import numpy as np
import pytest
from scipy.stats import norm

import snowdrip as sd
from snowdrip.melt import _melt_monthly
from cases import COLD_T_M, SIGMA, SW_SUMMER, TEMPERATE_T_M


def degree_day_melt(PDD, S, R, ddf_snow, ddf_ice, p_max):
    """Sequential degree-day melt from a PDD sum; snow until its pores are full, then ice."""
    PDD_snow = max(0.0, (p_max * S - R) / ddf_snow)
    return ddf_snow * min(PDD, PDD_snow) + ddf_ice * max(0.0, PDD - PDD_snow)


def test_degree_days_threshold():
    T = TEMPERATE_T_M
    x = T / SIGMA
    pdd = np.sum(sd.DAYS_IN_MONTH * (SIGMA * norm.pdf(x) + T * norm.cdf(x)))
    assert sd.degree_days(T, SIGMA) == pytest.approx(pdd, rel=1e-12)
    assert sd.degree_days(T, SIGMA, -3.0) == pytest.approx(sd.degree_days(T + 3.0, SIGMA), rel=1e-12)


@pytest.mark.parametrize("S, R", [(0.0, 0.0), (0.2, 0.05), (1.0, 0.04), (3.0, 0.5)])
@pytest.mark.parametrize("dT", [-6.0, 0.0, 4.0])
def test_melt_degree_day_form(S, R, dT):
    """With no shortwave term the melt is the degree-day scheme (pores full -> ice factor)."""
    T = TEMPERATE_T_M + dT
    p = sd.degree_day_params(sd.Params(sigma_w=SIGMA))
    expected = degree_day_melt(sd.degree_days(T, SIGMA), S, R, 0.004, 0.008, p.p_max)
    assert _melt_monthly(T, S, R, None, p) == pytest.approx(expected, rel=1e-12, abs=1e-15)


def test_melt_shortwave():
    T, S, R = COLD_T_M, 4.0, 0.04                                 # thick snow: its pores do not fill
    p = sd.Params(sigma_w=SIGMA)
    with pytest.raises(ValueError):
        _melt_monthly(T, S, R, None, p)                    # SW is required
    m0 = _melt_monthly(T, S, R, np.zeros(12), p)
    m1 = _melt_monthly(T, S, R, SW_SUMMER, p)
    m2 = _melt_monthly(T, S, R, 1.5 * SW_SUMMER, p)
    assert 0.0 < m0 < m1 < m2                                     # sunshine adds melt
    T_star = -(p.beta_SW_snow * SW_SUMMER + p.beta_0_snow) / p.beta_T_snow   # snow only (pores not full here)
    assert m1 == pytest.approx(p.beta_T_snow * sd.degree_days(T, SIGMA, T_star), rel=1e-12)
    # pores full: the rest melts as ice, with T*_ice = 0 (beta_SW_ice = 0)
    S2 = 0.2
    A = p.p_max * S2 - R
    Ms = p.beta_T_snow * sd.degree_days(TEMPERATE_T_M, SIGMA, -(p.beta_SW_snow * SW_SUMMER + p.beta_0_snow) / p.beta_T_snow)
    Mi = p.beta_T_ice * sd.degree_days(TEMPERATE_T_M, SIGMA)
    assert Ms > A > 0
    assert _melt_monthly(TEMPERATE_T_M, S2, R, SW_SUMMER, p) == pytest.approx(A + (1 - A / Ms) * Mi, rel=1e-12)


def test_vectorised():
    """Many climates in one call give the same results as one call per climate."""
    rng = np.random.default_rng(5)
    T = rng.normal(-5, 8, (40, 12)); P = rng.uniform(0.1, 3, (40, 12)); SW = rng.uniform(0, 350, (40, 12))
    for p in (sd.Params(sigma_w=SIGMA), sd.degree_day_params(sd.Params(sigma_w=SIGMA))):
        sw = SW if p.beta_SW_snow else None
        V = sd.surface_forcing(T, P, sw, p)
        for i in range(len(T)):
            one = sd.surface_forcing(T[i], P[i], None if sw is None else SW[i], p)
            for k in ("S", "R", "M", "T_sfc", "f_c", "T_mu"):
                assert V[k][i] == pytest.approx(one[k], rel=1e-12, abs=1e-15)
    assert isinstance(_melt_monthly(T[0], 1.0, 0.1, SW[0], sd.Params(sigma_w=SIGMA)), float)


def test_melt_daily():
    rng = np.random.default_rng(1)
    T = np.repeat(TEMPERATE_T_M, 30) + rng.normal(0, SIGMA, 360)
    SW = np.repeat(SW_SUMMER, 30)
    p = sd.Params(sigma_w=SIGMA)
    for S, R in [(3.0, 0.1), (0.3, 0.05)]:
        Ms = np.sum(np.maximum(p.beta_T_snow * T + p.beta_SW_snow * SW + p.beta_0_snow, 0))
        Mi = np.sum(np.maximum(p.beta_T_ice * T, 0))
        A = p.p_max * S - R
        exp = Ms if Ms <= A else A + (1 - A / Ms) * Mi
        assert sd.melt_daily(T, S, R, SW, p) == pytest.approx(exp, rel=1e-12)
    q = sd.degree_day_params(sd.Params(sigma_w=SIGMA))
    pdd = np.sum(np.maximum(T, 0))
    assert sd.melt_daily(T, 0.3, 0.05, p=q) == pytest.approx(degree_day_melt(pdd, 0.3, 0.05, 0.004, 0.008, q.p_max))
