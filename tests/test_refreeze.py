"""Regression tests and invariants of the firn model. Reference values were produced by the version used
for the paper; a change means the model changed."""
import numpy as np
import pytest

import snowdrip as sd
from cases import T_M, model

REFERENCE = {
    # (case, tau_w): expected values
    ("cold", 1.0): dict(S=0.9605229460133432, R=0.039477053986656824, M=0.1819662132292444,
                        T_sfc=-10.749549009739388, f_c=0.9405696064142766, W=0.2214432672159012,
                        runoff=0.0, bdot=1.0, T10m=-6.150568230293637, rho_deep=475.558332511416, w_T=0.0),
    ("cold", 7.0): dict(W=0.2214432672159012, runoff=0.0, bdot=1.0, T10m=-5.504529410629165,
                        rho_deep=475.558332511416, w_T=0.0),
    ("temperate", 1.0): dict(S=0.9250909480949735, R=0.07490905190502647, M=0.3735085859769446,
                             T_sfc=-9.255652618738498, f_c=0.8954233241579068, W=0.32572797933578757,
                             runoff=0.12268965854618352, bdot=0.8773103414538165,
                             T10m=-0.004788255541186961, rho_deep=560.5158865189485, w_T=1.0),
}


@pytest.mark.parametrize("key", list(REFERENCE))
def test_reference_cases(key):
    case, tau = key
    f, r = model(case, tau_w=tau)
    for name, expected in REFERENCE[key].items():
        got = f[name] if name in f else r[name]
        assert got == pytest.approx(expected, rel=1e-9, abs=1e-12), name


@pytest.mark.parametrize("case", ["cold", "temperate"])
def test_mass_conservation(case):
    f, r = model(case)
    assert r["W"] + r["runoff"] == pytest.approx(f["M"] + f["R"])
    assert r["bdot"] == pytest.approx(f["S"] - f["M"] + r["W"])
    assert 0.0 <= r["W"] <= sd.Params().p_max * f["S"] + 1e-12


@pytest.mark.parametrize("case", ["cold", "temperate"])
@pytest.mark.parametrize("params", [dict(tau_w=1.0), dict(tau_w=20.0), dict(f_pf=0.5)])
def test_heat_only_parameters(case, params):
    """tau_w and f_pf move latent heat only: retention, runoff and bdot are unchanged."""
    _, r0 = model(case)
    _, r1 = model(case, **params)
    for name in ("W", "runoff", "bdot"):
        assert r1[name] == pytest.approx(r0[name], rel=1e-12, abs=1e-15)


def test_density_profile():
    for case in T_M:
        _, r = model(case)
        rho, z = r["rho"], r["z"]
        assert rho[0] == pytest.approx(sd.RHO_SNOW)
        assert np.all(np.diff(rho) >= 0) and np.all(rho <= r["rho_deep"] + 1e-9)
        assert np.interp(r["ell_d"], z, rho) == pytest.approx(r["rho_deep"] - (r["rho_deep"] - sd.RHO_SNOW) / np.e, rel=1e-3)


def test_defaults():
    p = sd.Params()
    assert not hasattr(p, "gamma")
    assert p.tau_w == 7.0 and p.f_pf == 0.0
    assert p.k == 0.34 and p.rain_threshold == 1.0 and p.beta_SW_ice == 0.0
    assert p.beta_SW_snow > 0.0 and 0.0 < p.beta_T_snow < p.beta_T_ice
    q = sd.degree_day_params()
    assert (q.beta_T_snow, q.beta_SW_snow, q.beta_T_ice, q.beta_SW_ice) == (0.004, 0.0, 0.008, 0.0)


def test_temperature_bounds():
    for case in T_M:
        f, r = model(case)
        assert np.all(r["T"] <= 0.0)
        assert r["T10m"] >= f["T_sfc"] - 1e-12           # latent heat only warms


def test_below_ela():
    f, r = model("temperate")
    p = sd.degree_day_params()
    r = sd.firn(f["S"], 5.0 * f["M"], f["T_sfc"], f["f_c"], f["R"], p)   # much more melt than snow
    assert not r["firn"] and r["bdot"] < 0
    assert np.isnan(r["T10m"]) and r["rho_deep"] == sd.RHO_I
    for key in ("w_T", "E_C", "E_T", "phi", "deep_share"):
        assert np.isnan(r[key]), key
    assert r["temperate"] is False
    assert r["ell_d"] > 0 and r["H_snow"] > 0 and np.isfinite(r["ell_c"])
    assert r["W"] <= p.p_max * f["S"] + 1e-12
