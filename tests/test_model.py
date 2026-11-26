import dataclasses

import numpy as np
import pytest

import snowdrip as sd
from cases import P_M, SIGMA, SW_SUMMER, T_M


@pytest.mark.parametrize("case", ["cold", "temperate"])
def test_ela_functions(case):
    p = sd.degree_day_params(sd.Params(sigma_w=SIGMA))
    gamma = 6.5                                         # K/km
    z_ela = sd.ela_elevation(T_M[case], 0.0, gamma, P_M, p=p)
    T_ela = T_M[case] - gamma * z_ela / 1000.0
    assert abs(sd.net_balance(T_ela, P_M, p=p)) < 2e-3
    dT = sd.ela_offset(T_M[case], P_M, p=p)
    assert dT == pytest.approx(-gamma * z_ela / 1000.0, abs=0.02)


@pytest.mark.parametrize("case", ["cold", "temperate"])
def test_ela_functions_shortwave(case):
    p = sd.Params(sigma_w=SIGMA)
    gamma = 6.5
    z_ela = sd.ela_elevation(T_M[case], 0.0, gamma, P_M, SW_SUMMER, p)
    T_ela = T_M[case] - gamma * z_ela / 1000.0
    assert abs(sd.net_balance(T_ela, P_M, SW_SUMMER, p)) < 2e-3
    dT = sd.ela_offset(T_M[case], P_M, SW_SUMMER, p)
    assert dT == pytest.approx(-gamma * z_ela / 1000.0, abs=0.02)
    sw = lambda z: SW_SUMMER * (1 + z / 1e3)                     # shortwave changing with elevation
    z2 = sd.ela_elevation(T_M[case], 0.0, gamma, P_M, sw, p)
    assert z2 != pytest.approx(z_ela, abs=1.0)
    assert abs(sd.net_balance(T_M[case] - gamma * z2 / 1000.0, P_M, sw(z2), p)) < 2e-3


def test_run_is_the_two_steps():
    p = sd.Params(sigma_w=SIGMA)
    r = sd.run(T_M["temperate"], P_M, SW_SUMMER, p)
    f = sd.surface_forcing(T_M["temperate"], P_M, SW_SUMMER, p)
    o = sd.firn(f["S"], f["M"], f["T_sfc"], f["f_c"], f["R"], p)
    for key, value in (f | o).items():
        assert np.array_equal(r[key], value, equal_nan=True) if np.ndim(value) else r[key] == value or np.isnan(value), key


def test_keyword_overrides():
    T = T_M["temperate"]
    a = sd.run(T, P_M, SW_SUMMER, p_max=0.55, tau_w=1.0)
    b = sd.run(T, P_M, SW_SUMMER, dataclasses.replace(sd.Params(), p_max=0.55, tau_w=1.0))
    assert a["W"] == b["W"] and a["T10m"] == b["T10m"]
    assert a["W"] != sd.run(T, P_M, SW_SUMMER)["W"]
    with pytest.raises(TypeError):
        sd.run(T, P_M, SW_SUMMER, p_maks=0.55)
    with pytest.raises(dataclasses.FrozenInstanceError):
        sd.Params().p_max = 0.55
