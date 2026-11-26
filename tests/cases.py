"""Reference climates of the CFM benchmark pair (v1.0 paper): precipitation 1 m w.e./yr, the sweep's seasonal
cycle at +250 m (cold firn) and 0 m (temperate firn), daily weather sigma 3.834 K."""
import dataclasses

import numpy as np

import snowdrip as sd

SIGMA = 3.834057903
P_M = np.full(12, 1.0)
COLD_T_M = np.array([-18.520489, -16.476914, -12.925528, -8.67992, -4.96013, -2.756922,
                     -2.680714, -4.784882, -8.441309, -12.68806, -16.382494, -18.520489])
TEMPERATE_T_M = COLD_T_M + 1.625          # 250 m lower at 6.5 K/km
T_M = {"cold": COLD_T_M, "temperate": TEMPERATE_T_M}
SW_SUMMER = np.array([20, 60, 130, 210, 280, 310, 290, 220, 150, 80, 30, 10], float)    # W m-2


def model(case, **params):
    p = dataclasses.replace(sd.degree_day_params(), sigma_w=SIGMA, **params)
    f = sd.surface_forcing(T_M[case], P_M, p=p)
    return f, sd.firn(f["S"], f["M"], f["T_sfc"], f["f_c"], f["R"], p)
