# SnowDRIP

SnowDRIP is a steady-state temperature-index model with firn refreezing for glaciers. From a
climatological year it computes, for one elevation at a time, the melt, the refrozen meltwater, the
runoff, the net balance, the 10 m firn temperature and the density of the firn leaving the seasonal
active layer. Snow melt is the simplified energy balance of Giesen and Oerlemans (2012) over Gaussian
daily weather; refreezing follows from a mass and heat budget built to reproduce the time-mean state of
the Community Firn Model.

See [docs/user_guide.md](docs/user_guide.md) for inputs, outputs, forcing and parameters. The model is
described in the SnowDRIP v1.0 paper (in preparation for Geoscientific Model Development).

The name started as *Snow Degree-day Refreezing In Pores*.

## Install

    pip install .                 # numpy and scipy
    pip install ".[examples]"     # + xarray, netCDF4, matplotlib, jupyter for examples/

## Use

```python
import numpy as np
import snowdrip as sd

T_m = np.array([-16.9, -14.9, -11.3, -7.1, -3.3, -1.1, -1.1, -3.2, -6.8, -11.1, -14.8, -16.9])  # monthly means [C]
P_m = np.full(12, 1.0)                                                       # precipitation [m w.e./yr]
SW_m = np.array([20, 60, 130, 210, 280, 310, 290, 220, 150, 80, 30, 10.0])   # shortwave [W m-2]

r = sd.run(T_m, P_m, SW_m)
print(r["M"], r["W"], r["runoff"], r["bdot"], r["T10m"], r["rho_deep"])

r = sd.run(T_m, P_m, SW_m, tau_w=1.0)                                        # parameters as keywords

z_ela = sd.ela_elevation(T_m, z_ref=2000.0, lapse_rate=6.5, P_m=P_m, SW_m=SW_m)   # T_m given at 2000 m
r = sd.run(T_m, P_m, p=sd.degree_day_params())                               # degree-day melt, no shortwave
```

The melt coefficients are calibrated with ERA5 forcing prepared in a particular way (lapse rate,
precipitation factor, shortwave at elevation); see the user guide before using other forcing.

## Layout

| | |
|---|---|
| `src/snowdrip/model.py` | the full model (`run`), net balance and equilibrium-line altitude |
| `src/snowdrip/melt.py` | the climatological year at the surface, with temperature-index melt (`surface_forcing`) |
| `src/snowdrip/refreeze.py` | the firn model: refreezing, runoff, temperature and density (`firn`) |
| `src/snowdrip/params.py` | constants, forcing conventions and `Params` |
| `docs/user_guide.md` | how to use the model (start here) |
| `docs/model_internals.md` | derivations and the CFM benchmark |
| `examples/` | the ERA5 forcing (`era5_forcing.py`) and a glacier profile (`fog_glacier.ipynb`) |

## Tests

    pytest tests

## Citation

See CITATION.cff. Please cite the model paper and the archived version of the code you used.

## License

MIT (see LICENSE).
