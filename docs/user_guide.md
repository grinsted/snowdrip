# SnowDRIP user guide

SnowDRIP computes the time-mean surface mass balance and near-surface firn of a glacier at one elevation,
from a climatological year: melt, refreezing, runoff, net balance, the 10 m firn temperature and the density
of the firn leaving the seasonal active layer. It is steady state: it describes the mean state the firn would
have under a given climate, not its evolution through time. The derivations and the benchmark against the
Community Firn Model are in the SnowDRIP v1.0 paper (in preparation for Geoscientific Model Development).

```python
import snowdrip as sd
```

## Running the model

Give 12 monthly means at the elevation of interest:

```python
r = sd.run(T_m, P_m, SW_m)
```

| input | unit | |
|---|---|---|
| `T_m` | °C | monthly mean air temperature |
| `P_m` | m w.e./yr | monthly precipitation rate, already corrected to the glacier (see Forcing) |
| `SW_m` | W m⁻² | monthly mean downward shortwave radiation |

`run()` returns one dict with everything below. It does two steps, which can also be called separately:

```python
s = sd.surface_forcing(T_m, P_m, SW_m)                        # the climatological year at the surface
r = sd.firn(s["S"], s["M"], s["T_sfc"], s["f_c"], s["R"])     # the firn
```

`surface_forcing()` integrates over Gaussian daily weather (spread `sigma_w`) in each month and returns:

| output | unit | |
|---|---|---|
| `S`, `R` | m w.e./yr | snowfall and rain (rain above `p.rain_threshold`, 1 °C) |
| `M` | m w.e./yr | gross melt |
| `T_sfc` | °C | time-mean surface temperature, mean(min(T, 0)) |
| `f_c` | – | fraction of the year below 0 °C |
| `T_mu` | °C | annual mean air temperature |

Its inputs may carry leading dimensions, `(..., 12)`, to evaluate many climates in one call (e.g. every
elevation band of a glacier); the outputs then have shape `(...)`. `run()` and `firn()` take one climate at a
time.

With your own melt model or measured melt, call `firn(S, M, T_sfc, f_c, rain)` directly with annual totals.
For daily series (e.g. a weather station), `sd.melt_daily(T_d, S, R, SW_d, n_years=...)` gives the melt with
the same melt model.

**Parameters.** Every function takes the parameters as `p` (default `sd.Params()`) and as keyword overrides:

```python
r = sd.run(T_m, P_m, SW_m, tau_w=1.0, p_max=0.55)
p = sd.Params(tau_w=1.0, p_max=0.55)                          # the same, reusable
r = sd.run(T_m, P_m, SW_m, p)
```

When calling `surface_forcing()` and `firn()` separately, give both the same parameters: `p_max` is used by
both steps.

## Forcing

The melt coefficients in `Params` are calibrated with a particular forcing and only hold with it:

- **Temperature** moved from the reanalysis surface to the elevation of interest with a lapse rate of
  6.5 K/km in all months. The FoG mass-balance profiles constrain the product of the ice-melt coefficient
  and the lapse rate, so a different lapse rate needs a different `beta_T_ice`.
- **Precipitation** from ERA5 multiplied by 1.70, constant with elevation. Gridded products miss the snow
  that collects on glaciers.
- **Shortwave radiation** moved to the elevation of interest with the overlying air mass and cloud water.

The model itself knows none of this: it takes temperature, precipitation and shortwave at the glacier.
The lapse rate and precipitation factor are `LAPSE_RATE` and `PRECIP_FACTOR` in `examples/era5_forcing.py`.

`examples/era5_forcing.py` builds this forcing from the ERA5 climatologies
([doi:10.5281/zenodo.24079220](https://doi.org/10.5281/zenodo.24079220)); `examples/fog_glacier.ipynb`
uses it for one glacier. Other sources (weather stations, regional climate models, other reanalyses) work as
long as they give monthly temperature, precipitation and shortwave at the glacier; expect a different
precipitation correction, and check melt against observations before trusting absolute values.

Without shortwave radiation, use the degree-day form, `sd.run(T_m, P_m, p=sd.degree_day_params())` (4 and
8 mm w.e./K/day for snow and ice).

## Outputs of `firn()`

Also in the dict from `run()`. Fluxes in m w.e./yr, depths in m, temperatures in °C, densities in kg m⁻³.

| key | |
|---|---|
| `W` | retained (refrozen) meltwater and rain |
| `runoff` | M + R − W |
| `bdot` | net balance S − M + W |
| `ablation` | max(−bdot, 0) |
| `firn` | `True` above the equilibrium line (bdot > 0) |
| `T10m` | firn temperature at 10 m |
| `T_deep` | temperature below the active layer |
| `rho_deep` | density of the firn leaving the active layer |
| `z`, `T`, `rho` | depth grid (default 0–20 m) and the temperature and density profiles on it |
| `ell_d` | depth of the active layer (where the year's water refreezes; e-folding depth of `rho`) |
| `ell_c`, `ell_s` | depths reached by the winter cold and by a cold spell |
| `H_snow` | thickness of the year's snowfall |
| `temperate`, `w_T` | firn regime: `w_T` is the weight of the temperate state (0 cold, 1 temperate); `temperate` is `w_T >= 0.5` |
| `E_C`, `E_T` | water excess of the cold and the temperate state, which set `w_T` |
| `phi`, `deep_share`, `xi` | open porosity of the firn below, share of the water refrozen deep, refreezing per unit cold supply |

Below the equilibrium line (`firn` is `False`) the year's snow lies on ice: retention follows the seasonal
snow, `rho` and `rho_deep` are ice, and `T`, `T10m`, `T_deep` and the regime diagnostics are NaN.

The density profile does not include dry densification. To continue it below the active layer, use
`rho_deep`, `ell_d`, `T_deep` and `bdot` as the upper boundary of a steady-state densification model, e.g.
d rho/dz = max[(rho_deep − rho)/ell_d, Herron–Langway rate].

## Parameters

| parameter | default | change it? |
|---|---|---|
| `beta_T_snow`, `beta_SW_snow`, `beta_0_snow` | 2.85 mm/K/day, 0.0546 mm/day per W m⁻², −5.34 mm/day | no: calibrated on PROMICE / GC-Net melt |
| `beta_T_ice`, `beta_SW_ice` | 6.0 mm/K/day, 0 | only with different forcing (calibrated on FoG with the ERA5 forcing above) |
| `k` | 0.34 W/m/K | no: calibrated on the Community Firn Model |
| `sigma_w` | 3.05 K | rarely: daily spread of temperature about the monthly mean, fitted to the degree days at PROMICE / GC-Net stations; the snow coefficients belong to it |
| `tau_w` | 7 days | yes, to match the forcing: decorrelation time of the weather; 7 d for real weather, 1 d for white-noise daily forcing |
| `rain_threshold` | 1 °C | rarely |
| `p_max` | 0.62 | no: pore space of fresh snow, 1 − ρ_snow/ρ_ice |
| `f_pf` | 0 | optional: share of the water that bypasses the winter cold as preferential flow; 0 in v1.0 |

`tau_w` and `f_pf` move latent heat only: they change temperatures, not retention, runoff or net balance.

## Recipes

**A mass-balance profile.** Loop over elevations, moving the climatological year to each one:

```python
for z in elevations:
    T = T_ref - 6.5 * (z - z_ref) / 1000
    r = sd.run(T, P_m, SW_at(z))
```

**The equilibrium-line altitude.** `sd.ela_elevation(T_m, z_ref, lapse_rate, P_m, SW_m)`
returns the elevation where bdot = 0, with `T_m` given at `z_ref` and the lapse rate in K/km (positive).
`SW_m` may be 12 values or a function of elevation returning them. `sd.ela_offset(T_m, P_m, SW_m)` instead
returns the uniform temperature shift that puts a site at its equilibrium line. `sd.net_balance(T_m, P_m, SW_m)`
is bdot alone. All three take parameters like `run()`.

## How it works

In short; the derivations and the benchmark against CFM are in [model_internals.md](model_internals.md).

- **Melt.** Daily melt is max(0, β_T T + β_SW SW + β_0): degree days above a threshold that sunshine lowers
  below 0 °C (Giesen and Oerlemans, 2012). Snow melts until melt and rain fill the pores of the year's snow;
  after that the surface melts as ice.
- **Retention.** At most the pore space of the year's snowfall, p_max S, can hold meltwater; the rest runs
  off.
- **Refreezing and temperature.** Water refreezes in two reservoirs: a shallow one, cooled by the next cold
  spell, and a deep one at about half the snowpack thickness, reached only by the winter cold. Its latent heat
  warms the firn in a steady balance with conduction and burial.
- **Cold and temperate firn.** When there is more water than the winter cold can freeze, the firn becomes
  temperate (0 °C below the active layer) and retention is limited by the cold supply. Near the transition
  the two states are blended (`w_T`).
- **Density.** The firn leaving the active layer is denser by the refrozen share of the surviving mass.

## Limitations

- Steady state only; interannual variability enters only through the daily weather spread.
- The near-surface seasonal structure (the warm summer layer) is not represented; the model targets the
  firn below a few metres.
- Near the cold–temperate transition the result is a blend of two states, and the real firn may be in either.
- In strongly seasonal climates retention can drop by up to about 0.5 m w.e./yr across the equilibrium line,
  where the cold firn above gives way to snow on ice.
- In the wet percolation zone the default (no preferential flow) tends to be too cold; `f_pf` is the lever.

## Reference

The SnowDRIP v1.0 paper, in preparation for Geoscientific Model Development. See CITATION.cff for how to
cite the model and the code version.
