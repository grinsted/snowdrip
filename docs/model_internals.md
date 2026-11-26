# SnowDRIP model internals

Derivations and CFM benchmark behind the model, for readers who want the details; to use the model,
start with the [user guide](user_guide.md). SnowDRIP is a minimal steady-state model of the time-mean firn temperature
(T10m), density and meltwater retention in the accumulation zone of mountain glaciers, to be
applied globally. Benchmark: CFM with off-the-shelf physics (snowdrip-analysis: `lib/cfm_sweep.py`, `checks/cfm_benchmark.py`).

## Inputs and notation

Fluxes in m w.e./yr, depths in m of snow.

| symbol | meaning |
|---|---|
| S, R | snowfall and rain; P = S + R total precipitation |
| M | gross melt (temperature-index melt, below) |
| T_sfc | time-mean surface temperature mean(min(T, 0)) |
| f_c | cold-season fraction Prob(T < 0) (from the climatological year) |
| W | retained (refrozen) water; corresponds to RF in RACMO/MAR budgets |
| Q | runoff = M + R - W (RU) |
| bdot | net accumulation S - M + W |
| G(ell) | conductive cold supply across depth ell, as an equivalent accumulation |
| H_snow | snowpack thickness S rho_w / rho_snow |
| ell_c, ell_s, ell_d | winter cold-wave depth, weather-scale depth, deep-source depth |
| xi | c abs(T_sfc) / L, the Stefan number of the snow: refreezing per unit cold supply |

### Melt

Daily melt on surface s (snow, ice) is max(0, beta_T,s T + beta_SW,s SW + beta_0,s), for snow the simplified
energy balance of Giesen and Oerlemans (2012): degree days above the threshold
T*_s = -(beta_SW,s SW + beta_0,s) / beta_T,s, which the shortwave radiation lowers below 0 C, so that
snow melts on sunny days with sub-freezing air (half of the melt at KAN_U falls on such days). Over
the Gaussian daily weather of each month the melt potential is M_s = beta_T,s PDD(T*_s)
(`degree_days`, `surface_forcing`). The surface melts as snow until melt and rain fill the pores of
the year's snow, A = max(0, P_max S - R) of melt, and as ice after that:

    M = min(M_snow, A) + max(0, 1 - A / M_snow) M_ice

beta_SW,ice = beta_0,ice = 0 in v1.0 (ice melts by degree days: the FoG profiles do not support a shortwave
term for ice). With beta_SW = beta_0 = 0 for both surfaces this is the sequential degree-day scheme of Reeh (1991), switching
to ice when the snow's pores are saturated rather than when the snow is gone (`degree_day_params`,
DDF 4 and 8 mm/K/day; the daily form `cfm_sweep.melt_degree_day_daily` drives the CFM benchmark).

The shortwave radiation at elevation z (`shortwave_at_elevation`, examples/era5_forcing.py) comes from the reanalysis values at
its surface z_ref: SW(z) = SW_toa tau_ref^exp(-(z - z_ref)/8 km) [1 - a_ref F(z)], with
tau_ref = SW_clear/SW_toa, a_ref = 1 - SW/SW_clear, and F(z) the share of the column's cloud water
above z relative to above z_ref (`cloud_water_share`). `ERA5Forcing` (same file) assembles the default forcing, with
temperatures moved by LAPSE_RATE = 6.5 K/km (examples/era5_forcing.py) in all months (lapse_rate="era5" for the ERA5 monthly free-air field).
Coefficients (Params): beta_SW,snow = 0.0546 mm/day per W m-2, all the absorbed shortwave at the median albedo of
the PROMICE snow months (0.79); beta_T,snow = 2.85 mm/K/day and beta_0,snow = -5.34 mm/day from the PROMICE / GC-Net
energy-balance melt, forced by the stations' own monthly air temperature and downward shortwave (independent of the
ERA5 downscaling); beta_T,ice = 6.0 mm/K/day from the FoG mass-balance series 1951-2025 with the default ERA5 forcing
and the snow coefficients fixed (snowdrip-analysis: checks/melt_calibration.py). The series constrain
beta_T,ice x LAPSE_RATE (~39 mm w.e./d/km) rather than beta_T,ice alone.

## Steady heat budget

Annual-mean heat equation for the firn column, burial velocity w = rho_w bdot / rho, surface
held at T_sfc, no heat flux at depth, and an exponential latent-heat source
s(z) ~ W_i exp(-z/ell) from water refreezing around depth ell. The exact solution is

    T(z) - T_sfc = dT (1 - exp(-z/ell)),    dT = W_i (L/c + abs(T_sfc)) / (bdot + G(ell))

with

    G(ell) = k / ell * yr / (rho_w c)    [m w.e./yr]

Conduction is expressed as the accumulation rate whose advected cold content equals the
conductive heat loss to the surface across ell. bdot / G is a Peclet number; with ell of a few
metres G is a few m w.e./yr, so conduction dominates except at high accumulation. Sources at
different depths add.

## Depth scales and conductivity

* **Winter cold-wave depth** ell_c = sqrt(kappa_snow f_c yr / pi): the depth winter cold reaches
  into snow during the cold season. About 1.9 m in the CFM sweep (about 1.5 m in a maritime
  climate).
* **Weather-scale depth** ell_s = sqrt(kappa_snow tau_w / pi): shallow meltwater refreezes in the
  next cold spell; tau_w is the decorrelation time of the weather, a property of the forcing:
  about 7 days for real weather (the default; PROMICE summer daily anomalies about a smooth
  daily climatology, the same definition as sigma_w = 2.2 K), 1 day for the white-noise daily
  forcing of the CFM benchmark. About 0.27 m and 0.10 m. Heat only: retention and runoff do not
  depend on it.
* **Deep-source depth** ell_d = H_snow / 2: water that gets below the surface layer is spread
  through the annual snowpack.
* **h = H_snow / ell_c** compares the snowpack with the reach of the winter cold.

Both cold waves travel through fresh snow, so ell_c and ell_s use the conductivity of snow at
rho_snow = 350 kg m-3, k_snow = 0.29 W/m/K (Calonne et al., 2019, as in CFM), and
kappa_snow = k_snow / (rho_snow c). Buried heat at ell_d escapes through the older, denser column:
G(ell_d) uses k = 0.34 W/m/K, calibrated on CFM. It is close to the winter (Nov-Apr) series
(harmonic-mean) Calonne conductivity of the top ell_c of CFM's own firn (0.33; 0.43 in summer): the cold
enters in winter, through the winter snowpack. A conductivity from the density rho_deep of the firn
leaving the active layer, k(rho_deep), would degrade the CFM benchmark (T10m rms 0.58 -> 4.4 K).

## Mass

    W_C = min(M + R, P_max S),   P_max = 1 - rho_snow / rho_ice = 0.62
    bdot = S - M + W,   Q = M + R - W

Retention is limited by water, or by the pore space of the year's snowfall: if melt removes snow
from the top while its water refreezes in the pores below, the remaining snow is solid ice once
a fraction P_max of the snowfall has melted (Reeh's P_max). The limit uses the whole winter
snowpack S, not the surviving snow S - M: gross melt includes refrozen ice that melts again in
the same summer.

## Cold regime

A share s of the retained water is held beyond the reach of the winter cold and refreezes at
ell_d; the rest refreezes at ell_s, where its latent heat is mostly lost to the surface
(G(ell_s) with k_snow is about 16 m w.e./yr at tau_w = 7 d, 42 at 1 d). The cold wave decays as exp(-z/ell_c), so averaged
over the annual layer the share beyond its reach is 1 - (1 - exp(-h))/h. The part of it beyond
the hard cut-off max(0, 1 - 1/h) has to get into older firn, which is only possible through its
pores:

    s(h, phi) = (1 - phi) s_closed + phi s_open,    s_closed = max(0, 1 - 1/h),    s_open = 1 - (1 - exp(-h))/h
    phi = clip((rho_c - rho_deep_c) / (rho_c - rho_snow), 0, 1),   rho_c = 830 kg m-3

with rho_deep_c the firn density of the cold branch (below) and rho_c pore close-off (also CFM's
impermeable density). Rain falls on a mostly cold snowpack and percolates past the warm surface
layer: if h > 1 all retained rain refreezes at ell_d, otherwise it is split like melt.

Optional preferential flow (`Params.f_pf`, default 0 = v1.0, the CFM emulation): a fraction f_pf
of the remaining water bypasses the cold wave, s' = s + f_pf (1 - s). This moves latent heat
only (retention, runoff and bdot are unchanged). CFM's bucket scheme has no preferential flow,
and in Greenland's wet percolation zone v1.0 is too cold (see the paper notes); f_pf = 0.5
removes most of that bias at the station-driven PROMICE sites at the cost of the CFM agreement.

    dT_s = (1 - s) W L' / (bdot + G_s(ell_s)),   dT_d = s W L' / (bdot + G(ell_d)),   L' = L/c + abs(T_sfc)
    T_deep = min(T_sfc + dT_s + dT_d, 0)
    T(z) = T_deep - dT_s exp(-z/ell_s) - dT_d exp(-z/ell_d)

## Temperate regime

Wet firn at 0 C below the deep reservoir (depth ~ell_d) lets water through and stays temperate. Its time-mean
temperature is T_sfc exp(-z/ell_d). Retention is what offsets the cold supply of that profile,
plus the capillary (irreducible, Coleou and Lesaffre 1998) water carried down with the firn;
the rest drains:

    W_T = min(M + R, P_max S, xi (bdot_T + G(ell_d)) + capillary export)

solved by fixed-point iteration (bdot_T depends on W_T). Just below the ELA bdot_T can be negative
while the cold branch still accumulates; no firn is then carried down, so the burial term and the
capillary export use max(bdot_T, 0), and so does the water excess E_T below. In CFM the water held in temperate firn
is capillary (3-9 % of pore space), not a saturated aquifer: with ponding off, water above the
irreducible content that reaches an impermeable layer (>= 830 kg m-3, >= 0.4 m) runs off
immediately, so no slope term is needed.

## From cold to temperate

Only water that reaches the firn within ell_c (the reach of the winter cold) counts against its cold supply:

    f_a = [f_new + (1 - f_new) phi] * W_C / (M + R),   f_new = min(1, h)

This year's snow is accessible; older firn within ell_c in proportion to its open porosity; and
water that runs off over the saturated surface never gets there (which makes the wettest,
lowest sites warm later than slightly higher ones, as in CFM). The relative water excess of each state:

    E_C = f_a W_C / (xi (bdot_c + G(ell_c))) - 1     > 0: the cold state cannot hold
    E_T = f_a W_C / (xi (bdot_T + G(ell_d))) - 1     > 0: the temperate state sustains itself

Both positive: temperate; both negative: cold. Where they disagree the model has two states
(E_T > 0 > E_C: bistable, thick snowpacks, 0.1-0.7 K wide in T_mu, less than the interannual
variability) or none (E_C > 0 > E_T: the gap, P about 0.5-1.2 in the sweep climate, up to 2 K wide).
The result is a linear blend of the two branches,

    w_T = max(E_C, E_T) / (abs(E_C) + abs(E_T)),

0 on the edge beyond which only the cold state exists and 1 on the edge beyond which only the
temperate state exists. Temperature profile, retention, runoff and density are all blended with
w_T. No parameter is involved.

## Density

Firn leaving the active zone has

    rho_deep = rho_snow + (rho_ice - rho_snow) W / bdot,

linear in the refrozen share of the surviving mass. The density approaches rho_deep exponentially,

    rho(z) = rho_deep - (rho_deep - rho_snow) exp(-z / ell_d),   i.e.   d rho/dz = (rho_deep - rho) / ell_d,

with the depth of the deep source, ell_d, as e-folding depth: the active layer is where the
year's water ends up, and the densification rate is proportional to the remaining density deficit.
Against the CFM time-mean profiles (0 to ell_d + 1 m, 30 runs) this fits better than a linear ramp
to rho_deep at ell_d (rms 40 vs 57 kg m-3; with CFM's own rho_deep, i.e. shape only, 48 vs 59). Faster
approaches (e-folding 0.1-0.3 ell_d, or ell_s) and a two-stage step at ell_s fit worse (59-95).
No new parameter. (Earlier versions: a linear ramp to rho_deep at gamma H_snow, gamma = 0.6, then at
ell_d.) T10m, rho_deep and runoff do not depend on the profile shape.

Dry densification is not part of SnowDRIP: rho_deep, ell_d, T_deep and bdot are the upper
boundary conditions for a steady-state densification model. A natural coupling is the single
equation d rho/dz = max[(rho_deep - rho)/ell_d, Herron-Langway rate at T and bdot]: melt-driven
densification where it is faster, Herron-Langway in dry firn and as rho approaches rho_deep, with no
switching depth (in the mass study: `density_model.py`). `firn()` returns the profile without dry
densification.

## Below the ELA

Where the firn model gives bdot <= 0 no firn survives: seasonal snow lies on ice. `firn()` then
returns `firn = False` and

    W = min(M + R, P_max S, xi G(ell_d)),    bdot = S - M + W

the cold supply without burial, as the temperate branch at bdot = 0. bdot is monotonic in
temperature and continuous at the ELA where the firn just above it is temperate. Where it is cold, in
strongly seasonal climates, retention drops across the ELA by up to ~0.5 m w.e./yr (seasonal amplitude
18 K, P 1.5; continuous in most climates): above the ELA the cold reservoir is old firn, partly
shielded from runoff over ice lenses (the factor W/(M + R) in f_a), below it the ice surface that all
water touches. Letting the firn regimes run below the ELA removes the step but makes retention jump
back to the full pore space a few K below the ELA, because f_a shrinks down-glacier.

Density below the ELA is ice (`ablation_column_density`, the place for a future near-surface
low-density layer such as observed by Chandra et al.). Temperature is not modelled there: T, T10m and
T_deep are NaN, as are the firn-only diagnostics (w_T, E_C, E_T, phi, deep_share). Also returned for later use: H_snow, xi = c|T_sfc|/L and the ablation rate -bdot.

## From climate to ELA

`surface_forcing(T_m, P_m, SW_m)` gives the model inputs at one elevation from a
climatological year: 12 monthly mean temperatures, each with Gaussian daily weather noise sigma_w
(sigma_w = 3.05 K: the spread about the monthly mean that reproduces the degree days of the daily temperatures at 46 PROMICE / GC-Net stations; the raw summer weather residual is ~2.2 K), and 12 monthly precipitation rates. PDDs, T_sfc,
f_c and the rain fraction (above `rain_threshold`) are sums of monthly Gaussians, then the temperature-index melt. A
single Gaussian around the annual mean is not used: it fits the FoG profiles below the ELA worse (rms 0.61
vs 0.58 m w.e./yr). `net_balance` is bdot of the full model for such a year (where the cold branch
has no surviving firn, the cold branch's negative bdot; the sign is consistent but the value is not
a profile below the ELA).

`ela_offset(T_m, P_m)` is the uniform temperature shift of the year that gives bdot = 0;
`ela_elevation(T_m, z_ref, lapse_rate, P_m)` the elevation of bdot = 0 with monthly (or one) lapse
rates and precipitation constant with height. The cold branch alone retains at least as much as the
full model, so its root (cheap, no firn model) is an upper bound; the full-model root is bracketed
below it and found with Brent's method (about 1.5 ms). In rare climates bdot crosses zero more
than once in a narrow band (bistable regimes; one of 80 pre-2000 FoG glaciers, three crossings
within 60 m); one crossing is returned.

## Parameters

| parameter | value | role |
|---|---|---|
| beta_T,snow, beta_SW,snow, beta_0,snow | 2.85 mm/K/day, 0.0546 mm/day per W m-2, -5.34 mm/day | snow melt; PROMICE / GC-Net energy-balance melt |
| beta_T,ice, beta_SW,ice | 6.0 mm/K/day, 0 | ice melt; FoG mass-balance series, with LAPSE_RATE |
| sigma_w | 3.05 K | daily spread about the monthly mean (PROMICE / GC-Net degree days); the snow coefficients belong to it |
| LAPSE_RATE | 6.5 K/km | temperature lapse rate of the ERA5 forcing (examples/era5_forcing.py), all months; not a model parameter |
| DDF_snow, DDF_ice | 4, 8 mm w.e./K/day | degree-day form (`degree_day_params`), used in the CFM benchmark |
| rain threshold | 1 C | rain/snow partition (USACE 1956) |
| PRECIP_FACTOR | 1.70 | input correction in examples/era5_forcing.py (not a model parameter): global multiplier on gridded precipitation: glaciers sit where snow collects (wind drift, avalanches, sub-grid orography); median needed for the FoG winter balances. Precipitation is constant with elevation: a slope-dependent orographic gradient is a possible extension, left out to keep the model easy to apply |
| P_max | 1 - rho_snow/rho_ice = 0.62 | pore-space limit on retention (not fitted) |
| k | 0.34 W/m/K | conductivity of the buried-heat path, sets G; calibrated on CFM |
| tau_w | 7 days (1 day in the CFM benchmark) | weather decorrelation time, sets ell_s; property of the forcing, measured (PROMICE), not fitted |
| rho_c | 830 kg m-3 | pore close-off (not fitted) |
| f_pf | 0 | optional preferential-flow bypass of the cold wave; 0 in v1.0 (CFM emulation) |

## Comparison with CFM

CFM variant `iav_clampT_imp04_dt2_gC_rain_ddf2b`: interannual variability, surface temperature
clamped at 0 C, impermeable layers >= 0.4 m, 2-day steps, grid C, rain above 1 C, degree-day melt as SnowDRIP's.
Reference climate T_mu = -9 C at +0 m (lapse 6.5 K/km), sigma_T = 7 K (70 % seasonal). The model
is driven with CFM's own snowfall, rain, melt and T_sfc over the averaging window.

* Sweep (+0/+100/+250 m, P = 0.4/1.0/1.5/2.5/4.0, 15 runs): rms T10m 0.66 K, runoff 0.10 of
  M + R.
* 15 extra runs around the transition (9 through the gap, 6 at low snowfall): 0.49 K.
* All 30 runs: 0.58 K.

## Limitations and open questions

* The transition is tested in one climate only; a maritime climate (sigma_T about 4 K) is the
  obvious next check.
* Cold edge of the transition: +250 m, P = 1.5 is 1.4 K too cold and 0 m, P = 0.6 1.4 K too cold.
  CFM's transition is smoother than the blend there, plausibly from interannual variability.
* The temperate P = 4 cases are about 1.2 K too cold (CFM's 20 m domain may bias them).
* Temperate runoff at +100 m is too low (0-0.03 vs 0.18-0.22 of M + R in CFM): the model
  retains all the water its temperate cold supply and capillary export allow.
* The near-surface structure (warm summer layer) is not represented; the model targets T below
  a few metres.
