"""Physical constants, forcing conventions and model parameters."""

from dataclasses import dataclass, replace

RHO_WATER = 1000.0
RHO_I = 917.0
RHO_SNOW = 350.0
RHO_CO = 830.0              # pore close-off
CP = 2097.0                 # J/kg/K (as CFM)
LF = 3.335e5                # J/kg (as CFM)
SEC_PER_YEAR = 31557600.0

DDF_SNOW, DDF_ICE = 0.004, 0.008   # m w.e./K/day, degree-day form (CFM benchmark)


@dataclass(frozen=True)
class Params:
    p_max: float = 1.0 - RHO_SNOW / RHO_I   # maximum retention per unit snowfall: pore space of the snow
    k: float = 0.34                         # W/m/K, conductivity of the buried-heat path (calibrated on CFM)
    sigma_w: float = 3.05                   # K, daily spread of temperature about the monthly mean
    tau_w: float = 7.0                      # days, decorrelation time of the weather (1 d for white-noise forcing)
    beta_T_snow: float = 2.85e-3            # m w.e./K/day         (PROMICE energy-balance melt)
    beta_SW_snow: float = 0.0546e-3         # m w.e./day per W m-2 (absorbed shortwave at albedo 0.79)
    beta_0_snow: float = -5.34e-3           # m w.e./day           (Giesen and Oerlemans 2012)
    beta_T_ice: float = 6.0e-3              # m w.e./K/day         (FoG series 1951-2025, ERA5 at 6.5 K/km and P x 1.70)
    beta_SW_ice: float = 0.0                # m w.e./day per W m-2 (0 in v1.0)
    rain_threshold: float = 1.0             # C
    f_pf: float = 0.0                       # preferential-flow bypass to the deep reservoir (0 in v1.0)


DEFAULT = Params()


def resolve(p=None, overrides=None):
    """p (default Params()) with the keyword overrides applied."""
    p = DEFAULT if p is None else p
    return replace(p, **overrides) if overrides else p


def degree_day_params(p=None, ddf_snow=DDF_SNOW, ddf_ice=DDF_ICE):
    """Params with the melt in degree-day form: no shortwave and no constant term."""
    return replace(resolve(p), beta_T_snow=ddf_snow, beta_SW_snow=0.0, beta_0_snow=0.0,
                   beta_T_ice=ddf_ice, beta_SW_ice=0.0)
