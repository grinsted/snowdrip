"""
Climate forcing at a glacier from the ERA5 climatologies (doi:10.5281/zenodo.24079220), as used to calibrate
the SnowDRIP melt coefficients: temperature moved with LAPSE_RATE, precipitation times PRECIP_FACTOR, and
shortwave radiation moved to the elevation of interest with the overlying air mass and cloud water. Needs xarray and netCDF4.
"""

import tomllib
from pathlib import Path

import numpy as np
import xarray as xr

DATA = Path(__file__).resolve().parents[1] / "data"
LAPSE_RATE = 6.5            # K/km, all months; beta_T_ice is calibrated with it
PRECIP_FACTOR = 1.70        # multiplier on ERA5 precipitation (median of FoG winter balances)
G0 = 9.80665                # m s-2
H_SW = 8000.0               # m, scale height of the clear-sky optical depth



def data_folder(name):
    """Folder of the data set `name`: data/<name> in the repository, unless data/folders.local.toml
    (not committed; see data/folders.example.toml) points elsewhere."""
    local = DATA / "folders.local.toml"
    if local.exists():
        with open(local, "rb") as f:
            folders = tomllib.load(f)
        if name in folders:
            return DATA / Path(folders[name]).expanduser()
    return DATA / name


def cloud_water_above(z, z_lev, q_lev, p_lev):
    """Cloud water path [kg m-2] above height z [m], from the specific cloud water q_lev [kg/kg] on pressure
    levels p_lev [Pa] (decreasing upwards) at heights z_lev [m]. z_lev and q_lev may be (..., n_levels)."""
    z_lev, q_lev, p_lev = (np.asarray(a, float) for a in (z_lev, q_lev, p_lev))
    z = np.asarray(z, float)
    qm, dp = 0.5 * (q_lev[..., 1:] + q_lev[..., :-1]), p_lev[:-1] - p_lev[1:]
    frac = np.clip((z_lev[..., 1:] - z[..., None]) / (z_lev[..., 1:] - z_lev[..., :-1]), 0.0, 1.0)
    wp = np.sum(qm * dp / G0 * frac, axis=-1)
    below = q_lev[..., 0] * (p_lev[0] * np.exp((z_lev[..., 0] - z) / 8000.0) - p_lev[0]) / G0
    out = wp + np.where(z < z_lev[..., 0], below, 0.0)
    return float(out) if np.ndim(out) == 0 else out


def cloud_water_share(z, z_ref, z_lev, q_lev, p_lev):
    """F(z): cloud water above z relative to that above the reanalysis surface z_ref."""
    return cloud_water_above(z, z_lev, q_lev, p_lev) / np.maximum(cloud_water_above(z_ref, z_lev, q_lev, p_lev), 1e-9)


def shortwave_at_elevation(z, z_ref, SW, SW_clear, SW_toa, F, H=H_SW):
    """
    Downward shortwave [W m-2] at elevation z from the reanalysis all-sky SW, clear-sky SW_clear and
    top-of-atmosphere SW_toa at z_ref, and the cloud-water share F:
        SW(z) = SW_toa tau^exp(-(z - z_ref)/H) [1 - a F],  tau = SW_clear/SW_toa,  a = 1 - SW/SW_clear.
    """
    SW, SW_clear, SW_toa, F = (np.asarray(a, float) for a in (SW, SW_clear, SW_toa, F))
    tau = np.clip(SW_clear / np.maximum(SW_toa, 1e-6), 1e-3, 1.0)
    a = np.clip(1.0 - SW / np.maximum(SW_clear, 1e-6), 0.0, 1.0)
    return SW_toa * tau ** np.exp(-(z - z_ref) / H) * np.clip(1.0 - a * np.clip(F, 0.0, None), 0.0, None)


class ERA5Forcing:
    """
    The climatological year at the nearest ERA5 grid cell; forcing(z) returns T_m [C], P_m [m w.e./yr] and
    SW_m [W m-2] at elevation z [m]. lapse_rate="era5" uses the monthly free-air lapse rate instead.
    """
    LEVELS = (850, 700, 500, 300)
    _DATASETS = {}          # kept open: reopening a netCDF file while another handle is open can crash netCDF

    @classmethod
    def _open(cls, path):
        key = str(path)
        if key not in cls._DATASETS:
            cls._DATASETS[key] = xr.open_dataset(path)
        return cls._DATASETS[key]

    def __init__(self, folder, lat, lon, precip_factor=PRECIP_FACTOR, lapse_rate=LAPSE_RATE):
        d = Path(folder)
        sel = dict(lat=lat, lon=lon % 360, method="nearest")

        def get(f, v):
            return self._open(d / f)[v].sel(**sel).values

        self.z_ref = float(get("era5_orography.nc", "z"))
        self.T_ref = get("era5_t2m_climatology.nc", "t2m") - 273.15
        self.gamma = (get("era5_lapserate.nc", "gamma") if isinstance(lapse_rate, str) and lapse_rate == "era5"
                      else np.full(12, float(lapse_rate)))
        self.P_m = get("era5_tp_climatology.nc", "tp") * precip_factor
        self.SW, self.SW_clear, self.SW_toa = (get(f"era5_{v}_climatology.nc", v) / 86400.0
                                               for v in ("ssrd", "ssrdc", "tisr"))
        cs = self._open(d / "era5_cloudwater_climatology_1deg.nc").sel(**sel)
        lev = cs.level.values                                               # hPa
        self.q_lev = (cs.clwc + cs.ciwc).transpose("month", "level").values
        self.p_lev = lev * 100.0
        Z = np.array([get(f"era5_z{p}_climatology.nc", f"z{p}") for p in self.LEVELS])
        self.z_lev = self._level_heights(Z, lev)
        self._sw = {}

    def _level_heights(self, Z, lev):
        """Heights of the cloud-water levels per month: geopotential heights at 850-300 hPa interpolated in
        ln p, extrapolated below 850 hPa with the 850-700 hPa slope."""
        pk, lp = np.log(np.array(self.LEVELS, float)), np.log(lev)
        out = np.empty((12, len(lev)))
        for m in range(12):
            slope = (Z[1, m] - Z[0, m]) / (pk[1] - pk[0])
            out[m] = np.where(lp > pk[0], Z[0, m] + slope * (lp - pk[0]), np.interp(-lp, -pk, Z[:, m]))
        return out

    def shortwave(self, z):
        key = round(float(z), 3)
        if key not in self._sw:
            F = cloud_water_share(z, self.z_ref, self.z_lev, self.q_lev, self.p_lev)
            self._sw[key] = shortwave_at_elevation(z, self.z_ref, self.SW, self.SW_clear, self.SW_toa, F)
        return self._sw[key].copy()

    def __call__(self, z, dT=0.0):
        return dict(T_m=self.T_ref - self.gamma * (z - self.z_ref) / 1000.0 + dT, P_m=self.P_m.copy(),
                    SW_m=self.shortwave(z))
