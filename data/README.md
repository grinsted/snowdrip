# Input data

The fields used to drive SnowDRIP in the v1.0 paper are not kept in this repository. They are
archived on Zenodo: ERA5 monthly climatologies 1950-2025, lapse rates and invariant fields,
https://doi.org/10.5281/zenodo.24079220. `examples/era5_forcing.py` reads them from `data/era5_climatology/`:

- 2 m temperature and total precipitation, with the ERA5 orography as their reference height;
- the monthly free-air lapse rate at the annual-mean -5 C isotherm;
- the all-sky, clear-sky and top-of-atmosphere shortwave radiation, the geopotential heights at
  850-300 hPa and the cloud liquid and ice water on 20 pressure levels (1 deg), which move the
  shortwave radiation to the elevation of interest.

Any source of monthly mean temperature, precipitation and shortwave radiation at the glacier can be
used instead; the melt coefficients are calibrated with this forcing.

To keep the data elsewhere on your machine, copy `folders.example.toml` to `folders.local.toml` (not
committed) and set the path there.
