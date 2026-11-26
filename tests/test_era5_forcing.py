import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("xarray")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
import era5_forcing as ef  # noqa: E402


def test_shortwave_at_elevation():
    lev = np.array([1000.0, 850.0, 700.0, 500.0, 300.0])
    p_lev, z_lev = lev * 100.0, np.array([100.0, 1500.0, 3000.0, 5600.0, 9200.0])
    q = np.array([2e-4, 3e-4, 1e-4, 2e-5, 0.0])
    F = [ef.cloud_water_share(z, 500.0, z_lev, q, p_lev) for z in (500.0, 1500.0, 3000.0, 6000.0)]
    assert F[0] == pytest.approx(1.0) and np.all(np.diff(F) < 0) and F[-1] >= 0.0
    SW, SWc, toa = 180.0, 260.0, 450.0
    assert ef.shortwave_at_elevation(500.0, 500.0, SW, SWc, toa, 1.0) == pytest.approx(SW)
    sw = [ef.shortwave_at_elevation(z, 500.0, SW, SWc, toa, f) for z, f in zip((500.0, 1500.0, 3000.0, 6000.0), F)]
    assert np.all(np.diff(sw) > 0) and sw[-1] < toa
    # with no cloud left above, the clear-sky radiation at that height
    assert ef.shortwave_at_elevation(500.0, 500.0, SW, SWc, toa, 0.0) == pytest.approx(SWc)
