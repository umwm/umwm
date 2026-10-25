from datetime import datetime

import numpy as np
import pytest
import xarray as xr

import umwm
from umwm.forcing import rankine_vortex


def test_forcing_filename():
    assert umwm.forcing_filename(datetime(2026, 1, 2, 3, 4, 5)) == "umwmin_2026-01-02_03:04:05.nc"
    assert umwm.forcing_filename("2026-01-02 03:04:05") == "umwmin_2026-01-02_03:04:05.nc"


def test_forcing_files():
    config = umwm.Config()
    config.domain.stoptimestr = "2026-01-01 02:00:00"
    assert umwm.forcing_files(config) == []
    config.forcing.steady = True
    assert umwm.forcing_files(config) == []  # no forcing from file
    config.forcing.currents = True
    assert umwm.forcing_files(config) == ["umwmin.nc"]
    config.forcing.steady = False
    assert umwm.forcing_files(config) == [umwm.forcing_filename(time) for time in config.step_times(include_start=True)]
    assert len(umwm.forcing_files(config)) == 3


def test_write_forcing(tmp_path):
    ny, nx = 5, 7
    uw = np.arange(ny * nx, dtype=np.float64).reshape(ny, nx)
    path = umwm.write_forcing(tmp_path, "2026-01-01 06:00:00", uw=uw, vw=-uw, rhoa=1.25,
                              lon=np.zeros((ny, nx)), lat=np.ones((ny, nx)), x=np.arange(nx), y=np.arange(ny))
    assert path == tmp_path / "umwmin_2026-01-01_06:00:00.nc"
    with xr.open_dataset(path) as dataset:
        assert dataset.uw.dims == ("time", "y", "x")
        assert dataset.uw.dtype == np.float32
        assert dataset.sizes == {"time": 1, "y": ny, "x": nx}
        np.testing.assert_array_equal(dataset.uw.values[0], uw.astype(np.float32))
        np.testing.assert_array_equal(dataset.rhoa.values[0], np.full((ny, nx), 1.25, dtype=np.float32))
        assert dataset.time.values[0] == np.datetime64("2026-01-01T06:00:00")
        assert dataset.attrs["valid_time"] == "2026-01-01 06:00:00"
        assert "uc" not in dataset


def test_write_forcing_to_file_path(tmp_path):
    path = umwm.write_forcing(tmp_path / "custom.nc", "2026-01-01 00:00:00", uw=np.zeros((3, 4)), vw=np.zeros((3, 4)))
    assert path == tmp_path / "custom.nc" and path.is_file()


def test_write_forcing_errors(tmp_path):
    with pytest.raises(ValueError, match="at least one"):
        umwm.write_forcing(tmp_path, "2026-01-01 00:00:00")
    with pytest.raises(ValueError, match="2-D"):
        umwm.write_forcing(tmp_path, "2026-01-01 00:00:00", rhoa=1.2)
    with pytest.raises(ValueError, match="shape"):
        umwm.write_forcing(tmp_path, "2026-01-01 00:00:00", uw=np.zeros((3, 4)), vw=np.zeros((4, 3)))


def test_write_grid(tmp_path):
    lon, lat = np.meshgrid(np.linspace(-80, -79, 6), np.linspace(25, 26, 4))
    depth = np.full(lon.shape, 50.0)
    depth[0, :] = 0.0  # land
    path = umwm.write_grid(tmp_path, lon, lat)
    assert path.name == "umwm.grid"
    path = umwm.write_grid(tmp_path, lon, lat, depth=depth)
    assert path.name == "umwm.gridtopo"
    with xr.open_dataset(path) as dataset:
        assert dataset.lon.dims == ("y", "x")
        np.testing.assert_array_equal(dataset.z.values, -depth.astype(np.float32))


def test_rankine_vortex():
    x = np.linspace(-200e3, 200e3, 81, dtype=np.float32)
    x2d, y2d = np.meshgrid(x, x)
    uw, vw = rankine_vortex(x2d, y2d, vmax=30.0, rmw=50e3)
    speed = np.hypot(uw, vw)
    assert uw.dtype == np.float32
    assert speed.max() == pytest.approx(30.0)
    assert speed[40, 40] == 0.0
    # counter-clockwise: east of the center the wind blows north
    assert vw[40, 50] > 0 and abs(uw[40, 50]) < 1e-5
    uw_cw, vw_cw = rankine_vortex(x2d, y2d, vmax=30.0, rmw=50e3, clockwise=True)
    np.testing.assert_allclose(uw_cw, -uw)


def test_write_rankine_forcing(tmp_path):
    paths = umwm.write_rankine_forcing(tmp_path, "2026-01-01 00:00:00", 2, nx=11, ny=7, center_x=0.0, translation_x=10e3)
    assert [path.name for path in paths] == [f"umwmin_2026-01-01_0{hour}:00:00.nc" for hour in range(3)]
    with xr.open_dataset(paths[2]) as dataset:
        assert dataset.sizes == {"time": 1, "y": 7, "x": 11}
        assert dataset.attrs["vortex_center_x_m"] == 20e3
        assert set(dataset.data_vars) == {"lon", "lat", "uw", "vw", "uc", "vc", "rhoa", "rhow"}


def test_write_forcing_dataset(tmp_path):
    times = np.array(["2026-01-01T00:00:00", "2026-01-01T01:00:00"], dtype="datetime64[s]")
    ny, nx = 4, 6
    u10 = np.arange(2 * ny * nx, dtype=np.float32).reshape(2, ny, nx)
    dataset = xr.Dataset(
        {"u10": (("time", "y", "x"), u10), "v10": (("time", "y", "x"), -u10), "rhoa": (("time", "y", "x"), np.full_like(u10, 1.2)),
         "lon": (("y", "x"), np.zeros((ny, nx))), "lat": (("y", "x"), np.ones((ny, nx))), "extra": (("time",), [1, 2])},
        coords={"time": times},
    )
    paths = umwm.write_forcing_dataset(tmp_path, dataset, rename={"u10": "uw", "v10": "vw"})
    assert [path.name for path in paths] == ["umwmin_2026-01-01_00:00:00.nc", "umwmin_2026-01-01_01:00:00.nc"]
    with xr.open_dataset(paths[1]) as written:
        np.testing.assert_array_equal(written.uw.values[0], u10[1])
        np.testing.assert_array_equal(written.vw.values[0], -u10[1])
        assert set(written.data_vars) == {"lon", "lat", "uw", "vw", "rhoa"}
    # dimensions in (x, y) order are transposed
    transposed = dataset.transpose("time", "x", "y")
    assert umwm.write_forcing_dataset(tmp_path / "t", transposed, rename={"u10": "uw", "v10": "vw"})
    with xr.open_dataset(tmp_path / "t" / "umwmin_2026-01-01_00:00:00.nc") as written:
        assert written.uw.dims == ("time", "y", "x")
    with pytest.raises(ValueError, match="time dimension"):
        umwm.write_forcing_dataset(tmp_path, dataset.isel(time=0))
    with pytest.raises(ValueError, match="none of the forcing variables"):
        umwm.write_forcing_dataset(tmp_path, dataset[["extra"]])
