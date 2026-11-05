"""Tests of the data sources of the input files, with a small synthetic store
in the layout of the ARCO ERA5 store, and of the real ERA5 store when
UMWM_TEST_ERA5 is set."""

import math
import os
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

import umwm
from umwm.cli import main
from umwm.sources import SourceError, source_problems
from umwm.sources.regrid import Regridder, fill_missing, model_grid, open_boundary_cells
from umwm.sources.waves import regrid_partitions, standard_partitions, to_model_direction

START = datetime(2024, 1, 1)
HOURS = 7
LAND = (slice(10, 21), slice(0, 21))  # 70N to 50N, 0E to 40E: latitude and longitude indices


def synthetic_era5(path):
    """A store of all the ERA5 variables of the source on a 2-degree global grid
    for 7 hours from 2024-01-01: uniform fields, NaN on a block of land, with
    the wind speed growing by 1 m/s per hour."""
    from umwm.sources.era5 import VARIABLES

    lat = np.arange(90.0, -90.1, -2.0)
    lon = np.arange(0.0, 360.0, 2.0)
    shape = (HOURS, lat.size, lon.size)
    values = {
        "u10": 10.0 + np.arange(HOURS)[:, None, None], "v10": 0.0, "rhoa": 1.25, "fice": 0.0, "depth": 4000.0,
        "hs": math.sqrt(1 + 4 + 0.25), "tm": 8.0, "dir": 225.0, "spr": 0.6,
        "hs0": 1.0, "tm0": 5.0, "dir0": 270.0, "spr0": 0.5,
        "hs1": 2.0, "tm1": 12.0, "dir1": 180.0, "spr1": 0.3,
        "hs2": 0.0, "tm2": 0.0, "dir2": 0.0, "hs3": 0.0, "tm3": 0.0, "dir3": 0.0,
    }
    data_vars = {}
    for name, variable in VARIABLES.items():
        if variable in data_vars:
            continue
        field = np.broadcast_to(np.asarray(values[name], dtype=np.float32), shape).copy()
        if name not in ("u10", "v10"):
            field[(slice(None),) + LAND] = np.nan
        data_vars[variable] = (("time", "latitude", "longitude"), field)
    hours = (START - datetime(1900, 1, 1)) // timedelta(hours=1) + np.arange(HOURS, dtype=np.int64)
    dataset = xr.Dataset(data_vars, coords={
        "time": ("time", hours, {"units": "hours since 1900-01-01 00:00:00"}),
        "latitude": ("latitude", lat.astype(np.float32)),
        "longitude": ("longitude", lon.astype(np.float32))})
    dataset.to_zarr(path, consolidated=True, zarr_format=2)
    return path


@pytest.fixture(scope="module")
def era5_store(tmp_path_factory):
    pytest.importorskip("zarr")
    return synthetic_era5(tmp_path_factory.mktemp("era5") / "era5.zarr")


def era5_config(isglobal=False):
    """A configuration of a run from ERA5 on a 1-degree regional grid in the
    Indian Ocean, or a 10-degree global grid."""
    config = umwm.Config()
    config.domain.starttimestr = "2024-01-01 00:00:00"
    config.domain.stoptimestr = "2024-01-01 01:00:00"
    config.grid.topofromfile = True
    for name in ("winds", "air_density", "seaice", "waves"):
        setattr(config.forcing, name, True)
    config.forcing.source = "era5"
    if isglobal:
        config.domain.isglobal = True
        config.domain.mm, config.domain.nm = 36, 13
        config.grid.lon0, config.grid.lat0, config.grid.dlon, config.grid.dlat = 0.0, -60.0, 10.0, 10.0
        config.output.xpl, config.output.ypl = 10, 7
    else:
        config.domain.mm, config.domain.nm = 21, 11
        config.grid.lon0, config.grid.lat0, config.grid.dlon, config.grid.dlat = 70.0, -30.0, 1.0, 1.0
    config.validate()
    return config


# ---- regridding --------------------------------------------------------------

def test_regridder_periodic_and_descending():
    lon = np.arange(0.0, 360.0, 2.0)
    lat = np.arange(90.0, -90.1, -2.0)
    field = (np.cos(np.deg2rad(lon))[None, :] * np.ones((lat.size, 1)) + lat[:, None] / 100.0)
    points_lon = np.array([[359.0, -1.0, 181.0, 3.0]])
    points_lat = np.array([[0.5, -45.25, 10.0, 89.9]])
    regridder = Regridder(lon, lat, points_lon, points_lat)
    expected = np.cos(np.deg2rad(points_lon)) + points_lat / 100.0
    # bilinear between grid lines 2 degrees apart: within the curvature of the cosine
    assert np.allclose(regridder(field), expected, atol=2e-4)
    assert regridder.rows == slice(0, 69)


def test_regridder_missing_values():
    lon = np.arange(0.0, 360.0, 1.0)
    lat = np.arange(-10.0, 10.1, 1.0)
    field = np.ones((lat.size, lon.size))
    field[:, 1] = np.nan
    field[:, 2] = 3.0
    regridder = Regridder(lon, lat, np.array([0.5, 1.5, 1.25]), np.array([0.0, 0.0, 0.0]))
    # the missing corners are left out, with the weights renormalized
    assert np.allclose(regridder(field), [1.0, 3.0, 3.0])
    assert np.allclose(regridder.coverage(field), [0.5, 0.5, 0.25])
    field[:, 2] = np.nan
    assert np.isnan(regridder(field)[1])
    assert fill_missing(np.array([[1.0, np.nan, np.nan, np.nan]]), 7.0, iterations=2)[0].tolist() == [1, 1, 1, 7]


def test_model_grid(tmp_path):
    config = umwm.Config()
    with pytest.raises(ValueError, match="longitude-latitude"):
        model_grid(config)
    config.domain.mm, config.domain.nm = 5, 3
    config.grid.lon0, config.grid.lat0, config.grid.dlon, config.grid.dlat = 178.0, 10.0, 1.0, 0.5
    lon, lat = model_grid(config)
    assert lon[0].tolist() == [178.0, 179.0, -180.0, -179.0, -178.0] and lat[:, 0].tolist() == [10.0, 10.5, 11.0]
    config.grid.gridfromfile = True
    umwm.write_grid(tmp_path, lon + 0.5, lat)
    assert np.allclose(model_grid(config, tmp_path)[0], lon + 0.5)
    depth = np.full((3, 5), 100.0)
    assert not open_boundary_cells(config, depth)[[0, 2]].any() and open_boundary_cells(config, depth)[1, 1:4].all()
    config.domain.isglobal = True
    assert not open_boundary_cells(config, depth).any()


# ---- wave partitions ---------------------------------------------------------

def test_model_direction():
    # from north, from east, from south, from west
    assert np.allclose(to_model_direction([0.0, 90.0, 180.0, 270.0]), [-np.pi / 2, np.pi, np.pi / 2, 0.0])


def test_standard_partitions():
    fields = {"hs": np.array([3.0, 2.0, np.nan]), "tm": np.full(3, 8.0), "dir": np.full(3, 225.0),
              "spr": np.full(3, 0.6),
              "hs0": np.array([1.0, 2.0, np.nan]), "tm0": np.full(3, 5.0), "dir0": np.full(3, 270.0),
              "spr0": np.full(3, 0.5),
              "hs1": np.array([2.0, 1.0, np.nan]), "tm1": np.array([12.0, 0.0, np.nan]), "dir1": np.full(3, 180.0)}
    partitions = standard_partitions(fields, 2)
    assert len(partitions) == 3
    wind_sea, swell, residual = partitions
    # 1 + 4 < 9: the residual 4 m^2 has the parameters of the whole wave field
    assert np.allclose(residual.energy[0], 4.0) and residual.period[0] == 8.0
    assert np.isclose(residual.direction[0], np.deg2rad(45.0)) and residual.spread[0] == 0.6
    # the swell without a period is absent, and the wind sea is the whole wave field
    assert swell.energy[1] == 0.0 and np.isclose(wind_sea.energy[1], 4.0) and residual.energy[1] == 0.0
    assert np.isnan(wind_sea.energy[2]) and np.isnan(residual.energy[2])
    fields["hs"] = np.array([1.0, 2.0, np.nan])
    wind_sea, swell, residual = standard_partitions(fields, 2)
    # more energy in the partitions than in the whole wave field: scaled down
    assert np.isclose(wind_sea.energy[0] + swell.energy[0], 1.0) and residual.energy[0] == 0.0


def test_regrid_partitions():
    energy = np.zeros((2, 360))
    energy[:, 0] = 1.0
    energy[:, 1] = 9.0
    energy[:, 2] = np.nan
    period = np.where(energy > 4, 10.0, 6.0)
    direction = np.where(energy > 4, np.pi / 2, 0.0)
    from umwm.sources.waves import Partition

    variables = regrid_partitions(Regridder(np.arange(0.0, 360.0, 1.0), np.array([0.0, 1.0]),
                                            np.array([0.5, 2.0, 3.5]), np.zeros(3)),
                                  [Partition(energy, period, direction, np.full(energy.shape, 0.4))])
    # halfway between 1 and 9 m^2: the variance is interpolated, the period and
    # direction weighted by it
    assert np.isclose(variables["phs0"][0], math.sqrt(5.0))
    assert np.isclose(variables["ptm0"][0], (1 * 6 + 9 * 10) / 10)
    assert np.isclose(variables["pdir0"][0], math.atan2(9, 1)) and np.isclose(variables["pspr0"][0], 0.4)
    # no data at any corner
    assert variables["phs0"][1] == -1.0 and variables["ptm0"][1] == 0.0
    assert variables["phs0"][2] == 0.0


def test_write_forcing_partitions(tmp_path):
    path = umwm.write_forcing(tmp_path, "2026-01-01 00:00:00", phs0=np.array([[1.0, np.nan], [-3.0, 0.0]]),
                              ptm0=np.array([[5.0, np.nan], [5.0, 5.0]]), pdir0=0.0)
    with xr.open_dataset(path) as dataset:
        assert dataset.phs0.values[0].tolist() == [[1.0, -1.0], [-1.0, 0.0]]
        assert dataset.ptm0.values[0].tolist() == [[5.0, 0.0], [5.0, 5.0]]
        assert "_FillValue" not in dataset.phs0.encoding or np.isnan(dataset.phs0.encoding["_FillValue"]) is False
        assert dataset.pdir0.attrs["units"] == "rad"
    with pytest.raises(TypeError, match="phs10"):
        umwm.write_forcing(tmp_path, "2026-01-01 00:00:00", phs10=1.0, uw=np.zeros((2, 2)))


# ---- configuration -----------------------------------------------------------

def test_source_configuration():
    config = era5_config()
    assert config.source_name == "era5" and not source_problems(umwm.Config())
    config.forcing.source = "gfs-2050"
    with pytest.raises(umwm.ConfigError, match="unknown data source"):
        config.validate()
    config = era5_config()
    config.forcing.currents = True
    with pytest.raises(umwm.ConfigError, match="does not provide currents"):
        config.validate()
    config = era5_config()
    config.grid.dlon = config.grid.dlat = 0.0
    with pytest.raises(umwm.ConfigError, match="longitude-latitude grid"):
        config.validate()
    config = era5_config(isglobal=True)
    config.grid.dlon = 9.0
    with pytest.raises(umwm.ConfigError, match="mm \\* dlon = 360"):
        config.validate()
    # the new variables may be left out of a namelist, which the model defaults the same
    text = umwm.Config().to_namelist()
    assert "source" in text and "dlon" in text
    assert not umwm.config.unset_variables("\n".join(line for line in text.splitlines()
                                                     if not line.lstrip().startswith(("source", "dlon", "waves"))))


# ---- preparation of the input files ------------------------------------------

def test_prepare_regional(era5_store, tmp_path):
    from umwm.sources import prepare_inputs

    config = era5_config()
    config.domain.dtg = 1800
    paths = prepare_inputs(config, tmp_path, url=str(era5_store), log=lambda message: None)
    names = sorted(path.name for path in paths)
    assert names == ["umwm.gridtopo", "umwmin_2024-01-01_00:00:00.nc", "umwmin_2024-01-01_00:30:00.nc",
                     "umwmin_2024-01-01_01:00:00.nc"]
    with xr.open_dataset(tmp_path / "umwm.gridtopo") as grid:
        assert np.all(grid.z.values == -4000.0) and grid.lon.values[0, 0] == 70.0
    with xr.open_dataset(tmp_path / "umwmin_2024-01-01_00:30:00.nc") as forcing:
        values = {name: forcing[name].values[0] for name in forcing.data_vars}
    # the wind interpolated in time; the waves of the nearest time, at every
    # time for the open boundaries of the regional domain
    assert np.allclose(values["uw"], 10.5) and np.allclose(values["vw"], 0.0)
    assert np.allclose(values["rhoa"], 1.25) and np.allclose(values["fice"], 0.0)
    expected = {"phs0": 1.0, "ptm0": 5.0, "pdir0": 0.0, "pspr0": 0.5, "phs1": 2.0, "ptm1": 12.0,
                "pdir1": np.pi / 2, "pspr1": 0.3, "phs2": 0.0, "phs3": 0.0, "phs4": 0.5, "ptm4": 8.0,
                "pdir4": np.pi / 4, "pspr4": 0.6}
    for name, value in expected.items():
        assert np.allclose(values[name], value, atol=1e-6), name
    # in place: nothing to do, unless forced
    assert prepare_inputs(config, tmp_path, url=str(era5_store), log=lambda message: None) == []
    assert len(prepare_inputs(config, tmp_path, url=str(era5_store), force=True, log=lambda message: None)) == 4


def test_prepare_global(era5_store, tmp_path):
    from umwm.sources import prepare_inputs

    config = era5_config(isglobal=True)
    prepare_inputs(config, tmp_path, url=str(era5_store), log=lambda message: None)
    with xr.open_dataset(tmp_path / "umwm.gridtopo") as grid:
        land = grid.z.values >= 0  # the model is land where the bottom is not below sea level
        # 50N and 60N, 0E to 40E are in the block of land
        assert land[11:13, 0:5].all() and land.sum() == 10
    with xr.open_dataset(tmp_path / "umwmin_2024-01-01_00:00:00.nc") as first, \
            xr.open_dataset(tmp_path / "umwmin_2024-01-01_01:00:00.nc") as last:
        # the waves only at the start, without open boundaries
        assert "phs0" in first and "phs0" not in last
        assert np.all(first.phs0.values[0][land] == -1.0)
        assert np.allclose(last.uw.values, 11.0)


def test_prepare_errors(era5_store, tmp_path):
    from umwm.sources import prepare_inputs

    config = era5_config()
    config.domain.stoptimestr = "2024-01-01 12:00:00"
    with pytest.raises(SourceError, match="not covered"):
        prepare_inputs(config, tmp_path, url=str(era5_store), log=lambda message: None)
    with pytest.raises(SourceError, match="cannot open the ERA5 store"):
        prepare_inputs(era5_config(), tmp_path, url=str(tmp_path / "missing.zarr"), log=lambda message: None)
    with pytest.raises(SourceError, match="no data source"):
        prepare_inputs(umwm.Config(), tmp_path)


def test_prepare_for_another_grid(era5_store, tmp_path):
    """The input files that umwm prepared for another grid, such as those of
    another case in the same directory, are prepared again; one that it did
    not prepare stops the preparation, rather than being replaced."""
    from umwm.forcing import write_forcing
    from umwm.sources import prepare_inputs
    from umwm.sources.prepare import input_grid_problem, pending_inputs

    def quiet(message):
        pass

    prepare_inputs(era5_config(isglobal=True), tmp_path, url=str(era5_store), log=quiet)
    regional = era5_config()
    messages = []
    bathymetry, files = pending_inputs(regional, tmp_path, log=messages.append)
    assert bathymetry and [path.name for _, path in files] == ["umwmin_2024-01-01_00:00:00.nc",
                                                                "umwmin_2024-01-01_01:00:00.nc"]
    assert len(messages) == 1 and "3 input file(s)" in messages[0]
    assert "umwm.gridtopo: lon is 36 x 13, but the grid is 21 x 11 (mm x nm)" in messages[0]
    assert len(prepare_inputs(regional, tmp_path, url=str(era5_store), log=quiet)) == 3
    with xr.open_dataset(tmp_path / "umwm.gridtopo") as grid:
        assert grid.z.shape == (11, 21) and grid.lon.values[0, 0] == 70.0
    assert pending_inputs(regional, tmp_path) == (False, [])

    # another region of the same size, which only the coordinates tell
    assert input_grid_problem(tmp_path / "umwm.gridtopo", (11, 21), ((70.0, -30.0), (90.0, -20.0))) is None
    assert input_grid_problem(tmp_path / "umwm.gridtopo", (11, 21), ((-290.0, -30.0), (-270.0, -20.0))) is None
    assert input_grid_problem(tmp_path / "umwm.gridtopo", (11, 21), ((80.0, -30.0), (100.0, -20.0))) == (
        "its first and last cells are at (70.000, -30.000) and (90.000, -20.000), but those of the grid are at "
        "(80.000, -30.000) and (100.000, -20.000)")
    shifted = era5_config()
    shifted.grid.lon0 = 80.0
    bathymetry, files = pending_inputs(shifted, tmp_path)
    assert bathymetry and len(files) == 2

    # a file for another grid that umwm did not prepare is not replaced
    handmade = tmp_path / "handmade"
    write_forcing(handmade, START, uw=np.zeros((5, 5)), vw=np.zeros((5, 5)))
    with pytest.raises(SourceError, match="for another grid: uw is 5 x 5, but the grid is 21 x 11.*not replaced"):
        pending_inputs(regional, handmade)


def test_domain_from_another_grid(library, era5_store, tmp_path, monkeypatch):
    """A domain whose input directory has the input files of another grid,
    prepared by umwm, prepares them again; run --check tells so."""
    monkeypatch.setenv("UMWM_ERA5_URL", str(era5_store))
    from umwm.sources import prepare_inputs

    prepare_inputs(era5_config(isglobal=True), tmp_path / "run" / "input", log=lambda message: None)
    with umwm.Domain(era5_config(), tmp_path / "run", quiet=True) as domain:
        assert domain._inputs_to_prepare()
        domain.initialize()
        assert not domain._inputs_to_prepare()
        domain.diagnostics()
        assert np.allclose(domain.grid_field("swh")[1:-1, 1:-1], math.sqrt(5.25), rtol=1e-4)


def test_domain_from_era5(library, era5_store, tmp_path, monkeypatch):
    """A domain with ERA5 as its source prepares its input files and starts
    from the significant wave height of ERA5, with open boundaries."""
    monkeypatch.setenv("UMWM_ERA5_URL", str(era5_store))
    with umwm.Domain(era5_config(), tmp_path / "run", quiet=True) as domain:
        domain.initialize()
        assert len(domain.missing_inputs()) == 0
        domain.diagnostics()
        swh = domain.grid_field("swh")
        assert np.allclose(swh[1:-1, 1:-1], math.sqrt(5.25), rtol=1e-4)
        assert domain.config.forcing.source == "era5"
        domain.run()
        # the wind of the last source time step, just before the stop time
        assert np.isclose(domain.grid_field("wspd")[5, 10], 11.0, rtol=1e-3)


def test_cli_forcing_era5(era5_store, tmp_path, capsys):
    workdir = tmp_path / "run"
    config = era5_config()
    config.write(workdir / "namelists" / "main.nml")
    assert main(["run", str(workdir), "--check"]) == 0
    assert "3 missing input file(s) will be prepared from era5" in capsys.readouterr().out
    assert main(["forcing", "era5", str(workdir), "--url", str(era5_store)]) == 0
    assert "prepared 3 input files" in capsys.readouterr().out
    assert sorted(path.name for path in (workdir / "input").iterdir())[0] == "umwm.gridtopo"
    assert main(["forcing", "era5", str(workdir), "--url", str(era5_store)]) == 0
    assert "in place" in capsys.readouterr().out
    # the files of another grid, prepared by umwm, are prepared again
    other = era5_config()
    other.grid.lon0 = 80.0
    other.write(workdir / "namelists" / "main.nml")
    assert main(["run", str(workdir), "--check"]) == 0
    assert "3 input file(s) prepared for another grid will be prepared from era5" in capsys.readouterr().out
    assert main(["forcing", "era5", str(workdir), "-s", "forcing.currents=.true.", "--url", str(era5_store)]) == 1
    assert "does not provide currents" in capsys.readouterr().err


@pytest.mark.skipif(not os.environ.get("UMWM_TEST_ERA5"), reason="set UMWM_TEST_ERA5=1 to read the ARCO ERA5 store")
def test_real_era5(tmp_path):
    """Two hours of the Gulf of Mexico during hurricane Helene from the ARCO
    ERA5 store: the input files, with ERA5 waves at the open boundaries."""
    pytest.importorskip("zarr")
    from umwm.sources import prepare_inputs

    config = era5_config()
    config.domain.starttimestr, config.domain.stoptimestr = "2024-09-26 00:00:00", "2024-09-26 01:00:00"
    config.domain.mm, config.domain.nm = 79, 57
    config.grid.lon0, config.grid.lat0, config.grid.dlon, config.grid.dlat = -98.5, 17.5, 0.25, 0.25
    paths = prepare_inputs(config, tmp_path, log=lambda message: None)
    assert len(paths) == 3
    with xr.open_dataset(paths[1]) as forcing:
        assert 20 < float(np.hypot(forcing.uw, forcing.vw).max()) < 40
        assert 5 < float(forcing.phs0.max()) < 15
