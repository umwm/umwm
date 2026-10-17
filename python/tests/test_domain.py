"""Tests of the model driven from Python; they need the built library."""

from datetime import datetime

import numpy as np
import pytest

import umwm
from conftest import RANKINE_REFERENCE, matches_reference


@pytest.fixture
def short_config():
    config = umwm.Config()
    config.domain.stoptimestr = "2026-01-01 02:00:00"
    config.output.outrst = 1
    config.output.outspec = 1
    config.spectrum_points = [umwm.SpectrumPoint("center", x=11, y=6)]
    return config


def test_lifecycle_constant_wind(library, short_config, tmp_path):
    with umwm.Domain(short_config, tmp_path / "run") as domain:
        assert not domain.initialized
        assert domain.time == datetime(2026, 1, 1, 0)
        assert domain.rank == 0 and domain.size == 1
        assert (domain.workdir / "input").is_dir() and (domain.workdir / "output").is_dir()

        domain.initialize()
        assert domain.initialized and not domain.finished
        assert (domain.workdir / "namelists" / "main.nml").is_file()
        assert (domain.workdir / "namelists" / "spectrum.nml").read_text().split() == ["XY", "11", "6", "center"]
        sizes = domain.sizes
        assert sizes["mm"] == 21 and sizes["nm"] == 11 and sizes["im"] == 171 and sizes["imm"] == 231
        assert sizes["om"] == 37 and sizes["pm"] == 36 and sizes["lm"] == 42
        assert domain.shape == (11, 21)
        assert domain.attribute("dth") == pytest.approx(2 * np.pi / 36)

        times = list(domain.iterate())
        assert times == [datetime(2026, 1, 1, 1), datetime(2026, 1, 1, 2)]
        assert domain.finished and domain.time == times[-1]

        swh = domain.grid_field("swh")
        assert swh.shape == (11, 21) and swh.dtype == np.float32
        assert np.isnan(swh).sum() == 231 - 171  # land cells
        assert 0 < np.nanmin(swh) <= np.nanmax(swh) < 10

        # output files and in-memory state agree
        output = domain.output()
        assert output.sizes["time"] == 3 and output.sizes["z"] == 42
        np.testing.assert_allclose(output.swh.isel(time=-1).values, swh, rtol=1e-6, equal_nan=True)
        state = domain.state(ocean_stress=True)
        assert state.time.values[0] == np.datetime64("2026-01-01T02:00:00")
        for name in ("swh", "mwp", "cd", "taux_ocn", "epsx_ocn", "d_stokes"):
            np.testing.assert_allclose(state[name].values[0], output[name].isel(time=-1).values, rtol=1e-5, equal_nan=True)
        np.testing.assert_allclose(state.u_stokes.values[0], output.u_stokes.isel(time=-1).values, rtol=1e-5, equal_nan=True)
        assert state.swh.attrs["units"] == "m"

        # spectrum on the grid, its point output, and its integral
        spectrum = domain.spectrum()
        assert spectrum.dims == ("y", "x", "direction", "frequency")
        assert spectrum.shape == (11, 21, 36, 37)
        points = domain.spectrum_output()
        assert set(points) == {"center"} and points["center"].sizes["time"] == 3
        np.testing.assert_allclose(points["center"].F.isel(time=-1).values, spectrum.isel(y=5, x=10).values, rtol=1e-6)
        kdk = domain.remap(domain.field("kdk").T)  # (y, x, frequency)
        variance = (spectrum.values * kdk[:, :, None, :]).sum(axis=(2, 3)) * domain.attribute("dth")
        np.testing.assert_allclose(4 * np.sqrt(variance), swh, rtol=1e-4, equal_nan=True)

        assert [path.name for path in domain.restart_files()] == [
            "umwmrst_2026-01-01_01:00:00.nc", "umwmrst_2026-01-01_02:00:00.nc"]
        assert dict(domain.grid().sizes) == {"y": 11, "x": 21}
    assert not domain.initialized


def test_views_and_fields(library, short_config, tmp_path):
    with umwm.Domain(short_config, tmp_path / "run") as domain:
        domain.initialize()
        spectrum = domain.view("F")
        assert spectrum.shape == (37, 36, 172) and spectrum.flags.f_contiguous
        assert domain.bounds("F") == (1, 1, 0)
        assert domain.field("frequency").shape == (37,)
        assert domain.field("wavenumber").shape == (37, 171)
        assert domain.field("ii").dtype == np.int32 and domain.field("ii").shape == (21, 11)
        mask = domain.grid_field("seamask")
        assert mask.shape == (11, 21) and mask.sum() == 171
        assert np.all(mask[0, :] == 0) and np.all(mask[:, -1] == 0)
        depth = domain.grid_field("depth")
        assert np.all(depth[mask == 1] == 1000.0)
        with pytest.raises(KeyError):
            domain.field("no_such_field")
        with pytest.raises(KeyError):
            domain.grid_field("wavenumber")
        # the view reflects the model as it advances
        before = spectrum[:, :, 1:].sum()
        domain.advance()
        assert spectrum[:, :, 1:].sum() > before
        # cell-indexed arrays remap to the grid
        cells = domain.field("swh")
        grid = domain.remap(cells)
        np.testing.assert_allclose(grid, domain.grid_field("swh"), equal_nan=True)
        assert domain.remap(domain.field("wspd")).shape == (11, 21)
        # writing to a view of the constant forcing changes the model
        wind = domain.view("wspd")
        wind[:] = 25.0
        domain.advance()
        assert np.nanmax(domain.grid_field("wspd")) == pytest.approx(25.0)


def test_run_until_and_restart_run(library, short_config, tmp_path):
    with umwm.Domain(short_config, tmp_path / "run") as domain:
        domain.run(until="2026-01-01 01:00:00")
        assert domain.time == datetime(2026, 1, 1, 1) and not domain.finished
        domain.run()
        assert domain.finished
    restart = short_config.copy()
    restart.domain.restart = True
    restart.domain.starttimestr = "2026-01-01 01:00:00"
    with umwm.Domain(restart, tmp_path / "run") as domain:
        domain.run()
        assert domain.time == datetime(2026, 1, 1, 2)
        assert np.nanmax(domain.grid_field("swh")) > 0


def test_rankine_regression(rankine_domain):
    """The idealized hurricane case reproduces the reference values of tests/test_domain.f90."""
    domain = rankine_domain
    assert domain.missing_inputs() == []
    domain.run()
    mismatches = []
    for name, (mean, maximum, tolerance) in RANKINE_REFERENCE.items():
        values = domain.grid_field(name)
        if not matches_reference(values, mean, maximum, tolerance):
            mismatches.append(f"{name}: mean {np.nanmean(values):.7g} max {np.nanmax(values):.7g} "
                              f"expected {mean:.7g} {maximum:.7g}")
    assert not mismatches, "\n".join(mismatches)
    points = domain.spectrum_output()
    assert list(points) == ["rankine"] and points["rankine"].sizes["time"] == 2
    assert points["rankine"].F.dims == ("time", "direction", "frequency")


def test_invalid_configuration_raises(library, tmp_path):
    config = umwm.Config()
    config.domain.mm = 2
    with pytest.raises(umwm.ConfigError, match="mm and nm"):
        umwm.Domain(config, tmp_path / "bad").initialize()


def test_missing_forcing_raises(library, tmp_path):
    config = umwm.Config()
    config.forcing.winds = True
    config.domain.stoptimestr = "2026-01-01 01:00:00"
    domain = umwm.Domain(config, tmp_path / "run")
    assert [path.name for path in domain.missing_inputs()] == [
        "umwmin_2026-01-01_00:00:00.nc", "umwmin_2026-01-01_01:00:00.nc"]
    with pytest.raises(FileNotFoundError, match="umwmin_2026-01-01_01:00:00.nc"):
        domain.initialize()
    assert not domain.initialized


def test_namelist_file_is_preserved(library, rankine_config, repo_root, tmp_path):
    """A namelist the configuration was read from is not rewritten unless changed."""
    workdir = tmp_path / "run"
    (workdir / "namelists").mkdir(parents=True)
    source = repo_root / "tests" / "namelists" / "forcing_rankine.nml"
    target = workdir / "namelists" / "main.nml"
    target.write_text(source.read_text())
    domain = umwm.Domain(target, workdir)
    domain.write_namelists()
    assert target.read_text() == source.read_text()
    domain.config.domain.stoptimestr = "2026-01-01 02:00:00"
    domain.write_namelists()
    assert target.read_text() != source.read_text()
    assert umwm.Config.from_namelist(target).domain.stoptimestr == "2026-01-01 02:00:00"


def test_input_dir_link(library, rankine_config, tmp_path):
    inputs = tmp_path / "shared-inputs"
    inputs.mkdir()
    domain = umwm.Domain(rankine_config, tmp_path / "run", input_dir=inputs)
    assert domain.input_dir.is_symlink() and domain.input_dir.resolve() == inputs
    with pytest.raises(FileExistsError):
        umwm.Domain(rankine_config, tmp_path / "run", input_dir=tmp_path)
    with pytest.raises(FileNotFoundError):
        umwm.Domain(rankine_config, tmp_path / "other", input_dir=tmp_path / "nope")


def test_run_function(library, tmp_path):
    config = umwm.Config()
    config.domain.stoptimestr = "2026-01-01 01:00:00"
    config.output.stokes = False
    domain = umwm.run(config, tmp_path / "run")
    try:
        assert domain.finished
        assert "u_stokes" not in domain.output()
        assert "u_stokes" not in domain.state()
    finally:
        domain.finalize()
