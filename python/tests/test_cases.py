import numpy as np
import pytest

import umwm
from umwm.config import unset_variables

CASES = ["fetch-limited", "duration-limited", "jonswap", "wave-tank", "titan-ontario-lacus", "mars-jezero-lake",
         "seamount", "eddy-field"]
WITH_INPUTS = {"seamount": ["umwm.gridtopo"], "eddy-field": ["umwmin.nc"]}


def test_list_cases(repo_root, tmp_path, monkeypatch):
    """The cases found by the package outside of the repository, in the copy
    in the package or in the repository for an editable install, are those of
    the repository, and each is listed in cases.txt."""
    examples = repo_root / "examples"
    monkeypatch.chdir(tmp_path)
    cases = umwm.list_cases()
    assert [case.name for case in cases] == CASES
    assert (umwm.examples_dir() / "cases.txt").read_text() == (examples / "cases.txt").read_text()
    for case in cases:
        assert case.description
        for name in ("main.nml", "spectrum.nml"):
            assert (case.directory / name).read_text() == (examples / case.name / name).read_text()
        inputs = sorted(path.name for path in case.input_dir.iterdir()) if case.input_dir else []
        assert inputs == WITH_INPUTS.get(case.name, [])
        for name in inputs:
            assert (case.input_dir / name).read_bytes() == (examples / case.name / "input" / name).read_bytes()
        assert case.notes().startswith(f"UMWM example case: {case.name}\n")
    assert {path.parent.name for path in examples.glob("*/main.nml")} == set(CASES)


@pytest.mark.parametrize("name", CASES)
def test_case_namelist_is_complete_and_valid(name):
    """Every case sets every namelist variable, so that the model, whose
    defaults differ from those of the package, reads the same configuration."""
    case = umwm.get_case(name)
    unset = unset_variables(case.namelist.read_text(encoding="utf-8"))
    assert not unset, f"{name} does not set {unset}"
    config = case.config()
    config.validate()
    assert (config.forcing_from_file or config.grid.topofromfile) == (case.input_dir is not None)
    assert config.forcing.steady or not config.forcing_from_file
    assert not (config.grid.gridfromfile or config.domain.restart)
    assert config.output.outrst == 0
    assert config.spectrum_points and all(point.by_index for point in config.spectrum_points)


def test_planetary_parameters():
    titan = umwm.Config.from_case("titan-ontario-lacus")
    assert (titan.physics.g, titan.forcing_constant.rhoa0, titan.forcing_constant.rhow0) == (1.352, 5.564, 588.15)
    assert (titan.physics.nu_air, titan.physics.nu_water, titan.physics.sfct) == (1.156e-6, 8.084e-7, 0.032766)
    mars = umwm.Config.from_case("mars-jezero-lake")
    assert (mars.physics.g, mars.forcing_constant.rhoa0, mars.forcing_constant.rhow0) == (3.71, 0.9188, 998.21)
    assert (mars.physics.nu_air, mars.physics.nu_water, mars.physics.sfct) == (1.576e-5, 1.156e-6, 0.073836)
    assert umwm.Config.from_case("wave-tank").physics.fcut_fac == 10.0
    assert umwm.Config.from_case("jonswap").physics.fcut_fac == umwm.Config().physics.fcut_fac


def test_get_case(tmp_path):
    assert umwm.get_case("titan").name == "titan-ontario-lacus"
    assert umwm.Config.from_case("jon").domain.nm == 101
    with pytest.raises(KeyError, match="unknown case"):
        umwm.get_case("venus")
    (tmp_path / "cases.txt").write_text("# comment\nlake-a  first\n\nlake-b  second case\n")
    assert [(case.name, case.description) for case in umwm.list_cases(tmp_path)] == [
        ("lake-a", "first"), ("lake-b", "second case")]
    assert umwm.get_case("lake-b", tmp_path).directory == tmp_path / "lake-b"
    with pytest.raises(KeyError, match="more than one case"):
        umwm.get_case("lake", tmp_path)


def test_examples_dir_from_environment(tmp_path, monkeypatch):
    monkeypatch.setenv("UMWM_EXAMPLES", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="UMWM_EXAMPLES"):
        umwm.examples_dir()
    (tmp_path / "cases.txt").write_text("only  the only case\n")
    assert umwm.examples_dir() == tmp_path
    assert [case.name for case in umwm.list_cases()] == ["only"]
    # Without it, the examples directory of the current directory comes
    # first, as for the executable.
    monkeypatch.delenv("UMWM_EXAMPLES")
    monkeypatch.chdir(tmp_path)
    (tmp_path / "examples").mkdir()
    (tmp_path / "examples" / "cases.txt").write_text("local  a case of the current directory\n")
    assert umwm.examples_dir() == (tmp_path / "examples").resolve()
    assert [case.name for case in umwm.list_cases()] == ["local"]


def test_namelist_text_keeps_comments(tmp_path):
    case = umwm.get_case("mars")
    config = case.config()
    assert config.namelist_text(case.namelist) == case.namelist.read_text(encoding="utf-8")
    config.forcing_constant.wspd0 = 3.0
    text = config.namelist_text(case.namelist)
    assert "Jezero" not in text and umwm.Config.from_namelist_text(text).forcing_constant.wspd0 == 3.0
    assert umwm.Config().namelist_text() == umwm.Config().to_namelist()
    # A namelist that leaves out a variable, which the model would set to
    # zero, is written in full, with the value of the package.
    partial = tmp_path / "main.nml"
    lines = case.namelist.read_text(encoding="utf-8").splitlines()
    partial.write_text("\n".join(line for line in lines if not line.lstrip().startswith("sbf_fac")) + "\n")
    assert unset_variables(partial.read_text()) == ["physics.sbf_fac"]
    config = umwm.Config.from_namelist(partial)
    assert config.namelist_text(partial) == config.to_namelist()
    assert not unset_variables(config.to_namelist())


@pytest.mark.parametrize("name", sorted(WITH_INPUTS))
def test_case_inputs(library, tmp_path, name):
    """The input files of a case are complete and match its grid, and are
    copied into the run directory; the steady forcing of the eddy field is
    one file, umwmin.nc."""
    case = umwm.get_case(name)
    domain = case.domain(tmp_path / name)
    assert not domain.input_dir.is_symlink()
    for path in case.input_dir.iterdir():
        copied = domain.input_dir / path.name
        assert not copied.is_symlink() and copied.read_bytes() == path.read_bytes()
    assert domain.input_problems() == []
    # An input file that differs is kept, unless overwritten.
    changed = domain.input_dir / WITH_INPUTS[name][0]
    changed.write_bytes(b"changed")
    with pytest.raises(FileExistsError, match="differs from the input file"):
        case.domain(tmp_path / name)
    assert changed.read_bytes() == b"changed"
    assert case.copy_inputs(domain.input_dir, overwrite=True) == [changed]
    assert changed.read_bytes() == (case.input_dir / changed.name).read_bytes()
    assert case.copy_inputs(domain.input_dir) == []
    if domain.config.forcing_from_file:
        domain.config.forcing.steady = False
        assert len(domain.missing_inputs()) == len(domain.config.step_times())


@pytest.mark.parametrize("name", CASES)
def test_case_runs(library, tmp_path, name):
    """Each case runs for one global time step from its namelist, which is
    copied with its comments into the run directory, and writes its output."""
    case = umwm.get_case(name)
    with case.domain(tmp_path / name, quiet=True) as domain:
        domain.advance()
        swh = domain.grid_field("swh")
        assert np.all(np.isfinite(swh[1:-1, 1:-1])) and np.nanmin(swh) > 0
        assert domain.namelist_path.read_text(encoding="utf-8") == case.namelist.read_text(encoding="utf-8")
        assert set(domain.spectrum_output()) == {point.name for point in domain.config.spectrum_points}
        assert umwm.output_files(domain.workdir)


def test_refraction_gradients(library, tmp_path):
    """Refraction turns the waves over the seamount at the phase speed
    gradient at fixed frequency times cg/cp, the rate at which the energy
    crosses it, and through the eddy field at the velocity gradient of the
    current, both by centered differences, away from the land rows. The
    action per unit ln(frequency) is the variance spectrum times k/cg."""
    def centered(field, axis, spacing):
        return (np.roll(field, -1, axis) - np.roll(field, 1, axis)) / (2 * spacing)

    def assert_matches(actual, expected):
        rows = slice(2, -2)
        scale = np.nanmax(np.abs(expected[rows]))
        assert scale > 0
        np.testing.assert_allclose(actual[rows], expected[rows], rtol=1e-4, atol=1e-6 * scale)

    with umwm.get_case("seamount").domain(tmp_path / "seamount", quiet=True) as domain:
        domain.initialize()
        grid = domain.config.grid
        cp = domain.remap(domain.field("cp0").T, lbound=domain.bounds("cp0")[1])  # (y, x, frequency)
        cg = domain.remap(domain.field("cg0").T, lbound=domain.bounds("cg0")[1])
        assert_matches(domain.remap(domain.field("cgdlncp_dx").T), cg / cp * centered(cp, 1, grid.delx))
        assert_matches(domain.remap(domain.field("cgdlncp_dy").T), cg / cp * centered(cp, 0, grid.dely))
        kovcg = domain.remap(domain.field("kovcg").T, lbound=domain.bounds("kovcg")[1])
        np.testing.assert_allclose(kovcg, domain.remap(domain.field("wavenumber").T) / cg, rtol=1e-5)

    with umwm.get_case("eddy-field").domain(tmp_path / "eddy-field", quiet=True) as domain:
        domain.advance()  # the current gradients are evaluated in the time steps
        grid = domain.config.grid
        uc, vc = domain.grid_field("uc"), domain.grid_field("vc")
        assert_matches(domain.grid_field("du_dx"), centered(uc, 1, grid.delx))
        assert_matches(domain.grid_field("du_dy"), centered(uc, 0, grid.dely))
        assert_matches(domain.grid_field("dv_dx"), centered(vc, 1, grid.delx))
        assert_matches(domain.grid_field("dv_dy"), centered(vc, 0, grid.dely))


def test_wave_tank_cutoff(library, tmp_path):
    """In the tank, the raised prognostic cutoff lets the waves grow along the
    fetch; with the default cutoff the whole spectrum is diagnostic and the
    wave height is the same, and far too high, all along the tank."""
    heights = {}
    for fcut_fac in (10.0, 0.53):
        config = umwm.Config.from_case("wave-tank")
        config.physics.fcut_fac = fcut_fac
        config.domain.nm = 21  # a 4 m tank
        config.output.ypl = 11
        config.spectrum_points = []
        config.domain.stoptimestr = "2026-01-01 00:00:20"
        config.output.outgrid = config.output.outspec = 0
        config.output.stokes = False
        with umwm.Domain(config, tmp_path / f"tank{fcut_fac}", quiet=True) as domain:
            domain.run()
            heights[fcut_fac] = domain.grid_field("swh")[1:-1, 0]
    growing, uniform = heights[10.0], heights[0.53]
    assert np.all(np.diff(growing) > 0) and 0.004 < growing[0] < 0.006 and 0.015 < growing[-1] < 0.025
    assert np.ptp(uniform) < 0.02 * uniform.mean() and uniform.mean() > 3 * growing.max()
