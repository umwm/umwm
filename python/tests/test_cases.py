from dataclasses import fields

import numpy as np
import pytest

import umwm
from umwm import namelist
from umwm.config import GROUPS

CASES = ["fetch-limited", "duration-limited", "jonswap", "wave-tank", "titan-ontario-lacus", "mars-jezero-lake"]

# Variables that the model defaults to the same values as the package; the
# model defaults the others to zero.
OPTIONAL = {"domain.reftimestr", "physics.fcut_fac"}


def test_list_cases(repo_root):
    """The cases found by the package, in the repository or the copy in the
    package, are those of the repository, and each is listed in cases.txt."""
    examples = repo_root / "examples"
    cases = umwm.list_cases()
    assert [case.name for case in cases] == CASES
    assert (umwm.examples_dir() / "cases.txt").read_text() == (examples / "cases.txt").read_text()
    for case in cases:
        assert case.description
        for name in ("main.nml", "spectrum.nml"):
            assert (case.directory / name).read_text() == (examples / case.name / name).read_text()
        assert case.notes().startswith(f"UMWM example case: {case.name}\n")
    assert {path.parent.name for path in examples.glob("*/main.nml")} == set(CASES)


@pytest.mark.parametrize("name", CASES)
def test_case_namelist_is_complete_and_valid(name):
    """Every case sets every namelist variable, so that the model, whose
    defaults differ from those of the package, reads the same configuration."""
    case = umwm.get_case(name)
    groups = namelist.loads(case.namelist.read_text(encoding="utf-8"))
    for group_name, group in GROUPS.items():
        assert group_name in groups
        missing = {item.name for item in fields(group)
                   if f"{group_name}.{item.name}" not in OPTIONAL} - set(groups[group_name])
        assert not missing, f"{name}: &{group_name.upper()} does not set {sorted(missing)}"
    config = case.config()
    config.validate()
    assert not config.forcing_from_file
    assert not (config.grid.gridfromfile or config.grid.topofromfile or config.domain.restart)
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


def test_namelist_text_keeps_comments():
    case = umwm.get_case("mars")
    config = case.config()
    assert config.namelist_text(case.namelist) == case.namelist.read_text(encoding="utf-8")
    config.forcing_constant.wspd0 = 3.0
    text = config.namelist_text(case.namelist)
    assert "Jezero" not in text and umwm.Config.from_namelist_text(text).forcing_constant.wspd0 == 3.0
    assert umwm.Config().namelist_text() == umwm.Config().to_namelist()


@pytest.mark.parametrize("name", CASES)
def test_case_runs(library, tmp_path, name):
    """Each case runs for one global time step from its namelist, which is
    copied with its comments into the run directory, and writes its output."""
    case = umwm.get_case(name)
    with umwm.Domain(case.namelist, tmp_path / name, quiet=True) as domain:
        domain.advance()
        swh = domain.grid_field("swh")
        assert np.all(np.isfinite(swh[1:-1, 1:-1])) and np.nanmin(swh) > 0
        assert domain.namelist_path.read_text(encoding="utf-8") == case.namelist.read_text(encoding="utf-8")
        assert set(domain.spectrum_output()) == {point.name for point in domain.config.spectrum_points}
        assert umwm.output_files(domain.workdir)


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
