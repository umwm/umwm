from datetime import datetime, timedelta

import pytest

import umwm
from umwm.config import GROUPS, Config, ConfigError, SpectrumPoint, parse_value


def test_defaults_match_repository_namelist(repo_root):
    """The Python defaults are those of namelists/main.nml in the repository."""
    config = Config.from_namelist(repo_root / "namelists" / "main.nml")
    defaults = Config()
    for group in GROUPS:
        assert getattr(config, group) == getattr(defaults, group), group
    # spectrum.nml and exclude.nml next to main.nml are picked up
    assert [point.name for point in config.spectrum_points] == ["test_loc1", "test_loc2", "test_loc3"]
    assert config.exclude_points == [(5, 5)]


def test_rankine_namelist(repo_root):
    config = Config.from_namelist(repo_root / "tests" / "namelists" / "forcing_rankine.nml")
    assert config.forcing.winds and config.forcing.currents
    assert not config.forcing.seaice
    assert config.forcing_from_file
    assert config.stokes.depths == [0.1]
    assert config.spectrum_points == [SpectrumPoint("rankine", x=16, y=6)]
    assert config.stop_time - config.start_time == timedelta(hours=1)
    assert config.step_times() == [datetime(2026, 1, 1, 0), datetime(2026, 1, 1, 1)]


def test_missing_group_raises():
    with pytest.raises(ConfigError, match="&PHYSICS"):
        Config.from_namelist_text("&DOMAIN\n mm = 21\n/\n")


def test_set_and_get_with_conversion():
    config = Config()
    config.set("domain.mm", "31")
    config["forcing_constant.wspd0"] = "15"
    config.set("forcing.winds", ".true.")
    config.set("output.stokes", "F")
    config.set("stokes.depths", "0.5, 1.5 2.5")
    config.update({"physics.g": 3.71}, domain__dtg=1800)
    assert config.domain.mm == 31 and isinstance(config.domain.mm, int)
    assert config.forcing_constant.wspd0 == 15.0
    assert config.forcing.winds is True and config.output.stokes is False
    assert config.stokes.depths == [0.5, 1.5, 2.5]
    assert config["physics.g"] == 3.71 and config.domain.dtg == 1800.0
    assert dict(config.items())["domain.mm"] == 31


@pytest.mark.parametrize("key, value", [
    ("domain.nope", 1), ("nope.mm", 1), ("mm", 1), ("domain.mm", "many"), ("forcing.winds", "maybe"),
    ("domain.mm", 2.5),
])
def test_set_invalid_raises(key, value):
    with pytest.raises(ConfigError):
        Config().set(key, value)


def test_time_properties():
    config = Config()
    config.start_time = datetime(2027, 3, 4, 5, 6, 7)
    config.stop_time = "2027-03-04_11:06:07"
    assert config.domain.starttimestr == "2027-03-04 05:06:07"
    assert config.domain.stoptimestr == "2027-03-04 11:06:07"
    assert len(config.step_times(include_start=False)) == 6
    assert config.time_step == timedelta(hours=1)


def test_dict_and_json_round_trip():
    config = Config()
    config.domain.mm = 41
    config.spectrum_points = [SpectrumPoint("a", x=3, y=4), SpectrumPoint("b", x=5, y=6)]
    config.exclude_points = [(7, 8)]
    assert Config.from_dict(config.to_dict()) == config
    assert Config.from_json(config.to_json()) == config


def test_namelist_round_trip():
    config = Config()
    config.physics.nu_air = 1.1e-5
    config.domain.starttimestr = "2026-02-03 04:05:06"
    assert Config.from_namelist_text(config.to_namelist()) == config
    assert Config.from_namelist_text(config.to_namelist(comments=False)) == config


def test_write_files(tmp_path):
    config = Config()
    config.spectrum_points = [SpectrumPoint("p", lon=-80.5, lat=25.25)]
    config.exclude_points = [(5, 5), (6, 7)]
    path = config.write(tmp_path / "namelists" / "main.nml")
    assert path.is_file()
    assert (tmp_path / "namelists" / "spectrum.nml").read_text().splitlines() == ["LL", "-80.5 25.25 p"]
    assert (tmp_path / "namelists" / "exclude.nml").read_text() == "5 5\n6 7\n"
    assert Config.from_namelist(path) == config


def test_spectrum_points_validation():
    with pytest.raises(ConfigError):
        SpectrumPoint("p")
    with pytest.raises(ConfigError):
        SpectrumPoint("p", x=1, y=2, lon=3.0, lat=4.0)
    with pytest.raises(ConfigError):
        SpectrumPoint("two words", x=1, y=2)
    config = Config()
    config.spectrum_points = [SpectrumPoint("a", x=3, y=4), SpectrumPoint("b", lon=1.0, lat=2.0)]
    with pytest.raises(ConfigError, match="all by lon/lat"):
        config.spectrum_points_text()


@pytest.mark.parametrize("change, message", [
    ({"domain.mm": 2}, "mm and nm"),
    ({"domain.pm": 30}, "divisible by 4"),
    ({"domain.fmin": 3.0}, "fmin must be < fmax"),
    ({"domain.stoptimestr": "2025-01-01 00:00:00"}, "after starttimestr"),
    ({"domain.starttimestr": "yesterday"}, "YYYY-MM-DD"),
    ({"output.outgrid": 5}, "outgrid"),
    ({"output.xpl": 50}, "xpl and ypl"),
    ({"forcing_constant.fice_lth": 0.9}, "fice_lth"),
    ({"stokes.depths": [-1.0]}, "Stokes depth"),
])
def test_validate_errors(change, message):
    config = Config().update(change)
    with pytest.raises(ConfigError, match=message):
        config.validate()


def test_validate_ok():
    Config().validate()


def test_parse_value():
    assert parse_value("3") == 3
    assert parse_value("3.") == 3.0
    assert parse_value(".t.") is True
    assert parse_value("'a b'") == "a b"
    assert parse_value("1 2,3") == [1, 2, 3]


def test_package_exports():
    assert umwm.Config is Config
    assert isinstance(umwm.__version__, str)


def test_describe_and_diff(repo_root):
    config = Config()
    text = config.describe()
    assert text.startswith("&DOMAIN")
    assert "Gravitational acceleration [m/s^2]" in text
    assert "&STOKES" in config.describe("stokes") and "&DOMAIN" not in config.describe("stokes")
    with pytest.raises(ConfigError):
        config.describe("nope")
    other = Config.from_namelist(repo_root / "tests" / "namelists" / "forcing_rankine.nml")
    differences = config.diff(other)
    assert differences["domain.stoptimestr"] == ("2026-01-01 06:00:00", "2026-01-01 01:00:00")
    assert differences["forcing.winds"] == (False, True)
    assert "spectrum_points" in differences
    assert config.diff(config.copy()) == {}


def test_from_namelist_directory(tmp_path):
    workdir = tmp_path / "case"
    config = Config()
    config.domain.mm = 33
    config.write(workdir / "namelists" / "main.nml")
    assert Config.from_namelist(workdir).domain.mm == 33
    with pytest.raises(FileNotFoundError):
        Config.from_namelist(tmp_path)


@pytest.mark.parametrize("change, ok", [
    ({"output.xpl": 1}, False), ({"output.xpl": 21}, False), ({"output.ypl": 1}, False),
    ({"output.xpl": 2}, True), ({"output.xpl": 1, "domain.isglobal": True}, True),
])
def test_screen_point_must_be_interior(change, ok):
    config = Config().update(change)
    if ok:
        config.validate()
    else:
        with pytest.raises(ConfigError, match="interior"):
            config.validate()
