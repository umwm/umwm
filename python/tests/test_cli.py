import json

import numpy as np
import pytest

import umwm
from umwm.cli import main


def test_version(capsys):
    with pytest.raises(SystemExit) as info:
        main(["--version"])
    assert info.value.code == 0
    assert umwm.__version__ in capsys.readouterr().out


def test_config_defaults_and_overrides(capsys):
    assert main(["config", "-s", "domain.mm=31", "-s", "forcing.winds=.true."]) == 0
    text = capsys.readouterr().out
    config = umwm.Config.from_namelist_text(text)
    assert config.domain.mm == 31 and config.forcing.winds is True
    assert main(["config", "--json", "-s", "domain.mm=31"]) == 0
    assert json.loads(capsys.readouterr().out)["domain"]["mm"] == 31


def test_config_check(capsys):
    assert main(["config", "--check"]) == 0
    assert "valid" in capsys.readouterr().out
    assert main(["config", "--check", "-s", "domain.mm=2"]) == 1
    assert "mm and nm" in capsys.readouterr().err
    assert main(["config", "-s", "oops"]) == 1
    assert main(["config", "missing.nml"]) == 1


def test_config_write(tmp_path, capsys):
    path = tmp_path / "namelists" / "main.nml"
    assert main(["config", "-o", str(path), "-s", "domain.nm=15"]) == 0
    assert umwm.Config.from_namelist(path).domain.nm == 15


def test_init(tmp_path, capsys):
    config_dir = tmp_path / "case"
    assert main(["init", str(config_dir), "-s", "domain.stoptimestr=2026-01-01 01:00:00"]) == 0
    assert sorted(path.name for path in config_dir.iterdir()) == ["main.nml"]
    config = umwm.Config.from_namelist(config_dir)
    assert config.domain.stoptimestr == "2026-01-01 01:00:00"
    assert main(["init", str(config_dir)]) == 1
    assert "exists" in capsys.readouterr().err
    assert main(["init", str(config_dir), "--force", "--from", str(config_dir / "main.nml")]) == 0


def test_forcing_rankine(tmp_path, capsys):
    assert main(["forcing", "rankine", "-i", str(tmp_path / "input"), "--hours", "1", "--nx", "9", "--ny", "5"]) == 0
    out = capsys.readouterr().out
    assert "umwmin_2026-01-01_01:00:00.nc" in out
    assert len(list((tmp_path / "input").glob("umwmin_*.nc"))) == 2


def test_run_and_output(library, tmp_path, capsys):
    """umwm run CONFIG -o OUTPUT writes the output and the namelists, with the
    overrides, to OUTPUT, and leaves the configuration as it is."""
    config_dir, output = tmp_path / "case", tmp_path / "runs" / "one"
    assert main(["init", str(config_dir), "-s", "domain.stoptimestr=2026-01-01 01:00:00", "-s", "output.stokes=.false."]) == 0
    namelist = (config_dir / "main.nml").read_text()
    assert main(["run", str(config_dir), "-o", str(output), "--check"]) == 0
    assert main(["run", str(config_dir), "-o", str(output), "-s", "forcing_constant.wspd0=12"]) == 0
    assert "run finished at 2026-01-01 01:00:00" in capsys.readouterr().out
    assert (config_dir / "main.nml").read_text() == namelist
    assert umwm.Config.from_namelist(output / "umwmout_main.nml").forcing_constant.wspd0 == 12.0
    assert len(umwm.output_files(output)) == 2 and (output / "umwmout.grid").is_file()

    assert main(["output", str(output), "--stats", "-v", "swh", "wspd"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].split() == ["time", "swh", "wspd"]
    assert lines[-1].startswith("2026-01-01 01:00:00")
    assert float(lines[-1].split()[-1]) == pytest.approx(12.0)

    assert main(["output", str(output), "-v", "swh"]) == 0
    assert "swh" in capsys.readouterr().out
    assert main(["output", str(output), "-v", "nope"]) == 1
    assert main(["output", str(tmp_path / "nothing")]) == 1

    # the output directory defaults to the current directory, and the
    # configuration may be given by its main namelist file
    assert main(["run", str(config_dir / "main.nml"), "-o", str(tmp_path / "runs" / "two"), "-q"]) == 0
    assert len(umwm.output_files(tmp_path / "runs" / "two")) == 2


def test_run_in_the_configuration_directory(library, tmp_path, capsys, monkeypatch):
    """The output directory may be the configuration or input directory: the
    model writes only umwmout_*, umwmspc_* and umwmrst_* files there, and
    umwm clean removes only those."""
    monkeypatch.chdir(tmp_path)
    case = umwm.get_case("seamount")
    assert main(["init", "--case", "seamount", "lake"]) == 0
    config_dir = tmp_path / "lake"
    settings = ["-s", "domain.stoptimestr=2026-01-01 01:00:00", "-q"]
    files = {path: path.read_bytes() for path in config_dir.rglob("*") if path.is_file()}
    monkeypatch.chdir(config_dir)
    assert main(["run", ".", *settings]) == 0  # in the current directory
    assert main(["run", ".", "-o", "input", *settings]) == 0  # in the input directory
    assert all(path.read_bytes() == content for path, content in files.items())
    for output in (config_dir, config_dir / "input"):
        assert len(umwm.output_files(output)) == 2 and set(umwm.spectrum_output_files(output)) == {"upwind", "top", "lee", "away"}
        assert (output / "umwmout_main.nml").is_file() and (output / "umwmout_spectrum.nml").is_file()
    # run again from the copies of the namelists, in place
    assert main(["run", "umwmout_main.nml", "-i", "input", *settings]) == 0
    capsys.readouterr()
    assert main(["clean", "."]) == 0 and main(["clean", "input"]) == 0
    assert "removed" in capsys.readouterr().out
    assert sorted(path for path in config_dir.rglob("*") if path.is_file()) == sorted(files)
    assert (config_dir / "input" / "umwm.gridtopo").read_bytes() == (case.input_dir / "umwm.gridtopo").read_bytes()


def test_run_errors(library, tmp_path, capsys):
    assert main(["run", str(tmp_path / "missing")]) == 1
    assert "configuration not found" in capsys.readouterr().err
    (tmp_path / "empty").mkdir()
    assert main(["run", str(tmp_path / "empty")]) == 1
    assert "no main.nml in the configuration directory" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        main(["run"])  # the configuration is required
    config_dir = tmp_path / "case"
    assert main(["init", str(config_dir)]) == 0
    assert main(["run", str(config_dir), "-o", str(tmp_path / "bad"), "-s", "domain.mm=2", "--check"]) == 1
    assert main(["run", str(config_dir), "-o", str(tmp_path / "nofiles"), "-s", "forcing.winds=.true.", "--check"]) == 1
    assert "missing input files" in capsys.readouterr().err


def test_info(library, capsys):
    assert main(["info"]) == 0
    out = capsys.readouterr().out
    assert "library:" in out and "mpi:" in out


def test_relaunch_arguments():
    from umwm.cli import relaunch_arguments

    assert relaunch_arguments(["run", "-n", "4", "case", "-s", "domain.mm=5"]) == ["run", "case", "-s", "domain.mm=5"]
    assert relaunch_arguments(["run", "--processes=2", "-n8", "case"]) == ["run", "case"]
    assert relaunch_arguments(["run", "--processes", "3"]) == ["run"]


def test_run_processes_without_mpi_library(library, tmp_path, capsys, monkeypatch):
    if library.mpi:
        pytest.skip("the library is built with MPI")
    for name in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_RANK"):
        monkeypatch.delenv(name, raising=False)
    assert main(["init", str(tmp_path)]) == 0
    assert main(["run", "-n", "2", str(tmp_path)]) == 1
    assert "built without MPI" in capsys.readouterr().err


def test_run_log_and_clean(library, tmp_path, capsys):
    config_dir, output = tmp_path / "case", tmp_path / "output"
    assert main(["init", str(config_dir), "-s", "domain.stoptimestr=2026-01-01 01:00:00", "-s", "output.stokes=.false."]) == 0
    log = tmp_path / "run.log"
    assert main(["run", str(config_dir), "-o", str(output), "--log", str(log)]) == 0
    out = capsys.readouterr().out
    assert "run finished" in out and "umwm: solver" not in out
    assert "umwm: solver: integrating to time" in log.read_text()
    assert main(["run", str(config_dir), "-o", str(output), "--check"]) == 0
    assert f"output to {output}" in capsys.readouterr().out
    assert main(["clean", str(output), "--keep-restart"]) == 0
    assert "removed 4 file(s)" in capsys.readouterr().out  # 2 output files, the grid, and the namelist
    assert umwm.output_files(output) == [] and list(output.iterdir()) == []
    assert main(["clean", str(tmp_path / "nothing")]) == 1


def test_config_describe_and_diff(tmp_path, capsys):
    assert main(["config", "--describe", "physics"]) == 0
    out = capsys.readouterr().out
    assert out.startswith("&PHYSICS") and "Von Karman constant" in out
    other = tmp_path / "other.nml"
    assert main(["config", "-o", str(other), "-s", "domain.mm=99"]) == 0
    capsys.readouterr()
    assert main(["config", "--diff", str(other)]) == 0
    out = capsys.readouterr().out
    assert "domain.mm" in out and "99" in out
    assert main(["config", str(other), "--diff", str(other)]) == 0
    assert "do not differ" in capsys.readouterr().out


def test_output_plot(library, tmp_path, capsys):
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    config_dir, output = tmp_path / "case", tmp_path / "output"
    assert main(["init", str(config_dir), "-s", "domain.stoptimestr=2026-01-01 01:00:00"]) == 0
    assert main(["run", str(config_dir), "-o", str(output), "--quiet"]) == 0
    image = tmp_path / "swh.png"
    assert main(["output", str(output), "--plot", "swh", "--save", str(image)]) == 0
    assert image.is_file() and image.stat().st_size > 1000
    assert main(["output", str(output), "--plot", "u_stokes", "--level", "2", "--time", "0", "--save", str(image)]) == 0
    assert main(["output", str(output), "--plot", "nope", "--save", str(image)]) == 1


def test_cases(capsys):
    assert main(["cases"]) == 0
    out = capsys.readouterr().out
    assert all(name in out for name in ("fetch-limited", "wave-tank", "titan-ontario-lacus", "mars-jezero-lake"))
    assert main(["cases", "titan"]) == 0
    out = capsys.readouterr().out
    case = umwm.get_case("titan")
    assert out.startswith("UMWM example case: titan-ontario-lacus") and f"umwm run {case.directory} -o OUTPUT" in out
    assert main(["cases", "venus"]) == 1
    assert "unknown case 'venus'" in capsys.readouterr().err


def test_config_case(tmp_path, capsys):
    assert main(["config", "--case", "jonswap", "--check"]) == 0
    assert "valid" in capsys.readouterr().out
    assert main(["config", "--case", "mars", "--diff", "defaults"]) == 0
    out = capsys.readouterr().out
    assert "physics.g" in out and "3.71" in out and "9.80665" in out
    assert main(["config", "--case", "wave", "-o", str(tmp_path / "main.nml"), "-s", "forcing_constant.wspd0=5"]) == 0
    config = umwm.Config.from_namelist(tmp_path / "main.nml")
    assert config.physics.fcut_fac == 10.0 and config.forcing_constant.wspd0 == 5.0
    assert main(["config", str(tmp_path / "main.nml"), "--case", "jonswap"]) == 1
    assert "either a namelist or --case" in capsys.readouterr().err


def test_init_case(tmp_path, capsys, monkeypatch):
    monkeypatch.chdir(tmp_path)
    case = umwm.get_case("wave-tank")
    assert main(["init", "--case", "wave"]) == 0
    assert (tmp_path / "wave-tank" / "main.nml").read_text() == case.namelist.read_text()
    assert umwm.Config.from_namelist(tmp_path / "wave-tank").spectrum_points == case.config().spectrum_points
    assert main(["init", "runs/tank", "--case", "wave", "-s", "forcing_constant.wspd0=5"]) == 0
    config = umwm.Config.from_namelist(tmp_path / "runs" / "tank")
    assert config.forcing_constant.wspd0 == 5.0 and config.physics.fcut_fac == 10.0
    assert main(["init", "--case", "wave", "--from", str(case.namelist)]) == 1


def test_run_example_case(library, tmp_path, capsys, monkeypatch):
    """An example case runs from its directory, which the run leaves as it is;
    run takes no --case."""
    monkeypatch.chdir(tmp_path)
    case = umwm.get_case("duration-limited")
    before = {path: path.read_bytes() for path in case.directory.rglob("*") if path.is_file()}
    assert main(["run", str(case.directory), "-s", "domain.stoptimestr=2026-01-01 03:00:00", "-q"]) == 0
    out = capsys.readouterr().out
    assert f"running {case.directory}" in out and "run finished at 2026-01-01 03:00:00" in out
    assert len(umwm.output_files(tmp_path)) == 4  # in the current directory
    assert set(umwm.spectrum_output_files(tmp_path)) == {"sea"}
    assert {path: path.read_bytes() for path in case.directory.rglob("*") if path.is_file()} == before
    with pytest.raises(SystemExit):
        main(["run", "--case", "duration"])


def test_case_with_input_files(library, tmp_path, capsys, monkeypatch):
    monkeypatch.chdir(tmp_path)
    case = umwm.get_case("seamount")
    assert main(["init", "--case", "seamount"]) == 0
    copied = tmp_path / "seamount" / "input" / "umwm.gridtopo"
    assert copied.is_file() and not copied.is_symlink()
    assert copied.read_bytes() == (case.input_dir / "umwm.gridtopo").read_bytes()
    assert "copied umwm.gridtopo" in capsys.readouterr().out
    settings = ["-s", "domain.stoptimestr=2026-01-01 01:00:00", "-q"]
    eddies = umwm.get_case("eddy-field")
    assert main(["run", str(eddies.directory), "-o", "runs/eddies", *settings]) == 0
    assert len(umwm.output_files(tmp_path / "runs" / "eddies")) == 2
    # another input directory
    assert main(["run", "seamount", "-o", "runs/seamount", "-i", str(case.input_dir), "--check"]) == 0
    copied.write_bytes(b"edited")
    assert main(["run", "seamount", "-o", "runs/seamount", "--check", *settings]) == 1
    assert "cannot be read" in capsys.readouterr().err
    assert main(["init", "--case", "seamount", "--force"]) == 0
    assert copied.read_bytes() == (case.input_dir / "umwm.gridtopo").read_bytes()
    assert sorted(path.name for path in case.input_dir.iterdir()) == ["umwm.gridtopo"]
