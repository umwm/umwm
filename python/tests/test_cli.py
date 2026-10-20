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
    workdir = tmp_path / "case"
    assert main(["init", str(workdir), "-s", "domain.stoptimestr=2026-01-01 01:00:00"]) == 0
    for name in ("namelists", "input", "output", "restart"):
        assert (workdir / name).is_dir()
    config = umwm.Config.from_namelist(workdir / "namelists" / "main.nml")
    assert config.domain.stoptimestr == "2026-01-01 01:00:00"
    assert main(["init", str(workdir)]) == 1
    assert "exists" in capsys.readouterr().err
    assert main(["init", str(workdir), "--force", "--from", str(workdir / "namelists" / "main.nml")]) == 0


def test_forcing_rankine(tmp_path, capsys):
    assert main(["forcing", "rankine", "-w", str(tmp_path), "--hours", "1", "--nx", "9", "--ny", "5"]) == 0
    out = capsys.readouterr().out
    assert "umwmin_2026-01-01_01:00:00.nc" in out
    assert len(list((tmp_path / "input").glob("umwmin_*.nc"))) == 2


def test_run_and_output(library, tmp_path, capsys):
    workdir = tmp_path / "case"
    assert main(["init", str(workdir), "-s", "domain.stoptimestr=2026-01-01 01:00:00", "-s", "output.stokes=.false."]) == 0
    assert main(["run", "-w", str(workdir), "--check"]) == 0
    assert main(["run", "-w", str(workdir), "-s", "forcing_constant.wspd0=12"]) == 0
    assert "run finished at 2026-01-01 01:00:00" in capsys.readouterr().out
    assert umwm.Config.from_namelist(workdir / "namelists" / "main.nml").forcing_constant.wspd0 == 12.0
    assert len(umwm.output_files(workdir)) == 2

    assert main(["output", str(workdir), "--stats", "-v", "swh", "wspd"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].split() == ["time", "swh", "wspd"]
    assert lines[-1].startswith("2026-01-01 01:00:00")
    assert float(lines[-1].split()[-1]) == pytest.approx(12.0)

    assert main(["output", str(workdir), "-v", "swh"]) == 0
    assert "swh" in capsys.readouterr().out
    assert main(["output", str(workdir), "-v", "nope"]) == 1
    assert main(["output", str(tmp_path / "nothing")]) == 1


def test_run_errors(library, tmp_path, capsys):
    assert main(["run", "-w", str(tmp_path / "empty")]) == 1
    assert "umwm init" in capsys.readouterr().err
    assert main(["run", "-w", str(tmp_path / "bad"), "--defaults", "-s", "domain.mm=2", "--check"]) == 1
    assert main(["run", "-w", str(tmp_path / "nofiles"), "--defaults", "-s", "forcing.winds=.true.", "--check"]) == 1
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
    assert main(["run", "-n", "2", str(tmp_path), "--defaults"]) == 1
    assert "built without MPI" in capsys.readouterr().err


def test_run_directory_target_log_and_clean(library, tmp_path, capsys):
    workdir = tmp_path / "case"
    assert main(["init", str(workdir), "-s", "domain.stoptimestr=2026-01-01 01:00:00", "-s", "output.stokes=.false."]) == 0
    log = tmp_path / "run.log"
    assert main(["run", str(workdir), "--log", str(log)]) == 0
    out = capsys.readouterr().out
    assert "run finished" in out and "umwm: solver" not in out
    assert "umwm: solver: integrating to time" in log.read_text()
    assert main(["run", str(workdir), "--check"]) == 0
    assert "UMWM domain in" in capsys.readouterr().out
    assert main(["clean", str(workdir), "--keep-restart"]) == 0
    assert "removed 3 file(s)" in capsys.readouterr().out
    assert umwm.output_files(workdir) == []
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
    workdir = tmp_path / "case"
    assert main(["init", str(workdir), "-s", "domain.stoptimestr=2026-01-01 01:00:00"]) == 0
    assert main(["run", str(workdir), "--quiet"]) == 0
    image = tmp_path / "swh.png"
    assert main(["output", str(workdir), "--plot", "swh", "--save", str(image)]) == 0
    assert image.is_file() and image.stat().st_size > 1000
    assert main(["output", str(workdir), "--plot", "u_stokes", "--level", "2", "--time", "0", "--save", str(image)]) == 0
    assert main(["output", str(workdir), "--plot", "nope", "--save", str(image)]) == 1
