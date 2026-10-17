"""The ``umwm`` command-line tool.

Subcommands::

    umwm init [WORKDIR]            create a run directory with default namelists
    umwm config [NAMELIST]         show, check, or write a configuration
    umwm run [NAMELIST]            run the model in a run directory
    umwm forcing rankine           write idealized hurricane forcing files
    umwm output [SOURCE]           summarize the output files
    umwm info                      show the library and package information

Each subcommand accepts ``--help``. ``umwm run`` works per process under
``mpiexec`` when the library is built with MPI.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import Config, ConfigError, format_time
from .namelist import NamelistError

__all__ = ["main", "build_parser"]


def _add_set_option(parser):
    parser.add_argument(
        "-s", "--set", metavar="GROUP.NAME=VALUE", action="append", default=[], dest="settings",
        help="override a namelist variable, e.g. -s domain.stoptimestr='2026-01-01 12:00:00'",
    )


def _apply_settings(config: Config, settings) -> Config:
    for setting in settings:
        if "=" not in setting:
            raise ConfigError(f"expected GROUP.NAME=VALUE, got {setting!r}")
        key, value = setting.split("=", 1)
        config.set(key.strip(), value.strip())
    return config


def _load_config(namelist, settings, defaults_if_missing=False) -> Config:
    if namelist is None:
        config = Config()
    else:
        path = Path(namelist)
        if not path.is_file():
            if defaults_if_missing:
                config = Config()
            else:
                raise FileNotFoundError(f"namelist file not found: {path}")
        else:
            config = Config.from_namelist(path)
    return _apply_settings(config, settings)


# ---- subcommands -------------------------------------------------------------

def command_init(args) -> int:
    workdir = Path(args.workdir).resolve()
    namelist = workdir / "namelists" / "main.nml"
    if namelist.exists() and not args.force:
        raise FileExistsError(f"{namelist} exists; use --force to overwrite it")
    config = _load_config(args.from_namelist, args.settings)
    for name in ("namelists", "input", "output", "restart"):
        (workdir / name).mkdir(parents=True, exist_ok=True)
    config.write(namelist)
    print(f"umwm: wrote {namelist}")
    print(f"umwm: run directory {workdir} is ready; put input files in {workdir / 'input'} and run 'umwm run -w {workdir}'")
    return 0


def command_config(args) -> int:
    config = _load_config(args.namelist, args.settings)
    if args.check:
        config.validate()
        print("umwm: configuration is valid")
    text = config.to_json() + "\n" if args.json else config.to_namelist()
    if args.output:
        path = Path(args.output)
        if args.json:
            path.write_text(text, encoding="utf-8")
        else:
            config.write(path)
        print(f"umwm: wrote {path}")
    elif not args.check or args.show:
        sys.stdout.write(text)
    return 0


def command_run(args) -> int:
    from .domain import Domain

    workdir = Path(args.workdir).resolve()
    namelist = Path(args.namelist) if args.namelist else workdir / "namelists" / "main.nml"
    if not namelist.is_file():
        if not args.defaults:
            raise FileNotFoundError(
                f"namelist file not found: {namelist}; create one with 'umwm init', "
                "or pass --defaults to run the default configuration")
        config = Config()
    else:
        config = Config.from_namelist(namelist)
    _apply_settings(config, args.settings)
    if args.check:
        config.validate()
        domain = Domain(config, workdir, input_dir=args.input_dir)
        domain.check_inputs()
        if domain.rank == 0:
            print("umwm: configuration and input files are in place")
        return 0
    domain = Domain(config, workdir, input_dir=args.input_dir)
    domain.run(until=args.until)
    if domain.rank == 0:
        print(f"umwm: run finished at {domain.time:%Y-%m-%d %H:%M:%S}; output in {domain.output_dir}")
    domain.finalize()
    return 0


def command_forcing_rankine(args) -> int:
    from .forcing import write_rankine_forcing

    if args.output_dir:
        directory = Path(args.output_dir)
    else:
        directory = Path(args.workdir) / "input"
    paths = write_rankine_forcing(
        directory, args.start_time, args.hours, nx=args.nx, ny=args.ny, dx=args.dx, dy=args.dy,
        vmax=args.vmax, rmw=args.rmw, center_x=args.center_x, center_y=args.center_y,
        translation_x=args.translation_x, translation_y=args.translation_y, clockwise=args.clockwise,
        uc=args.uc, vc=args.vc, rhoa=args.rhoa, rhow=args.rhow, center_lon=args.center_lon,
        center_lat=args.center_lat,
    )
    for path in paths:
        print(f"umwm: wrote {path}")
    print("umwm: use with forcing.winds, currents, air_density and water_density set to .true.")
    return 0


def command_output(args) -> int:
    from . import output as _output

    source = args.source or "."
    if args.spectrum:
        for name, dataset in _output.open_spectrum_output(source).items():
            print(f"=== spectrum point {name} ===")
            print(dataset)
        return 0
    dataset = _output.open_output(source)
    if args.variables:
        missing = [name for name in args.variables if name not in dataset]
        if missing:
            raise KeyError(f"variables not in the output: {', '.join(missing)}")
    if args.stats:
        import numpy as np

        names = args.variables or ["swh", "mwp", "mwd", "cd", "ust"]
        names = [name for name in names if name in dataset]
        print(f"{'time':20s}" + "".join(f"{name:>36s}" for name in names))
        print(f"{'':20s}" + "".join(f"{'min':>12s}{'mean':>12s}{'max':>12s}" for _ in names))
        for index, time in enumerate(dataset["time"].values):
            label = str(np.datetime_as_string(time, unit="s")).replace("T", " ")
            row = f"{label:20s}"
            for name in names:
                values = dataset[name].isel(time=index).values
                row += f"{np.nanmin(values):12.5g}{np.nanmean(values):12.5g}{np.nanmax(values):12.5g}"
            print(row)
        return 0
    print(dataset[args.variables] if args.variables else dataset)
    return 0


def command_info(args) -> int:
    import numpy
    import xarray

    from ._library import LibraryNotFoundError, find_library, library

    print(f"umwm {__version__}")
    print(f"python {sys.version.split()[0]} at {sys.executable}")
    print(f"numpy {numpy.__version__}, xarray {xarray.__version__}")
    try:
        path = find_library()
    except LibraryNotFoundError as error:
        print(f"library: not found ({error})")
        return 1
    lib = library(path)
    print(f"library: {path}")
    print(f"mpi: {'yes' if lib.mpi else 'no'}")
    lib.env_init()
    print(f"denormals flushed to zero: {'yes' if lib.denormals_flushed else 'no'}")
    return 0


# ---- parser ------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="umwm",
        description="University of Miami Wave Model: configure and run the model, and inspect its output.",
    )
    parser.add_argument("--version", action="version", version=f"umwm {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    subparsers.required = True

    init = subparsers.add_parser("init", help="create a run directory with the default namelists",
                                 description="Create a run directory with the namelists, input, output and "
                                             "restart subdirectories and a main.nml namelist.")
    init.add_argument("workdir", nargs="?", default=".", help="the run directory (default: current directory)")
    init.add_argument("--from", dest="from_namelist", metavar="NAMELIST", help="start from this namelist file instead of the defaults")
    _add_set_option(init)
    init.add_argument("-f", "--force", action="store_true", help="overwrite an existing main.nml")
    init.set_defaults(func=command_init)

    config = subparsers.add_parser("config", help="show, check, or write a configuration",
                                   description="Print the configuration as a namelist (or JSON), after applying "
                                               "overrides; without a namelist file, the defaults are shown.")
    config.add_argument("namelist", nargs="?", help="a main.nml file (default: the built-in defaults)")
    _add_set_option(config)
    config.add_argument("-o", "--output", metavar="FILE", help="write the configuration to this file instead of printing it")
    config.add_argument("--json", action="store_true", help="use JSON instead of the namelist format")
    config.add_argument("--check", action="store_true", help="validate the configuration")
    config.add_argument("--show", action="store_true", help="with --check, also print the configuration")
    config.set_defaults(func=command_config)

    run = subparsers.add_parser("run", help="run the model",
                                description="Run the model in a run directory, like the umwm executable, "
                                            "reading namelists/main.nml there unless a namelist is given. "
                                            "Under mpiexec, each process runs its tile.")
    run.add_argument("namelist", nargs="?", help="the main.nml file (default: WORKDIR/namelists/main.nml)")
    run.add_argument("-w", "--workdir", default=".", help="the run directory with input, output and restart (default: current directory)")
    run.add_argument("-i", "--input-dir", metavar="DIR", help="link this directory as WORKDIR/input if that does not exist")
    _add_set_option(run)
    run.add_argument("--until", metavar="TIME", help="stop at this time (YYYY-MM-DD HH:MM:SS) instead of the configured stop time")
    run.add_argument("--defaults", action="store_true", help="use the default configuration if there is no namelist file")
    run.add_argument("--check", action="store_true", help="only validate the configuration and check the input files")
    run.set_defaults(func=command_run)

    forcing = subparsers.add_parser("forcing", help="write input forcing files")
    forcing_commands = forcing.add_subparsers(dest="forcing_command", metavar="KIND")
    forcing_commands.required = True
    rankine = forcing_commands.add_parser("rankine", help="idealized, translating Rankine hurricane vortex",
                                          description="Write hourly forcing files for a translating Rankine vortex on a "
                                                      "Cartesian grid centered at the origin.")
    rankine.add_argument("-w", "--workdir", default=".", help="the run directory; files go to its input subdirectory")
    rankine.add_argument("-o", "--output-dir", metavar="DIR", help="write the files to this directory instead")
    rankine.add_argument("--start-time", default="2026-01-01 00:00:00", help="time of the first file (default: %(default)s)")
    rankine.add_argument("--hours", type=int, default=6, help="last hourly file to write, inclusive (default: %(default)s)")
    rankine.add_argument("--nx", type=int, default=51, help="grid cells in x (default: %(default)s)")
    rankine.add_argument("--ny", type=int, default=21, help="grid cells in y (default: %(default)s)")
    rankine.add_argument("--dx", type=float, default=10000.0, help="grid spacing in x in m (default: %(default)s)")
    rankine.add_argument("--dy", type=float, default=10000.0, help="grid spacing in y in m (default: %(default)s)")
    rankine.add_argument("--vmax", type=float, default=30.0, help="maximum wind speed in m/s (default: %(default)s)")
    rankine.add_argument("--rmw", type=float, default=50000.0, help="radius of maximum wind in m (default: %(default)s)")
    rankine.add_argument("--center-x", type=float, default=-150000.0, help="initial vortex center x in m (default: %(default)s)")
    rankine.add_argument("--center-y", type=float, default=0.0, help="initial vortex center y in m (default: %(default)s)")
    rankine.add_argument("--translation-x", type=float, default=50000.0, help="translation in x per hour in m (default: %(default)s)")
    rankine.add_argument("--translation-y", type=float, default=0.0, help="translation in y per hour in m (default: %(default)s)")
    rankine.add_argument("--clockwise", action="store_true", help="rotate clockwise instead of counter-clockwise")
    rankine.add_argument("--uc", type=float, default=0.0, help="uniform current x-component in m/s (default: %(default)s)")
    rankine.add_argument("--vc", type=float, default=0.0, help="uniform current y-component in m/s (default: %(default)s)")
    rankine.add_argument("--rhoa", type=float, default=1.2, help="air density in kg/m^3 (default: %(default)s)")
    rankine.add_argument("--rhow", type=float, default=1030.0, help="water density in kg/m^3 (default: %(default)s)")
    rankine.add_argument("--center-lon", type=float, default=0.0, help="longitude of the grid center (default: %(default)s)")
    rankine.add_argument("--center-lat", type=float, default=0.0, help="latitude of the grid center (default: %(default)s)")
    rankine.set_defaults(func=command_forcing_rankine)

    output = subparsers.add_parser("output", help="summarize the output files",
                                   description="Open the gridded output of a run directory, an output directory, "
                                               "a file, or a glob with xarray and print it.")
    output.add_argument("source", nargs="?", help="run directory, output directory, file, or glob (default: current directory)")
    output.add_argument("-v", "--variables", nargs="+", metavar="VAR", help="only these variables")
    output.add_argument("--stats", action="store_true", help="print the min, mean and max of the variables at each time")
    output.add_argument("--spectrum", action="store_true", help="show the spectrum output instead")
    output.set_defaults(func=command_output)

    info = subparsers.add_parser("info", help="show the library and package information")
    info.set_defaults(func=command_info)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args) or 0)
    except (ConfigError, NamelistError, FileNotFoundError, FileExistsError, KeyError, ValueError, ImportError, RuntimeError) as error:
        message = error.args[0] if isinstance(error, KeyError) and error.args else error
        print(f"umwm: error: {message}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
