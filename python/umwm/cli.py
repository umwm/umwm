"""The ``umwm`` command-line tool.

Subcommands::

    umwm init [WORKDIR]            create a run directory with default namelists
    umwm config [NAMELIST]         show, check, describe, compare, or write a configuration
    umwm run [NAMELIST|WORKDIR]    run the model in a run directory
    umwm cases [NAME]              list the example cases, or describe one
    umwm forcing rankine           write idealized hurricane forcing files
    umwm output [SOURCE]           summarize or plot the output files
    umwm clean [WORKDIR]           remove the output and restart files of a run
    umwm info                      show the library and package information

Each subcommand accepts ``--help``. ``umwm run`` works per process under
``mpiexec`` when the library is built with MPI, and ``umwm run -n N``
launches it under ``mpiexec`` itself. ``init``, ``config`` and ``run`` take
``--case NAME`` to start from an example case instead of a namelist file.
"""

from __future__ import annotations

import argparse
import os
import shlex
import sys
from pathlib import Path

from . import __version__
from .config import Config, ConfigError, namelist_path
from .namelist import NamelistError

__all__ = ["main", "build_parser", "relaunch_arguments"]

MPI_ENVIRONMENT = ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_RANK", "MPI_LOCALNRANKS",
                   "MV2_COMM_WORLD_SIZE", "SLURM_STEP_NUM_TASKS")
"""Environment variables set by MPI launchers, used to tell that a run is already under one."""


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


def _add_case_option(parser, action):
    parser.add_argument("-c", "--case", metavar="NAME",
                        help=f"{action} the example case NAME (see 'umwm cases'); a unique prefix of the name is enough")


def _source_namelist(namelist, case):
    """The namelist file given by a path or by an example case, or None for the defaults."""
    if case is None:
        return namelist
    if namelist is not None:
        raise ValueError("give either a namelist or --case, not both")
    from .cases import get_case

    return get_case(case).namelist


def _load_config(namelist, settings, case=None) -> Config:
    namelist = _source_namelist(namelist, case)
    config = Config() if namelist is None else Config.from_namelist(namelist)
    return _apply_settings(config, settings)


def relaunch_arguments(argv) -> list[str]:
    """The command-line arguments without the ``-n``/``--processes`` option, for
    relaunching the command under ``mpiexec``."""
    result = []
    skip = False
    for argument in argv:
        if skip:
            skip = False
            continue
        if argument in ("-n", "--processes"):
            skip = True
            continue
        if argument.startswith("--processes=") or (argument.startswith("-n") and argument[2:].isdigit()):
            continue
        result.append(argument)
    return result


def _relaunch_under_mpi(processes: int, argv) -> None:
    """Replace this process by ``mpiexec -n PROCESSES python -m umwm.cli ARGS``."""
    from ._library import library

    if not library().mpi:
        raise RuntimeError("the UMWM library was built without MPI; rebuild the package with "
                           "UMWM_MPI=yes, or build the library with 'make shared MPI=yes' and set UMWM_LIBRARY")
    launcher = shlex.split(os.environ.get("UMWM_MPIEXEC", "mpiexec"))
    command = [*launcher, "-n", str(processes), sys.executable, "-m", "umwm.cli", *relaunch_arguments(argv)]
    print("umwm: " + " ".join(shlex.quote(argument) for argument in command))
    sys.stdout.flush()
    os.execvp(command[0], command)


# ---- subcommands -------------------------------------------------------------

def command_init(args) -> int:
    source = _source_namelist(args.from_namelist, args.case)
    if args.workdir is None:
        args.workdir = Path(source).parent.name if args.case else "."
    workdir = Path(args.workdir).resolve()
    namelist = workdir / "namelists" / "main.nml"
    if namelist.exists() and not args.force:
        raise FileExistsError(f"{namelist} exists; use --force to overwrite it")
    config = _load_config(source, args.settings)
    for name in ("namelists", "input", "output", "restart"):
        (workdir / name).mkdir(parents=True, exist_ok=True)
    config.write(namelist, source=source)
    print(f"umwm: wrote {namelist}")
    print(f"umwm: run directory {workdir} is ready; put input files in {workdir / 'input'} and run 'umwm run {workdir}'")
    return 0


def command_config(args) -> int:
    config = _load_config(args.namelist, args.settings, args.case)
    if args.check:
        config.validate()
        print("umwm: configuration is valid")
    if args.diff is not None:
        other = Config() if args.diff == "defaults" else Config.from_namelist(args.diff)
        differences = config.diff(other)
        this = args.namelist or args.case or "defaults"
        if not differences:
            print(f"umwm: {this} and {args.diff} do not differ")
        else:
            print(f"{'variable':<28}{this:<28}{args.diff}")
            for key, (value, other_value) in differences.items():
                print(f"{key:<28}{str(value):<28}{other_value}")
        return 0
    if args.describe is not None:
        print(config.describe(args.describe or None))
        return 0
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
    if args.processes and args.processes > 1 and not any(name in os.environ for name in MPI_ENVIRONMENT):
        _relaunch_under_mpi(args.processes, args.argv)

    from .domain import Domain

    if args.case:
        return _run_case(args)
    workdir = Path(args.workdir).resolve() if args.workdir else None
    if args.target and Path(args.target).is_dir():
        if workdir is None:
            workdir = Path(args.target).resolve()
        namelist = workdir / "namelists" / "main.nml"
    elif args.target:
        namelist = Path(args.target)
    else:
        namelist = (workdir or Path.cwd()) / "namelists" / "main.nml"
    workdir = workdir or Path.cwd()
    if not namelist.is_file():
        if not args.defaults:
            raise FileNotFoundError(
                f"namelist file not found: {namelist}; create one with 'umwm init', "
                "or pass --defaults to run the default configuration")
        config = Config()
    else:
        config = Config.from_namelist(namelist)
    _apply_settings(config, args.settings)
    domain = Domain(config, workdir, input_dir=args.input_dir, log=args.log, quiet=args.quiet)
    return _run_domain(domain, args)


def _run_case(args) -> int:
    """Run an example case in the run directory given as the target or with
    --workdir, by default a directory named after the case."""
    from .cases import get_case
    from .domain import Domain

    case = get_case(args.case)
    if args.target and args.workdir:
        raise ValueError("with --case, give the run directory either as the target or with --workdir")
    workdir = Path(args.target or args.workdir or case.name)
    if workdir.is_file():
        raise ValueError(f"with --case, the target is the run directory, not a file: {workdir}")
    domain = Domain(case.namelist, workdir.resolve(), input_dir=args.input_dir, log=args.log, quiet=args.quiet)
    _apply_settings(domain.config, args.settings)
    if domain.rank == 0 and not args.check:
        print(f"umwm: running example case {case.name} in {domain.workdir}")
    return _run_domain(domain, args)


def _run_domain(domain, args) -> int:
    if args.check:
        domain.config.validate()
        domain.check_inputs()
        if domain.rank == 0:
            print(domain.summary())
            print("umwm: configuration and input files are in place")
        return 0
    domain.run(until=args.until)
    if domain.rank == 0:
        print(f"umwm: run finished at {domain.time:%Y-%m-%d %H:%M:%S}; output in {domain.output_dir}")
    domain.finalize()
    return 0


def command_cases(args) -> int:
    from .cases import examples_dir, get_case, list_cases

    if args.name:
        case = get_case(args.name)
        print(case.notes())
        print(f"\numwm: the namelists of {case.name} are in {case.directory}")
        print(f"umwm: run it with 'umwm run --case {case.name}', or create a run directory "
              f"from it with 'umwm init --case {case.name}'")
        return 0
    cases = list_cases()
    width = max(len(case.name) for case in cases)
    for case in cases:
        print(f"{case.name:<{width}}  {case.description}")
    print(f"\numwm: example cases in {examples_dir()}; 'umwm cases NAME' describes one, and "
          "'umwm run --case NAME' runs it")
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
    if args.plot:
        from .plotting import _pyplot, plot_field

        if args.plot not in dataset:
            raise KeyError(f"variable not in the output: {args.plot}")
        grid = None
        try:
            grid = _output.open_grid(source)
        except (FileNotFoundError, OSError, ValueError):
            pass
        plot_field(dataset, args.plot, time=args.time, level=args.level, grid=grid)
        plt = _pyplot()
        if args.save:
            plt.savefig(args.save, dpi=150, bbox_inches="tight")
            print(f"umwm: wrote {args.save}")
        else:
            plt.show()
        return 0
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


def command_clean(args) -> int:
    from .domain import Domain

    workdir = Path(args.workdir).resolve()
    if not (workdir / "output").is_dir() and not (workdir / "restart").is_dir():
        raise FileNotFoundError(f"{workdir} has no output or restart directory")
    config = Config.from_namelist(workdir) if (workdir / "namelists" / "main.nml").is_file() else Config()
    domain = Domain(config, workdir)
    removed = domain.clean(output=not args.keep_output, restart=not args.keep_restart)
    print(f"umwm: removed {len(removed)} file(s) from {workdir}")
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
    init.add_argument("workdir", nargs="?",
                      help="the run directory (default: current directory, or the name of the case with --case)")
    init.add_argument("--from", dest="from_namelist", metavar="NAMELIST", help="start from this namelist file instead of the defaults")
    _add_case_option(init, "start from")
    _add_set_option(init)
    init.add_argument("-f", "--force", action="store_true", help="overwrite an existing main.nml")
    init.set_defaults(func=command_init)

    config = subparsers.add_parser("config", help="show, check, describe, compare, or write a configuration",
                                   description="Print the configuration as a namelist (or JSON), after applying "
                                               "overrides; without a namelist file, the defaults are shown.")
    config.add_argument("namelist", nargs="?", help="a main.nml file or a run directory (default: the built-in defaults)")
    _add_case_option(config, "use")
    _add_set_option(config)
    config.add_argument("-o", "--output", metavar="FILE", help="write the configuration to this file instead of printing it")
    config.add_argument("--json", action="store_true", help="use JSON instead of the namelist format")
    config.add_argument("--check", action="store_true", help="validate the configuration")
    config.add_argument("--show", action="store_true", help="with --check, also print the configuration")
    config.add_argument("--describe", nargs="?", const="", metavar="GROUP",
                        help="list the variables with their descriptions and units, optionally of one group")
    config.add_argument("--diff", metavar="NAMELIST",
                        help="show the variables that differ from this namelist file or run directory, or 'defaults'")
    config.set_defaults(func=command_config)

    run = subparsers.add_parser("run", help="run the model",
                                description="Run the model in a run directory, like the umwm executable, "
                                            "reading namelists/main.nml there unless a namelist is given. "
                                            "Under mpiexec, each process runs its tile.")
    run.add_argument("target", nargs="?", metavar="NAMELIST|WORKDIR",
                     help="the main.nml file, or the run directory (default: WORKDIR/namelists/main.nml)")
    run.add_argument("-w", "--workdir", help="the run directory with input, output and restart (default: current directory)")
    _add_case_option(run, "run")
    run.add_argument("-i", "--input-dir", metavar="DIR", help="link this directory as WORKDIR/input if that does not exist")
    _add_set_option(run)
    run.add_argument("-n", "--processes", type=int, metavar="N",
                     help="launch N processes under mpiexec (or $UMWM_MPIEXEC); the library must be built with MPI")
    run.add_argument("--until", metavar="TIME", help="stop at this time (YYYY-MM-DD HH:MM:SS) instead of the configured stop time")
    run.add_argument("--log", metavar="FILE", help="append the output of the model to this file instead of the terminal")
    run.add_argument("-q", "--quiet", action="store_true", help="discard the output of the model")
    run.add_argument("--defaults", action="store_true", help="use the default configuration if there is no namelist file")
    run.add_argument("--check", action="store_true", help="only validate the configuration and check the input files")
    run.epilog = ("With --case, the run directory is the target or WORKDIR, by default a directory named after "
                  "the case in the current directory.")
    run.set_defaults(func=command_run)

    cases = subparsers.add_parser("cases", help="list the example cases, or describe one",
                                  description="List the example cases, idealized configurations without input "
                                              "files, or describe one; run one with 'umwm run --case NAME'.")
    cases.add_argument("name", nargs="?", help="the case to describe; a unique prefix of its name is enough")
    cases.set_defaults(func=command_cases)

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

    output = subparsers.add_parser("output", help="summarize or plot the output files",
                                   description="Open the gridded output of a run directory, an output directory, "
                                               "a file, or a glob with xarray and print, summarize, or plot it.")
    output.add_argument("source", nargs="?", help="run directory, output directory, file, or glob (default: current directory)")
    output.add_argument("-v", "--variables", nargs="+", metavar="VAR", help="only these variables")
    output.add_argument("--stats", action="store_true", help="print the min, mean and max of the variables at each time")
    output.add_argument("--spectrum", action="store_true", help="show the spectrum output instead")
    output.add_argument("--plot", metavar="VAR", help="draw this variable with matplotlib (pip install umwm[plot])")
    output.add_argument("--time", type=int, default=-1, help="time index to plot (default: the last)")
    output.add_argument("--level", type=int, default=0, help="depth index of the Stokes drift to plot (default: 0)")
    output.add_argument("--save", metavar="FILE", help="save the plot to this file instead of showing it")
    output.set_defaults(func=command_output)

    clean = subparsers.add_parser("clean", help="remove the output and restart files of a run",
                                  description="Remove the files the model wrote to the output and restart "
                                              "directories of a run directory, and nothing else.")
    clean.add_argument("workdir", nargs="?", default=".", help="the run directory (default: current directory)")
    clean.add_argument("--keep-output", action="store_true", help="keep the gridded and spectrum output")
    clean.add_argument("--keep-restart", action="store_true", help="keep the restart files")
    clean.set_defaults(func=command_clean)

    info = subparsers.add_parser("info", help="show the library and package information")
    info.set_defaults(func=command_info)

    return parser


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(argv)
    args.argv = argv
    try:
        return int(args.func(args) or 0)
    except (ConfigError, NamelistError, FileNotFoundError, FileExistsError, KeyError, ValueError, ImportError,
            RuntimeError) as error:
        message = error.args[0] if isinstance(error, KeyError) and error.args else error
        print(f"umwm: error: {message}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
