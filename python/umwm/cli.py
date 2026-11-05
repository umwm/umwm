"""The ``umwm`` command-line tool.

Subcommands::

    umwm init [DIR]                 create a configuration directory
    umwm config [CONFIG]            show, check, describe, compare, or write a configuration
    umwm run CONFIG [-o OUTPUT]     run the model, like the umwm executable
    umwm cases [NAME]               list the example cases, or describe one
    umwm forcing rankine            write idealized hurricane forcing files
    umwm forcing era5 CONFIG        prepare the input files of a configuration from ERA5
    umwm output [OUTPUT]            summarize or plot the output files
    umwm clean [OUTPUT]             remove the output and restart files of a run
    umwm info                       show the library and package information

A configuration ``CONFIG`` is a directory with the ``main.nml`` namelist, the
``spectrum.nml`` points of the spectrum output and, with ``fillLakes``, the
``exclude.nml`` cells, or a main namelist file whose directory holds the
others. ``umwm run CONFIG`` reads the input files from ``CONFIG/input`` (or
``-i INPUT``) and writes the output and restart files, with copies of the
namelists, to the current directory (or ``-o OUTPUT``), like the ``umwm``
executable. Each subcommand accepts ``--help``. ``umwm run`` works per
process under ``mpiexec`` when the library is built with MPI, and ``umwm run
-n N`` launches it under ``mpiexec`` itself. ``init`` and ``config`` take
``--case NAME`` to start from an example case.
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


def _source_case(namelist, case):
    """The example case given with --case, or None."""
    if case is None:
        return None
    if namelist is not None:
        raise ValueError("give either a namelist or --case, not both")
    from .cases import get_case

    return get_case(case)


def _source_namelist(namelist, case):
    """The namelist file given by a path or by an example case, or None for the defaults."""
    case = _source_case(namelist, case)
    return namelist if case is None else case.namelist


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
    case = _source_case(args.from_namelist, args.case)
    source = args.from_namelist if case is None else case.namelist
    if args.directory is None:
        args.directory = "." if case is None else case.name
    directory = Path(args.directory).resolve()
    namelist = directory / "main.nml"
    if namelist.exists() and not args.force:
        raise FileExistsError(f"{namelist} exists; use --force to overwrite it")
    config = _load_config(source, args.settings)
    config.write(namelist, source=None if source is None else namelist_path(source))
    print(f"umwm: wrote {namelist}")
    if case is not None:
        for path in case.copy_inputs(directory / "input", overwrite=args.force):
            print(f"umwm: copied {path.name} into {path.parent}")
    if config.source_name is not None:
        print(f"umwm: configuration {directory} is ready; run it with 'umwm run {directory} -o OUTPUT', which "
              f"prepares its input files in {directory / 'input'} from {config.source_name}, or prepare them now "
              f"with 'umwm forcing {config.source_name} {directory}'")
    else:
        print(f"umwm: configuration {directory} is ready; its input files go in {directory / 'input'}; run it "
              f"with 'umwm run {directory} -o OUTPUT'")
    return 0


def command_config(args) -> int:
    config = _load_config(args.namelist, args.settings, args.case)
    if args.check:
        config.validate()
        print("umwm: configuration is valid")
    if args.diff is not None:
        other = Config() if args.diff == "defaults" and not Path(args.diff).exists() else Config.from_namelist(args.diff)
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

    domain = Domain(args.config, args.output, input_dir=args.input, log=args.log, quiet=args.quiet)
    _apply_settings(domain.config, args.settings)
    return _run_domain(domain, args)


def _run_domain(domain, args) -> int:
    if args.check:
        domain.config.validate()
        source = domain.config.source_name
        missing = domain.missing_inputs() if source is not None else []
        stale = []
        if source is not None:
            from .sources.prepare import pending_inputs

            bathymetry, files = pending_inputs(domain.config, domain.input_dir)
            pending = ([domain.input_dir / "umwm.gridtopo"] if bathymetry else []) + [path for _, path in files]
            stale = [path for path in pending if path.is_file()]
        if not missing and not stale:
            domain.check_inputs()
        if domain.rank == 0:
            print(domain.summary())
            if missing or stale:
                files = " and ".join(([f"{len(missing)} missing input file(s)"] if missing else [])
                                     + ([f"{len(stale)} input file(s) prepared for another grid"] if stale else []))
                print(f"umwm: configuration is valid; {files} will be prepared from {source} when the model "
                      f"runs, or now with 'umwm forcing {source}'")
            else:
                print("umwm: configuration and input files are in place")
        return 0
    if domain.rank == 0:
        print(f"umwm: running {domain.config_dir or 'the configuration'} with input from {domain.input_dir} and "
              f"output to {domain.output_dir}")
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
        print(f"\numwm: the configuration of {case.name} is {case.directory}")
        print(f"umwm: run it with 'umwm run {case.directory} -o OUTPUT', or copy it into a directory to edit "
              f"with 'umwm init --case {case.name} DIR'")
        return 0
    cases = list_cases()
    width = max(len(case.name) for case in cases)
    for case in cases:
        print(f"{case.name:<{width}}  {case.description}")
    print(f"\numwm: example cases in {examples_dir()}; 'umwm cases NAME' describes one, "
          f"'umwm run {examples_dir()}/NAME -o OUTPUT' runs it, and 'umwm init --case NAME DIR' copies it")
    return 0


def command_forcing_rankine(args) -> int:
    from .forcing import write_rankine_forcing

    paths = write_rankine_forcing(
        Path(args.input), args.start_time, args.hours, nx=args.nx, ny=args.ny, dx=args.dx, dy=args.dy,
        vmax=args.vmax, rmw=args.rmw, center_x=args.center_x, center_y=args.center_y,
        translation_x=args.translation_x, translation_y=args.translation_y, clockwise=args.clockwise,
        uc=args.uc, vc=args.vc, rhoa=args.rhoa, rhow=args.rhow, center_lon=args.center_lon,
        center_lat=args.center_lat,
    )
    for path in paths:
        print(f"umwm: wrote {path}")
    print("umwm: use with forcing.winds, currents, air_density and water_density set to .true.")
    return 0


def command_forcing_source(args) -> int:
    """Prepare the input files of a configuration from a data source."""
    from .sources import prepare_inputs

    namelist = namelist_path(args.config)
    config = _load_config(namelist, args.settings)
    input_dir = Path(args.input) if args.input else namelist.parent / "input"
    options = {"workers": args.workers}
    if args.url:
        options["url"] = args.url
    paths = prepare_inputs(config, input_dir, source=args.source, force=args.force, **options)
    if not paths and not config.forcing_from_file and not config.grid.topofromfile:
        print(f"umwm: {namelist} reads no input files: set forcing.winds, air_density, seaice or waves, "
              "or grid.topofromfile")
    return 0


def command_output(args) -> int:
    from . import output as _output

    source = args.output or "."
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
    from .domain import clean_output

    directory = Path(args.output).resolve()
    if not directory.is_dir():
        raise FileNotFoundError(f"output directory not found: {directory}")
    removed = clean_output(directory, output=not args.keep_output, restart=not args.keep_restart)
    print(f"umwm: removed {len(removed)} file(s) from {directory}")
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

    init = subparsers.add_parser("init", help="create a configuration directory",
                                 description="Create a configuration directory with a main.nml namelist, from "
                                             "the defaults, a namelist file, or an example case, whose "
                                             "spectrum.nml and input files are copied too.")
    init.add_argument("directory", nargs="?", metavar="DIR",
                      help="the configuration directory (default: current directory, or the name of the case "
                           "with --case)")
    init.add_argument("--from", dest="from_namelist", metavar="CONFIG",
                      help="start from this configuration directory or namelist file instead of the defaults")
    _add_case_option(init, "start from")
    _add_set_option(init)
    init.add_argument("-f", "--force", action="store_true",
                      help="overwrite an existing main.nml, and with --case, input files that differ from those of the case")
    init.set_defaults(func=command_init)

    config = subparsers.add_parser("config", help="show, check, describe, compare, or write a configuration",
                                   description="Print the configuration as a namelist (or JSON), after applying "
                                               "overrides; without a namelist file, the defaults are shown.")
    config.add_argument("namelist", nargs="?", metavar="CONFIG",
                        help="a configuration directory or main namelist file (default: the built-in defaults)")
    _add_case_option(config, "use")
    _add_set_option(config)
    config.add_argument("-o", "--output", metavar="FILE", help="write the configuration to this file instead of printing it")
    config.add_argument("--json", action="store_true", help="use JSON instead of the namelist format")
    config.add_argument("--check", action="store_true", help="validate the configuration")
    config.add_argument("--show", action="store_true", help="with --check, also print the configuration")
    config.add_argument("--describe", nargs="?", const="", metavar="GROUP",
                        help="list the variables with their descriptions and units, optionally of one group")
    config.add_argument("--diff", metavar="NAMELIST",
                        help="show the variables that differ from this configuration directory or namelist file, or "
                             "from the defaults with 'defaults' (unless a file or directory of that name exists)")
    config.set_defaults(func=command_config)

    run = subparsers.add_parser("run", help="run the model",
                                description="Run the model with the configuration CONFIG, like the umwm "
                                            "executable: read the input files from CONFIG/input, or INPUT, and "
                                            "write the output and restart files, with copies of the namelists "
                                            "(umwmout_main.nml and the others), to the current directory, or "
                                            "OUTPUT, created if needed; a restart reads its restart file from "
                                            "there too. Under mpiexec, each process runs its tile.")
    run.add_argument("config", metavar="CONFIG",
                     help="the configuration directory, with main.nml and spectrum.nml, or its main namelist file")
    run.add_argument("-o", "--output", metavar="DIR",
                     help="the directory of the output and restart files (default: current directory); it may be "
                          "the configuration or input directory")
    run.add_argument("-i", "--input", metavar="DIR", help="the directory of the input files (default: CONFIG/input)")
    _add_set_option(run)
    run.add_argument("-n", "--processes", type=int, metavar="N",
                     help="launch N processes under mpiexec (or $UMWM_MPIEXEC); the library must be built with MPI")
    run.add_argument("--until", metavar="TIME", help="stop at this time (YYYY-MM-DD HH:MM:SS) instead of the configured stop time")
    run.add_argument("--log", metavar="FILE", help="append the output of the model to this file instead of the terminal")
    run.add_argument("-q", "--quiet", action="store_true", help="discard the output of the model")
    run.add_argument("--check", action="store_true", help="only validate the configuration and check the input files")
    run.epilog = ("Overrides of -s apply to the run only: they are in the copy of the namelist in the output "
                  "directory, umwmout_main.nml, and not written to the configuration.")
    run.set_defaults(func=command_run)

    cases = subparsers.add_parser("cases", help="list the example cases, or describe one",
                                  description="List the example cases, or describe one with its configuration "
                                              "directory; run one with 'umwm run DIR -o OUTPUT', or copy it "
                                              "into a directory to edit with 'umwm init --case NAME DIR'.")
    cases.add_argument("name", nargs="?", help="the case to describe; a unique prefix of its name is enough")
    cases.set_defaults(func=command_cases)

    forcing = subparsers.add_parser("forcing", help="write input forcing files")
    forcing_commands = forcing.add_subparsers(dest="forcing_command", metavar="KIND")
    forcing_commands.required = True
    rankine = forcing_commands.add_parser("rankine", help="idealized, translating Rankine hurricane vortex",
                                          description="Write hourly forcing files for a translating Rankine vortex on a "
                                                      "Cartesian grid centered at the origin.")
    rankine.add_argument("-i", "--input", metavar="DIR", default="input",
                         help="the input directory to write the files to (default: %(default)s)")
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

    from .sources import SOURCES

    for info in SOURCES.values():
        source = forcing_commands.add_parser(
            info.name, help=f"prepare the input files of a run from {info.name}",
            description=f"Prepare the input files of a configuration that are missing from its input "
                        f"directory, or that it prepared for another grid, from "
                        f"{info.description}: those that the namelist reads (with grid.topofromfile, the "
                        "bathymetry, and the forcing of forcing.winds, air_density, seaice and waves), on its "
                        "grid (grid.lon0, lat0, dlon and dlat, or the grid file of grid.gridfromfile), at its "
                        "times. The model does it itself when forcing.source names the data source. Needs "
                        f"pip install 'umwm[{info.extra}]'.")
        source.add_argument("config", metavar="CONFIG",
                            help="the configuration directory, or its main namelist file")
        source.add_argument("-i", "--input", metavar="DIR",
                            help="write the input files to this directory (default: CONFIG/input)")
        _add_set_option(source)
        source.add_argument("-f", "--force", action="store_true", help="replace the input files that exist")
        source.add_argument("--url", help=f"the location of the dataset (default: that of the {info.name} source)")
        source.add_argument("--workers", type=int, default=16,
                            help="concurrent reads from the dataset (default: %(default)s)")
        source.set_defaults(func=command_forcing_source, source=info.name)

    output = subparsers.add_parser("output", help="summarize or plot the output files",
                                   description="Open the gridded output of an output directory, a file, or a "
                                               "glob with xarray and print, summarize, or plot it.")
    output.add_argument("output", nargs="?", metavar="OUTPUT",
                        help="output directory, file, or glob (default: current directory)")
    output.add_argument("-v", "--variables", nargs="+", metavar="VAR", help="only these variables")
    output.add_argument("--stats", action="store_true", help="print the min, mean and max of the variables at each time")
    output.add_argument("--spectrum", action="store_true", help="show the spectrum output instead")
    output.add_argument("--plot", metavar="VAR", help="draw this variable with matplotlib (pip install umwm[plot])")
    output.add_argument("--time", type=int, default=-1, help="time index to plot (default: the last)")
    output.add_argument("--level", type=int, default=0, help="depth index of the Stokes drift to plot (default: 0)")
    output.add_argument("--save", metavar="FILE", help="save the plot to this file instead of showing it")
    output.set_defaults(func=command_output)

    clean = subparsers.add_parser("clean", help="remove the output and restart files of a run",
                                  description="Remove the output and restart files the model writes to an "
                                              "output directory (umwmout_*.nc, umwmout.grid, umwmspc_*.nc and "
                                              "umwmrst_*.nc), and nothing else: the copies of the namelists, "
                                              "umwmout_main.nml and the others, are kept.")
    clean.add_argument("output", nargs="?", default=".", metavar="OUTPUT",
                       help="the output directory (default: current directory)")
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
            RuntimeError, OSError) as error:
        message = error.args[0] if isinstance(error, KeyError) and error.args else error
        print(f"umwm: error: {message}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
