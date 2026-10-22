# UMWM Python interface

The `umwm` Python package configures and runs the University of Miami Wave
Model from Python and gives access to the model state and output with NumPy
and xarray, and the `umwm` command-line tool wraps the package. The time
integration runs in the compiled Fortran library (`libumwm.so`, built from the
same sources and with the same compiler flags as the `umwm` executable); Python
drives it through a small C API (`src/umwm_capi.f90`) and views the model
memory without copying, so there is no performance cost compared to the
executable.

## Installation

The package is built from the UMWM repository and needs the same system
dependencies as the model: a Fortran compiler, netCDF-Fortran (`nf-config` on
the path, or see below), `make`, and optionally MPI.

```
pip install .
```

The build runs `make shared`, which compiles the Fortran sources with `-fPIC`
into `libumwm.so` and places the library in the package. The build
configuration of the Makefile applies, with the variables passed as
`UMWM_`-prefixed environment variables:

```
UMWM_MPI=yes UMWM_FC=mpifort pip install .            # MPI-enabled library
UMWM_FC=ifx UMWM_FCFLAGS="-O3 -xHost" pip install .   # another compiler
UMWM_NF_CONFIG=/path/to/nf-config pip install .      # a specific netCDF
UMWM_NETCDF=/prefix pip install .                     # or its install prefix
```

Alternatively, build the library yourself and point the package to it, at
build time (the library is then copied into the package) or at run time:

```
make shared MPI=yes                         # builds src/shared/libumwm.so
UMWM_LIBRARY=src/shared/libumwm.so pip install .
# or, with any installation of the package:
export UMWM_LIBRARY=/path/to/libumwm.so
```

`pip install -e .` installs the package in editable mode, with the library
in `python/umwm/`. `umwm info` shows which library is in use and whether it
was built with MPI.

### With uv

[uv](https://docs.astral.sh/uv/) manages the environment, the dependencies,
and the build from the repository; `uv.lock` pins the dependencies.

```
uv sync                                   # .venv with umwm, its dependencies, pytest and matplotlib
uv run umwm info
uv run umwm init runs/case1
uv run umwm run runs/case1
uv run pytest                             # the Python test suite
uv run python python/examples/rankine.py
uv build                                  # sdist and wheel in dist/
uv tool install .                         # the umwm command outside the project
```

`uv sync` builds the library like `pip install .` does and installs the
package in editable mode. The build options are the same `UMWM_` variables,
and uv rebuilds the library when they change, and when the Fortran sources
or the build files change, so after editing a `.f90` file `uv sync` is all
that is needed:

```
UMWM_MPI=yes UMWM_FC=mpifort uv sync      # MPI-enabled library
uv run umwm run -n 4 runs/case1           # or: uv run mpiexec -n 4 umwm run runs/case1
uv run mpiexec -n 2 python python/tests/mpi_smoke.py
uv sync                                   # back to the serial library
```

The scripts in `tools/` and the Fortran test suite use Python too; with
uv's environment, run them as `uv run python script.py` and
`make test PYTHON="uv run python"` (or `PYTHON="uv run python" fpm test
--runner ./tests/fpm-mpi-test-runner.sh`).

## Running the model

A `Domain` wraps one model domain in a run directory with the layout the
model expects: `namelists/`, `input/`, `output/` and `restart/`. The
configuration is written to `namelists/main.nml` (and `spectrum.nml` and
`exclude.nml` when points are configured), and the model reads its input
files from `input/` and writes to `output/` and `restart/`, exactly as the
`umwm` executable does when run in that directory.

```python
import umwm

config = umwm.Config()                       # the defaults of namelists/main.nml
config.domain.mm, config.domain.nm = 51, 21  # namelist variables, by group
config.domain.stoptimestr = "2026-01-01 12:00:00"
config.forcing_constant.wspd0 = 15.0
config.set("physics.sds_fac", 40.0)          # or by dotted key

with umwm.Domain(config, workdir="runs/case1") as domain:
    domain.run()                             # to the stop time, in Fortran
    swh = domain.grid_field("swh")           # NumPy (ny, nx), NaN on land
    output = domain.output()                 # xarray.Dataset of output/umwmout_*.nc
```

`umwm.run(config, workdir)` does the same and returns the domain. A
configuration can also be a path to an existing `main.nml`, or a run
directory (`umwm.Domain("runs/case1")` reads its `namelists/main.nml` and
runs there); a file the configuration was read from is left untouched
unless the configuration is changed, so its comments are preserved.

The model prints its progress to the standard output, through the Fortran
runtime, so Python's `sys.stdout` redirection does not catch it (and in a
Jupyter kernel it goes to the kernel's terminal). `Domain(..., log="run.log")`
appends it to a file instead, and `Domain(..., quiet=True)` discards it;
in parallel runs the processes other than the root log to files named with
their rank. `domain.summary()` describes the run in a few lines.

The model advances in global time steps (`domain.dtg` seconds, the forcing
and output interval), each an inner loop of dynamic source time steps in
Fortran. The steps can be driven one at a time:

```python
domain = umwm.Domain(config, "runs/case1")
domain.initialize()                          # grids, forcing, initial state
for time in domain.iterate():                # one global time step each
    print(time, domain.grid_field("swh").max())
domain.run(until="2026-01-01 06:00:00")      # or to a given time
domain.finalize()
```

`Domain.initialize` validates the configuration (`Config.validate`, the
checks of the model's own validation) and checks the input files the
configuration needs: that they exist (`Domain.missing_inputs`), and that
they hold the required variables with the shape of the grid
(`Domain.input_problems`), raising `ConfigError`, `FileNotFoundError` or
`ValueError` before the model starts; the model's own stop on an invalid
configuration is turned into a `ConfigError` as well.

Output is written at the configured intervals, and on demand at the
current model time: `domain.write_output()` (gridded), `domain.write_spectrum_output()`
and `domain.write_restart()`, so a run can be saved at any global time
step; `domain.clean()` removes the files the model wrote to the output and
restart directories.

### Input files

Forcing is read, per global time step, from `input/umwmin_<time>.nc` when
`config.forcing.winds` (wind components `uw`, `vw`), `currents` (`uc`,
`vc`), `air_density` (`rhoa`), `water_density` (`rhow`) or `seaice` (`fice`)
is set. `write_forcing` writes such a file from arrays of shape `(ny, nx)`,
and `write_forcing_dataset` one file per time of an xarray dataset with
dimensions `(time, y, x)`, for example from a reanalysis interpolated to the
model grid. With `config.forcing.steady`, the forcing is read once, from
`input/umwmin.nc`, and held constant; `write_forcing(..., steady=True)`
writes that file, and `umwm.forcing_files(config)` lists the files the model
reads:

```python
for time in config.step_times():
    domain.write_forcing(time, uw=uw, vw=vw)          # into domain.input_dir
umwm.write_forcing("runs/case1/input", time, uw=uw, vw=vw, lon=lon, lat=lat)
umwm.write_forcing_dataset("runs/case1/input", winds, rename={"u10": "uw", "v10": "vw"})
```

With `config.grid.gridfromfile` the longitudes and latitudes are read from
`input/umwm.grid` (or `umwm.gridtopo`), and with `config.grid.topofromfile`
the water depth from `input/umwm.gridtopo`; `write_grid(path, lon, lat,
depth=None)` writes them. `Domain(config, workdir, input_dir=...)` links an
existing directory of input files as `workdir/input`.

`write_rankine_forcing` writes hourly files for an idealized, translating
hurricane vortex, the test case of the model; see `python/examples/rankine.py`.

### Example cases

The [example cases](../examples/README.md) of the repository are idealized
configurations with constant forcing: fetch-limited and duration-limited
growth, JONSWAP, a laboratory wave tank, lakes on Titan and Mars, and waves
over a seamount and through an eddy field, which bring their input files
(a bathymetry and steady currents) in their `input/` directory. The package
includes them. Like the `umwm` executable, it reads instead the `examples`
directory of the current directory, such as the repository, if there is
one, or the copy of the `examples` directory that `UMWM_EXAMPLES` points to.

```python
umwm.list_cases()                            # [Case(name="fetch-limited", ...), ...]
case = umwm.get_case("titan")                # by name, or a unique prefix of it
print(case.description, case.notes())        # one line, and the comments of its namelist
config = umwm.Config.from_case("jonswap")    # with its spectrum output points
domain = umwm.run(config, "runs/jonswap")
domain = umwm.run(case.namelist, "runs/titan")   # keeps the comments in runs/titan/namelists
with umwm.get_case("seamount").domain("runs/seamount") as domain:   # with a copy of its input files
    domain.run()
```

### Model state

Every field of the domain is accessible by name (see `umwm.FIELDS` and the
module `umwm.fields`); the names are those of the output files where they
exist (`swh`, `mwp`, `cd`, `ust`, `taux_form`, `u_stokes`, `F`, `Sin`,
...), with the Fortran component names as aliases.

- `domain.field(name)` is a NumPy array in the Fortran index order: cell
  fields over the cells of this process, spectral fields with the frequency
  and direction leading, such as the spectrum `F` of shape `(om, pm, cells)`,
  grid fields `(mm, nm)`. `domain.view(name)` (or `field(name, copy=False)`)
  views the Fortran memory without copying: it reflects the model as it
  advances, assignments change the model state (for example the wind speed
  `wspd` of a constant-forcing run), and it must not be used after
  `finalize`. `domain.bounds(name)` gives the Fortran lower bounds, such as
  the first cell index.
- `domain.grid_field(name, level=0)` is a field on the grid, shape
  `(ny, nx)` with NaN on land, like a variable of the gridded output; cell
  fields are gathered from all processes to the root process. `level`
  selects the depth of the Stokes drift components.
- `domain.state(variables=None, ocean_stress=False)` is an
  `xarray.Dataset` of the gridded state at the current time, like a gridded
  output file (dimensions `time, y, x`, Stokes drift with `z`), with the
  descriptions and units of the fields. The ocean-side fluxes are
  evaluated by the model at gridded output times only; `ocean_stress=True`
  updates them first.
- `domain.spectrum(name="F")` is an `xarray.DataArray` of the spectrum (or a
  spectral source term) on the grid, dimensions
  `(y, x, direction, frequency)`; `domain.spectrum_at(ix, iy)` is the
  spectrum at one grid point (0-based indices), `(direction, frequency)`.
- `domain.plot(name)` draws a field on the grid with matplotlib
  (`pip install umwm[plot]`); `umwm.plot_field(dataset, name)` does the same
  for an output dataset.
- `domain.remap(values)` remaps any array indexed by cell to the grid.
- `domain.sizes` holds the grid and spectrum sizes, `domain.shape` the grid
  shape `(ny, nx)`, `domain.time` the model time, `domain.finished` whether
  the stop time is reached, and `domain.attribute("dts")` the current source
  time step (also `sumt`, `dtamin`, `dtr`, `dth`, `dlnf`).

### Output files

The output files are opened with xarray and the netCDF4 engine, from a
domain or from paths:

```python
output = domain.output()                     # gridded: (time, y, x), Stokes (time, z, y, x)
spectra = domain.spectrum_output()           # {point name: Dataset (time, direction, frequency)}
restarts = domain.restart_files()            # paths; umwm.open_restart(path)
grid = domain.grid()                         # output/umwmout.grid

umwm.open_output("runs/case1")               # a run or output directory, a file, a glob, or paths
umwm.open_spectrum_output("runs/case1/output")
```

Gridded output files, one per output time, are concatenated along `time`,
lazily when [dask](https://www.dask.org/) is installed (`pip install
"umwm[lazy]"`), otherwise in memory.

### Parallel runs

With a library built with MPI, scripts run under `mpiexec`, each process
holding a tile of the cells (`domain.sizes["istart"]` to `"iend"`):

```
mpiexec -n 4 python script.py
mpiexec -n 4 umwm run -w runs/case1
```

`Domain` initialization, `advance`, `run`, `grid_field` and `state` are
collective and must be called by all processes; `field` and `view` are local
to the tile. Gathered fields are complete on the root process only; other
processes receive NaN fields. `domain.rank`, `domain.size` and
`domain.barrier()` (or `umwm.rank()`, `umwm.size()`, `umwm.barrier()`) help
writing input files from the root process only, as `Domain.write_forcing`
does. MPI is initialized on first use unless it already is, so the package
works alongside mpi4py.

### Performance

The model loop runs in Fortran, compiled from the same sources with the
same optimization flags as the executable, and Python only calls into it
once per global time step, so a run from Python takes the same time as the
executable and produces the same output. One detail is handled by the
library: with fast-math optimization (the default flags) the compiler
flushes denormal floating-point results to zero at the startup of an
executable but not of a shared library, and the wave spectrum underflows
in much of its range, where denormal arithmetic is many times slower; the
library therefore enables flushing to zero when it initializes, for the
whole process. `umwm info` reports it, and `umwm.flush_denormals(False)`
restores gradual underflow.

## Configuration

`Config` mirrors `main.nml`: one dataclass per namelist group, with the
namelist variable names, and the defaults of `namelists/main.nml`.

| Group | Class | Contents |
| --- | --- | --- |
| `config.domain` | `DomainConfig` | `isglobal`, `mm`, `nm`, `om`, `pm`, `fmin`, `fmax`, `fprog`, `starttimestr`, `stoptimestr`, `reftimestr`, `dtg`, `restart` |
| `config.physics` | `PhysicsConfig` | `g`, `nu_air`, `nu_water`, `sfct`, `kappa`, `z`, `gustiness`, `dmin`, `explim`, `sin_fac`, `sin_diss1`, `sin_diss2`, `sds_fac`, `sds_power`, `mss_fac`, `snl_fac`, `sdt_fac`, `sbf_fac`, `sbp_fac`, `fcut_fac` |
| `config.grid` | `GridConfig` | `gridfromfile`, `delx`, `dely`, `topofromfile`, `dpt`, `fillestuaries`, `filllakes` |
| `config.forcing` | `ForcingConfig` | `winds`, `currents`, `air_density`, `water_density`, `seaice` |
| `config.forcing_constant` | `ForcingConstantConfig` | `wspd0`, `wdir0`, `uc0`, `vc0`, `rhoa0`, `rhow0`, `fice0`, `fice_lth`, `fice_uth` |
| `config.output` | `OutputConfig` | `outgrid`, `outspec`, `outrst`, `xpl`, `ypl`, `stokes` |
| `config.stokes` | `StokesConfig` | `depths` |

`config.describe()` lists the variables with their descriptions and units,
and `config.diff(other)` the variables that differ between two
configurations.

`config.spectrum_points` lists `SpectrumPoint`s for the spectrum output
(`SpectrumPoint("name", x=16, y=6)` by 1-based grid indices, or
`SpectrumPoint("name", lon=-80.5, lat=25.0)`), written to `spectrum.nml`;
`config.exclude_points` lists `(x, y)` cells of enclosed seas to exclude,
written to `exclude.nml`, for `grid.filllakes`.

```python
config = umwm.Config.from_namelist("namelists/main.nml")   # also reads spectrum.nml, exclude.nml next to it
config = umwm.Config.from_case("wave-tank")                  # an example case
config["domain.dtg"] = 1800                                  # dotted keys, with conversion from strings
config.update({"output.outgrid": 1}, forcing__winds=True)
config.start_time, config.stop_time                          # datetimes of the time strings
config.step_times()                                          # the global time steps
config.validate()                                            # raises ConfigError
config.write("runs/case1/namelists/main.nml")                # main.nml, spectrum.nml, exclude.nml
config.to_namelist(), config.to_dict(), config.to_json()
umwm.Config.from_dict(...), umwm.Config.from_json(...)
```

The namelist files are read and written by `umwm.namelist`, which handles
the namelist syntax the model uses (groups, scalars, lists, comments).

## Command-line tool

```
umwm init [WORKDIR] [-s GROUP.NAME=VALUE ...] [--from NAMELIST | --case NAME] [-f]
umwm config [NAMELIST | --case NAME] [-s GROUP.NAME=VALUE ...] [-o FILE] [--json] [--check] [--describe [GROUP]] [--diff NAMELIST]
umwm run [NAMELIST|WORKDIR] [-w WORKDIR] [-c NAME] [-i INPUT_DIR] [-s GROUP.NAME=VALUE ...] [-n N] [--until TIME] [--log FILE] [-q] [--check]
umwm cases [NAME]
umwm forcing rankine [-w WORKDIR] [--hours N] [--nx N --ny N --dx M --dy M] [--vmax V --rmw R ...]
umwm output [SOURCE] [-v VAR ...] [--stats] [--spectrum] [--plot VAR [--time T] [--level L] [--save FILE]]
umwm clean [WORKDIR] [--keep-output] [--keep-restart]
umwm info
```

A typical session: create a run directory with default namelists, adjust
them, generate or copy the input files, run, and look at the output.

```
umwm init runs/case1 -s domain.stoptimestr="2026-01-01 12:00:00" -s forcing.winds=.true.
umwm forcing rankine -w runs/case1 --hours 12 --nx 21 --ny 11
umwm run runs/case1 -s forcing.currents=.true. --log runs/case1/run.log
umwm output runs/case1 --stats -v swh mwp cd
umwm output runs/case1 --plot swh --save swh.png
umwm clean runs/case1
```

The example cases run with `--case NAME` (or a unique prefix of the name),
by default in a run directory named after the case, into which the input
files of the case are copied, and `umwm cases` lists them; `umwm init --case
NAME` creates a run directory from a case to edit, which `umwm run --case`
does not overwrite unless given `--force`:

```
umwm cases                                       # list the example cases
umwm cases titan                                 # describe titan-ontario-lacus
umwm run --case jonswap                          # in ./jonswap
umwm run --case jonswap runs/j15 -s forcing_constant.wspd0=15
umwm output jonswap --stats
umwm init --case wave-tank runs/tank             # then edit and 'umwm run runs/tank'
umwm config --case mars --diff defaults          # what the case changes
```

`umwm run` reads `namelists/main.nml` of the run directory (the current
directory by default) like the `umwm` executable does; `-s` overrides are
written back to the namelist so that the run is reproducible. `--check`
validates the configuration and the input files without running and prints
a summary of the run; `--log` and `--quiet` send the output of the model to
a file or discard it. Under `mpiexec`, `umwm run` runs the model in
parallel, and `umwm run -n 4` launches `mpiexec -n 4` itself (the launcher
command can be set with `UMWM_MPIEXEC`). `umwm config --describe` documents
the namelist variables, and `umwm config A --diff B` compares two
configurations.

## Tests

```
python -m pytest python/tests                         # serial
mpiexec -n 2 python python/tests/mpi_smoke.py         # with an MPI library
uv run pytest                                         # with uv
```

The tests include the idealized hurricane case of the Fortran test suite,
run from Python and compared with the same reference values.
