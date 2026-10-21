"""The model domain: configure, run, and inspect a UMWM simulation from Python.

A :class:`Domain` wraps one ``domain_type`` instance of the Fortran model.
It manages a run directory with the layout the model expects (``namelists/``,
``input/``, ``output/`` and ``restart/``), writes the namelists from its
:class:`~umwm.Config`, and drives the model through the C API of the shared
library, so that the time integration runs entirely in compiled Fortran. The
model state is accessible as NumPy arrays viewing the Fortran memory, as
fields gathered to the horizontal grid, and as xarray datasets, and the
output files as xarray datasets.
"""

from __future__ import annotations

import contextlib
import os
from datetime import datetime
from pathlib import Path
from typing import Iterator

import numpy as np
import xarray as xr

from . import fields as _fields
from . import output as _output
from ._library import library
from ._redirect import redirect_stdout_fd
from .config import Config, ConfigError, namelist_path, parse_time
from .forcing import FORCING_VARIABLES, forcing_filename, write_forcing, write_grid

__all__ = ["Domain", "run"]

INT_ATTRIBUTES = ("mm", "nm", "im", "imm", "istart", "iend", "iistart", "iiend", "om", "pm", "lm")


class Domain:
    """A UMWM model domain in a run directory.

    Parameters
    ----------
    config
        A :class:`~umwm.Config`, the path of a ``main.nml`` namelist file, a
        run directory (whose ``namelists/main.nml`` is read, and which
        becomes the ``workdir`` unless one is given), or ``None`` for the
        default configuration.
    workdir
        The run directory; the current directory by default. It is created
        with the ``namelists``, ``input``, ``output`` and ``restart``
        subdirectories if they do not exist.
    input_dir
        An existing directory with the input files (forcing and grid files),
        linked as ``workdir/input`` if that does not exist yet.
    log
        A file to which the output of the model (and anything else written
        to the standard output while the model runs) is appended instead of
        the terminal; in parallel runs, the processes other than the root
        write to files named with their rank.
    quiet
        Discard the output of the model (unless ``log`` is given).

    The model is set up by :meth:`initialize`, which writes the namelists
    and reads the inputs, and integrated by :meth:`run` or :meth:`advance`;
    both initialize the domain if needed. :meth:`finalize` releases the
    model memory, and is also called on leaving a ``with`` block.

    >>> with Domain(config, workdir="runs/case1") as domain:
    ...     domain.run()
    ...     swh = domain.grid_field("swh")
    ...     output = domain.output()
    """

    def __init__(self, config=None, workdir=None, *, input_dir=None, log=None, quiet=False):
        self._lib = library()
        self._handle = None
        self._source_namelist = None

        if config is None:
            self.config = Config()
        elif isinstance(config, Config):
            self.config = config
        else:
            source = Path(config)
            if source.is_dir() and workdir is None:
                workdir = source
            self._source_namelist = namelist_path(source).resolve()
            self.config = Config.from_namelist(self._source_namelist)

        self.workdir = Path(workdir or Path.cwd()).resolve()
        self.namelist_dir = self.workdir / "namelists"
        self.namelist_path = self.namelist_dir / "main.nml"
        self.input_dir = self.workdir / "input"
        self.output_dir = self.workdir / "output"
        self.restart_dir = self.workdir / "restart"

        self.quiet = bool(quiet)
        self.log = None
        if log is not None:
            path = Path(log)
            if self.rank > 0:
                path = path.with_name(f"{path.stem}.rank{self.rank}{path.suffix}")
            self.log = path

        if self.rank == 0:
            for directory in (self.namelist_dir, self.output_dir, self.restart_dir):
                directory.mkdir(parents=True, exist_ok=True)
            if input_dir is not None:
                self._link_input_dir(Path(input_dir).resolve())
            self.input_dir.mkdir(parents=True, exist_ok=True)
        self.barrier()

    def _link_input_dir(self, input_dir: Path) -> None:
        if not input_dir.is_dir():
            raise FileNotFoundError(f"input directory not found: {input_dir}")
        if self.input_dir.exists() or self.input_dir.is_symlink():
            if self.input_dir.resolve() != input_dir:
                raise FileExistsError(
                    f"{self.input_dir} exists and is not {input_dir}; remove it to link the input directory")
            return
        self.input_dir.symlink_to(input_dir, target_is_directory=True)

    # ---- environment --------------------------------------------------------

    @property
    def rank(self) -> int:
        """Rank of this process, 0 in serial mode."""
        return self._lib.rank

    @property
    def size(self) -> int:
        """Number of processes, 1 in serial mode."""
        return self._lib.size

    @property
    def mpi(self) -> bool:
        """Whether the model library was built with MPI."""
        return self._lib.mpi

    def barrier(self) -> None:
        """Wait for all processes; a no-op in serial mode."""
        self._lib.barrier()

    @contextlib.contextmanager
    def _in_workdir(self):
        """Run the model with the run directory as the current directory,
        where the model expects its input, output, and namelist files."""
        previous = os.getcwd()
        os.chdir(self.workdir)
        try:
            yield
        finally:
            os.chdir(previous)

    @contextlib.contextmanager
    def _model_output(self):
        """Send the output of the model to the log file, or discard it when quiet."""
        if self.log is None and not self.quiet:
            yield
        else:
            with redirect_stdout_fd(self.log, flush=self._lib.flush):
                yield

    # ---- setup --------------------------------------------------------------

    @property
    def initialized(self) -> bool:
        return self._handle is not None

    def _require_handle(self) -> int:
        if self._handle is None:
            raise RuntimeError("the domain is not initialized; call initialize() first")
        return self._handle

    def write_namelists(self) -> Path:
        """Write ``namelists/main.nml``, and ``spectrum.nml`` and ``exclude.nml``
        when points are configured, from the configuration.

        A namelist file the configuration was read from is left untouched
        when the configuration still matches it, so that its comments are
        preserved, and is copied with them into another run directory. Only
        the root process writes; all processes wait for it.
        """
        if self.rank == 0:
            if not self._matches_file(self.namelist_path, lambda text: Config.from_namelist_text(text) == self._bare_config()):
                self.namelist_path.write_text(self.config.namelist_text(self._source_namelist), encoding="utf-8")
            if self.config.spectrum_points:
                path = self.namelist_dir / "spectrum.nml"
                text = self.config.spectrum_points_text()
                if not (path.is_file() and path.read_text(encoding="utf-8").split() == text.split()):
                    path.write_text(text, encoding="utf-8")
            if self.config.exclude_points:
                path = self.namelist_dir / "exclude.nml"
                text = self.config.exclude_points_text()
                if not (path.is_file() and path.read_text(encoding="utf-8").split() == text.split()):
                    path.write_text(text, encoding="utf-8")
        self.barrier()
        return self.namelist_path

    def _bare_config(self) -> Config:
        """The configuration without the points held in separate files."""
        bare = self.config.copy()
        bare.spectrum_points = []
        bare.exclude_points = []
        return bare

    @staticmethod
    def _matches_file(path: Path, matches) -> bool:
        if not path.is_file():
            return False
        try:
            return bool(matches(path.read_text(encoding="utf-8")))
        except Exception:
            return False

    def _required_forcing_variables(self) -> list[str]:
        forcing = self.config.forcing
        names = []
        if forcing.winds:
            names += ["uw", "vw"]
        if forcing.currents:
            names += ["uc", "vc"]
        if forcing.seaice:
            names.append("fice")
        if forcing.air_density:
            names.append("rhoa")
        if forcing.water_density:
            names.append("rhow")
        return names

    def missing_inputs(self) -> list[Path]:
        """The input files the configuration requires that do not exist."""
        config = self.config
        missing = []
        gridtopo = self.input_dir / "umwm.gridtopo"
        if config.grid.topofromfile and not gridtopo.is_file():
            missing.append(gridtopo)
        elif config.grid.gridfromfile and not (gridtopo.is_file() or (self.input_dir / "umwm.grid").is_file()):
            missing.append(self.input_dir / "umwm.grid")
        if config.grid.filllakes and not (self.namelist_dir / "exclude.nml").is_file() and not config.exclude_points:
            missing.append(self.namelist_dir / "exclude.nml")
        if config.forcing_from_file:
            for time in config.forcing_times():
                path = self.input_dir / forcing_filename(time)
                if not path.is_file():
                    missing.append(path)
        if config.domain.restart:
            path = self.restart_dir / f"umwmrst_{config.start_time:%Y-%m-%d_%H:%M:%S}.nc"
            if not path.is_file():
                missing.append(path)
        return missing

    def input_problems(self) -> list[str]:
        """Problems with the input files the configuration requires: files that
        are missing, and files without the required variables or with a
        shape that does not match the grid (``mm`` by ``nm``)."""
        from netCDF4 import Dataset

        config = self.config
        mm, nm = config.domain.mm, config.domain.nm
        problems = [f"missing {path}" for path in self.missing_inputs()]

        def check(path: Path, variables: list[str], minimum_rank: int) -> None:
            if not path.is_file():
                return
            try:
                with Dataset(path) as dataset:
                    for name in variables:
                        if name not in dataset.variables:
                            problems.append(f"{path}: variable {name!r} not found")
                            continue
                        shape = tuple(dataset.variables[name].shape)
                        if len(shape) < minimum_rank or shape[-2:] != (nm, mm):
                            expected = "(..., " if minimum_rank > 2 else "("
                            problems.append(f"{path}: variable {name!r} has shape {shape}, expected "
                                            f"{expected}{nm}, {mm}) for nm = {nm} and mm = {mm}")
            except OSError as error:
                problems.append(f"{path}: cannot be read ({error})")

        gridtopo = self.input_dir / "umwm.gridtopo"
        if config.grid.topofromfile:
            check(gridtopo, ["lon", "lat", "z"] if config.grid.gridfromfile else ["z"], 2)
        elif config.grid.gridfromfile:
            check(gridtopo if gridtopo.is_file() else self.input_dir / "umwm.grid", ["lon", "lat"], 2)
        variables = self._required_forcing_variables()
        if variables:
            for time in config.forcing_times():
                check(self.input_dir / forcing_filename(time), variables, 3)
        return problems

    def check_inputs(self) -> None:
        """Raise ``FileNotFoundError`` if required input files are missing, or
        ``ValueError`` if they lack the required variables or have the wrong
        shape; see :meth:`input_problems`."""
        problems = self.input_problems()
        if not problems:
            return
        listing = "\n  ".join(problems[:10])
        more = f"\n  ... and {len(problems) - 10} more" if len(problems) > 10 else ""
        if all(problem.startswith("missing ") for problem in problems):
            listing = "\n  ".join(problem[len("missing "):] for problem in problems[:10])
            raise FileNotFoundError(f"missing input files:\n  {listing}{more}")
        raise ValueError(f"input file problems:\n  {listing}{more}")

    def initialize(self) -> "Domain":
        """Validate the configuration, check the inputs, write the namelists,
        and initialize the model: construct the grids, load the forcing at
        the start time, and compute the initial wave state."""
        if self._handle is not None:
            raise RuntimeError("the domain is already initialized")
        self.config.validate()
        self.check_inputs()
        self.write_namelists()
        namelist = os.path.relpath(self.namelist_path, self.workdir)
        with self._in_workdir(), self._model_output():
            handle = self._lib.domain_new()
            stat = self._lib.domain_initialize(handle, namelist)
            if stat != 0:
                self._lib.domain_delete(handle)
                raise ConfigError(f"the model rejected the configuration in {self.namelist_path}; see the messages above")
        self._handle = handle
        return self

    # ---- time stepping ------------------------------------------------------

    @property
    def time(self) -> datetime:
        """The current model time (the start time before initialization)."""
        if self._handle is None:
            return self.config.start_time
        return parse_time(self._lib.domain_time(self._handle))

    @property
    def finished(self) -> bool:
        """Whether the model has reached the stop time."""
        return self._handle is not None and self._lib.domain_finished(self._handle)

    def advance(self) -> datetime:
        """Integrate over one global time step (``domain.dtg`` seconds):
        update the forcing, integrate the wave spectrum, and write output
        due at the new time. Returns the new model time."""
        if self._handle is None:
            self.initialize()
        with self._in_workdir(), self._model_output():
            self._lib.domain_advance(self._handle)
        return self.time

    def iterate(self) -> Iterator[datetime]:
        """Advance step by step, yielding the model time after each step."""
        if self._handle is None:
            self.initialize()
        while not self.finished:
            yield self.advance()

    def run(self, until=None) -> "Domain":
        """Integrate from the current time to the stop time, or to ``until``
        (a datetime or time string) if that comes first."""
        if self._handle is None:
            self.initialize()
        if until is None:
            with self._in_workdir(), self._model_output():
                self._lib.domain_run(self._handle)
            return self
        until = until if isinstance(until, datetime) else parse_time(until)
        while not self.finished and self.time < until:
            self.advance()
        return self

    def diagnostics(self) -> None:
        """Update the integrated wave diagnostics (``swh``, ``mwp``, ...) from the current spectrum."""
        self._lib.domain_diagnostics(self._require_handle())

    def ocean_stress(self) -> None:
        """Update the ocean-side momentum and energy fluxes and the form drag
        components, which the model otherwise evaluates at gridded output times."""
        self._lib.domain_ocean_stress(self._require_handle())

    # ---- output on demand ---------------------------------------------------

    def write_output(self) -> Path:
        """Write the gridded output file for the current model time, as the
        model does at its gridded output times (updating the ocean stress
        first), and return its path. Collective in parallel runs."""
        handle = self._require_handle()
        with self._in_workdir(), self._model_output():
            self._lib.domain_write_output(handle)
        return self.output_dir / f"umwmout_{self.time:%Y-%m-%d_%H:%M:%S}.nc"

    def write_spectrum_output(self) -> dict[str, list[Path]]:
        """Append the spectra at the configured points at the current model
        time to the spectrum output files, and return the files by point name."""
        handle = self._require_handle()
        with self._in_workdir(), self._model_output():
            self._lib.domain_write_spectrum(handle)
        return _output.spectrum_output_files(self.output_dir)

    def write_restart(self) -> Path:
        """Write a restart file for the current model time and return its
        path. A later run restarts from it with ``config.domain.restart``
        set and the start time set to this time. Collective in parallel runs."""
        handle = self._require_handle()
        with self._in_workdir(), self._model_output():
            self._lib.domain_write_restart(handle)
        return self.restart_dir / f"umwmrst_{self.time:%Y-%m-%d_%H:%M:%S}.nc"

    # ---- teardown -----------------------------------------------------------

    def finalize(self) -> None:
        """Release the model memory. Views of fields become invalid."""
        if self._handle is not None:
            handle, self._handle = self._handle, None
            self._lib.domain_finalize(handle)
            self._lib.domain_delete(handle)

    close = finalize

    def clean(self, output: bool = True, restart: bool = True) -> list[Path]:
        """Remove the files the model wrote to the output and restart
        directories (and nothing else), and return their paths."""
        removed = []
        if self.rank == 0:
            patterns = []
            if output:
                patterns += [(self.output_dir, "umwmout_*.nc"), (self.output_dir, "umwmspc_*.nc"),
                             (self.output_dir, "umwmout.grid")]
            if restart:
                patterns.append((self.restart_dir, "umwmrst_*.nc"))
            for directory, pattern in patterns:
                for path in sorted(directory.glob(pattern)):
                    path.unlink()
                    removed.append(path)
        self.barrier()
        return removed

    def __enter__(self) -> "Domain":
        return self

    def __exit__(self, *exc_info) -> None:
        self.finalize()

    def __del__(self):
        try:
            self.finalize()
        except Exception:
            pass

    def __repr__(self) -> str:
        state = f"time={self.time:%Y-%m-%d %H:%M:%S}" if self.initialized else "not initialized"
        return f"Domain(workdir={str(self.workdir)!r}, {state})"

    def summary(self) -> str:
        """A short description of the domain: grid, spectrum, time span,
        forcing sources, output, and processes."""
        config = self.config
        domain, grid, forcing, constant, output = (
            config.domain, config.grid, config.forcing, config.forcing_constant, config.output)
        lines = [f"UMWM domain in {self.workdir}"]
        cells = f"{domain.mm} x {domain.nm} grid points"
        if self.initialized:
            sizes = self.sizes
            cells += f", {sizes['im']} sea cells and {sizes['imm'] - sizes['im']} land"
        geometry = "geometry from input/umwm.grid" if grid.gridfromfile else f"{grid.delx:g} x {grid.dely:g} m spacing"
        depth = "bathymetry from input/umwm.gridtopo" if grid.topofromfile else f"{grid.dpt:g} m depth"
        lines.append(f"  grid: {cells}, {geometry}, {depth}" + (", global" if domain.isglobal else ""))
        lines.append(f"  spectrum: {domain.om} frequencies from {domain.fmin:g} to {domain.fmax:g} Hz "
                     f"(prognostic to {domain.fprog:g} Hz), {domain.pm} directions")
        time = f"  time: {domain.starttimestr} to {domain.stoptimestr}, global step {domain.dtg:g} s"
        if self.initialized:
            time += f", now at {self.time:%Y-%m-%d %H:%M:%S}"
        if domain.restart:
            time += ", from a restart file"
        lines.append(time)
        sources = [name.replace("_", " ") for name in ("winds", "currents", "air_density", "water_density", "seaice")
                   if getattr(forcing, name)]
        constants = []
        if not forcing.winds:
            constants.append(f"wind {constant.wspd0:g} m/s at {constant.wdir0:g} rad")
        if not forcing.currents:
            constants.append(f"current ({constant.uc0:g}, {constant.vc0:g}) m/s")
        if not forcing.air_density:
            constants.append(f"air density {constant.rhoa0:g} kg/m^3")
        if not forcing.water_density:
            constants.append(f"water density {constant.rhow0:g} kg/m^3")
        if not forcing.seaice:
            constants.append(f"sea ice fraction {constant.fice0:g}")
        if sources and forcing.steady:
            sources[-1] += " (steady)"
        lines.append("  forcing: " + ("from input files: " + ", ".join(sources) + "; " if sources else "")
                     + "constant " + ", ".join(constants))

        def interval(value):
            return "off" if value == 0 else "every global step" if value == -1 else f"every {value} h"

        lines.append(f"  output: gridded {interval(output.outgrid)}, spectrum {interval(output.outspec)} at "
                     f"{len(config.spectrum_points)} point(s), restart {interval(output.outrst)}, Stokes drift "
                     + (f"at {len(config.stokes.depths)} depths" if output.stokes else "off"))
        lines.append(f"  processes: {self.size}, library {'with' if self.mpi else 'without'} MPI")
        return "\n".join(lines)

    # ---- fields -------------------------------------------------------------

    @property
    def sizes(self) -> dict[str, int]:
        """Grid and spectrum sizes: ``mm``, ``nm`` (grid), ``im`` (sea cells),
        ``imm`` (all cells), ``istart``/``iend`` (cells of this process),
        ``iistart``/``iiend`` (with halo), ``om``, ``pm`` (spectrum), ``lm`` (Stokes depths)."""
        handle = self._require_handle()
        return {name: self._lib.domain_int(handle, name) for name in INT_ATTRIBUTES}

    @property
    def shape(self) -> tuple[int, int]:
        """The grid shape ``(ny, nx)`` of gridded fields."""
        handle = self._require_handle()
        return self._lib.domain_int(handle, "nm"), self._lib.domain_int(handle, "mm")

    def attribute(self, name: str) -> float:
        """A scalar attribute of the domain: ``dts`` (current source time step), ``sumt``,
        ``dtamin`` (advection stability limit), ``dtr``, ``dtr_max``, ``dth``, ``dlnf``."""
        return self._lib.domain_real(self._require_handle(), name)

    def field(self, name: str, copy: bool = True) -> np.ndarray:
        """A field of the model state as a NumPy array in the Fortran index
        order (column-major), see :mod:`umwm.fields` for the names.

        Cell fields are indexed by the cell index of this process; the
        spectrum ``F`` is ``(frequency, direction, cell)`` including the
        halo cells; grid fields are ``(x, y)``. With ``copy=False`` the array
        views the Fortran memory without copying: it reflects the model as it
        advances, writes modify the model state, and it must not be used
        after :meth:`finalize`.
        """
        view, _ = self._lib.domain_field(self._require_handle(), name)
        return view.copy(order="F") if copy else view

    def view(self, name: str) -> np.ndarray:
        """A zero-copy view of a field; see :meth:`field`."""
        return self.field(name, copy=False)

    def bounds(self, name: str) -> tuple[int, ...]:
        """The Fortran lower bounds of the dimensions of a field, such as the
        first cell index of the cell dimension."""
        _, lbounds = self._lib.domain_field(self._require_handle(), name)
        return lbounds

    def grid_field(self, name: str, level: int = 0) -> np.ndarray:
        """A field on the horizontal grid as an array of shape ``(ny, nx)``.

        Cell fields are gathered from all processes to the root process and
        remapped to the grid, with NaN on land; on other processes the array
        is all NaN. Grid fields (``lon``, ``depth``, ``seamask``, ...) are
        returned as they are. For the Stokes drift components, ``level``
        selects the depth level (0-based).
        """
        handle = self._require_handle()
        nm, mm = self.shape
        view, _ = self._lib.domain_field(handle, name)
        if view.ndim == 2 and view.shape == (mm, nm):
            return np.array(view.T)
        values = np.full((nm, mm), np.nan, dtype=np.float32)
        self._lib.domain_gather(handle, name, level + 1, values)
        return values

    def remap(self, values: np.ndarray, lbound: int | None = None) -> np.ndarray:
        """Remap an array indexed by cell along its first dimension to the grid,
        returning an array of shape ``(ny, nx, ...)`` with NaN on land.

        ``lbound`` is the Fortran index of the first cell; by default it is
        1 for arrays over all cells and the first cell of this process
        otherwise. In parallel runs, only arrays over all cells can be
        remapped, since a process holds only its own cells.
        """
        handle = self._require_handle()
        sizes = self.sizes
        values = np.asarray(values)
        if lbound is None:
            lbound = 1 if values.shape[0] in (sizes["imm"], sizes["im"]) else sizes["istart"]
        first, last = lbound, lbound + values.shape[0] - 1
        if self.size > 1 and (first > 1 or last < sizes["im"]):
            raise ValueError("only arrays over all cells can be remapped in a parallel run; use grid_field")
        ii = self._lib.domain_field(handle, "ii")[0].T  # (ny, nx)
        sea = (ii >= max(first, 1)) & (ii <= min(last, sizes["im"]))
        result = np.full(ii.shape + values.shape[1:], np.nan, dtype=np.result_type(values.dtype, np.float32))
        result[sea] = values[ii[sea] - lbound]
        return result

    # ---- xarray -------------------------------------------------------------

    def _grid_coords(self) -> dict:
        return {
            "lon": (("y", "x"), self.grid_field("lon"), _fields.attrs("lon")),
            "lat": (("y", "x"), self.grid_field("lat"), _fields.attrs("lat")),
        }

    def state(self, variables=None, ocean_stress: bool = False) -> xr.Dataset:
        """The gridded model state at the current time as a dataset, like a
        gridded output file: the fields in ``variables`` (by default all
        gridded output fields) with dimensions ``(time, y, x)``, the Stokes
        drift with a ``z`` dimension when enabled, and the grid.

        With ``ocean_stress=True``, the ocean-side fluxes and the form drag
        components are updated first, as the model does before writing
        gridded output. In parallel runs, the fields are gathered to the
        root process; other processes get NaN fields.
        """
        handle = self._require_handle()
        if ocean_stress:
            self.ocean_stress()
        names = list(variables) if variables is not None else list(_fields.GRID_OUTPUT_FIELDS)
        time = np.array([np.datetime64(self.time, "s")])
        data_vars = {}
        for name in names:
            data_vars[name] = (("time", "y", "x"), self.grid_field(name)[np.newaxis], _fields.attrs(name))
        for name in ("seamask", "depth"):
            data_vars[name] = (("y", "x"), self.grid_field(name), _fields.attrs(name))
        coords = {"time": ("time", time, {"long_name": "time"})}
        coords.update(self._grid_coords())
        if variables is None and self.config.output.stokes and self._lib.domain_int(handle, "lm") > 0:
            depths = self.field("z")
            coords["z"] = ("z", depths, {"long_name": "depth", "units": "m"})
            for name in ("u_stokes", "v_stokes"):
                levels = np.stack([self.grid_field(name, level) for level in range(len(depths))])
                data_vars[name] = (("time", "z", "y", "x"), levels[np.newaxis], _fields.attrs(name))
            data_vars["d_stokes"] = (("time", "y", "x"), self.grid_field("d_stokes")[np.newaxis], _fields.attrs("d_stokes"))
        return xr.Dataset(data_vars, coords=coords, attrs={"title": "UMWM model state", "workdir": str(self.workdir)})

    def _spectral_coords(self) -> dict:
        return {
            "frequency": ("frequency", self.field("frequency"), _fields.attrs("frequency")),
            "direction": ("direction", self.field("direction"), _fields.attrs("direction")),
        }

    def spectrum(self, name: str = "F") -> xr.DataArray:
        """The wave spectrum ``F`` (or a spectral source term: ``Sin``, ``Sds``,
        ``Snl``) on the grid as a data array with dimensions
        ``(y, x, direction, frequency)``, NaN on land.

        In parallel runs, a process holds only its own cells, and the array
        has dimensions ``(cell, direction, frequency)`` instead, with the
        Fortran cell index as coordinate.
        """
        handle = self._require_handle()
        values, lbounds = self._lib.domain_field(handle, name)
        if values.ndim != 3:
            raise ValueError(f"{name!r} is not a spectral field over the cells")
        sizes = self.sizes
        coords = self._spectral_coords()
        columns = values.T  # (cell, direction, frequency), C-order view
        if self.size == 1:
            first = lbounds[2]
            ii = self._lib.domain_field(handle, "ii")[0].T
            sea = (ii >= 1) & (ii <= sizes["im"])
            data = np.full(ii.shape + columns.shape[1:], np.nan, dtype=columns.dtype)
            data[sea] = columns[ii[sea] - first]
            coords.update(self._grid_coords())
            return xr.DataArray(data, dims=("y", "x", "direction", "frequency"), coords=coords,
                                name=name, attrs=_fields.attrs(name))
        cells = np.arange(lbounds[2], lbounds[2] + values.shape[2])
        coords["cell"] = ("cell", cells, {"long_name": "cell index (Fortran, 1-based)"})
        return xr.DataArray(np.array(columns), dims=("cell", "direction", "frequency"), coords=coords,
                            name=name, attrs=_fields.attrs(name))

    def spectrum_at(self, ix: int, iy: int, name: str = "F") -> xr.DataArray:
        """The spectrum ``F`` (or a spectral source term) at one grid point as
        a data array ``(direction, frequency)``, without building the whole
        grid. ``ix`` and ``iy`` are 0-based array indices (namelist spectrum
        points are 1-based). In parallel runs the point must belong to the
        cells of this process.
        """
        handle = self._require_handle()
        sizes = self.sizes
        if not (0 <= ix < sizes["mm"] and 0 <= iy < sizes["nm"]):
            raise IndexError(f"grid point ({ix}, {iy}) is outside the {sizes['mm']} x {sizes['nm']} grid")
        cell = int(self._lib.domain_field(handle, "ii")[0][ix, iy])
        if cell > sizes["im"]:
            raise ValueError(f"grid point ({ix}, {iy}) is land")
        values, lbounds = self._lib.domain_field(handle, name)
        if values.ndim != 3:
            raise ValueError(f"{name!r} is not a spectral field over the cells")
        first, last = lbounds[2], lbounds[2] + values.shape[2] - 1
        if not first <= cell <= last:
            raise ValueError(f"grid point ({ix}, {iy}) belongs to another process")
        coords = self._spectral_coords()
        coords["x"] = ix
        coords["y"] = iy
        for coordinate in ("lon", "lat"):
            coords[coordinate] = ((), float(self._lib.domain_field(handle, coordinate)[0][ix, iy]), _fields.attrs(coordinate))
        return xr.DataArray(np.array(values[:, :, cell - first].T), dims=("direction", "frequency"),
                            coords=coords, name=name, attrs=_fields.attrs(name))

    # ---- files --------------------------------------------------------------

    def output_files(self) -> list[Path]:
        """The gridded output files written so far."""
        return _output.output_files(self.output_dir)

    def output(self, **kwargs) -> xr.Dataset:
        """The gridded output files as one dataset along ``time``; see :func:`umwm.open_output`."""
        return _output.open_output(self.output_dir, **kwargs)

    def spectrum_output(self, **kwargs) -> dict[str, xr.Dataset]:
        """The spectrum output as a dataset per point; see :func:`umwm.open_spectrum_output`."""
        return _output.open_spectrum_output(self.output_dir, **kwargs)

    def restart_files(self) -> list[Path]:
        """The restart files written so far."""
        return _output.restart_files(self.restart_dir)

    def grid(self, **kwargs) -> xr.Dataset:
        """The grid definition written by the model to ``output/umwmout.grid``."""
        return _output.open_grid(self.output_dir, **kwargs)

    def write_forcing(self, time, **fields) -> Path:
        """Write a forcing input file for ``time`` into the input directory;
        see :func:`umwm.write_forcing`. Only the root process writes."""
        if self.rank == 0:
            path = write_forcing(self.input_dir, time, **fields)
        else:
            path = self.input_dir / forcing_filename(time)
        self.barrier()
        return path

    def write_grid(self, lon, lat, depth=None, **kwargs) -> Path:
        """Write the grid input file into the input directory; see :func:`umwm.write_grid`.
        Only the root process writes."""
        if self.rank == 0:
            path = write_grid(self.input_dir, lon, lat, depth, **kwargs)
        else:
            path = self.input_dir / ("umwm.grid" if depth is None else "umwm.gridtopo")
        self.barrier()
        return path

    # ---- plotting -----------------------------------------------------------

    def plot(self, name: str = "swh", level: int = 0, **kwargs):
        """Draw a field of the current state on the grid with matplotlib; see
        :func:`umwm.plotting.plot_field`. Returns the matplotlib axes."""
        from .plotting import plot_field

        data = self.grid_field(name, level)
        dataset = xr.Dataset({name: (("y", "x"), data, _fields.attrs(name))}, coords=self._grid_coords(),
                             attrs={"time": f"{self.time:%Y-%m-%d %H:%M:%S}"})
        grid = self.grid() if (self.output_dir / "umwmout.grid").is_file() else None
        return plot_field(dataset, name, grid=grid, **kwargs)


def run(config=None, workdir=None, **kwargs) -> Domain:
    """Run the model to the stop time and return the domain, still initialized
    so that its state and output can be inspected; see :class:`Domain`."""
    domain = Domain(config, workdir, **kwargs)
    domain.run()
    return domain
