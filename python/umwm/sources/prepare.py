"""Preparation of the input files of a run from a data source."""

from __future__ import annotations

import time as _time
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

from .. import __version__
from ..forcing import forcing_filename, write_forcing, write_grid
from . import SourceError, get_source, source_problems
from .base import FORCING_FIELDS, wave_fields
from .regrid import Regridder, fill_missing, model_grid, open_boundary_cells
from .waves import regrid_partitions, standard_partitions

__all__ = ["prepare_inputs", "pending_inputs", "input_grid_problem", "MEMORY_LIMIT", "CORNER_TOLERANCE"]

MEMORY_LIMIT = 512 * 2**20
"""The bytes of source fields read at a time."""

CORNER_TOLERANCE = 1e-3
"""The largest difference [degrees] of the longitude or latitude of the first
or last cell of an input file from that of the grid of the model, for the file
to be for the grid, as in the model."""

PREPARED = "prepared by umwm"
"""The start of the history attribute of the input files that umwm prepared."""

_GRID_FIELDS = ("lon", "lat", "z", "uw", "vw", "uc", "vc", "fice", "rhoa", "rhow") + tuple(f"phs{p}" for p in range(10))


def input_grid_problem(path, shape, corners=None) -> str | None:
    """Why the input file at ``path`` is not for a grid of ``shape`` ``(nm,
    mm)``, or None if it is, as far as these tell: its fields are of another
    shape, or, if it has ``lon`` and ``lat`` and ``corners`` are given, those
    of its first and last cells are not ``corners``, the ``(lon, lat)`` of the
    first and last cells of the grid. The model checks its input files in the
    same way."""
    from netCDF4 import Dataset

    with Dataset(path) as dataset:
        variables = dataset.variables
        for name in _GRID_FIELDS:
            if name in variables and variables[name].ndim >= 2:
                found = tuple(variables[name].shape[-2:])
                if found != tuple(shape):
                    return (f"{name} is {found[1]} x {found[0]}, but the grid is {shape[1]} x {shape[0]} "
                            "(mm x nm)")
        if corners is None or not all(name in variables and variables[name].ndim >= 2 for name in ("lon", "lat")):
            return None
        lon, lat = variables["lon"], variables["lat"]
        index = (0,) * (lon.ndim - 2)
        found = [(float(lon[index + (0, 0)]), float(lat[index + (0, 0)])),
                 (float(lon[index + (-1, -1)]), float(lat[index + (-1, -1)]))]
    if all(abs((flon - glon + 180.0) % 360.0 - 180.0) <= CORNER_TOLERANCE and abs(flat - glat) <= CORNER_TOLERANCE
           for (flon, flat), (glon, glat) in zip(found, corners)):
        return None

    def text(point):
        return f"({point[0]:.3f}, {point[1]:.3f})"

    return (f"its first and last cells are at {text(found[0])} and {text(found[1])}, but those of the grid "
            f"are at {text(corners[0])} and {text(corners[1])}")


def _prepared_by_umwm(path) -> bool:
    """Whether umwm prepared the input file at ``path``, as its history tells."""
    from netCDF4 import Dataset

    with Dataset(path) as dataset:
        return str(getattr(dataset, "history", "")).startswith(PREPARED)


def pending_inputs(config, input_dir, force: bool = False, log=None) -> tuple[bool, list[tuple[datetime, Path]]]:
    """The input files of a configuration to prepare in ``input_dir``: whether
    the bathymetry ``umwm.gridtopo`` (with ``grid.topofromfile``) is to be
    prepared, and the forcing files with their times. They are those that do
    not exist, and those for another grid (see :func:`input_grid_problem`)
    that umwm prepared, which are prepared again, as ``log`` tells; or all
    of them with ``force``. A file for another grid that umwm did not
    prepare raises :class:`SourceError`, rather than being replaced."""
    input_dir = Path(input_dir)
    grid = {}
    stale = []

    def expected():
        # the shape of the grid, and its first and last cells if it is a lon-lat grid
        if not grid:
            try:
                lon, lat = model_grid(config, input_dir)
            except (ValueError, OSError):
                grid.update(shape=(config.domain.nm, config.domain.mm), corners=None)
            else:
                grid.update(shape=lon.shape, corners=((float(lon[0, 0]), float(lat[0, 0])),
                                                      (float(lon[-1, -1]), float(lat[-1, -1]))))
        return grid["shape"], grid["corners"]

    def pending(path):
        if force or not path.is_file():
            return True
        problem = input_grid_problem(path, *expected())
        if problem is None:
            return False
        if not _prepared_by_umwm(path):
            raise SourceError(f"the input file {path} is for another grid: {problem}. As umwm did not prepare it, "
                              "it is not replaced: remove it, or prepare the input files in another directory")
        stale.append((path, problem))
        return True

    bathymetry = config.grid.topofromfile and pending(input_dir / "umwm.gridtopo")
    files = []
    if config.forcing_from_file:
        times = [config.start_time] if config.forcing.steady else config.step_times(include_start=True)
        for time in times:
            path = input_dir / forcing_filename(time, config.forcing.steady)
            if pending(path):
                files.append((time, path))
    if stale and log is not None:
        path, problem = stale[0]
        log(f"umwm: {len(stale)} input file(s) in {input_dir} were prepared for another grid ({path.name}: "
            f"{problem}); preparing them again")
    return bathymetry, files


def prepare_inputs(config, input_dir, source: str | None = None, *, force: bool = False, log=print,
                   **options) -> list[Path]:
    """Prepare the input files of a configuration that are missing from
    ``input_dir``, or that are for another grid, or all of them with
    ``force`` (see :func:`pending_inputs`), from its data source
    (``forcing.source``), or from ``source``; ``options`` are those of the
    source, such as the ``url`` of its dataset. Returns the paths of the
    files written.

    The files are on the grid of the model (see
    :func:`umwm.sources.regrid.model_grid`): the bathymetry
    ``umwm.gridtopo``, with ``grid.topofromfile``, where the model is sea if
    at least half of the weight of its bilinear interpolation is on sea in
    the source; and the forcing files of the start time and of the end of
    every global time step, with the fields that the configuration reads
    from file: the wind (``forcing.winds``), the air density
    (``air_density``), the sea ice fraction (``seaice``), and the wave
    partitions (``waves``), at the start time, for the initial spectrum,
    and at every time if the domain has open boundaries, regional with sea
    at its edges. The fields are interpolated in time between the times of
    the source, and the wave partitions are those of the nearest time.
    """
    name = (source or config.source_name or "").strip().lower()
    if not name or name == "none":
        raise SourceError("the configuration has no data source; set forcing.source, such as 'era5'")
    problems = source_problems(config, name)
    if problems:
        raise SourceError("; ".join(problems))
    if config.forcing.currents or config.forcing.water_density:
        raise SourceError(f"the data source {name} does not provide currents or water density")
    input_dir = Path(input_dir)
    bathymetry, pending = pending_inputs(config, input_dir, force, log=log)
    if not bathymetry and not pending:
        log(f"umwm: {name}: the input files are in place in {input_dir}")
        return []
    input_dir.mkdir(parents=True, exist_ok=True)
    lon, lat = model_grid(config, input_dir)
    written = []
    started = _time.monotonic()

    with get_source(name, **options) as dataset:
        log(f"umwm: {name}: preparing the input files in {input_dir} from {dataset.description()}")
        regridder = Regridder(*dataset.grid(), lon, lat)
        first, last = dataset.coverage()
        step = dataset.time_step
        attrs = {"source": dataset.description(), "history": f"prepared by umwm {__version__} on "
                 f"{datetime.now():%Y-%m-%d %H:%M:%S}"}

        def native(time: datetime) -> tuple[datetime, datetime, float]:
            # the times of the source before and after a time, and the weight of the latter
            before = first + (time - first) // step * step
            weight = (time - before) / step
            return before, before + step if weight > 0 else before, weight

        needed = [config.start_time] + [time for time, _ in pending]
        for time in needed:
            before, after, _ = native(time)
            if before < first or after > last:
                raise SourceError(f"{time:%Y-%m-%d %H:%M} is not covered by {name}, which covers "
                                  f"{first:%Y-%m-%d %H:%M} to {last:%Y-%m-%d %H:%M}")

        # bathymetry, written or read, which tells the open boundaries
        gridtopo = input_dir / "umwm.gridtopo"
        if bathymetry:
            time = native(config.start_time)[0]
            source_depth = dataset.read(["depth"], [time], rows=regridder.rows)[time]["depth"]
            sea = regridder.coverage(source_depth) >= 0.5
            depth = np.where(sea, fill_missing(regridder(source_depth), config.physics.dmin), 0.0)
            written.append(write_grid(input_dir, lon, lat, depth, attrs=attrs))
            log(f"umwm: {name}: wrote {written[-1]}: {np.count_nonzero(sea)} sea cells of {sea.size}")
        elif config.grid.topofromfile:
            from netCDF4 import Dataset

            with Dataset(gridtopo) as file:
                depth = -np.asarray(file.variables["z"][:], dtype=np.float64)
        else:
            depth = np.full(lon.shape, config.grid.dpt)

        forcing_names = [field for switch, fields in FORCING_FIELDS.items() if getattr(config.forcing, switch)
                         for field in fields]
        waves = config.forcing.waves
        wave_names = wave_fields(dataset.partitions, dataset.period, dataset.total) if waves else []
        boundaries = waves and bool(open_boundary_cells(config, depth).any())
        if waves:
            log(f"umwm: {name}: wave partitions at " + ("every time, for the open boundaries" if boundaries
                                                         else "the start time, for the initial spectrum"))

        # the fields of each file, from the times of the source it needs
        plan = []
        for time, path in pending:
            before, after, weight = native(time)
            # the model reads a file at every time; with no other field, it has the waves
            with_waves = waves and (boundaries or time == config.start_time or not forcing_names)
            nearest = after if weight >= 0.5 else before
            reads = {before: set(forcing_names), after: set(forcing_names)} if forcing_names else {}
            if with_waves:
                reads.setdefault(nearest, set()).update(wave_names)
            plan.append((time, path, before, after, weight, nearest if with_waves else None, reads))

        row_bytes = 4 * (regridder.rows.stop - regridder.rows.start) * dataset.grid()[0].size
        per_time = max(1, len(forcing_names) + len(wave_names)) * row_bytes
        block = max(1, int(MEMORY_LIMIT // per_time))
        for start in range(0, len(plan), block):
            chunk = plan[start:start + block]
            reads = {}
            for *_, file_reads in chunk:
                for time, names in file_reads.items():
                    reads.setdefault(time, set()).update(names)
            data = {}
            for names in {frozenset(names) for names in reads.values()}:
                times = sorted(time for time, time_names in reads.items() if frozenset(time_names) == names)
                data.update(dataset.read(sorted(names), times, rows=regridder.rows))

            for time, path, before, after, weight, nearest, _ in chunk:
                fields = {}
                if forcing_names:
                    def at_time(field):
                        values = data[before][field]
                        return values if weight == 0 else (1 - weight) * values + weight * data[after][field]

                    if config.forcing.winds:
                        fields["uw"] = fill_missing(regridder(at_time("u10")), 0.0)
                        fields["vw"] = fill_missing(regridder(at_time("v10")), 0.0)
                    if config.forcing.air_density:
                        fields["rhoa"] = fill_missing(regridder(at_time("rhoa")), config.forcing_constant.rhoa0)
                    if config.forcing.seaice:
                        fields["fice"] = np.clip(fill_missing(regridder(at_time("fice")), 0.0), 0.0, 1.0)
                if nearest is not None:
                    partitions = standard_partitions(data[nearest], dataset.partitions, dataset.period,
                                                     dataset.total)
                    fields.update(regrid_partitions(regridder, partitions))
                if not fields:
                    continue
                written.append(write_forcing(path, time, lon=lon, lat=lat, attrs=attrs, **fields))
                log(f"umwm: {name}: wrote {path.name}" + (" with wave partitions" if nearest is not None else ""))

    log(f"umwm: {name}: prepared {len(written)} input files in {_time.monotonic() - started:.0f} s")
    return written
