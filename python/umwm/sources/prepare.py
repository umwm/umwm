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

__all__ = ["prepare_inputs", "pending_inputs", "MEMORY_LIMIT"]

MEMORY_LIMIT = 512 * 2**20
"""The bytes of source fields read at a time."""


def pending_inputs(config, input_dir, force: bool = False) -> tuple[bool, list[tuple[datetime, Path]]]:
    """The input files of a configuration to prepare in ``input_dir``: whether
    the bathymetry ``umwm.gridtopo`` (with ``grid.topofromfile``), and the
    forcing files with their times, that do not exist, or all with ``force``."""
    input_dir = Path(input_dir)
    bathymetry = config.grid.topofromfile and (force or not (input_dir / "umwm.gridtopo").is_file())
    if not config.forcing_from_file:
        return bathymetry, []
    if config.forcing.steady:
        times = [config.start_time]
    else:
        times = config.step_times(include_start=True)
    files = [(time, input_dir / forcing_filename(time, config.forcing.steady)) for time in times]
    return bathymetry, [(time, path) for time, path in files if force or not path.is_file()]


def prepare_inputs(config, input_dir, source: str | None = None, *, force: bool = False, log=print,
                   **options) -> list[Path]:
    """Prepare the input files of a configuration that are missing from
    ``input_dir``, or all of them with ``force``, from its data source
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
    bathymetry, pending = pending_inputs(config, input_dir, force)
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
