"""Writing of UMWM input files: forcing fields and grids.

The model reads, for each global time step, a file ``input/umwmin_<time>.nc``
with the forcing fields on the native grid, or once a file ``input/umwmin.nc``
of steady forcing (``forcing.steady``), and optionally the grid geometry
from ``input/umwm.grid`` or ``input/umwm.gridtopo``. These functions write
such files from NumPy arrays of shape ``(ny, nx)``.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr

from .config import format_time, parse_time

__all__ = ["forcing_filename", "forcing_files", "write_forcing", "write_forcing_dataset", "write_grid",
           "rankine_vortex", "write_rankine_forcing", "STEADY_FORCING_FILE"]

FORCING_VARIABLES = {
    "uw": ("wind x-component", "m s-1"),
    "vw": ("wind y-component", "m s-1"),
    "uc": ("current x-component", "m s-1"),
    "vc": ("current y-component", "m s-1"),
    "rhoa": ("air density", "kg m-3"),
    "rhow": ("water density", "kg m-3"),
    "fice": ("sea ice fraction", "1"),
}
"""Forcing variables of the input files, with description and units."""


STEADY_FORCING_FILE = "umwmin.nc"
"""The input file of steady forcing (``forcing.steady``), read once."""


def forcing_filename(time, steady: bool = False) -> str:
    """The input file name for a forcing time: ``umwmin_YYYY-MM-DD_HH:MM:SS.nc``,
    or ``umwmin.nc`` for steady forcing (``forcing.steady``), whatever the time."""
    if steady:
        return STEADY_FORCING_FILE
    time = time if isinstance(time, datetime) else parse_time(time)
    return f"umwmin_{time:%Y-%m-%d_%H:%M:%S}.nc"


def forcing_files(config) -> list[str]:
    """The names of the forcing input files that the model reads for a
    configuration: one for the start time and for the end of each global time
    step, or the steady forcing file; none if it reads no forcing from file."""
    if not config.forcing_from_file:
        return []
    if config.forcing.steady:
        return [STEADY_FORCING_FILE]
    return [forcing_filename(time) for time in config.step_times(include_start=True)]


def _as_field(name, value, shape, dtype=np.float32):
    array = np.asarray(value, dtype=dtype)
    if array.ndim == 0:
        array = np.full(shape, array, dtype=dtype)
    if array.shape != shape:
        raise ValueError(f"{name} has shape {array.shape}, expected {shape} (ny, nx)")
    return array


def write_forcing(path, time, *, uw=None, vw=None, uc=None, vc=None, rhoa=None, rhow=None,
                  fice=None, lon=None, lat=None, x=None, y=None, attrs=None, steady=False) -> Path:
    """Write a forcing input file for one time.

    ``path`` is the input directory, in which case the file is named after
    the time, or ``umwmin.nc`` for steady forcing (``steady``, for
    ``forcing.steady``), or the path of the file to write. Fields are arrays of shape
    ``(ny, nx)`` or scalars; the fields that the configuration reads from
    file (``forcing.winds``: ``uw`` and ``vw``; ``forcing.currents``: ``uc``
    and ``vc``; ``forcing.air_density``: ``rhoa``; ``forcing.water_density``:
    ``rhow``; ``forcing.seaice``: ``fice``) must be given. ``lon`` and
    ``lat`` (2-D) and ``x`` and ``y`` (1-D) are optional coordinates written
    for reference. Returns the path of the file written.
    """
    time = time if isinstance(time, datetime) else parse_time(time)
    given = {name: value for name, value in dict(uw=uw, vw=vw, uc=uc, vc=vc, rhoa=rhoa, rhow=rhow, fice=fice).items()
             if value is not None}
    if not given:
        raise ValueError("at least one forcing field must be given")
    shapes = {np.shape(value) for value in given.values() if np.ndim(value) == 2}
    if len(shapes) != 1:
        raise ValueError("the forcing fields must share one shape (ny, nx)" if shapes else
                         "at least one forcing field must be a 2-D array of shape (ny, nx)")
    shape = shapes.pop()
    ny, nx = shape

    data_vars = {}
    if lon is not None and lat is not None:
        data_vars["lon"] = (("y", "x"), _as_field("lon", lon, shape), {"long_name": "longitude", "units": "degrees_east"})
        data_vars["lat"] = (("y", "x"), _as_field("lat", lat, shape), {"long_name": "latitude", "units": "degrees_north"})
    for name, value in given.items():
        description, units = FORCING_VARIABLES[name]
        data_vars[name] = (("time", "y", "x"), _as_field(name, value, shape)[np.newaxis], {"long_name": description, "units": units})

    coords = {"time": ("time", np.array([np.datetime64(time, "s")]), {"long_name": "time"})}
    if x is not None:
        coords["x"] = ("x", np.asarray(x, dtype=np.float32).reshape(nx), {"long_name": "x coordinate", "units": "m"})
    if y is not None:
        coords["y"] = ("y", np.asarray(y, dtype=np.float32).reshape(ny), {"long_name": "y coordinate", "units": "m"})

    dataset = xr.Dataset(data_vars=data_vars, coords=coords, attrs=dict(attrs or {}, valid_time=format_time(time)))
    path = Path(path)
    if path.is_dir() or path.suffix != ".nc":
        path = path / forcing_filename(time, steady)
    path.parent.mkdir(parents=True, exist_ok=True)
    encoding = {name: {"dtype": "float32"} for name in dataset.data_vars}
    dataset.to_netcdf(path, engine="netcdf4", format="NETCDF4", encoding=encoding, unlimited_dims=("time",))
    return path


def _to_datetime(value) -> datetime:
    """A datetime from a numpy datetime64, a cftime object, a string, or a datetime."""
    if isinstance(value, datetime):
        return value
    if isinstance(value, np.datetime64):
        return np.datetime64(value, "s").astype(datetime)
    if hasattr(value, "strftime"):
        return parse_time(value.strftime("%Y-%m-%d %H:%M:%S"))
    return parse_time(str(value))


def write_forcing_dataset(path, dataset, *, rename=None) -> list[Path]:
    """Write one forcing file per time of an xarray dataset into the directory ``path``.

    The dataset holds the forcing variables (``uw``, ``vw``, ``uc``, ``vc``,
    ``rhoa``, ``rhow``, ``fice``; others are ignored) with dimensions
    ``(time, y, x)``, and optionally ``lon`` and ``lat``. ``rename`` maps the
    names in the dataset to these, for example ``{"u10": "uw", "v10": "vw"}``.
    Returns the paths of the files written, in time order.
    """
    rename = {old: new for old, new in dict(rename or {}).items() if old in dataset}
    dataset = dataset.rename(rename) if rename else dataset
    if "time" not in dataset.dims:
        raise ValueError("the dataset needs a time dimension")
    names = [name for name in FORCING_VARIABLES if name in dataset]
    if not names:
        raise ValueError("the dataset has none of the forcing variables " + ", ".join(FORCING_VARIABLES))
    if {"y", "x"} <= set(dataset.dims):
        dataset = dataset.transpose(..., "y", "x")
    static = {}
    for coordinate in ("lon", "lat"):
        if coordinate in dataset:
            values = dataset[coordinate]
            static[coordinate] = (values.isel(time=0) if "time" in values.dims else values).values
    paths = []
    for index in range(dataset.sizes["time"]):
        time = _to_datetime(dataset["time"].values[index])
        fields = {name: dataset[name].isel(time=index).values for name in names}
        paths.append(write_forcing(path, time, **fields, **static))
    return paths


def write_grid(path, lon, lat, depth=None, *, attrs=None) -> Path:
    """Write the grid input file.

    ``lon`` and ``lat`` are arrays of shape ``(ny, nx)``. Without ``depth``,
    ``umwm.grid`` is written (for ``grid.gridfromfile``); with ``depth``, the
    water depth in meters, positive down, ``umwm.gridtopo`` is written (for
    ``grid.topofromfile``), with land where the depth is not positive.
    ``path`` is the input directory or the path of the file to write.
    Returns the path of the file written.
    """
    lon = np.asarray(lon, dtype=np.float32)
    lat = np.asarray(lat, dtype=np.float32)
    if lon.ndim != 2 or lon.shape != lat.shape:
        raise ValueError("lon and lat must be 2-D arrays of the same shape (ny, nx)")
    data_vars = {
        "lon": (("y", "x"), lon, {"long_name": "longitude", "units": "degrees_east"}),
        "lat": (("y", "x"), lat, {"long_name": "latitude", "units": "degrees_north"}),
    }
    name = "umwm.grid"
    if depth is not None:
        depth = _as_field("depth", depth, lon.shape)
        data_vars["z"] = (("y", "x"), -depth, {"long_name": "bottom elevation (negative below sea level)", "units": "m"})
        name = "umwm.gridtopo"
    dataset = xr.Dataset(data_vars=data_vars, attrs=dict(attrs or {}))
    path = Path(path)
    if path.is_dir() or path.suffix not in (".grid", ".gridtopo"):
        path = path / name
    path.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_netcdf(path, engine="netcdf4", format="NETCDF4")
    return path


# ---- idealized forcing -------------------------------------------------------

def rankine_vortex(x2d, y2d, vmax, rmw, clockwise=False):
    """Wind components of a Rankine vortex centered at the origin of ``x2d``, ``y2d``.

    The tangential wind speed grows linearly to ``vmax`` at the radius of
    maximum wind ``rmw`` and decays as ``1/r`` outside. Returns ``(uw, vw)``
    as float32 arrays.
    """
    if rmw <= 0:
        raise ValueError("rmw must be > 0")
    radius = np.hypot(x2d, y2d)
    wspd = np.empty_like(radius, dtype=np.float64)
    inside = radius <= rmw
    outside = ~inside
    wspd[inside] = vmax * radius[inside] / rmw
    wspd[outside] = vmax * rmw / radius[outside]
    wspd[radius == 0.0] = 0.0

    uw = np.zeros_like(radius, dtype=np.float64)
    vw = np.zeros_like(radius, dtype=np.float64)
    nonzero = radius > 0.0
    if clockwise:
        uw[nonzero] = wspd[nonzero] * y2d[nonzero] / radius[nonzero]
        vw[nonzero] = -wspd[nonzero] * x2d[nonzero] / radius[nonzero]
    else:
        uw[nonzero] = -wspd[nonzero] * y2d[nonzero] / radius[nonzero]
        vw[nonzero] = wspd[nonzero] * x2d[nonzero] / radius[nonzero]
    return uw.astype(np.float32), vw.astype(np.float32)


def write_rankine_forcing(path, start_time, hours, *, nx=51, ny=21, dx=10000.0, dy=10000.0,
                          vmax=30.0, rmw=50000.0, center_x=-150000.0, center_y=0.0,
                          translation_x=50000.0, translation_y=0.0, clockwise=False,
                          uc=0.0, vc=0.0, rhoa=1.2, rhow=1030.0,
                          center_lon=0.0, center_lat=0.0) -> list[Path]:
    """Write hourly forcing files for an idealized, translating Rankine hurricane vortex.

    The grid is Cartesian with ``nx`` by ``ny`` cells spaced ``dx`` by ``dy``
    meters and centered at the origin. The vortex center starts at
    ``(center_x, center_y)`` and moves by ``(translation_x, translation_y)``
    meters per hour; files are written for hours 0 to ``hours`` inclusive
    into the directory ``path``. Currents and densities are uniform.
    Returns the paths of the files written.
    """
    if nx <= 0 or ny <= 0:
        raise ValueError("nx and ny must be > 0")
    if dx <= 0 or dy <= 0:
        raise ValueError("dx and dy must be > 0")
    if hours < 0:
        raise ValueError("hours must be >= 0")
    start_time = start_time if isinstance(start_time, datetime) else parse_time(start_time)

    x = (np.arange(nx, dtype=np.float32) - (nx - 1) / 2) * dx
    y = (np.arange(ny, dtype=np.float32) - (ny - 1) / 2) * dy
    x2d, y2d = np.meshgrid(x, y)
    meters_per_degree = 111_320.0
    coslat = np.cos(np.deg2rad(center_lat))
    if abs(coslat) < 1.0e-6:
        raise ValueError("center_lat is too close to a pole")
    lat = (center_lat + y2d / meters_per_degree).astype(np.float32)
    lon = (center_lon + x2d / (meters_per_degree * coslat)).astype(np.float32)

    paths = []
    for hour in range(hours + 1):
        time = start_time + timedelta(hours=hour)
        vortex_x = center_x + hour * translation_x
        vortex_y = center_y + hour * translation_y
        uw, vw = rankine_vortex(x2d - vortex_x, y2d - vortex_y, vmax, rmw, clockwise)
        attrs = {
            "title": "Idealized moving Rankine hurricane vortex forcing",
            "rankine_vmax_m_s": vmax, "rankine_rmw_m": rmw, "grid_dx_m": dx, "grid_dy_m": dy,
            "vortex_center_x_m": vortex_x, "vortex_center_y_m": vortex_y,
            "rotation": "clockwise" if clockwise else "counter-clockwise",
        }
        paths.append(write_forcing(path, time, uw=uw, vw=vw, uc=uc, vc=vc, rhoa=rhoa, rhow=rhow,
                                   lon=lon, lat=lat, x=x, y=y, attrs=attrs))
    return paths
