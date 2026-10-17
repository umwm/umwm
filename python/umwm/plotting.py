"""Quick-look plots of gridded fields with matplotlib, an optional dependency
(``pip install umwm[plot]``)."""

from __future__ import annotations

import numpy as np
import xarray as xr

from . import fields as _fields

__all__ = ["plot_field"]

DIRECTION_FIELDS = {"wdir", "mwd", "dwd"}
"""Fields in radians on the full circle, drawn with a cyclic colormap."""


def _pyplot():
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("matplotlib is required for plotting; install it with "
                          "'pip install matplotlib' or 'pip install umwm[plot]'") from None
    return plt


def _axes(dataset, grid, shape):
    """The coordinate arrays and axis labels: grid coordinates in km from the
    grid file, longitude and latitude when they are set, or cell indices."""
    ny, nx = shape
    if grid is not None and "xx" in grid and "yy" in grid:
        return (np.asarray(grid["xx"].values) / 1000.0, np.asarray(grid["yy"].values) / 1000.0,
                "x [km]", "y [km]", "equal")
    if "lon" in dataset and "lat" in dataset:
        lon, lat = dataset["lon"], dataset["lat"]
        if "time" in lon.dims:
            lon, lat = lon.isel(time=0), lat.isel(time=0)
        lon, lat = np.asarray(lon.values), np.asarray(lat.values)
        if lon.shape == shape and (np.any(lon != 0) or np.any(lat != 0)):
            return lon, lat, "longitude [degrees east]", "latitude [degrees north]", "auto"
    return np.arange(nx), np.arange(ny), "x [cell]", "y [cell]", "equal"


def plot_field(dataset, name, *, time=-1, level=0, ax=None, cmap=None, grid=None, title=None, **kwargs):
    """Draw a gridded field of an output or state dataset as a color map.

    ``time`` selects the time index when the field has a time dimension,
    ``level`` the depth index of the Stokes drift. ``grid`` is the dataset
    of ``umwmout.grid`` (see :func:`umwm.open_grid`), used for the axes in
    kilometers of a Cartesian grid; otherwise longitude and latitude are used
    when they are set, or cell indices. Direction fields use a cyclic
    colormap, signed fields a diverging one centered on zero, and others a
    sequential one; ``cmap`` and keyword arguments go to ``pcolormesh``.
    Returns the matplotlib axes.
    """
    plt = _pyplot()
    field = dataset[name]
    if "time" in field.dims:
        field = field.isel(time=time)
    if "z" in field.dims:
        field = field.isel(z=level)
    values = np.asarray(field.values, dtype=float)
    if values.ndim != 2:
        raise ValueError(f"{name!r} is not a two-dimensional field; its dimensions are {field.dims}")
    x, y, xlabel, ylabel, aspect = _axes(dataset, grid, values.shape)

    options = dict(shading="auto")
    finite = values[np.isfinite(values)]
    if name in DIRECTION_FIELDS:
        options.update(cmap="twilight", vmin=-np.pi, vmax=np.pi)
    elif finite.size and finite.min() < 0 < finite.max():
        limit = float(np.max(np.abs(finite)))
        options.update(cmap="RdBu_r", vmin=-limit, vmax=limit)
    else:
        options.update(cmap="viridis")
    if cmap is not None:
        options["cmap"] = cmap
    options.update(kwargs)

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4.5))
    mesh = ax.pcolormesh(x, y, values, **options)
    description, units = _fields.FIELDS.get(name, (field.attrs.get("long_name", name), field.attrs.get("units")))
    label = f"{description} [{units}]" if units else description
    ax.figure.colorbar(mesh, ax=ax, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_aspect(aspect)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    if title is None:
        stamp = field.coords["time"].values if "time" in field.coords else dataset.attrs.get("time")
        if stamp is not None:
            stamp = str(np.datetime_as_string(stamp, unit="s")).replace("T", " ") if isinstance(stamp, np.datetime64) else str(stamp)
        title = f"{description}" + (f" at {stamp}" if stamp else "")
    ax.set_title(title, fontsize=10)
    return ax
