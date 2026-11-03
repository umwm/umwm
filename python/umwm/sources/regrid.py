"""The grid of the model, and interpolation to it from the grid of a dataset."""

from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = ["model_grid", "Regridder", "fill_missing", "open_boundary_cells"]


def model_grid(config, input_dir=None) -> tuple[np.ndarray, np.ndarray]:
    """The longitudes and latitudes of the cells of the model, arrays of shape
    ``(nm, mm)``: the regular lon-lat grid of ``grid.lon0``, ``lat0``, ``dlon``
    and ``dlat``, with longitudes in [-180, 180) as the model makes it, or the
    grid of ``input_dir/umwm.gridtopo`` (or ``umwm.grid``) with
    ``grid.gridfromfile``."""
    domain, grid = config.domain, config.grid
    if grid.gridfromfile:
        from netCDF4 import Dataset

        directory = Path(input_dir) if input_dir is not None else Path("input")
        for name in ("umwm.gridtopo", "umwm.grid"):
            path = directory / name
            if path.is_file():
                with Dataset(path) as dataset:
                    lon = np.asarray(dataset.variables["lon"][:], dtype=np.float64)
                    lat = np.asarray(dataset.variables["lat"][:], dtype=np.float64)
                if lon.shape != (domain.nm, domain.mm):
                    raise ValueError(f"{path}: the grid has shape {lon.shape}, expected "
                                     f"({domain.nm}, {domain.mm}) for nm = {domain.nm} and mm = {domain.mm}")
                return lon, lat
        raise FileNotFoundError(f"grid.gridfromfile is set, but {directory / 'umwm.grid'} is not found")
    if grid.dlon > 0 and grid.dlat > 0:
        lon = np.mod(grid.lon0 + np.arange(domain.mm) * grid.dlon + 180.0, 360.0) - 180.0
        lat = grid.lat0 + np.arange(domain.nm) * grid.dlat
        return np.broadcast_to(lon, (domain.nm, domain.mm)).copy(), \
            np.broadcast_to(lat[:, np.newaxis], (domain.nm, domain.mm)).copy()
    raise ValueError("the grid is not a longitude-latitude grid: set grid.lon0, grid.lat0, grid.dlon and "
                     "grid.dlat, or read it from a file with grid.gridfromfile")


class Regridder:
    """Bilinear interpolation from a regular longitude-latitude grid to points.

    ``source_lon`` (increasing, periodic if it covers the whole circle) and
    ``source_lat`` (increasing or decreasing) are the coordinates of the
    dataset, and ``lon`` and ``lat`` those of the points, arrays of any
    shape. Missing values (NaN) of a field are left out of the interpolation,
    with the weights of the other corners renormalized; a point whose four
    corners are all missing is missing. Points beyond the first or last
    latitude take the values of that latitude.
    """

    def __init__(self, source_lon, source_lat, lon, lat):
        source_lon = np.asarray(source_lon, dtype=np.float64)
        source_lat = np.asarray(source_lat, dtype=np.float64)
        lon = np.asarray(lon, dtype=np.float64)
        lat = np.asarray(lat, dtype=np.float64)
        if lon.shape != lat.shape:
            raise ValueError("lon and lat must have the same shape")
        if source_lon.ndim != 1 or source_lat.ndim != 1 or source_lon.size < 2 or source_lat.size < 2:
            raise ValueError("the source coordinates must be 1-D with at least 2 values")
        if np.any(np.diff(source_lon) <= 0):
            raise ValueError("the source longitudes must increase")
        self.shape = lon.shape
        nx, ny = source_lon.size, source_lat.size

        # longitude: offsets from the first source longitude, wrapped
        offsets = source_lon - source_lon[0]
        last_step = offsets[-1] - offsets[-2]
        periodic = abs(offsets[-1] + last_step - 360.0) < 1e-3 * last_step
        x = np.mod(lon.ravel() - source_lon[0], 360.0)
        i0 = np.clip(np.searchsorted(offsets, x, side="right") - 1, 0, nx - 1)
        i1 = i0 + 1
        span = np.where(i1 < nx, offsets[np.minimum(i1, nx - 1)] - offsets[i0], 360.0 - offsets[-1])
        if periodic:
            i1 = np.where(i1 < nx, i1, 0)
            wx = (x - offsets[i0]) / span
        else:
            inside = (i1 < nx) & (x <= offsets[-1])
            i1 = np.where(inside, i1, i0)
            wx = np.where(inside, (x - offsets[i0]) / np.where(inside, span, 1.0), 0.0)
        wx = np.clip(wx, 0.0, 1.0)

        # latitude: in increasing order, mapped back to the source index
        descending = source_lat[0] > source_lat[-1]
        ascending = source_lat[::-1] if descending else source_lat
        y = np.clip(lat.ravel(), ascending[0], ascending[-1])
        j0 = np.clip(np.searchsorted(ascending, y, side="right") - 1, 0, ny - 2)
        wy = np.clip((y - ascending[j0]) / (ascending[j0 + 1] - ascending[j0]), 0.0, 1.0)
        j1 = j0 + 1
        if descending:
            j0, j1 = ny - 1 - j0, ny - 1 - j1

        self.rows = slice(int(min(j0.min(), j1.min())), int(max(j0.max(), j1.max())) + 1)
        """The latitude index range of the source that the points need."""
        start = self.rows.start
        self._j = np.stack([j0, j0, j1, j1]) - start
        self._i = np.stack([i0, i1, i0, i1])
        self._w = np.stack([(1 - wx) * (1 - wy), wx * (1 - wy), (1 - wx) * wy, wx * wy])

    def corners(self, field) -> tuple[np.ndarray, np.ndarray]:
        """The values of a field, given on the rows :attr:`rows` of the source
        grid or on the whole grid, at the four corners of each point, and their
        bilinear weights, arrays of shape ``(4, points)``."""
        field = np.asarray(field)
        if field.shape[0] != self.rows.stop - self.rows.start:
            field = field[self.rows]
        return field[self._j, self._i], self._w

    def coverage(self, field) -> np.ndarray:
        """The sum of the weights of the corners where the field is not missing,
        between 0 and 1, of the shape of the points."""
        values, weights = self.corners(field)
        return (weights * np.isfinite(values)).sum(axis=0).reshape(self.shape)

    def __call__(self, field) -> np.ndarray:
        """The field interpolated to the points, NaN where it is missing at all corners."""
        values, weights = self.corners(field)
        valid = np.isfinite(values)
        weights = weights * valid
        total = weights.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            result = (weights * np.where(valid, values, 0.0)).sum(axis=0) / total
        result[total <= 0] = np.nan
        return result.reshape(self.shape)


def fill_missing(field, default: float, iterations: int = 4) -> np.ndarray:
    """The field with its missing values (NaN) filled with the mean of their
    neighbors that are not missing, repeatedly, and then with ``default``."""
    field = np.array(field, dtype=np.float64)
    for _ in range(iterations):
        missing = ~np.isfinite(field)
        if not missing.any():
            break
        padded = np.pad(field, 1, mode="edge")
        neighbors = np.stack([padded[:-2, 1:-1], padded[2:, 1:-1], padded[1:-1, :-2], padded[1:-1, 2:]])
        valid = np.isfinite(neighbors)
        count = valid.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(valid, neighbors, 0.0).sum(axis=0) / count
        field = np.where(missing & (count > 0), mean, field)
    return np.where(np.isfinite(field), field, default)


def open_boundary_cells(config, depth) -> np.ndarray:
    """Whether each cell of the grid, an array of shape ``(nm, mm)``, is an
    open boundary cell of the model: a sea cell (positive depth) next to the
    edges of a regional domain, which are land. A global domain has none."""
    nm, mm = config.domain.nm, config.domain.mm
    ring = np.zeros((nm, mm), dtype=bool)
    if config.domain.isglobal:
        return ring
    ring[1, 1:mm - 1] = ring[nm - 2, 1:mm - 1] = True
    ring[1:nm - 1, 1] = ring[1:nm - 1, mm - 2] = True
    return ring & (np.asarray(depth) > 0)
