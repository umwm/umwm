"""Access to UMWM output files with xarray.

The model writes gridded output to ``output/umwmout_<time>.nc``, one file per
output time, spectrum output to ``output/umwmspc_<point>_<start time>.nc``,
one file per point with a time series, and restart files to
``restart/umwmrst_<time>.nc``. All are opened with the netCDF4 engine.
"""

from __future__ import annotations

import glob
import re
from pathlib import Path

import xarray as xr

__all__ = ["output_files", "open_output", "spectrum_output_files", "open_spectrum_output",
           "restart_files", "open_restart", "open_grid"]

_GRID_OUTPUT = "umwmout_*.nc"
_SPECTRUM_OUTPUT = "umwmspc_*.nc"
_RESTART = "umwmrst_*.nc"
_SPECTRUM_NAME = re.compile(r"umwmspc_(?P<point>.+)_(?P<time>\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2})\.nc$")


def _directory(path, subdirectory):
    """The output (or restart) directory of a run directory, or the directory itself."""
    path = Path(path)
    if path.is_dir() and (path / subdirectory).is_dir() and not any(path.glob("umwm*.nc")):
        return path / subdirectory
    return path


def _files(source, pattern, subdirectory) -> list[Path]:
    if isinstance(source, (str, Path)):
        path = Path(source)
        if path.is_dir():
            return sorted(_directory(path, subdirectory).glob(pattern))
        if path.is_file():
            return [path]
        files = sorted(Path(match) for match in glob.glob(str(path)))
        if not files:
            raise FileNotFoundError(f"no files match {source}")
        return files
    return [Path(item) for item in source]


def output_files(source) -> list[Path]:
    """The gridded output files of a run or output directory, a glob, or a list of files."""
    return _files(source, _GRID_OUTPUT, "output")


def restart_files(source) -> list[Path]:
    """The restart files of a run or restart directory, a glob, or a list of files."""
    return _files(source, _RESTART, "restart")


def spectrum_output_files(source) -> dict[str, list[Path]]:
    """The spectrum output files of a run or output directory, by point name."""
    files: dict[str, list[Path]] = {}
    for path in _files(source, _SPECTRUM_OUTPUT, "output"):
        match = _SPECTRUM_NAME.match(path.name)
        name = match.group("point") if match else path.stem
        files.setdefault(name, []).append(path)
    return files


def _open_along_time(files, **kwargs) -> xr.Dataset:
    """Open files and concatenate them along time, lazily with dask when available."""
    if not files:
        raise FileNotFoundError("no output files found")
    kwargs.setdefault("engine", "netcdf4")
    try:
        import dask  # noqa: F401
    except ImportError:
        datasets = [xr.open_dataset(path, **kwargs) for path in files]
        if len(datasets) == 1:
            return datasets[0]
        return xr.concat(datasets, dim="time", data_vars="minimal", coords="minimal", compat="override")
    return xr.open_mfdataset(files, combine="nested", concat_dim="time", data_vars="minimal",
                             coords="minimal", compat="override", **kwargs)


def open_output(source, **kwargs) -> xr.Dataset:
    """Open the gridded output as one dataset along ``time``.

    ``source`` is a run directory, an output directory, a glob pattern, a
    single file, or a list of files. Keyword arguments are passed to xarray.
    Variables are dimensioned ``(time, y, x)``, Stokes drift components
    ``(time, z, y, x)``; land cells are NaN.
    """
    return _open_along_time(output_files(source), **kwargs)


def open_spectrum_output(source, **kwargs) -> dict[str, xr.Dataset]:
    """Open the spectrum output as a dataset per point name, along ``time``.

    The spectrum ``F`` and the source terms are dimensioned
    ``(time, direction, frequency)``.
    """
    return {name: _open_along_time(files, **kwargs) for name, files in spectrum_output_files(source).items()}


def open_restart(path, **kwargs) -> xr.Dataset:
    """Open a restart file; the spectrum ``F`` is dimensioned ``(x, th, f)`` over the sea cells."""
    kwargs.setdefault("engine", "netcdf4")
    return xr.open_dataset(path, **kwargs)


def open_grid(source, **kwargs) -> xr.Dataset:
    """Open the grid definition file ``output/umwmout.grid`` of a run or output directory."""
    path = Path(source)
    if path.is_dir():
        path = _directory(path, "output") / "umwmout.grid"
    kwargs.setdefault("engine", "netcdf4")
    return xr.open_dataset(path, **kwargs)
