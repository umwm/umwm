"""Data sources of the input files: reanalyses and forecasts, such as ERA5.

The model reads its input files from the input directory of a run: the grid
and bathymetry ``umwm.gridtopo``, and per global time step the forcing file
``umwmin_<time>.nc`` with the wind, air density, sea ice and, with
``forcing.waves``, the wave partitions from which the model builds its
initial spectrum and the spectrum at the open boundaries of a regional
domain. A data source prepares these files on the grid of the model from a
dataset: with ``forcing.source = 'era5'`` in the namelist, the missing input
files are prepared from ERA5 when the model runs, by ``umwm run``, the
``umwm`` executable, or :meth:`umwm.Domain.initialize`, or beforehand by
``umwm forcing era5`` or :func:`prepare_inputs`.

A source reads its fields on its own grid in the standard names and
conventions of :mod:`umwm.sources.base`; regridding, the partitions of the
wave field, and the files are common to all sources (:mod:`umwm.sources.prepare`),
so that a new source, such as NCEP GFS with its WAVEWATCH III wave fields,
only maps the variables of its dataset to the standard names. The sources
are listed in :data:`SOURCES`, imported only when used.
"""

from __future__ import annotations

import importlib
from dataclasses import dataclass

__all__ = ["SourceInfo", "SOURCES", "available_sources", "source_info", "get_source", "source_problems",
           "prepare_inputs", "SourceError"]


class SourceError(RuntimeError):
    """Raised when the input files cannot be prepared from a data source."""


@dataclass(frozen=True)
class SourceInfo:
    """A data source: its name in ``forcing.source``, a description, the input
    files it provides (``"grid"`` for the bathymetry, and the switches of the
    ``FORCING`` namelist it can serve), the class that reads it, and the
    optional dependencies it needs, as a pip extra of the package."""

    name: str
    description: str
    provides: frozenset
    target: str
    extra: str
    dependencies: tuple = ()

    def load(self):
        """The class of the source, which imports its dependencies."""
        module, _, name = self.target.partition(":")
        try:
            return getattr(importlib.import_module(module), name)
        except ImportError as error:
            raise SourceError(f"the {self.name} data source needs {', '.join(self.dependencies)}: "
                              f"pip install 'umwm[{self.extra}]' (or pip install "
                              f"{' '.join(self.dependencies)}); {error}") from error


SOURCES = {
    "era5": SourceInfo(
        name="era5",
        description="ECMWF ERA5 reanalysis, hourly from 1940, from the analysis-ready, cloud-optimized "
                    "(ARCO) Zarr store of Google Cloud Public Datasets",
        provides=frozenset({"grid", "winds", "air_density", "seaice", "waves"}),
        target="umwm.sources.era5:ERA5",
        extra="era5",
        dependencies=("zarr", "fsspec", "aiohttp"),
    ),
}
"""The data sources by name."""


def available_sources() -> list[str]:
    """The names of the data sources."""
    return list(SOURCES)


def source_info(name: str) -> SourceInfo:
    """The data source of the given name."""
    try:
        return SOURCES[name.strip().lower()]
    except KeyError:
        raise SourceError(f"unknown data source {name!r}; the data sources are: "
                          f"{', '.join(SOURCES)}") from None


def get_source(name: str, **options):
    """An open data source of the given name; ``options`` are those of its class,
    such as the ``url`` of its dataset."""
    return source_info(name).load()(**options)


def source_problems(config, name: str | None = None) -> list[str]:
    """The reasons why the input files of a configuration cannot be prepared
    from its data source (``forcing.source``, or ``name``): an unknown source,
    a grid that is not a longitude-latitude grid, or input files that the source
    does not provide."""
    name = name or config.source_name
    if name is None:
        return []
    if name not in SOURCES:
        return [f"unknown data source {name!r} in forcing.source; the data sources are: {', '.join(SOURCES)}"]
    info = SOURCES[name]
    problems = []
    if not (config.grid.gridfromfile or config.lonlat_grid):
        problems.append(f"the data source {name} needs a longitude-latitude grid: set grid.lon0, grid.lat0, "
                        "grid.dlon and grid.dlat, or read the grid from a file with grid.gridfromfile")
    missing = [switch for switch in ("winds", "currents", "air_density", "water_density", "seaice", "waves")
               if getattr(config.forcing, switch) and switch not in info.provides]
    if missing:
        problems.append(f"the data source {name} does not provide {', '.join(missing)}; set "
                        + " and ".join(f"forcing.{switch}" for switch in missing)
                        + " to .false. to use the constant values of forcing_constant")
    return problems


def prepare_inputs(config, input_dir, source: str | None = None, **options):
    """Prepare the missing input files of a configuration from a data source;
    see :func:`umwm.sources.prepare.prepare_inputs`."""
    from .prepare import prepare_inputs as prepare

    return prepare(config, input_dir, source=source, **options)
