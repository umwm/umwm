"""The interface of a data source.

A source reads fields of its dataset on its regular longitude-latitude grid,
in the standard names, units and conventions below, at the times of its
dataset. The preparation of the input files (:mod:`umwm.sources.prepare`)
interpolates them to the grid and times of the model, and writes them.

Standard fields, each an array of shape ``(latitude, longitude)``, with NaN
where the dataset has no data, such as over land:

========== ===================================================================
``u10``    eastward wind at 10 m [m/s]
``v10``    northward wind at 10 m [m/s]
``rhoa``   air density [kg/m^3]
``fice``   sea ice fraction [0-1]
``depth``  water depth [m], positive
``hs``     significant wave height of the whole wave field [m]
``tm``     its mean period Tm-1,0 = m-1/m0 [s]
``dir``    its mean direction, where the waves come from, clockwise from
           north [degrees], the convention of most datasets
``spr``    its directional spread sqrt(2 (1 - M1)) [rad]
``hs<p>``  significant wave height of wave partition p [m], 0 where the
           partition is absent; partition 0 is the wind sea, the others the
           swell, as many as :attr:`Source.partitions`
``tm<p>``  its mean period Tm-1,0 [s], or its peak period ``tp<p>`` if the
           :attr:`Source.period` of the source is ``"tp"``
``dir<p>`` its mean direction, as ``dir`` [degrees]
``spr<p>`` its directional spread [rad], optional
========== ===================================================================

The wind, air density, sea ice, and wave fields serve the switches
``winds``, ``air_density``, ``seaice`` and ``waves`` of the ``FORCING``
namelist; ``depth`` serves the bathymetry of ``grid.topofromfile``. A
source does not need to provide them all (see
:data:`umwm.sources.SOURCES`).
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np

__all__ = ["Source", "FORCING_FIELDS", "wave_fields"]

FORCING_FIELDS = {
    "winds": ("u10", "v10"),
    "air_density": ("rhoa",),
    "seaice": ("fice",),
}
"""The standard fields of the forcing switches, other than the waves."""


def wave_fields(partitions: int, period: str = "tm", total: bool = True) -> list[str]:
    """The standard names of the wave fields of a source with that many
    partitions and partition period (``"tm"`` or ``"tp"``), with the fields
    of the whole wave field if ``total``."""
    names = ["hs", "tm", "dir", "spr"] if total else []
    for p in range(partitions):
        names += [f"hs{p}", f"{period}{p}", f"dir{p}", f"spr{p}"]
    return names


class Source:
    """A data source of the input files of the model.

    Subclasses set the attributes and implement :meth:`grid`,
    :meth:`coverage` and :meth:`read`; the source may be used as a context
    manager, which closes it.
    """

    name: str = ""
    """The name of the source, as in ``forcing.source``."""

    time_step: timedelta = timedelta(hours=1)
    """The time step of the dataset; fields at other times are interpolated."""

    partitions: int = 0
    """The number of wave partitions: the wind sea (0) and the swell."""

    period: str = "tm"
    """The period of the partitions: mean period Tm-1,0 ``"tm"``, or peak period ``"tp"``."""

    total: bool = True
    """Whether the source provides the fields of the whole wave field, ``hs``,
    ``tm``, ``dir`` and ``spr``, to which the partitions are adjusted."""

    def description(self) -> str:
        """The dataset, for the attributes of the input files."""
        return self.name

    def grid(self) -> tuple[np.ndarray, np.ndarray]:
        """The longitudes (increasing, in degrees east) and the latitudes (increasing
        or decreasing, in degrees north) of the regular grid of the dataset."""
        raise NotImplementedError

    def coverage(self) -> tuple[datetime, datetime]:
        """The first and last times of the dataset."""
        raise NotImplementedError

    def read(self, names, times, rows=None) -> dict:
        """The standard fields ``names`` at the dataset times ``times``, as
        ``{time: {name: array}}``, with the arrays restricted to the latitude
        index range ``rows`` (a slice) if given."""
        raise NotImplementedError

    def close(self) -> None:
        """Release the resources of the source."""

    def __enter__(self):
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()
