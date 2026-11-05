"""ERA5 from the analysis-ready, cloud-optimized (ARCO) store of Google.

ERA5 (Hersbach et al., 2020) is the hourly reanalysis of ECMWF, from 1940,
with the wave model of the IFS, which partitions its wave spectrum into the
wind sea and up to three swell systems. Google Cloud Public Datasets serves
it as a Zarr store on a regular 0.25 degree grid (Carver and Merose, 2023),
read anonymously over HTTPS; each variable is a chunk per hour of the whole
globe, of 1.5 to 3.5 MB.

The source provides the wind at 10 m, the air density over the oceans, the
sea ice fraction, the bathymetry of the wave model, and the wave field: its
significant wave height, mean period Tm-1,0, mean direction and directional
spread, and the partitions: the wind sea, and the first, second and third
swell, with their significant wave height, mean period Tm-1,0 and mean
direction, and the directional spread of the wind sea and of the swell at
their spectral peak. The mean periods and directions are defined in the
documentation of the wave parameters of ECMWF (Bidlot, 2020).

The location of the store is given by the ``url`` of :class:`ERA5`, the
environment variable ``UMWM_ERA5_URL``, or else :data:`ERA5_URL`; it may be
a local copy of the store, or of the variables and times of a run.
"""

from __future__ import annotations

import os
import re
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

from . import SourceError
from .base import Source

__all__ = ["ERA5", "ERA5_URL", "URL_VARIABLE", "VARIABLES"]

ERA5_URL = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
"""The ARCO ERA5 store of Google Cloud Public Datasets, hourly at 0.25 degrees on 37 levels."""

URL_VARIABLE = "UMWM_ERA5_URL"
"""The environment variable that sets another location of the store."""

VARIABLES = {
    "u10": "10m_u_component_of_wind",
    "v10": "10m_v_component_of_wind",
    "rhoa": "air_density_over_the_oceans",
    "fice": "sea_ice_cover",
    "depth": "model_bathymetry",
    "hs": "significant_height_of_combined_wind_waves_and_swell",
    "tm": "mean_wave_period",
    "dir": "mean_wave_direction",
    "spr": "wave_spectral_directional_width",
    "hs0": "significant_height_of_wind_waves",
    "tm0": "mean_period_of_wind_waves",
    "dir0": "mean_direction_of_wind_waves",
    "spr0": "wave_spectral_directional_width_for_wind_waves",
}
"""The ERA5 variables of the standard fields (see :mod:`umwm.sources.base`)."""

for _p, _ordinal in ((1, "first"), (2, "second"), (3, "third")):
    VARIABLES[f"hs{_p}"] = f"significant_wave_height_of_{_ordinal}_swell_partition"
    VARIABLES[f"tm{_p}"] = f"mean_wave_period_of_{_ordinal}_swell_partition"
    VARIABLES[f"dir{_p}"] = f"mean_wave_direction_of_{_ordinal}_swell_partition"
    VARIABLES[f"spr{_p}"] = "wave_spectral_directional_width_for_swell"

DEFINED_EVERYWHERE = frozenset({VARIABLES["u10"], VARIABLES["v10"]})
"""The variables without missing values, which tell the times without data."""

_TIME_UNITS = re.compile(r"^\s*hours\s+since\s+(\d{4}-\d{2}-\d{2})(?:[ T](\d{1,2}:\d{2}(?::\d{2})?))?")


class ERA5(Source):
    """The ERA5 reanalysis in the ARCO Zarr store at ``url`` (by default
    ``$UMWM_ERA5_URL`` or :data:`ERA5_URL`), read with ``workers`` concurrent
    requests. ``storage_options`` are passed to fsspec."""

    name = "era5"
    time_step = timedelta(hours=1)
    partitions = 4
    period = "tm"
    total = True

    def __init__(self, url=None, workers: int = 16, storage_options=None):
        import zarr

        self.url = str(url or os.environ.get(URL_VARIABLE) or ERA5_URL)
        self.workers = max(1, int(workers))
        options = dict(storage_options or {})
        if self.url.startswith(("http://", "https://")):
            # through the proxy of the environment, if there is one
            options.setdefault("client_kwargs", {"trust_env": True})
        elif self.url.startswith("gs://"):
            options.setdefault("token", "anon")
        try:
            if "://" in self.url:
                self._group = zarr.open_consolidated(self.url, mode="r", storage_options=options)
            else:
                if not Path(self.url).exists():
                    raise FileNotFoundError(self.url)
                self._group = zarr.open_consolidated(self.url, mode="r")
        except Exception as error:  # noqa: BLE001 - any failure to reach the store
            raise SourceError(f"cannot open the ERA5 store {self.url}: {type(error).__name__}: {error}") from error

        time = self._group["time"]
        match = _TIME_UNITS.match(str(time.attrs.get("units", "")))
        if not match:
            raise SourceError(f"{self.url}: the time units {time.attrs.get('units')!r} are not hours since a date")
        self._origin = datetime.fromisoformat(match.group(1) + (" " + match.group(2) if match.group(2) else ""))
        values = np.asarray(time[:2], dtype=np.int64)
        if values.size == 2 and values[1] - values[0] != 1:
            raise SourceError(f"{self.url}: the times are not hourly")
        self._first = int(values[0])
        self._size = int(time.shape[0])
        self._lon = np.asarray(self._group["longitude"][:], dtype=np.float64)
        self._lat = np.asarray(self._group["latitude"][:], dtype=np.float64)

    def description(self) -> str:
        return f"ERA5 reanalysis of ECMWF, from {self.url}"

    def grid(self):
        return self._lon, self._lat

    def coverage(self):
        """The first and last times with data: those of the attributes of the
        store, the last including the preliminary ERA5T, or those of its
        time axis, which may extend beyond the data."""
        attrs = self._group.attrs
        first = self._time(0)
        last = self._time(self._size - 1)
        start = attrs.get("valid_time_start")
        stop = str(attrs.get("valid_time_stop_era5t") or attrs.get("valid_time_stop") or "")
        if start:
            first = max(first, datetime.fromisoformat(str(start)))
        if stop:
            # a date is the whole day
            last = min(last, datetime.fromisoformat(stop) + (timedelta(hours=23) if len(stop) <= 10 else timedelta()))
        return first, last

    def _time(self, index: int) -> datetime:
        return self._origin + timedelta(hours=self._first + index)

    def _index(self, time: datetime) -> int:
        hours = (time - self._origin) / timedelta(hours=1)
        if hours != int(hours):
            raise SourceError(f"ERA5 is hourly; {time} is not on the hour")
        index = int(hours) - self._first
        if not 0 <= index < self._size:
            raise SourceError(f"{time} is not in the ERA5 store {self.url}")
        return index

    def read(self, names, times, rows=None):
        names = list(names)
        unknown = [name for name in names if name not in VARIABLES]
        if unknown:
            raise SourceError(f"ERA5 does not provide {', '.join(unknown)}")
        variables = sorted({VARIABLES[name] for name in names})
        rows = rows or slice(None)
        requests = [(variable, time) for time in times for variable in variables]

        def fetch(request):
            variable, time = request
            values = np.asarray(self._group[variable][self._index(time), rows, :], dtype=np.float32)
            if variable in DEFINED_EVERYWHERE and not np.isfinite(values).any():
                first, last = self.coverage()
                raise SourceError(f"no ERA5 data of {variable} at {time:%Y-%m-%d %H:%M} in {self.url}, "
                                  f"which covers {first:%Y-%m-%d %H:%M} to {last:%Y-%m-%d %H:%M}")
            return values

        with ThreadPoolExecutor(max_workers=self.workers) as executor:
            arrays = dict(zip(requests, executor.map(fetch, requests)))
        return {time: {name: arrays[(VARIABLES[name], time)] for name in names} for time in times}
