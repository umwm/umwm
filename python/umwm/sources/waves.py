"""Wave partitions: from the integrated parameters of a dataset to the
partition fields of the input files of the model.

The model rebuilds its spectrum from the partitions of its input files
(``forcing.waves``): for each partition p, the significant wave height
``phs<p>``, the mean period Tm-1,0 ``ptm<p>`` (or the peak period
``ptp<p>``), the mean direction ``pdir<p>`` in the convention of the model,
the direction of propagation counterclockwise from east, in radians, and the
directional spread ``pspr<p>``, in radians. ``phs<p>`` is -1 at cells
without wave data, such as land or sea ice in the dataset.

The partitions of a dataset (the wind sea and the swell systems) do not
always add up to the significant wave height of the whole wave field, as the
datasets interpolate each of them to their grid; where the whole wave field
has more energy than the partitions, the difference is added as a residual
partition with the mean period, direction and spread of the whole wave
field, and where it has less, the partitions are scaled down to it, so that
the model starts from the significant wave height of the dataset.
"""

from __future__ import annotations

import numpy as np

__all__ = ["Partition", "standard_partitions", "regrid_partitions", "to_model_direction", "PARTITION_VARIABLES"]

PARTITION_VARIABLES = {
    "phs": ("significant wave height of partition {p}; -1 where there are no wave data", "m"),
    "ptm": ("mean period Tm-1,0 of partition {p}", "s"),
    "ptp": ("peak period of partition {p}", "s"),
    "pdir": ("mean direction of propagation of partition {p}, counterclockwise from east", "rad"),
    "pspr": ("directional spread sqrt(2 (1 - M1)) of partition {p}", "rad"),
    "pgam": ("JONSWAP peak enhancement factor of partition {p}", "1"),
}
"""The partition variables of the input files by prefix, with their description and units."""


class Partition:
    """One wave partition on a grid: its ``energy`` Hs^2, ``period``, mean
    ``direction`` in the convention of the model, and ``spread``; arrays of
    the same shape, the energy 0 where the partition is absent and NaN where
    there are no wave data. ``kind`` is the kind of period: mean period
    Tm-1,0 ``"tm"`` or peak period ``"tp"``."""

    def __init__(self, energy, period, direction, spread, kind="tm"):
        self.energy = energy
        self.period = period
        self.direction = direction
        self.spread = spread
        self.kind = kind


def to_model_direction(degrees):
    """The direction of propagation counterclockwise from east in radians, in
    (-pi, pi], of a direction where the waves come from, clockwise from north,
    in degrees."""
    radians = np.deg2rad(270.0 - np.asarray(degrees, dtype=np.float64))
    return np.angle(np.exp(1j * radians))


def standard_partitions(fields, partitions: int, period: str = "tm", total: bool = True) -> list[Partition]:
    """The partitions of the standard wave fields of a source on its grid (see
    :mod:`umwm.sources.base`), with the residual partition of the whole wave
    field if ``total``, and the partitions scaled down where they have more
    energy than it."""

    def get(name):
        return np.asarray(fields[name], dtype=np.float64)

    result = []
    for p in range(partitions):
        hs = get(f"hs{p}")
        t = get(f"{period}{p}")
        direction = get(f"dir{p}")
        spread = get(f"spr{p}") if f"spr{p}" in fields else np.zeros_like(hs)
        present = np.isfinite(hs) & (hs > 0) & np.isfinite(t) & (t > 0) & np.isfinite(direction)
        energy = np.where(present, hs**2, np.where(np.isfinite(hs), 0.0, np.nan))
        result.append(Partition(energy, np.where(present, t, 0.0), np.where(present, to_model_direction(direction), 0.0),
                                np.where(present & np.isfinite(spread), spread, 0.0), kind=period))

    if total and "hs" in fields:
        hs = get("hs")
        target = np.where(np.isfinite(hs), hs**2, np.nan)
        partial = np.zeros_like(target)
        for partition in result:
            partial += np.nan_to_num(partition.energy)
        with np.errstate(invalid="ignore", divide="ignore"):
            scale = np.where((partial > target) & (partial > 0), target / partial, 1.0)
        for partition in result:
            # partitions missing where the whole wave field has data are absent
            partition.energy = np.where(np.isfinite(target), np.nan_to_num(partition.energy) * scale, np.nan)
        residual = np.where(np.isfinite(target), np.maximum(target - partial, 0.0), np.nan)
        t, direction = get("tm"), get("dir")
        spread = get("spr") if "spr" in fields else np.zeros_like(t)
        present = (np.nan_to_num(residual) > 0) & np.isfinite(t) & (t > 0) & np.isfinite(direction)
        result.append(Partition(np.where(present | ~np.isfinite(residual), residual, 0.0),
                                np.where(present, t, 0.0), np.where(present, to_model_direction(direction), 0.0),
                                np.where(present & np.isfinite(spread), spread, 0.0)))
    return result


def regrid_partitions(regridder, partitions) -> dict:
    """The partition variables of the input files on the points of the
    regridder, from partitions on the grid of the source: the variance is
    interpolated, and the period, direction (as a vector) and spread weighted
    by the variance at each corner. ``phs<p>`` is -1 where the four corners
    have no wave data."""
    variables = {}
    for p, partition in enumerate(partitions):
        energy, weights = regridder.corners(partition.energy)
        valid = np.isfinite(energy)
        weights = weights * valid
        coverage = weights.sum(axis=0)
        energy = np.where(valid, energy, 0.0)
        weighted = weights * energy
        total = weighted.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            hs = np.sqrt(total / coverage)
            mean_period = (weighted * np.nan_to_num(regridder.corners(partition.period)[0])).sum(axis=0) / total
            spread = (weighted * np.nan_to_num(regridder.corners(partition.spread)[0])).sum(axis=0) / total
        direction = regridder.corners(partition.direction)[0]
        x = (weighted * np.cos(np.nan_to_num(direction))).sum(axis=0)
        y = (weighted * np.sin(np.nan_to_num(direction))).sum(axis=0)
        present = total > 0
        shape = regridder.shape
        variables[f"phs{p}"] = np.where(coverage > 0, np.where(present, hs, 0.0), -1.0).reshape(shape)
        variables[f"p{partition.kind}{p}"] = np.where(present, mean_period, 0.0).reshape(shape)
        variables[f"pdir{p}"] = np.where(present, np.arctan2(y, x), 0.0).reshape(shape)
        variables[f"pspr{p}"] = np.where(present, spread, 0.0).reshape(shape)
    return variables
