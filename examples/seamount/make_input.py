#!/usr/bin/env python
"""Write the bathymetry of the seamount case: input/umwm.gridtopo.

A Gaussian seamount rises from the deep floor at the center of the domain,
on the grid of main.nml next to this script:

    depth = DEPTH - (DEPTH - TOP) * exp(-r**2 / (2 * SIGMA**2))

The model reads the bottom elevation z (negative below sea level) from the
file. Change the parameters below and run the script again to vary the
seamount:

    python make_input.py
"""

from pathlib import Path

import numpy as np
import xarray as xr

import umwm

DEPTH = 1000.0   # depth of the floor around the seamount [m]
TOP = 20.0       # depth of the top of the seamount [m]
SIGMA = 20e3     # width of the seamount, the standard deviation of the Gaussian [m]

here = Path(__file__).resolve().parent
config = umwm.Config.from_namelist(here / "main.nml")
domain = config.domain
x = np.arange(domain.mm) * config.grid.delx
y = np.arange(domain.nm) * config.grid.dely
xc, yc = x[domain.mm // 2], y[domain.nm // 2]
r2 = (x[np.newaxis, :] - xc) ** 2 + (y[:, np.newaxis] - yc) ** 2
depth = DEPTH - (DEPTH - TOP) * np.exp(-r2 / (2 * SIGMA**2))

dataset = xr.Dataset(
    {"z": (("y", "x"), (-depth).astype(np.float32),
           {"long_name": "bottom elevation (negative below sea level)", "units": "m"})},
    coords={"x": ("x", x.astype(np.float32), {"units": "m"}), "y": ("y", y.astype(np.float32), {"units": "m"})},
    attrs={"description": f"Gaussian seamount at x = {xc / 1e3:g} km, y = {yc / 1e3:g} km, from {DEPTH:g} m "
                          f"to {TOP:g} m depth, sigma = {SIGMA / 1e3:g} km"},
)
path = here / "input" / "umwm.gridtopo"
path.parent.mkdir(exist_ok=True)
dataset.to_netcdf(path, engine="netcdf4", format="NETCDF4", encoding={"z": {"zlib": True}})
print(f"wrote {path}: depth {depth.min():.0f} to {depth.max():.0f} m")
