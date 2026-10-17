#!/usr/bin/env python
"""Run an idealized hurricane case from Python and summarize the output.

A translating Rankine vortex forces a 510 by 210 km Cartesian domain for six
hours. The forcing files are generated with the package, the model is run
step by step with the maximum wave height printed after each hour, and the
gridded output is then opened with xarray.

    python rankine.py [WORKDIR]
"""

import sys
from pathlib import Path

import numpy as np

import umwm

workdir = Path(sys.argv[1] if len(sys.argv) > 1 else "rankine")

config = umwm.Config()
config.domain.mm, config.domain.nm = 51, 21
config.domain.stoptimestr = "2026-01-01 06:00:00"
config.forcing.winds = True
config.forcing.currents = True
config.forcing_constant.rhoa0 = 1.2
config.output.outgrid = 1
config.output.outspec = 1
config.output.xpl, config.output.ypl = 26, 11
config.stokes.depths = [0.1, 1.0, 5.0, 10.0]
config.spectrum_points = [umwm.SpectrumPoint("center", x=26, y=11)]

domain = umwm.Domain(config, workdir)
umwm.write_rankine_forcing(domain.input_dir, config.start_time, hours=6, nx=51, ny=21,
                           vmax=40.0, rmw=40e3, center_x=-150e3, translation_x=50e3, uc=0.0, vc=0.0)

with domain:
    for time in domain.iterate():
        swh = domain.grid_field("swh")
        print(f"{time:%Y-%m-%d %H:%M}: maximum significant wave height {np.nanmax(swh):.2f} m, "
              f"source time step {domain.attribute('dts'):.1f} s")

    output = domain.output()
    print(output)
    print("\nMaximum significant wave height by hour:")
    print(output.swh.max(dim=("x", "y")).to_pandas())

    spectrum = domain.spectrum_output()["center"]
    peak = spectrum.F.isel(time=-1).argmax(dim=("direction", "frequency"))
    print(f"\nPeak of the spectrum at the center at {str(spectrum.time.values[-1])[:19]}: "
          f"{float(spectrum.frequency[peak['frequency']]):.3f} Hz, "
          f"{np.degrees(float(spectrum.direction[peak['direction']])):.0f} degrees")
