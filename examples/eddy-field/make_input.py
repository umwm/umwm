#!/usr/bin/env python
"""Write the steady currents of the eddy-field case: input/umwmin_<start time>.nc.

A synthetic field of mesoscale eddies, cyclones and anticyclones of
Gaussian shape, placed at random on the grid of main.nml next to this
script, periodic in x like the domain. Each eddy has the streamfunction

    psi = V * R * exp(1/2) * exp(-r**2 / (2 * R**2)),

whose azimuthal velocity is largest, V, at the radius R. The radii are
log-normally distributed around RADIUS, the speeds V uniformly between
SPEEDS, the senses of rotation random, and the eddies are kept apart by
SEPARATION times the sum of their radii, with their centers at least 1.5
radii beyond COAST from the coast at y = 0. The currents u = -dpsi/dy and
v = dpsi/dx are non-divergent. The model reads the currents once and holds them constant
(forcing.steady). Change the parameters below and run the script again to
vary the eddy field:

    python make_input.py
"""

from pathlib import Path

import numpy as np

import umwm

SEED = 5             # seed of the random eddy field
NUMBER = 14          # number of eddies
RADIUS = 25e3        # median radius of the largest speed [m]
SPREAD = 0.3         # standard deviation of the logarithm of the radius
SPEEDS = (0.1, 0.5)  # range of the largest speeds [m/s]
SEPARATION = 1.2     # smallest distance between eddy centers, in sums of radii
COAST = 60e3         # distance from the coast at y = 0 where the eddy field begins [m]


def eddy_field(x, y, rng):
    """The currents (uc, vc) of a random field of Gaussian eddies on the grid x, y,
    periodic in x, and the eddies as (x, y, radius, speed) tuples."""
    X, Y = np.meshgrid(x, y)
    dx = x[1] - x[0]
    width = x.size * dx
    psi = np.zeros(X.shape)
    eddies = []
    while len(eddies) < NUMBER:
        r = float(np.clip(RADIUS * np.exp(rng.normal(0, SPREAD)), RADIUS / 2, 2 * RADIUS))
        cx, cy = rng.uniform(0, width), rng.uniform(COAST + 1.5 * r, y[-1] - r)
        if any(np.hypot(min(abs(cx - ex), width - abs(cx - ex)), cy - ey) < SEPARATION * (r + er)
               for ex, ey, er, _ in eddies):
            continue
        speed = rng.uniform(*SPEEDS) * rng.choice([-1, 1])  # positive: clockwise
        eddies.append((cx, cy, r, speed))
        for shift in (-width, 0.0, width):  # the eddy and its periodic images
            psi += speed * r * np.exp(0.5) * np.exp(-((X - cx - shift) ** 2 + (Y - cy) ** 2) / (2 * r**2))
    uc = -np.gradient(psi, dx, axis=0)
    vc = (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) / (2 * dx)
    return uc, vc, eddies


here = Path(__file__).resolve().parent
config = umwm.Config.from_namelist(here / "main.nml")
x = np.arange(config.domain.mm) * config.grid.delx
y = np.arange(config.domain.nm) * config.grid.dely
uc, vc, eddies = eddy_field(x, y, np.random.default_rng(SEED))
speed = np.hypot(uc, vc)
path = umwm.write_forcing(here / "input", config.start_time, uc=uc, vc=vc, x=x, y=y,
                          attrs={"description": f"{len(eddies)} Gaussian eddies, seed {SEED}"})
print(f"wrote {path}: {len(eddies)} eddies ({sum(e[3] < 0 for e in eddies)} counterclockwise), "
      f"radii {min(e[2] for e in eddies) / 1e3:.0f} to {max(e[2] for e in eddies) / 1e3:.0f} km, "
      f"largest speed {speed.max():.2f} m/s, rms speed {np.sqrt(np.mean(speed[y > COAST] ** 2)):.2f} m/s")
