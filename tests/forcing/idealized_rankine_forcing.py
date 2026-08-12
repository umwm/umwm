#!/usr/bin/env python

"""Generate idealized UMWM forcing files with a Rankine hurricane vortex."""

import argparse
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create hourly idealized UMWM forcing files."
    )
    parser.add_argument(
        "--output-dir",
        default="input",
        help="directory for generated umwmin_*.nc files",
    )
    parser.add_argument(
        "--start-time",
        default="2026-01-01T00:00:00",
        help="first forcing time, e.g. 2026-01-01T00:00:00",
    )
    parser.add_argument(
        "--last-hour",
        type=int,
        default=6,
        help="last hourly file offset to write, inclusive",
    )
    parser.add_argument("--nx", type=int, default=51, help="number of x/lon cells")
    parser.add_argument("--ny", type=int, default=21, help="number of y/lat cells")
    parser.add_argument("--dx-m", type=float, default=10000.0, help="x grid spacing")
    parser.add_argument("--dy-m", type=float, default=10000.0, help="y grid spacing")
    parser.add_argument(
        "--vmax",
        type=float,
        default=30.0,
        help="maximum tangential wind speed in m/s",
    )
    parser.add_argument(
        "--rmw-m",
        type=float,
        default=50000.0,
        help="radius of maximum wind in meters",
    )
    parser.add_argument("--center-lon", type=float, default=0.0)
    parser.add_argument("--center-lat", type=float, default=0.0)
    parser.add_argument(
        "--initial-center-x-m",
        type=float,
        default=-150000.0,
        help="vortex center x-position at the first forcing time",
    )
    parser.add_argument(
        "--initial-center-y-m",
        type=float,
        default=0.0,
        help="vortex center y-position at the first forcing time",
    )
    parser.add_argument(
        "--translation-x-m-per-hour",
        type=float,
        default=50000.0,
        help="eastward vortex translation distance per hour",
    )
    parser.add_argument(
        "--translation-y-m-per-hour",
        type=float,
        default=0.0,
        help="northward vortex translation distance per hour",
    )
    parser.add_argument(
        "--clockwise",
        action="store_true",
        help="make the vortex rotate clockwise instead of counter-clockwise",
    )
    parser.add_argument("--uc", type=float, default=0.0, help="current x-component")
    parser.add_argument("--vc", type=float, default=0.0, help="current y-component")
    parser.add_argument("--rhoa", type=float, default=1.2, help="air density")
    parser.add_argument("--rhow", type=float, default=1030.0, help="water density")
    return parser.parse_args()


def centered_axis(num_cells, spacing_m):
    return (np.arange(num_cells, dtype=np.float32) - (num_cells - 1) / 2) * spacing_m


def lon_lat_from_xy(x2d, y2d, center_lon, center_lat):
    meters_per_degree_lat = 111_320.0
    coslat = np.cos(np.deg2rad(center_lat))
    if abs(coslat) < 1.0e-6:
        raise ValueError("center-lat is too close to a pole for lon metadata")

    lat = center_lat + y2d / meters_per_degree_lat
    lon = center_lon + x2d / (meters_per_degree_lat * coslat)
    return lon.astype(np.float32), lat.astype(np.float32)


def rankine_vortex(x2d, y2d, vmax, rmw_m, clockwise=False):
    radius = np.hypot(x2d, y2d)
    if rmw_m <= 0:
        raise ValueError("rmw-m must be > 0")

    wspd = np.empty_like(radius, dtype=np.float64)
    inside = radius <= rmw_m
    outside = ~inside
    wspd[inside] = vmax * radius[inside] / rmw_m
    wspd[outside] = vmax * rmw_m / radius[outside]
    wspd[radius == 0.0] = 0.0

    uw = np.zeros_like(radius, dtype=np.float64)
    vw = np.zeros_like(radius, dtype=np.float64)
    nonzero = radius > 0.0
    if clockwise:
        uw[nonzero] = wspd[nonzero] * y2d[nonzero] / radius[nonzero]
        vw[nonzero] = -wspd[nonzero] * x2d[nonzero] / radius[nonzero]
    else:
        uw[nonzero] = -wspd[nonzero] * y2d[nonzero] / radius[nonzero]
        vw[nonzero] = wspd[nonzero] * x2d[nonzero] / radius[nonzero]

    return (
        uw.astype(np.float32),
        vw.astype(np.float32),
    )


def add_time_axis(field):
    return field[np.newaxis, :, :]


def vortex_center(args, hour):
    return (
        args.initial_center_x_m + hour * args.translation_x_m_per_hour,
        args.initial_center_y_m + hour * args.translation_y_m_per_hour,
    )


def make_dataset(args, valid_time, hour):
    x = centered_axis(args.nx, args.dx_m)
    y = centered_axis(args.ny, args.dy_m)
    x2d, y2d = np.meshgrid(x, y)
    lon, lat = lon_lat_from_xy(x2d, y2d, args.center_lon, args.center_lat)
    vortex_x, vortex_y = vortex_center(args, hour)
    uw, vw = rankine_vortex(
        x2d - vortex_x, y2d - vortex_y, args.vmax, args.rmw_m, args.clockwise
    )

    shape = (args.ny, args.nx)
    ds = xr.Dataset(
        data_vars={
            "lon": (("y", "x"), lon),
            "lat": (("y", "x"), lat),
            "uw": (("time", "y", "x"), add_time_axis(uw)),
            "vw": (("time", "y", "x"), add_time_axis(vw)),
            "uc": (
                ("time", "y", "x"),
                add_time_axis(np.full(shape, args.uc, dtype=np.float32)),
            ),
            "vc": (
                ("time", "y", "x"),
                add_time_axis(np.full(shape, args.vc, dtype=np.float32)),
            ),
            "rhoa": (
                ("time", "y", "x"),
                add_time_axis(np.full(shape, args.rhoa, dtype=np.float32)),
            ),
            "rhow": (
                ("time", "y", "x"),
                add_time_axis(np.full(shape, args.rhow, dtype=np.float32)),
            ),
        },
        coords={
            "time": ("time", np.array([np.datetime64(valid_time, "s")])),
            "x": ("x", x.astype(np.float32)),
            "y": ("y", y.astype(np.float32)),
        },
        attrs={
            "title": "Idealized moving Rankine hurricane vortex forcing",
            "rankine_vmax_m_s": args.vmax,
            "rankine_rmw_m": args.rmw_m,
            "grid_dx_m": args.dx_m,
            "grid_dy_m": args.dy_m,
            "vortex_center_x_m": vortex_x,
            "vortex_center_y_m": vortex_y,
            "rotation": "clockwise" if args.clockwise else "counter-clockwise",
        },
    )

    ds["time"].attrs.update({"long_name": "time"})
    ds["x"].attrs.update({"long_name": "x coordinate", "units": "m"})
    ds["y"].attrs.update({"long_name": "y coordinate", "units": "m"})
    ds["lon"].attrs.update({"long_name": "longitude", "units": "degrees_east"})
    ds["lat"].attrs.update({"long_name": "latitude", "units": "degrees_north"})
    ds["uw"].attrs.update({"long_name": "wind x-component", "units": "m s-1"})
    ds["vw"].attrs.update({"long_name": "wind y-component", "units": "m s-1"})
    ds["uc"].attrs.update({"long_name": "current x-component", "units": "m s-1"})
    ds["vc"].attrs.update({"long_name": "current y-component", "units": "m s-1"})
    ds["rhoa"].attrs.update({"long_name": "air density", "units": "kg m-3"})
    ds["rhow"].attrs.update({"long_name": "water density", "units": "kg m-3"})

    return ds


def write_hourly_files(args):
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    start_time = datetime.fromisoformat(args.start_time.replace("_", "T"))

    for hour in range(args.last_hour + 1):
        valid_time = start_time + timedelta(hours=hour)
        outfile = output_dir / f"umwmin_{valid_time:%Y-%m-%d_%H:%M:%S}.nc"
        ds = make_dataset(args, valid_time, hour)
        encoding = {name: {"dtype": "float32"} for name in ds.data_vars}
        hour_ds = ds.assign_attrs(valid_time=valid_time.strftime("%Y-%m-%d %H:%M:%S"))
        hour_ds.to_netcdf(
            outfile,
            engine="netcdf4",
            format="NETCDF4",
            encoding=encoding,
            unlimited_dims=("time",),
        )
        print(f"Writing {outfile}")


def main():
    args = parse_args()
    if args.nx <= 0 or args.ny <= 0:
        raise ValueError("nx and ny must be > 0")
    if args.dx_m <= 0 or args.dy_m <= 0:
        raise ValueError("dx-m and dy-m must be > 0")
    if args.last_hour < 0:
        raise ValueError("last-hour must be >= 0")

    write_hourly_files(args)


if __name__ == "__main__":
    main()
