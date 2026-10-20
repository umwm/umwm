import numpy as np
import pytest
import xarray as xr

from umwm.plotting import plot_field


@pytest.fixture
def dataset():
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    ny, nx = 4, 6
    time = np.array(["2026-01-01T00:00:00", "2026-01-01T01:00:00"], dtype="datetime64[s]")
    swh = np.random.default_rng(0).random((2, ny, nx)).astype(np.float32)
    swh[:, 0, :] = np.nan
    return xr.Dataset(
        {"swh": (("time", "y", "x"), swh, {"long_name": "significant wave height", "units": "m"}),
         "taux_form": (("time", "y", "x"), swh - 0.5, {"units": "n/m^2"}),
         "lon": (("time", "y", "x"), np.zeros((2, ny, nx))), "lat": (("time", "y", "x"), np.zeros((2, ny, nx)))},
        coords={"time": time},
    )


def test_plot_field_sequential_and_diverging(dataset):
    ax = plot_field(dataset, "swh")
    assert ax.get_title() == "significant wave height at 2026-01-01 01:00:00"
    assert ax.get_xlabel() == "x [cell]"
    assert ax.collections[0].get_cmap().name == "viridis"
    ax = plot_field(dataset, "taux_form", time=0, cmap="Blues")
    assert ax.collections[0].get_cmap().name == "Blues"
    ax = plot_field(dataset, "taux_form")
    assert ax.collections[0].get_cmap().name == "RdBu_r"
    low, high = ax.collections[0].get_clim()
    assert low == -high


def test_plot_field_axes_from_grid_and_lonlat(dataset):
    grid = xr.Dataset({"xx": (("y", "x"), np.tile(np.arange(6) * 1000.0, (4, 1))),
                       "yy": (("y", "x"), np.tile(np.arange(4)[:, None] * 1000.0, (1, 6)))})
    ax = plot_field(dataset, "swh", grid=grid, title="custom")
    assert ax.get_xlabel() == "x [km]" and ax.get_title() == "custom"
    geographic = dataset.assign(lon=(("time", "y", "x"), np.tile(np.linspace(-80, -79, 6), (2, 4, 1))),
                                lat=(("time", "y", "x"), np.tile(np.linspace(25, 26, 4)[:, None], (2, 1, 6))))
    ax = plot_field(geographic, "swh")
    assert ax.get_xlabel().startswith("longitude")
    with pytest.raises(ValueError):
        plot_field(dataset, "swh", time=slice(None))
