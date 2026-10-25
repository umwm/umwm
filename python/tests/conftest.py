"""Shared fixtures of the UMWM Python test suite."""

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

# Reference domain mean and maximum of output fields of the one-hour idealized
# hurricane case, and the tolerances relative to the larger of the two, from
# tests/test_domain.f90.
RANKINE_REFERENCE = {
    "swh": (3.811700821, 6.315617085, 1e-4),
    "mwp": (5.526701450, 6.882987022, 1e-4),
    "mss": (5.120731890e-02, 6.678894162e-02, 1e-4),
    "cd": (4.989048932e-03, 4.808749259e-01, 1e-3),
    "ust": (9.779363275e-01, 1.610171556, 1e-4),
    "dwp": (7.718040943, 10.06773663, 1e-4),
    "taux_form": (-2.226837166e-02, 2.233541965, 1e-4),
    "taux_ocn": (-2.457904629e-02, 2.165972948, 1e-4),
    "epsx_ocn": (-3.290932626e-02, 4.102433681, 1e-4),
    "epsy_atm": (1.479409933e-01, 1.229280758e+01, 1e-4),
    "u_stokes": (-5.562203936e-03, 3.102710545e-01, 1e-4),
    "v_stokes": (-5.914379377e-03, 3.947870731e-01, 1e-4),
}

# Parameters of the forcing of that case, as in tests/Makefile.
RANKINE_FORCING = dict(nx=21, ny=11, dx=10000.0, dy=10000.0, vmax=30.0, rmw=50000.0,
                       center_x=0.0, center_y=0.0, translation_x=10000.0, translation_y=0.0,
                       uc=0.25, vc=-0.5, rhoa=1.3, rhow=1025.0)


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def library():
    """The UMWM library, or a skip when it is not built."""
    import umwm

    try:
        return umwm.library()
    except umwm.LibraryNotFoundError as error:
        pytest.skip(f"UMWM library not available: {error}")


@pytest.fixture
def rankine_config(repo_root):
    """The configuration of the idealized hurricane case of the Fortran test suite."""
    import umwm

    return umwm.Config.from_namelist(repo_root / "tests" / "namelists" / "forcing_rankine.nml")


@pytest.fixture
def rankine_domain(library, rankine_config, tmp_path):
    """A domain of the idealized hurricane case with its forcing files, not yet initialized."""
    import umwm

    domain = umwm.Domain(rankine_config, tmp_path / "rankine")
    umwm.write_rankine_forcing(domain.input_dir, rankine_config.start_time, 1, **RANKINE_FORCING)
    yield domain
    domain.finalize()


def matches_reference(values, mean, maximum, tolerance) -> bool:
    """Whether the mean and maximum over the sea cells of a gridded field match
    the reference within the tolerance relative to the larger of the two."""
    import numpy as np

    scale = max(abs(mean), abs(maximum))
    return (abs(np.nanmean(values) - mean) <= tolerance * scale
            and abs(np.nanmax(values) - maximum) <= tolerance * scale)
