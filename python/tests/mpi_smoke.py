"""Smoke test of the Python interface under MPI; run with mpiexec, for example::

    mpiexec -n 2 python python/tests/mpi_smoke.py [WORKDIR]

The idealized hurricane case of the Fortran test suite is run on all processes
and the gathered output fields are compared with the reference values on the
root process. The jonswap example case, periodic in x, is run for an hour, and
its wave field must be uniform along x, as in a serial run. The library must
be built with MPI (``MPI=yes``).
"""

import sys
from pathlib import Path

import numpy as np

import umwm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from conftest import RANKINE_FORCING, RANKINE_REFERENCE, REPO_ROOT, matches_reference  # noqa: E402


def main() -> int:
    workdir = Path(sys.argv[1] if len(sys.argv) > 1 else "mpi_smoke").resolve()
    config = umwm.Config.from_namelist(REPO_ROOT / "tests" / "namelists" / "forcing_rankine.nml")
    domain = umwm.Domain(config, workdir)
    if domain.rank == 0:
        print(f"mpi_smoke: mpi={domain.mpi} size={domain.size} workdir={workdir}")
        umwm.write_rankine_forcing(domain.input_dir, config.start_time, 1, **RANKINE_FORCING)
    domain.barrier()
    domain.run()
    local = domain.field("swh")
    print(f"mpi_smoke: rank {domain.rank} holds cells {domain.sizes['istart']}..{domain.sizes['iend']}, "
          f"local max swh {local.max():.4f}")
    failures = []
    for name, (mean, maximum, tolerance) in RANKINE_REFERENCE.items():
        values = domain.grid_field(name)  # collective; only rank 0 receives the field
        if domain.rank == 0 and not matches_reference(values, mean, maximum, tolerance):
            failures.append(f"{name}: mean {np.nanmean(values):.7g} max {np.nanmax(values):.7g} "
                            f"expected {mean:.7g} {maximum:.7g}")
    state = domain.state()
    if domain.rank == 0:
        output = domain.output()
        if not np.allclose(state.swh.values[0], output.swh.isel(time=-1).values, rtol=1e-6, equal_nan=True):
            failures.append("state and output swh differ")
        if domain.size > 1 and np.unique(domain.grid_field("nproc")).size < 2:
            failures.append("the grid is not partitioned")
    domain.finalize()
    failures += periodic_case(workdir / "jonswap")
    if domain.rank == 0:
        if failures:
            print("mpi_smoke: FAILED\n  " + "\n  ".join(failures))
            return 1
        print(f"mpi_smoke: OK on {domain.size} process(es)")
    return 0


def periodic_case(workdir: Path) -> list[str]:
    """Run the jonswap example case, periodic in x, for an hour, and check that
    its wave field is uniform along x on the root process."""
    config = umwm.Config.from_case("jonswap")
    config.domain.stoptimestr = "2026-01-01 01:00:00"
    config.output.outspec = 0
    config.output.stokes = False
    with umwm.Domain(config, workdir, quiet=True) as domain:
        domain.run()
        swh = domain.grid_field("swh")
        if domain.rank > 0:
            return []
    spread = np.nanmax(np.abs(swh - swh[:, :1]))
    if not np.isfinite(spread) or spread > 1e-5 * np.nanmax(swh):
        return [f"jonswap: swh is not uniform along the periodic x direction (spread {spread:.4g} m)"]
    return []


if __name__ == "__main__":
    sys.exit(main())
