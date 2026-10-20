"""University of Miami Wave Model (UMWM) Python interface.

Configure, run, and inspect the wave model from Python::

    import umwm

    config = umwm.Config()
    config.domain.stoptimestr = "2026-01-01 03:00:00"
    config.forcing_constant.wspd0 = 15.0

    with umwm.Domain(config, workdir="runs/case1") as domain:
        domain.run()
        swh = domain.grid_field("swh")   # NumPy array (ny, nx)
        output = domain.output()          # xarray.Dataset of the output files

The model runs in the compiled Fortran library; Python drives it through
its C API and views its memory without copying.
"""

from ._library import LibraryNotFoundError, barrier, find_library, flush_denormals, library, mpi, rank, size
from .config import (
    Config,
    ConfigError,
    DomainConfig,
    ForcingConfig,
    ForcingConstantConfig,
    GridConfig,
    OutputConfig,
    PhysicsConfig,
    SpectrumPoint,
    StokesConfig,
)
from .domain import Domain, run
from .fields import FIELDS
from .forcing import forcing_filename, write_forcing, write_forcing_dataset, write_grid, write_rankine_forcing
from .plotting import plot_field
from .output import (
    open_grid,
    open_output,
    open_restart,
    open_spectrum_output,
    output_files,
    restart_files,
    spectrum_output_files,
)

__version__ = "26.10"

__all__ = [
    "Config", "ConfigError", "Domain", "run", "SpectrumPoint",
    "DomainConfig", "PhysicsConfig", "GridConfig", "ForcingConfig", "ForcingConstantConfig",
    "OutputConfig", "StokesConfig", "FIELDS",
    "write_forcing", "write_forcing_dataset", "write_grid", "write_rankine_forcing", "forcing_filename",
    "plot_field",
    "open_output", "open_spectrum_output", "open_restart", "open_grid",
    "output_files", "restart_files", "spectrum_output_files",
    "library", "find_library", "LibraryNotFoundError", "rank", "size", "mpi", "barrier", "flush_denormals",
    "__version__",
]
