"""Loading of the UMWM shared library and its C API.

The library ``libumwm.so`` is built from the Fortran sources by the package
build (``pip install .``) and placed next to this module. It can also be
pointed to with the ``UMWM_LIBRARY`` environment variable, for example when
the library is built with ``make shared`` or with a different compiler or
MPI library.
"""

from __future__ import annotations

import atexit
import ctypes
import os
import sys
from ctypes import POINTER, byref, c_char_p, c_double, c_float, c_int, c_void_p
from pathlib import Path

import numpy as np

__all__ = ["Library", "LibraryNotFoundError", "library", "find_library", "rank", "size", "mpi", "barrier",
           "flush_denormals"]

LIBRARY_NAME = "libumwm.so"

# element type codes of umwm_domain_field
TYPE_REAL, TYPE_INTEGER, TYPE_LOGICAL = 1, 2, 3

_CTYPES = {
    (TYPE_REAL, 4): ctypes.c_float,
    (TYPE_REAL, 8): ctypes.c_double,
    (TYPE_INTEGER, 4): ctypes.c_int32,
    (TYPE_INTEGER, 8): ctypes.c_int64,
    (TYPE_LOGICAL, 4): ctypes.c_int32,
    (TYPE_LOGICAL, 8): ctypes.c_int64,
}


class LibraryNotFoundError(ImportError):
    """Raised when the UMWM shared library cannot be found."""


def find_library() -> Path:
    """Locate the UMWM shared library.

    ``UMWM_LIBRARY`` takes precedence; otherwise the library next to the
    package, placed there by the package build, is used.
    """
    override = os.environ.get("UMWM_LIBRARY")
    if override:
        path = Path(override).expanduser()
        if not path.is_file():
            raise LibraryNotFoundError(f"UMWM_LIBRARY points to a missing file: {path}")
        return path
    path = Path(__file__).resolve().parent / LIBRARY_NAME
    if path.is_file():
        return path
    raise LibraryNotFoundError(
        f"the UMWM shared library {LIBRARY_NAME} was not found in {path.parent}; "
        "build the package from the UMWM repository with 'pip install .' "
        "(or 'pip install -e .'), or build the library with 'make shared' and "
        "set UMWM_LIBRARY to its path"
    )


class Library:
    """The UMWM shared library and the Fortran parallel environment.

    The environment is initialized on first use (``MPI_Init`` when the library
    was built with MPI, unless MPI is already initialized, for example by
    mpi4py) and finalized when the interpreter exits.
    """

    def __init__(self, path=None):
        self.path = Path(path) if path else find_library()
        # RTLD_GLOBAL: MPI libraries load plugins that need the MPI symbols
        # to be globally visible.
        self._lib = ctypes.CDLL(str(self.path), mode=ctypes.RTLD_GLOBAL)
        self._bind()
        self._started = False

    def _bind(self) -> None:
        lib = self._lib
        lib.umwm_env_init.argtypes = []
        lib.umwm_env_init.restype = None
        lib.umwm_env_stop.argtypes = []
        lib.umwm_env_stop.restype = None
        lib.umwm_env_barrier.argtypes = []
        lib.umwm_env_barrier.restype = None
        lib.umwm_env_flush.argtypes = []
        lib.umwm_env_flush.restype = None
        lib.umwm_env_flush_denormals.argtypes = [c_int]
        lib.umwm_env_flush_denormals.restype = None
        lib.umwm_env_denormals_flushed.argtypes = []
        lib.umwm_env_denormals_flushed.restype = c_int
        for name in ("umwm_env_rank", "umwm_env_size", "umwm_env_mpi"):
            getattr(lib, name).argtypes = []
            getattr(lib, name).restype = c_int
        lib.umwm_domain_new.argtypes = []
        lib.umwm_domain_new.restype = c_void_p
        lib.umwm_domain_initialize.argtypes = [c_void_p, c_char_p]
        lib.umwm_domain_initialize.restype = c_int
        for name in ("umwm_domain_advance", "umwm_domain_run", "umwm_domain_diagnostics",
                     "umwm_domain_ocean_stress", "umwm_domain_write_output", "umwm_domain_write_spectrum",
                     "umwm_domain_write_restart", "umwm_domain_finalize", "umwm_domain_delete"):
            getattr(lib, name).argtypes = [c_void_p]
            getattr(lib, name).restype = None
        lib.umwm_domain_finished.argtypes = [c_void_p]
        lib.umwm_domain_finished.restype = c_int
        lib.umwm_domain_time.argtypes = [c_void_p, c_char_p]
        lib.umwm_domain_time.restype = None
        lib.umwm_domain_int.argtypes = [c_void_p, c_char_p, POINTER(c_int)]
        lib.umwm_domain_int.restype = c_int
        lib.umwm_domain_real.argtypes = [c_void_p, c_char_p, POINTER(c_double)]
        lib.umwm_domain_real.restype = c_int
        lib.umwm_domain_field.argtypes = [
            c_void_p, c_char_p, POINTER(c_void_p), POINTER(c_int), POINTER(c_int),
            POINTER(c_int), c_int * 3, c_int * 3,
        ]
        lib.umwm_domain_field.restype = c_int
        lib.umwm_domain_gather.argtypes = [c_void_p, c_char_p, c_int, POINTER(c_float)]
        lib.umwm_domain_gather.restype = c_int

    # ---- environment --------------------------------------------------------

    @property
    def mpi(self) -> bool:
        """Whether the library was built with MPI."""
        return bool(self._lib.umwm_env_mpi())

    def env_init(self) -> None:
        """Initialize the parallel environment, once."""
        if not self._started:
            self._lib.umwm_env_init()
            self._started = True
            atexit.register(self.env_stop)

    def env_stop(self) -> None:
        """Finalize the parallel environment, once."""
        if self._started:
            self._started = False
            self._lib.umwm_env_stop()

    @property
    def rank(self) -> int:
        """Rank of this process, 0 in serial mode."""
        self.env_init()
        return self._lib.umwm_env_rank()

    @property
    def size(self) -> int:
        """Number of processes, 1 in serial mode."""
        self.env_init()
        return self._lib.umwm_env_size()

    def barrier(self) -> None:
        """Wait for all processes; a no-op in serial mode."""
        self.env_init()
        self._lib.umwm_env_barrier()

    def flush(self) -> None:
        """Flush the standard output and error of the Fortran runtime."""
        self._lib.umwm_env_flush()

    def flush_denormals(self, enable: bool = True) -> None:
        """Enable or disable flushing denormal floating-point results to zero
        in this process.

        A library compiled with fast math (the default flags) enables it on
        initialization, as the compiler does at the startup of the executable
        but not of a shared library; the wave spectrum underflows in much of
        its range, where denormal arithmetic is very slow. Disabling it
        restores gradual underflow, for the whole process (NumPy included).
        """
        self._lib.umwm_env_flush_denormals(1 if enable else 0)

    @property
    def denormals_flushed(self) -> bool:
        """Whether denormal results are flushed to zero in this process."""
        return bool(self._lib.umwm_env_denormals_flushed())

    # ---- domain lifecycle ---------------------------------------------------

    def domain_new(self) -> int:
        self.env_init()
        return self._lib.umwm_domain_new()

    def domain_initialize(self, handle: int, path: str) -> int:
        sys.stdout.flush()
        return self._lib.umwm_domain_initialize(handle, os.fsencode(path))

    def domain_advance(self, handle: int) -> None:
        sys.stdout.flush()
        self._lib.umwm_domain_advance(handle)

    def domain_run(self, handle: int) -> None:
        sys.stdout.flush()
        self._lib.umwm_domain_run(handle)

    def domain_finished(self, handle: int) -> bool:
        return bool(self._lib.umwm_domain_finished(handle))

    def domain_time(self, handle: int) -> str:
        buffer = ctypes.create_string_buffer(20)
        self._lib.umwm_domain_time(handle, buffer)
        return buffer.value.decode()

    def domain_diagnostics(self, handle: int) -> None:
        self._lib.umwm_domain_diagnostics(handle)

    def domain_ocean_stress(self, handle: int) -> None:
        self._lib.umwm_domain_ocean_stress(handle)

    def domain_write_output(self, handle: int) -> None:
        sys.stdout.flush()
        self._lib.umwm_domain_write_output(handle)

    def domain_write_spectrum(self, handle: int) -> None:
        sys.stdout.flush()
        self._lib.umwm_domain_write_spectrum(handle)

    def domain_write_restart(self, handle: int) -> None:
        sys.stdout.flush()
        self._lib.umwm_domain_write_restart(handle)

    def domain_finalize(self, handle: int) -> None:
        self._lib.umwm_domain_finalize(handle)

    def domain_delete(self, handle: int) -> None:
        self._lib.umwm_domain_delete(handle)

    # ---- domain attributes and fields ---------------------------------------

    def domain_int(self, handle: int, name: str) -> int:
        value = c_int()
        if self._lib.umwm_domain_int(handle, name.encode(), byref(value)) != 0:
            raise KeyError(f"unknown integer attribute {name!r}")
        return value.value

    def domain_real(self, handle: int, name: str) -> float:
        value = c_double()
        if self._lib.umwm_domain_real(handle, name.encode(), byref(value)) != 0:
            raise KeyError(f"unknown real attribute {name!r}")
        return value.value

    def domain_field(self, handle: int, name: str) -> tuple[np.ndarray, tuple[int, ...]]:
        """A zero-copy view of a domain field in Fortran (column-major) order,
        with the Fortran lower bounds of its dimensions."""
        pointer = c_void_p()
        dtype, itemsize, rank = c_int(), c_int(), c_int()
        shape, lbounds = (c_int * 3)(), (c_int * 3)()
        stat = self._lib.umwm_domain_field(
            handle, name.encode(), byref(pointer), byref(dtype), byref(itemsize), byref(rank), shape, lbounds)
        if stat == 1:
            raise KeyError(f"unknown field {name!r}")
        if stat == 2:
            raise ValueError(f"field {name!r} is not allocated")
        try:
            ctype = _CTYPES[(dtype.value, itemsize.value)]
        except KeyError:
            raise TypeError(f"unsupported element type of field {name!r}") from None
        dims = tuple(shape[i] for i in range(rank.value))
        bounds = tuple(lbounds[i] for i in range(rank.value))
        if pointer.value is None or 0 in dims:
            return np.empty(dims, dtype=np.ctypeslib.as_ctypes_type(ctype), order="F"), bounds
        flat = ctypes.cast(pointer, POINTER(ctype))
        # A C-order view of the reversed shape, transposed, is the Fortran array.
        view = np.ctypeslib.as_array(flat, shape=dims[::-1]).T
        return view, bounds

    def domain_gather(self, handle: int, name: str, level: int, values: np.ndarray) -> None:
        """Gather a cell field to the root process into ``values``, a C-contiguous
        float32 array of shape (nm, mm); land cells are set to NaN."""
        if values.dtype != np.float32 or not values.flags.c_contiguous:
            raise TypeError("values must be a C-contiguous float32 array")
        stat = self._lib.umwm_domain_gather(handle, name.encode(), level, values.ctypes.data_as(POINTER(c_float)))
        if stat != 0:
            raise KeyError(f"field {name!r} cannot be gathered to the grid")


_library: Library | None = None


def library(path=None) -> Library:
    """The loaded UMWM library (loaded on first call)."""
    global _library
    if _library is None:
        _library = Library(path)
    elif path is not None and Path(path) != _library.path:
        raise RuntimeError(f"the UMWM library is already loaded from {_library.path}")
    return _library


def rank() -> int:
    """Rank of this process, 0 in serial mode."""
    return library().rank


def size() -> int:
    """Number of processes, 1 in serial mode."""
    return library().size


def mpi() -> bool:
    """Whether the UMWM library was built with MPI."""
    return library().mpi


def barrier() -> None:
    """Wait for all processes; a no-op in serial mode."""
    library().barrier()


def flush_denormals(enable: bool = True) -> None:
    """Enable or disable flushing denormal results to zero; see :meth:`Library.flush_denormals`."""
    library().flush_denormals(enable)
