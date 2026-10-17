"""Build of the UMWM shared library for the Python package.

The library is built from the Fortran sources with the UMWM Makefile
(``make shared``), so the same build configuration applies: ``FC``,
``FCFLAGS``, ``MPI=yes``, ``NF_CONFIG``, ``NETCDF``, ``NETCDF_FFLAGS`` and
``NETCDF_FLIBS`` are read from the environment, for example::

    pip install .                        # serial, gfortran, nf-config
    MPI=yes FC=mpifort pip install .     # parallel

A prebuilt library can be used instead by setting ``UMWM_LIBRARY`` to its
path, in which case it is copied into the package.
"""

import os
import shutil
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

try:
    from setuptools.command.bdist_wheel import bdist_wheel
except ImportError:  # setuptools < 70.1
    from wheel.bdist_wheel import bdist_wheel

ROOT = Path(__file__).resolve().parent
LIBRARY = "libumwm.so"
MAKE_VARIABLES = ("MPI", "FC", "FCFLAGS", "CPPFLAGS", "LDFLAGS", "LDLIBS", "NF_CONFIG",
                  "NETCDF", "NETCDF_FFLAGS", "NETCDF_FLIBS")


class BuildSharedLibrary(build_ext):
    """Build libumwm.so with the Makefile and place it in the package."""

    def get_ext_filename(self, ext_name):
        # The library is loaded with ctypes, not imported, so it keeps a plain
        # name without the Python ABI tag.
        return os.path.join(*ext_name.split(".")[:-1], LIBRARY)

    def build_extension(self, extension):
        destination = Path(self.get_ext_fullpath(extension.name))
        destination.parent.mkdir(parents=True, exist_ok=True)

        prebuilt = os.environ.get("UMWM_LIBRARY")
        if prebuilt:
            source = Path(prebuilt).expanduser()
            if not source.is_file():
                raise FileNotFoundError(f"UMWM_LIBRARY points to a missing file: {source}")
            self.announce(f"copying prebuilt library {source}", level=3)
            shutil.copy2(source, destination)
            return

        build_dir = (Path(self.build_temp).resolve() / "umwm-shared")
        command = [os.environ.get("MAKE", "make"), "--directory", str(ROOT / "src"), "shared",
                   f"SHARED_DIR={build_dir}", f"SHARED_LIB={LIBRARY}"]
        for name in MAKE_VARIABLES:
            value = os.environ.get(f"UMWM_{name}")
            if value is not None:
                command.append(f"{name}={value}")
        self.spawn(command)
        shutil.copy2(build_dir / LIBRARY, destination)


class PlatformWheel(bdist_wheel):
    """Tag the wheel as specific to the platform but not to the Python version:
    the library is loaded with ctypes and does not use the Python ABI."""

    def get_tag(self):
        _, _, platform = super().get_tag()
        return "py3", "none", platform


setup(
    ext_modules=[Extension("umwm.libumwm", sources=[])],
    cmdclass={"build_ext": BuildSharedLibrary, "bdist_wheel": PlatformWheel},
)
