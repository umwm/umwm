"""Build of the UMWM shared library for the Python package.

The library is built from the Fortran sources with the UMWM Makefile
(``make shared``), so the same build configuration applies: ``FC``,
``FCFLAGS``, ``MPI=yes``, ``GPU=yes``, ``GPU_ARCH``, ``NF_CONFIG``, ``NETCDF``,
``NETCDF_FFLAGS``, ``NETCDF_FLIBS`` and the others of ``MAKE_VARIABLES`` are
read from the environment with the prefix ``UMWM_``, for example::

    pip install .                                # serial, gfortran, nf-config
    UMWM_MPI=yes UMWM_FC=mpifort pip install .   # parallel
    UMWM_GPU=yes pip install .                   # NVIDIA GPU, nvfortran

A prebuilt library can be used instead by setting ``UMWM_LIBRARY`` to its
path, in which case it is copied into the package. The files of the example
cases in ``examples/`` that ``MANIFEST.in`` includes in the source
distribution are copied into the package as ``umwm/examples``.
"""

import os
import shutil
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py

try:
    from setuptools.command.bdist_wheel import bdist_wheel
except ImportError:  # setuptools < 70.1
    from wheel.bdist_wheel import bdist_wheel

ROOT = Path(__file__).resolve().parent
LIBRARY = "libumwm.so"
MAKE_VARIABLES = ("MPI", "GPU", "GPU_ARCH", "GPU_FLAGS", "ACC_KERNEL_FLAGS", "FC", "FCFLAGS",
                  "CPPFLAGS", "LDFLAGS", "LDLIBS", "NF_CONFIG", "NETCDF", "NETCDF_FFLAGS",
                  "NETCDF_FLIBS")


def example_patterns():
    """The patterns, relative to examples/, of the files of the example cases
    that MANIFEST.in includes, so that the wheel holds those of the sdist."""
    patterns = []
    for line in (ROOT / "MANIFEST.in").read_text(encoding="utf-8").splitlines():
        words = line.split()
        if words[:1] == ["include"]:
            patterns += [word[len("examples/"):] for word in words[1:] if word.startswith("examples/")]
    return patterns


EXAMPLES = example_patterns()


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


class BuildPy(build_py):
    """Copy the example cases into the package. An editable install finds
    them in the repository instead."""

    def run(self):
        super().run()
        if getattr(self, "editable_mode", False):
            return
        source = ROOT / "examples"
        destination = Path(self.build_lib) / "umwm" / "examples"
        for pattern in EXAMPLES:
            for path in sorted(source.glob(pattern)):
                target = destination / path.relative_to(source)
                self.mkpath(str(target.parent))
                self.copy_file(str(path), str(target))


class PlatformWheel(bdist_wheel):
    """Tag the wheel as specific to the platform but not to the Python version:
    the library is loaded with ctypes and does not use the Python ABI."""

    def get_tag(self):
        _, _, platform = super().get_tag()
        return "py3", "none", platform


setup(
    ext_modules=[Extension("umwm.libumwm", sources=[])],
    cmdclass={"build_ext": BuildSharedLibrary, "build_py": BuildPy, "bdist_wheel": PlatformWheel},
)
