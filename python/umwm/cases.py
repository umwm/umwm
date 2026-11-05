"""The example cases: idealized configurations of the model, with constant
wind and other parameters, ready to run. Most need no input files; the
others bring theirs, such as a bathymetry or steady currents.

The cases are in the ``examples`` directory of the UMWM repository, of which
the package includes a copy: one configuration directory per case, with its
``main.nml`` namelist (whose leading comments describe the case), the
``spectrum.nml`` points of its spectrum output, and the input files of the
cases that have them in ``input/``, listed with a one-line description in
``examples/cases.txt``. A case runs like any configuration directory, with
``umwm run DIR -o OUTPUT`` or the executable ``umwm DIR -o OUTPUT``; ``umwm
cases`` lists them with their directories, and ``umwm init --case NAME DIR``
copies one into a directory to edit.

>>> import umwm
>>> case = umwm.get_case("jonswap")
>>> domain = umwm.run(case.directory, output_dir="runs/jonswap")
"""

from __future__ import annotations

import filecmp
import os
import shutil
from dataclasses import dataclass
from pathlib import Path

from .config import Config

__all__ = ["Case", "examples_dir", "list_cases", "get_case"]

INDEX = "cases.txt"
"""The list of the cases in the examples directory."""


def examples_dir() -> Path:
    """The directory of the example cases: ``$UMWM_EXAMPLES`` if set, else the
    ``examples`` directory of the current directory, such as the repository,
    as for the ``umwm`` executable, else the copy included in the package, or
    the ``examples`` directory of the repository for a package installed from
    it in editable mode."""
    variable = os.environ.get("UMWM_EXAMPLES")
    if variable:
        path = Path(variable).expanduser()
        if not (path / INDEX).is_file():
            raise FileNotFoundError(f"no {INDEX} in UMWM_EXAMPLES={variable}")
        return path
    package = Path(__file__).resolve().parent
    for path in (Path("examples"), package / "examples", package.parents[1] / "examples"):
        if (path / INDEX).is_file():
            return path.resolve()
    raise FileNotFoundError("the example cases are not found; set UMWM_EXAMPLES to the examples "
                            "directory of the UMWM repository")


@dataclass(frozen=True)
class Case:
    """An example case: its name, one-line description, and directory."""

    name: str
    description: str
    directory: Path

    @property
    def namelist(self) -> Path:
        """The ``main.nml`` namelist of the case."""
        return self.directory / "main.nml"

    @property
    def input_dir(self) -> Path | None:
        """The directory of the input files of the case, or None if it has none."""
        path = self.directory / "input"
        return path if path.is_dir() else None

    def config(self) -> Config:
        """The configuration of the case, with its spectrum output points."""
        return Config.from_namelist(self.namelist)

    def domain(self, output_dir=None, **kwargs):
        """A :class:`~umwm.Domain` of the case, writing its output to
        ``output_dir`` (by default a directory named after the case) and
        reading the input files of the case from its ``input`` directory,
        unless ``kwargs`` give another ``input_dir``; ``kwargs`` are passed
        to the domain."""
        from .domain import Domain

        return Domain(self.namelist, Path(output_dir or self.name), **kwargs)

    def copy_inputs(self, directory, overwrite=False) -> list[Path]:
        """Copy the input files of the case into ``directory``, an input
        directory, and return the paths of the files copied. A file there
        that is the same as that of the case is kept, and one that differs
        raises ``FileExistsError`` unless ``overwrite``."""
        Path(directory).mkdir(parents=True, exist_ok=True)
        return self._copy(self._inputs_to_copy(directory, overwrite))

    def _inputs_to_copy(self, directory, overwrite) -> list[tuple[Path, Path]]:
        """The (source, target) pairs of the input files to copy into directory."""
        pending = []
        for source in sorted(self.input_dir.iterdir()) if self.input_dir else []:
            target = Path(directory) / source.name
            if target.is_file() and filecmp.cmp(source, target, shallow=False):
                continue
            if (target.exists() or target.is_symlink()) and not overwrite:
                raise FileExistsError(f"{target} differs from the input file {source.name} of the case "
                                      f"{self.name}; remove it, or replace it with overwrite=True (umwm init --force)")
            pending.append((source, target))
        return pending

    @staticmethod
    def _copy(pending) -> list[Path]:
        for source, target in pending:
            if target.exists() or target.is_symlink():
                target.unlink()  # not to write through a link
            shutil.copy2(source, target)
        return [target for _, target in pending]

    def notes(self) -> str:
        """The description of the case in the leading comments of its namelist."""
        lines = []
        for line in self.namelist.read_text(encoding="utf-8").splitlines():
            if not line.startswith("!"):
                break
            lines.append(line[2:] if line.startswith("! ") else line[1:])
        return "\n".join(lines).strip("\n")


def list_cases(directory=None) -> list[Case]:
    """The example cases, in the order of ``cases.txt``, from ``directory`` or
    else :func:`examples_dir`."""
    directory = Path(directory) if directory is not None else examples_dir()
    cases = []
    for line in (directory / INDEX).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        name, _, description = line.partition(" ")
        cases.append(Case(name, description.strip(), directory / name))
    return cases


def get_case(name: str, directory=None) -> Case:
    """The example case of the given name, or else the only case whose name
    starts with it, such as ``"titan"``. Raises ``KeyError`` otherwise."""
    cases = list_cases(directory)
    for case in cases:
        if case.name == name:
            return case
    matches = [case for case in cases if case.name.startswith(name)]
    if len(matches) == 1:
        return matches[0]
    if matches:
        raise KeyError(f"{name!r} matches more than one case: {', '.join(case.name for case in matches)}")
    raise KeyError(f"unknown case {name!r}; the example cases are: {', '.join(case.name for case in cases)}")
