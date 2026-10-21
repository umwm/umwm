"""Model configuration mirroring the UMWM namelists.

:class:`Config` holds one dataclass per namelist group of ``main.nml`` with
the variable names of the namelists, so that the Fortran documentation and
the ``--set group.name=value`` syntax of the command-line tool apply directly.
The defaults are those of ``namelists/main.nml`` in the UMWM repository. The
spectrum output points of ``spectrum.nml`` and the cells of enclosed seas to
exclude of ``exclude.nml`` are held alongside the namelist groups.
"""

from __future__ import annotations

import copy
import dataclasses
import json
from collections import OrderedDict
from dataclasses import dataclass, field, fields
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Iterator

from . import namelist
from .namelist import NamelistError

__all__ = [
    "Config",
    "ConfigError",
    "DomainConfig",
    "PhysicsConfig",
    "GridConfig",
    "ForcingConfig",
    "ForcingConstantConfig",
    "OutputConfig",
    "StokesConfig",
    "SpectrumPoint",
    "TIME_FORMAT",
    "parse_value",
    "namelist_path",
]

TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
"""Format of the time strings of the ``DOMAIN`` namelist."""

OUTPUT_INTERVALS = (-1, 0, 1, 2, 3, 4, 6, 8, 12, 24)
"""Allowed values of the output intervals in hours; -1 means every global time step."""


class ConfigError(ValueError):
    """Raised for an invalid configuration."""


def _doc(description: str, units: str | None = None) -> dict:
    """Field metadata: a description and optional units."""
    return {"description": description, "units": units}


@dataclass
class DomainConfig:
    """The ``DOMAIN`` namelist: grid and spectrum sizes and the time span."""

    isglobal: bool = field(default=False, metadata=_doc("Global (periodic in x) or regional domain"))
    mm: int = field(default=21, metadata=_doc("Domain size in x"))
    nm: int = field(default=11, metadata=_doc("Domain size in y"))
    om: int = field(default=37, metadata=_doc("Number of frequency bins"))
    pm: int = field(default=36, metadata=_doc("Number of directions"))
    fmin: float = field(default=0.0313, metadata=_doc("Lowest frequency bin", "Hz"))
    fmax: float = field(default=2.0, metadata=_doc("Highest frequency bin", "Hz"))
    fprog: float = field(default=2.0, metadata=_doc("Highest prognostic frequency bin", "Hz"))
    starttimestr: str = field(default="2026-01-01 00:00:00", metadata=_doc("Simulation start time"))
    stoptimestr: str = field(default="2026-01-01 06:00:00", metadata=_doc("Simulation end time"))
    reftimestr: str = field(default="1970-01-01 00:00:00", metadata=_doc("Reference time of the output time axis"))
    dtg: float = field(default=3600.0, metadata=_doc("Global (forcing and output) time step", "s"))
    restart: bool = field(default=False, metadata=_doc("Restart from file"))


@dataclass
class PhysicsConfig:
    """The ``PHYSICS`` namelist: physical constants and source function coefficients."""

    g: float = field(default=9.80665, metadata=_doc("Gravitational acceleration", "m/s^2"))
    nu_air: float = field(default=1.56e-5, metadata=_doc("Kinematic viscosity of air", "m^2/s"))
    nu_water: float = field(default=0.90e-6, metadata=_doc("Kinematic viscosity of water", "m^2/s"))
    sfct: float = field(default=0.07, metadata=_doc("Surface tension", "N/m"))
    kappa: float = field(default=0.4, metadata=_doc("Von Karman constant"))
    z: float = field(default=10.0, metadata=_doc("Height of the input wind speed", "m"))
    gustiness: float = field(default=0.0, metadata=_doc("Random wind gustiness factor (0 to 0.2)"))
    dmin: float = field(default=10.0, metadata=_doc("Depth limiter", "m"))
    explim: float = field(default=0.9, metadata=_doc("Exponent limiter (0.69 ~ 100% growth)"))
    sin_fac: float = field(default=0.11, metadata=_doc("Input factor from following winds"))
    sin_diss1: float = field(default=0.10, metadata=_doc("Damping factor from opposing winds"))
    sin_diss2: float = field(default=0.001, metadata=_doc("Damping factor from swell overrunning wind"))
    sds_fac: float = field(default=42.0, metadata=_doc("Breaking dissipation factor"))
    sds_power: float = field(default=2.4, metadata=_doc("Saturation spectrum power"))
    mss_fac: float = field(default=360.0, metadata=_doc("Mean-square-slope adjustment to Sds"))
    snl_fac: float = field(default=5.0, metadata=_doc("Wave energy downshifting factor"))
    sdt_fac: float = field(default=0.002, metadata=_doc("Dissipation due to turbulence factor"))
    sbf_fac: float = field(default=0.003, metadata=_doc("Bottom friction coefficient", "m/s"))
    sbp_fac: float = field(default=0.003, metadata=_doc("Bottom percolation coefficient", "m/s"))
    fcut_fac: float = field(default=0.53, metadata=_doc("Prognostic cutoff frequency factor, fc = fcut_fac*g/U"))


@dataclass
class GridConfig:
    """The ``GRID`` namelist: grid geometry and bathymetry sources."""

    gridfromfile: bool = field(default=False, metadata=_doc("Read lon/lat fields from input/umwm.grid or umwm.gridtopo"))
    delx: float = field(default=10000.0, metadata=_doc("Grid spacing in x if not from file", "m"))
    dely: float = field(default=10000.0, metadata=_doc("Grid spacing in y if not from file", "m"))
    topofromfile: bool = field(default=False, metadata=_doc("Read bathymetry from input/umwm.gridtopo"))
    dpt: float = field(default=1000.0, metadata=_doc("Constant water depth if not from file", "m"))
    fillestuaries: bool = field(default=False, metadata=_doc("Fill cells with 3 land neighbours"))
    filllakes: bool = field(default=False, metadata=_doc("Fill the enclosed seas listed in exclude.nml"))


@dataclass
class ForcingConfig:
    """The ``FORCING`` namelist: which forcing fields are read from input files."""

    winds: bool = field(default=False, metadata=_doc("Wind input from file"))
    currents: bool = field(default=False, metadata=_doc("Currents input from file"))
    air_density: bool = field(default=False, metadata=_doc("Air density input from file"))
    water_density: bool = field(default=False, metadata=_doc("Water density input from file"))
    seaice: bool = field(default=False, metadata=_doc("Sea ice fraction input from file"))


@dataclass
class ForcingConstantConfig:
    """The ``FORCING_CONSTANT`` namelist: forcing values used when not read from file."""

    wspd0: float = field(default=10.0, metadata=_doc("Wind speed", "m/s"))
    wdir0: float = field(default=0.0, metadata=_doc("Wind direction (mathematical convention)", "rad"))
    uc0: float = field(default=0.0, metadata=_doc("x-component of ocean current", "m/s"))
    vc0: float = field(default=0.0, metadata=_doc("y-component of ocean current", "m/s"))
    rhoa0: float = field(default=1.2, metadata=_doc("Air density", "kg/m^3"))
    rhow0: float = field(default=1030.0, metadata=_doc("Water density", "kg/m^3"))
    fice0: float = field(default=0.0, metadata=_doc("Sea ice fraction"))
    fice_lth: float = field(default=0.30, metadata=_doc("Sea ice fraction, lower threshold for attenuation"))
    fice_uth: float = field(default=0.75, metadata=_doc("Sea ice fraction, upper threshold for attenuation"))


@dataclass
class OutputConfig:
    """The ``OUTPUT`` namelist: output intervals and screen diagnostics."""

    outgrid: int = field(default=1, metadata=_doc("Gridded output interval (0 none, -1 every step)", "hours"))
    outspec: int = field(default=0, metadata=_doc("Spectrum output interval (0 none, -1 every step)", "hours"))
    outrst: int = field(default=6, metadata=_doc("Restart output interval (0 none, -1 every step)", "hours"))
    xpl: int = field(default=20, metadata=_doc("Grid cell in x for the screen diagnostics"))
    ypl: int = field(default=6, metadata=_doc("Grid cell in y for the screen diagnostics"))
    stokes: bool = field(default=True, metadata=_doc("Compute and output Stokes drift"))


@dataclass
class StokesConfig:
    """The ``STOKES`` namelist: depths of the Stokes drift output."""

    depths: list[float] = field(
        default_factory=lambda: [
            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
            2.4, 2.8, 3.2, 3.6, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0,
            16.0, 18.0, 20.0, 24.0, 28.0, 32.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
        ],
        metadata=_doc("Depths of the Stokes drift output, positive down", "m"),
    )


@dataclass
class SpectrumPoint:
    """A spectrum output point, given by grid indices (1-based) or by longitude and latitude."""

    name: str
    x: int | None = None
    y: int | None = None
    lon: float | None = None
    lat: float | None = None

    def __post_init__(self):
        by_index = self.x is not None and self.y is not None
        by_coord = self.lon is not None and self.lat is not None
        if by_index == by_coord:
            raise ConfigError(f"spectrum point {self.name!r} needs either x and y or lon and lat")
        if not self.name or any(char.isspace() for char in self.name):
            raise ConfigError(f"spectrum point name {self.name!r} must be a word without spaces")

    @property
    def by_index(self) -> bool:
        return self.x is not None


GROUPS = OrderedDict(
    domain=DomainConfig,
    physics=PhysicsConfig,
    grid=GridConfig,
    forcing=ForcingConfig,
    forcing_constant=ForcingConstantConfig,
    output=OutputConfig,
    stokes=StokesConfig,
)
"""Namelist group names of ``main.nml`` in file order and their dataclasses."""


@dataclass
class Config:
    """The complete UMWM configuration.

    Attributes are the namelist groups of ``main.nml`` (``domain``,
    ``physics``, ``grid``, ``forcing``, ``forcing_constant``, ``output`` and
    ``stokes``), plus the spectrum output points and the cells of enclosed
    seas to exclude from the domain.

    >>> config = Config()
    >>> config.domain.mm = 51
    >>> config.set("forcing_constant.wspd0", 15)
    >>> config.write("main.nml")
    """

    domain: DomainConfig = field(default_factory=DomainConfig)
    physics: PhysicsConfig = field(default_factory=PhysicsConfig)
    grid: GridConfig = field(default_factory=GridConfig)
    forcing: ForcingConfig = field(default_factory=ForcingConfig)
    forcing_constant: ForcingConstantConfig = field(default_factory=ForcingConstantConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    stokes: StokesConfig = field(default_factory=StokesConfig)
    spectrum_points: list[SpectrumPoint] = field(default_factory=list)
    exclude_points: list[tuple[int, int]] = field(default_factory=list)

    # ---- construction -------------------------------------------------------

    @classmethod
    def from_namelist(cls, path, spectrum_points=None, exclude_points=None) -> "Config":
        """Read the configuration from a ``main.nml`` file.

        ``path`` may also be a run directory, whose ``namelists/main.nml``
        is read. ``spectrum_points`` and ``exclude_points`` are optional paths
        to a ``spectrum.nml`` and an ``exclude.nml`` file; by default, files
        of those names next to the main namelist are read if they exist.
        """
        path = namelist_path(path)
        config = cls.from_namelist_text(path.read_text(encoding="utf-8"))
        if spectrum_points is None and (path.parent / "spectrum.nml").is_file():
            spectrum_points = path.parent / "spectrum.nml"
        if exclude_points is None and (path.parent / "exclude.nml").is_file():
            exclude_points = path.parent / "exclude.nml"
        if spectrum_points is not None:
            config.spectrum_points = read_spectrum_points(spectrum_points)
        if exclude_points is not None:
            config.exclude_points = read_exclude_points(exclude_points)
        return config

    @classmethod
    def from_case(cls, name: str) -> "Config":
        """The configuration of an example case (see :func:`umwm.list_cases`),
        given by its name or a unique prefix of it, such as ``"titan"``."""
        from .cases import get_case

        return cls.from_namelist(get_case(name).namelist)

    @classmethod
    def from_namelist_text(cls, text: str) -> "Config":
        """Read the configuration from the text of a ``main.nml`` file."""
        groups = namelist.loads(text)
        config = cls()
        missing = [name for name in GROUPS if name not in groups]
        if missing:
            raise ConfigError("namelist groups missing: " + ", ".join(f"&{name.upper()}" for name in missing))
        for group_name, variables in groups.items():
            if group_name not in GROUPS:
                continue  # extra groups are ignored, as the model does
            for name, value in variables.items():
                config.set(f"{group_name}.{name}", value)
        return config

    @classmethod
    def from_dict(cls, data: dict) -> "Config":
        """Build the configuration from nested dictionaries, as returned by :meth:`to_dict`."""
        config = cls()
        for group_name, variables in data.items():
            if group_name == "spectrum_points":
                config.spectrum_points = [SpectrumPoint(**point) if isinstance(point, dict) else point for point in variables]
            elif group_name == "exclude_points":
                config.exclude_points = [tuple(point) for point in variables]
            elif group_name in GROUPS:
                for name, value in variables.items():
                    config.set(f"{group_name}.{name}", value)
            else:
                raise ConfigError(f"unknown configuration group {group_name!r}")
        return config

    @classmethod
    def from_json(cls, text: str) -> "Config":
        """Build the configuration from a JSON document, as returned by :meth:`to_json`."""
        return cls.from_dict(json.loads(text))

    def copy(self) -> "Config":
        """A deep copy of the configuration."""
        return copy.deepcopy(self)

    # ---- access -------------------------------------------------------------

    def get(self, key: str) -> Any:
        """The value of a variable given by a dotted ``group.name`` key."""
        group_name, name = self._split_key(key)
        return getattr(getattr(self, group_name), name)

    def set(self, key: str, value: Any) -> None:
        """Set a variable given by a dotted ``group.name`` key.

        The value is converted to the type of the variable; strings are
        parsed as namelist values, so values from the command line work.
        """
        group_name, name = self._split_key(key)
        group = getattr(self, group_name)
        for item in fields(group):
            if item.name == name:
                setattr(group, name, _convert(value, item.type, key))
                return
        raise ConfigError(f"unknown variable {name!r} in namelist group &{group_name.upper()}")

    def update(self, values=None, **kwargs) -> "Config":
        """Set several variables from a mapping of dotted keys, or keyword arguments
        with ``__`` in place of the dot, and return the configuration."""
        for key, value in dict(values or {}, **kwargs).items():
            self.set(key.replace("__", "."), value)
        return self

    def _split_key(self, key: str) -> tuple[str, str]:
        if "." not in key:
            raise ConfigError(f"configuration key {key!r} must be of the form group.name")
        group_name, name = key.lower().split(".", 1)
        if group_name not in GROUPS:
            raise ConfigError(f"unknown namelist group {group_name!r}; expected one of {', '.join(GROUPS)}")
        return group_name, name

    def __getitem__(self, key: str) -> Any:
        return self.get(key)

    def __setitem__(self, key: str, value: Any) -> None:
        self.set(key, value)

    def items(self) -> Iterator[tuple[str, Any]]:
        """Iterate over ``(group.name, value)`` pairs of all namelist variables."""
        for group_name in GROUPS:
            group = getattr(self, group_name)
            for item in fields(group):
                yield f"{group_name}.{item.name}", getattr(group, item.name)

    def diff(self, other: "Config") -> dict[str, tuple[Any, Any]]:
        """The variables that differ from ``other``, as a mapping of dotted
        keys (and ``spectrum_points``, ``exclude_points``) to ``(this, other)``."""
        differences = {}
        for (key, value), (_, other_value) in zip(self.items(), other.items()):
            if value != other_value:
                differences[key] = (value, other_value)
        for name in ("spectrum_points", "exclude_points"):
            if getattr(self, name) != getattr(other, name):
                differences[name] = (getattr(self, name), getattr(other, name))
        return differences

    def describe(self, group: str | None = None) -> str:
        """The variables with their values, descriptions and units, as text;
        optionally only those of one namelist ``group``."""
        lines = []
        for group_name in GROUPS:
            if group is not None and group_name != group.lower():
                continue
            lines.append(f"&{group_name.upper()}")
            for item in fields(getattr(self, group_name)):
                value = namelist.format_value(getattr(getattr(self, group_name), item.name))
                if len(value) > 40:
                    value = value[:37] + "..."
                description = item.metadata.get("description", "")
                units = item.metadata.get("units")
                lines.append(f"  {item.name:<14} = {value:<40} {description}" + (f" [{units}]" if units else ""))
        if group is not None and not lines:
            raise ConfigError(f"unknown namelist group {group!r}; expected one of {', '.join(GROUPS)}")
        return "\n".join(lines)

    # ---- derived quantities -------------------------------------------------

    @property
    def start_time(self) -> datetime:
        return parse_time(self.domain.starttimestr)

    @start_time.setter
    def start_time(self, value) -> None:
        self.domain.starttimestr = format_time(value)

    @property
    def stop_time(self) -> datetime:
        return parse_time(self.domain.stoptimestr)

    @stop_time.setter
    def stop_time(self, value) -> None:
        self.domain.stoptimestr = format_time(value)

    @property
    def reference_time(self) -> datetime:
        return parse_time(self.domain.reftimestr)

    @property
    def time_step(self) -> timedelta:
        """The global time step as a timedelta."""
        return timedelta(seconds=round(self.domain.dtg))

    def step_times(self, include_start: bool = True) -> list[datetime]:
        """The model times at the end of each global time step, optionally with the start time."""
        times = [self.start_time] if include_start else []
        time = self.start_time
        while time < self.stop_time:
            time = time + self.time_step
            times.append(time)
        return times

    @property
    def forcing_from_file(self) -> bool:
        """Whether any forcing field is read from the input files."""
        forcing = self.forcing
        return any([forcing.winds, forcing.currents, forcing.air_density, forcing.water_density, forcing.seaice])

    # ---- validation ---------------------------------------------------------

    def validate(self) -> None:
        """Check the configuration for errors the model would reject at startup.

        This mirrors the checks of ``config_type % validate`` in the Fortran
        source, so that an invalid configuration can be reported as an
        exception before the model is started.
        """
        errors = []
        domain, physics, grid, forcing, constant, output = (
            self.domain, self.physics, self.grid, self.forcing, self.forcing_constant, self.output)
        if domain.mm < 3 or domain.nm < 3:
            errors.append("mm and nm must be > 2")
        if domain.om < 3:
            errors.append("om must be >= 3")
        if domain.pm % 4 != 0:
            errors.append("pm must be divisible by 4")
        if domain.fmin <= 0 or domain.fmax <= 0 or domain.fprog <= 0:
            errors.append("fmin, fmax and fprog must be > 0")
        if domain.fmin >= domain.fmax:
            errors.append("fmin must be < fmax")
        if domain.dtg <= 0:
            errors.append("dtg must be > 0")
        for name in ("starttimestr", "stoptimestr", "reftimestr"):
            try:
                parse_time(getattr(domain, name))
            except ValueError:
                errors.append(f"{name} must be YYYY-MM-DD HH:MM:SS")
        try:
            if self.stop_time <= self.start_time:
                errors.append("stoptimestr must be after starttimestr")
        except ValueError:
            pass
        if physics.z <= 0:
            errors.append("z must be > 0")
        if physics.gustiness < 0:
            errors.append("gustiness must be positive")
        if physics.dmin <= 0:
            errors.append("dmin must be > 0")
        if physics.fcut_fac <= 0:
            errors.append("fcut_fac must be > 0")
        if not grid.gridfromfile and (grid.delx <= 0 or grid.dely <= 0):
            errors.append("delx and dely must be > 0")
        if not grid.topofromfile and grid.dpt <= 0:
            errors.append("dpt must be > 0")
        if constant.rhoa0 <= 0 or constant.rhow0 <= 0:
            errors.append("rhoa0 and rhow0 must be > 0")
        if not all(0 <= value <= 1 for value in (constant.fice0, constant.fice_lth, constant.fice_uth)):
            errors.append("sea-ice fractions must be in [0, 1]")
        if constant.fice_lth > constant.fice_uth:
            errors.append("fice_lth must be <= fice_uth")
        for name in ("outgrid", "outspec", "outrst"):
            if getattr(output, name) not in OUTPUT_INTERVALS:
                errors.append(f"{name} must be one of {', '.join(map(str, OUTPUT_INTERVALS))}")
        x_interior = (1 <= output.xpl <= domain.mm) if domain.isglobal else (2 <= output.xpl <= domain.mm - 1)
        if not (x_interior and 2 <= output.ypl <= domain.nm - 1):
            errors.append("xpl and ypl must be interior grid points (the boundary rows, and the "
                          "boundary columns of a regional domain, are land)")
        if not any(depth > 0 for depth in self.stokes.depths):
            errors.append("at least one positive Stokes depth is required")
        if self.spectrum_points:
            if len({point.by_index for point in self.spectrum_points}) > 1:
                errors.append("spectrum points must all be given by indices or all by lon/lat")
            for point in self.spectrum_points:
                if point.by_index and not (2 <= point.x <= domain.mm - 1 and 2 <= point.y <= domain.nm - 1):
                    errors.append(f"spectrum point {point.name!r} is outside the interior of the grid")
        if errors:
            raise ConfigError("invalid configuration: " + "; ".join(errors))

    # ---- output -------------------------------------------------------------

    def to_dict(self) -> dict:
        """The configuration as nested dictionaries."""
        data = OrderedDict()
        for group_name in GROUPS:
            data[group_name] = OrderedDict(dataclasses.asdict(getattr(self, group_name)))
        data["spectrum_points"] = [
            {key: value for key, value in dataclasses.asdict(point).items() if value is not None}
            for point in self.spectrum_points
        ]
        data["exclude_points"] = [list(point) for point in self.exclude_points]
        return data

    def to_json(self, indent: int = 2) -> str:
        """The configuration as a JSON document."""
        return json.dumps(self.to_dict(), indent=indent)

    def to_namelist(self, comments: bool = True) -> str:
        """The text of the ``main.nml`` file for this configuration."""
        groups = OrderedDict()
        notes = {}
        for group_name in GROUPS:
            group = getattr(self, group_name)
            groups[group_name] = OrderedDict()
            for item in fields(group):
                value = getattr(group, item.name)
                groups[group_name][item.name] = value
                if comments:
                    description = item.metadata.get("description", "")
                    units = item.metadata.get("units")
                    notes[(group_name, item.name)] = f"{description} [{units}]" if units else description
        return namelist.dumps(groups, comments=notes)

    def namelist_text(self, source=None) -> str:
        """The text of the ``main.nml`` file for this configuration: the text of
        the namelist file ``source``, with its comments, if it holds this
        configuration, and :meth:`to_namelist` otherwise."""
        if source is not None:
            try:
                text = Path(source).read_text(encoding="utf-8")
                if list(Config.from_namelist_text(text).items()) == list(self.items()):
                    return text
            except (OSError, ConfigError, NamelistError):
                pass
        return self.to_namelist()

    def spectrum_points_text(self) -> str:
        """The text of the ``spectrum.nml`` file listing the spectrum output points."""
        if not self.spectrum_points:
            return ""
        by_index = {point.by_index for point in self.spectrum_points}
        if len(by_index) > 1:
            raise ConfigError("spectrum points must all be given by indices or all by lon/lat")
        lines = ["XY" if by_index.pop() else "LL"]
        for point in self.spectrum_points:
            if point.by_index:
                lines.append(f"{point.x} {point.y} {point.name}")
            else:
                lines.append(f"{point.lon!r} {point.lat!r} {point.name}")
        return "\n".join(lines) + "\n"

    def exclude_points_text(self) -> str:
        """The text of the ``exclude.nml`` file listing cells of enclosed seas to exclude."""
        return "".join(f"{m} {n}\n" for m, n in self.exclude_points)

    def write(self, path, spectrum_points: bool = True, exclude_points: bool = True, source=None) -> Path:
        """Write ``main.nml`` to ``path`` and, next to it, ``spectrum.nml`` and
        ``exclude.nml`` when points are configured. Returns the path. The
        namelist file ``source``, such as the one the configuration was read
        from, is copied with its comments if it holds this configuration."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.namelist_text(source), encoding="utf-8")
        if spectrum_points and self.spectrum_points:
            (path.parent / "spectrum.nml").write_text(self.spectrum_points_text(), encoding="utf-8")
        if exclude_points and self.exclude_points:
            (path.parent / "exclude.nml").write_text(self.exclude_points_text(), encoding="utf-8")
        return path

    def __str__(self) -> str:
        return self.to_namelist()


# ---- helpers ----------------------------------------------------------------

def parse_time(text: str) -> datetime:
    """Parse a ``YYYY-MM-DD HH:MM:SS`` time string, accepting ``_`` or ``T`` as the separator."""
    text = str(text).strip().replace("_", " ").replace("T", " ")
    return datetime.strptime(text, TIME_FORMAT)


def format_time(value) -> str:
    """Format a datetime, or a string accepted by :func:`parse_time`, as ``YYYY-MM-DD HH:MM:SS``."""
    if isinstance(value, datetime):
        return value.strftime(TIME_FORMAT)
    return parse_time(value).strftime(TIME_FORMAT)


def parse_value(text: str) -> Any:
    """Parse a value given on the command line as a namelist value or a list of them."""
    values = namelist.parse_values(text)
    if len(values) == 1:
        return values[0]
    if not values:
        return ""
    return values


def _convert(value: Any, kind: Any, key: str) -> Any:
    """Convert ``value`` to the annotated type ``kind`` of a configuration field."""
    kind = str(kind)
    if isinstance(value, str) and not kind.startswith("str"):
        value = parse_value(value)
    try:
        if kind.startswith("bool"):
            if isinstance(value, str):
                parsed = namelist.parse_scalar(value)
                if not isinstance(parsed, bool):
                    raise ValueError(value)
                return parsed
            return bool(value)
        if kind.startswith("int"):
            if isinstance(value, float) and not value.is_integer():
                raise ValueError(value)
            return int(value)
        if kind.startswith("float"):
            return float(value)
        if kind.startswith("str"):
            if isinstance(value, datetime):
                return format_time(value)
            return str(value)
        if kind.startswith("list"):
            if not isinstance(value, (list, tuple)):
                value = [value]
            return [float(item) for item in value]
    except (TypeError, ValueError) as error:
        raise ConfigError(f"invalid value {value!r} for {key} (expected {kind})") from error
    raise ConfigError(f"unsupported type {kind} for {key}")


def namelist_path(path) -> Path:
    """The ``main.nml`` file of ``path``: the file itself, or the
    ``namelists/main.nml`` (or ``main.nml``) of a directory."""
    path = Path(path)
    if path.is_dir():
        for candidate in (path / "namelists" / "main.nml", path / "main.nml"):
            if candidate.is_file():
                return candidate
        raise FileNotFoundError(f"no namelists/main.nml in {path}")
    if not path.is_file():
        raise FileNotFoundError(f"namelist file not found: {path}")
    return path


def read_spectrum_points(path) -> list[SpectrumPoint]:
    """Read the spectrum output points from a ``spectrum.nml`` file."""
    lines = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines()]
    lines = [line for line in lines if line]
    if not lines:
        return []
    mode = lines[0].upper()
    if mode not in ("XY", "LL"):
        raise ConfigError(f"{path}: the first line must be XY or LL, not {lines[0]!r}")
    points = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 3:
            raise ConfigError(f"{path}: cannot read a point from {line!r}")
        if mode == "XY":
            points.append(SpectrumPoint(name=parts[2], x=int(parts[0]), y=int(parts[1])))
        else:
            points.append(SpectrumPoint(name=parts[2], lon=float(parts[0]), lat=float(parts[1])))
    return points


def read_exclude_points(path) -> list[tuple[int, int]]:
    """Read the cells of enclosed seas to exclude from an ``exclude.nml`` file."""
    points = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) >= 2:
            points.append((int(parts[0]), int(parts[1])))
    return points
