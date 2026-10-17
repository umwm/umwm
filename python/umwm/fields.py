"""Names, descriptions, and units of the domain fields accessible from Python.

The names are those of the gridded and spectrum output files of the model
where they exist; the Fortran component names are accepted as aliases by
:meth:`umwm.Domain.field`.
"""

from __future__ import annotations

from collections import OrderedDict

__all__ = ["GRID_OUTPUT_FIELDS", "STOKES_FIELDS", "SPECTRAL_FIELDS", "KINEMATICS_FIELDS",
           "GRID_FIELDS", "FORCING_FIELDS", "FIELDS"]

GRID_OUTPUT_FIELDS = OrderedDict([
    # name: (description, units), in the order of the gridded output file
    ("wspd", ("wind speed", "m/s")),
    ("wdir", ("wind direction", "rad")),
    ("uc", ("ocean current, x-component", "m/s")),
    ("vc", ("ocean current, y-component", "m/s")),
    ("rhoa", ("air density", "kg/m^3")),
    ("rhow", ("water density", "kg/m^3")),
    ("fice", ("seaice fraction", "non-dimensional")),
    ("psim", ("universal stability function for momentum", "non-dimensional")),
    ("momx", ("momentum, x-component", "kgm/s")),
    ("momy", ("momentum, y-component", "kgm/s")),
    ("cgmxx", ("cg*momentum, xx-component", "kgm^2/s^2")),
    ("cgmxy", ("cg*momentum, xy-component", "kgm^2/s^2")),
    ("cgmyy", ("cg*momentum, yy-component", "kgm^2/s^2")),
    ("shelt", ("sheltering coefficient", "non-dimensional")),
    ("epsx_atm", ("wave energy growth flux, x-component", "kg/s^3")),
    ("epsy_atm", ("wave energy growth flux, y-component", "kg/s^3")),
    ("epsx_ocn", ("wave energy dissipation flux, x-component", "kg/s^3")),
    ("epsy_ocn", ("wave energy dissipation flux, y-component", "kg/s^3")),
    ("taux_form", ("form drag, x-component", "n/m^2")),
    ("tauy_form", ("form drag, y-component", "n/m^2")),
    ("taux_form_1", ("form drag, part 1, x-component", "n/m^2")),
    ("tauy_form_1", ("form drag, part 1, y-component", "n/m^2")),
    ("taux_form_2", ("form drag, part 2, x-component", "n/m^2")),
    ("tauy_form_2", ("form drag, part 2, y-component", "n/m^2")),
    ("taux_form_3", ("form drag, part 3, x-component", "n/m^2")),
    ("tauy_form_3", ("form drag, part 3, y-component", "n/m^2")),
    ("taux_skin", ("skin drag, x-component", "n/m^2")),
    ("tauy_skin", ("skin drag, y-component", "n/m^2")),
    ("taux_diag", ("diagnostic form drag, x-component", "n/m^2")),
    ("tauy_diag", ("diagnostic form drag, y-component", "n/m^2")),
    ("taux_ocn", ("momentum flux from breaking waves to ocean top, x-component", "n/m^2")),
    ("tauy_ocn", ("momentum flux from breaking waves to ocean top, y-component", "n/m^2")),
    ("taux_bot", ("momentum flux from waves to ocean bottom, x-component", "n/m^2")),
    ("tauy_bot", ("momentum flux from waves to ocean bottom, y-component", "n/m^2")),
    ("taux_snl", ("momentum flux due to snl, x-component", "n/m^2")),
    ("tauy_snl", ("momentum flux due to snl, y-component", "n/m^2")),
    ("tailatmx", ("atmosphere tail stress part, x-component", "n/m^2")),
    ("tailatmy", ("atmosphere tail stress part, y-component", "n/m^2")),
    ("tailocnx", ("ocean tail stress part, x-component", "n/m^2")),
    ("tailocny", ("ocean tail stress part, y-component", "n/m^2")),
    ("cd", ("drag coefficient of air", "non-dimensional")),
    ("ust", ("friction velocity of air", "m/s")),
    ("swh", ("significant wave height", "m")),
    ("mss", ("mean-squared slope", "non-dimensional")),
    ("mwp", ("mean wave period", "s")),
    ("mwl", ("mean wavelength", "m")),
    ("mwd", ("mean wave direction", "rad")),
    ("dwp", ("dominant wave period", "s")),
    ("dwl", ("dominant wavelength", "m")),
    ("dwd", ("dominant wave direction", "rad")),
    ("dcp0", ("dominant phase speed, intrinsic", "m/s")),
    ("dcg0", ("dominant group speed, intrinsic", "m/s")),
    ("dcp", ("dominant phase speed", "m/s")),
    ("dcg", ("dominant group speed", "m/s")),
    ("physics_time_step", ("Physics time step", "s")),
])
"""Gridded output fields, one value per sea cell, gathered to the grid by :meth:`umwm.Domain.grid_field`."""

OCEAN_STRESS_FIELDS = (
    "taux_ocn", "tauy_ocn", "taux_bot", "tauy_bot", "taux_snl", "tauy_snl", "tailocnx", "tailocny",
    "epsx_atm", "epsy_atm", "epsx_ocn", "epsy_ocn",
    "taux_form_1", "tauy_form_1", "taux_form_2", "tauy_form_2", "taux_form_3", "tauy_form_3",
)
"""Fields evaluated by the ocean stress routine, which the model runs at gridded output times."""

STOKES_FIELDS = OrderedDict([
    ("u_stokes", ("stokes drift x-component", "m/s")),
    ("v_stokes", ("stokes drift y-component", "m/s")),
    ("d_stokes", ("stokes drift e-folding depth", "m")),
])
"""Stokes drift fields; the velocity components have a depth dimension."""

SPECTRAL_FIELDS = OrderedDict([
    ("F", ("surface elevation variance spectrum", "m^4/rad")),
    ("Sin", ("wind input/export rate", "1/s")),
    ("Sds", ("wave breaking decay rate (spilling part; multiply by cothkd)", "1/s")),
    ("Snl", ("nonlinear wave-wave interaction absolute source/sink tendency", "m^4/rad/s")),
    ("Sdt", ("turbulent dissipation rate", "1/s")),
    ("Sdv", ("viscous dissipation rate", "1/s")),
    ("Sbf", ("combined bottom friction and bottom percolation rate", "1/s")),
    ("Sice", ("sea ice attenuation rate", "1/s")),
])
"""Spectral fields on the cells, with frequency (and direction) as leading dimensions."""

KINEMATICS_FIELDS = OrderedDict([
    ("frequency", ("frequency", "Hz")),
    ("direction", ("direction", "rad")),
    ("wavenumber", ("wavenumber", "rad/m")),
    ("cp0", ("intrinsic phase speed", "m/s")),
    ("cg0", ("intrinsic group speed", "m/s")),
    ("kdk", ("spectral weight k dk", "rad^2/m^2")),
    ("k3dk", ("spectral weight k^3 dk", "rad^4/m^4")),
    ("cothkd", ("plunging breaker factor coth(0.2 k d)", "non-dimensional")),
    ("oc", ("highest prognostic frequency bin", "index")),
])
"""Spectral coordinates and wave kinematics."""

GRID_FIELDS = OrderedDict([
    ("lon", ("longitude", "degrees east")),
    ("lat", ("latitude", "degrees north")),
    ("x", ("x coordinate", "m")),
    ("y", ("y coordinate", "m")),
    ("depth", ("ocean depth", "m")),
    ("dx", ("grid spacing in x", "m")),
    ("dy", ("grid spacing in y", "m")),
    ("area", ("cell area", "m^2")),
    ("curvature", ("grid rotation", "rad")),
    ("seamask", ("seamask (1 sea, 0 land)", "non-dimensional")),
    ("nproc", ("process owning the cell", "non-dimensional")),
    ("ii", ("cell index of each grid point (1-based)", "index")),
    ("mi", ("x index of each cell (1-based)", "index")),
    ("ni", ("y index of each cell (1-based)", "index")),
])
"""Grid geometry on the (x, y) grid and the index maps between grid points and cells."""

FORCING_FIELDS = OrderedDict([
    ("wspd", ("wind speed", "m/s")),
    ("wdir", ("wind direction", "rad")),
    ("uc", ("ocean current, x-component", "m/s")),
    ("vc", ("ocean current, y-component", "m/s")),
    ("rhoa", ("air density", "kg/m^3")),
    ("rhow", ("water density", "kg/m^3")),
    ("rhorat", ("air to water density ratio", "non-dimensional")),
    ("fice", ("seaice fraction", "non-dimensional")),
    ("uw", ("wind x-component on the grid", "m/s")),
    ("vw", ("wind y-component on the grid", "m/s")),
])
"""Forcing fields on the cells (all grid points, including land) and the native-grid wind snapshots."""

FIELDS = OrderedDict()
for _group in (GRID_OUTPUT_FIELDS, STOKES_FIELDS, SPECTRAL_FIELDS, KINEMATICS_FIELDS, GRID_FIELDS, FORCING_FIELDS):
    for _name, _meta in _group.items():
        FIELDS.setdefault(_name, _meta)
"""All documented fields: name to (description, units)."""


def attrs(name: str) -> dict:
    """The xarray attributes (description, long_name, units) of a field."""
    description, units = FIELDS.get(name, (name, None))
    result = {"description": description, "long_name": description}
    if units:
        result["units"] = units
    return result
