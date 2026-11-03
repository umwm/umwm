# UMWM example cases

A gallery of cases, ready to run: idealized cases, where the wind and the
other parameters are constant, and most cases need no input files, and two
real cases from the ERA5 reanalysis, which prepare their input files from
it when they run. Each case is a directory with its `main.nml` namelist,
whose leading comments describe the case, and the `spectrum.nml` points of
its spectrum output; `cases.txt` lists the cases. The `seamount` and
`eddy-field` cases bring their input files, a bathymetry and steady
currents, in their `input/` directory, with the `make_input.py` script that
writes them.

| Case | What it shows | Grid | Run time |
| --- | --- | --- | ---: |
| [`fetch-limited`](fetch-limited/main.nml) | Growth with fetch to full development, 20 m/s wind | 1000 km at 10 km, 48 h | 25 s |
| [`duration-limited`](duration-limited/main.nml) | Growth with time of a uniform sea, 20 m/s wind | one sea cell, 48 h | 1 s |
| [`jonswap`](jonswap/main.nml) | Fetch-limited growth over the fetches of JONSWAP, 10 m/s wind | 100 km at 1 km, 12 h | 15 s |
| [`wave-tank`](wave-tank/main.nml) | Waves in a laboratory wind-wave tank, 10 m/s wind | 10 m at 20 cm, 2 min | 15 s |
| [`titan-ontario-lacus`](titan-ontario-lacus/main.nml) | Methane-ethane waves on a lake of Titan, 1.5 m/s wind | 230 x 75 km at 5 km, 24 h | 60 s |
| [`mars-jezero-lake`](mars-jezero-lake/main.nml) | Waves on the ancient lake of Jezero crater, Mars, 5 m/s wind | 40 x 40 km at 1 km, 12 h | 20 s |
| [`seamount`](seamount/main.nml) | Waves over a Gaussian seamount, 20 m/s wind | 100 x 160 km at 2.5 km, 12 h | 45 s |
| [`eddy-field`](eddy-field/main.nml) | Waves through steady mesoscale eddies, 20 m/s wind | 200 x 300 km at 4 km, 15 h | 70 s |
| [`era5-helene`](era5-helene/main.nml) | Hurricane Helene in the Gulf of Mexico, from ERA5, with open boundaries | 79 x 57 at 0.25 degrees, 24 h | 75 s |
| [`era5-global`](era5-global/main.nml) | The global ocean, from ERA5 | 360 x 157 at 1 degree, 24 h | 15 min |

The run times are approximate, for a serial run on one core; those of the
ERA5 cases include the preparation of their input files, about 15 s.

## Running a case

With the `umwm` executable, `--case NAME` runs a case in the current
directory, which needs an `output/` directory for the output files, as the
UMWM directory has. A unique prefix of the name is enough:

```
./umwm --list-cases
./umwm --case jonswap
./umwm --case titan                  # titan-ontario-lacus
mpiexec -n 2 ./umwm --case jonswap   # in parallel
```

Each case writes its output to `output/` of the current directory, over
the output of any earlier run, so run each case in its own directory to
keep their output apart. The executable finds the cases in `examples/` of
the current directory or of the directory of the executable, also when it
is started through `PATH` or a symbolic link, or in `$UMWM_EXAMPLES`:

```
mkdir -p runs/jonswap/output && cd runs/jonswap
../../umwm --case jonswap
```

`./umwm path/to/main.nml` runs any other namelist, such as a copy of a case
that you changed; its `spectrum.nml` is read from the same directory. A case
reads its input files from its own `input/` directory, while other runs read
them from `input/` of the current directory.

In parallel, the processes share the sea cells and the results are the same
as those of a serial run. The cases are small, so a few processes are
enough: `fetch-limited`, `jonswap` and `wave-tank`, four grid columns wide,
run on up to four processes, and `duration-limited`, a single sea cell, only
in serial.

With the [Python command-line tool](../python/README.md), `umwm run --case
NAME` runs a case in a run directory, by default a directory named after the
case, with `-s` to change any namelist variable:

```
umwm cases                                       # list the cases
umwm cases titan                                 # describe one
umwm run --case jonswap                          # in ./jonswap
umwm run --case jonswap runs/j15 -s forcing_constant.wspd0=15
umwm output jonswap --stats                      # summarize the output
umwm output jonswap --plot swh --save swh.png
umwm init --case mars runs/mars                  # a run directory to edit
umwm init --case seamount runs/sm                # with a copy of the input files
umwm config --case titan --diff defaults         # what the case changes
```

`umwm run --case` and `umwm init --case` copy the input files of a case into
the `input/` directory of the run directory. `umwm run --case` does not
replace a namelist or input files there that differ from those of the case,
such as those of a run directory created with `umwm init --case` and
edited, unless given `--force`; run such a directory with `umwm run DIR`.
From Python, `umwm.Config.from_case(name)` reads the configuration of a
case, `umwm.list_cases()` lists them, and `case.domain(workdir)` is a
domain of a case with its input files:

```python
import umwm

config = umwm.Config.from_case("jonswap")
config.forcing_constant.wspd0 = 15.0
with umwm.Domain(config, workdir="runs/jonswap15") as domain:
    domain.run()
    output = domain.output()
    spectra = domain.spectrum_output()      # {"fetch010km": Dataset, ...}
```

The ERA5 cases need the Python package, with the dependencies of ERA5
(`pip install '.[era5]'` in the UMWM repository), and network access to the
ARCO ERA5 store of Google; the `umwm` executable prepares their input files
with it, in `input/` of the current directory (see [ERA5 cases](#era5-cases)).

## Fetch-limited and duration-limited growth

The ocean cases (`fetch-limited`, `duration-limited` and `jonswap`) model
frequencies from 0.04 to 14 Hz in 52 bins, the frequency spacing of the
default 37 bins from 0.0313 to 2 Hz, up to the smallest phase speed of sea
water at 13.7 Hz (a 1.65 cm wavelength), in the gravity-capillary range. The
longer high-frequency tail barely changes the wave heights, but lowers the
drag coefficient at 20 m/s by about 5%, to 2.0 × 10<sup>-3</sup>.

**`fetch-limited`**: a steady, uniform 20 m/s wind blows offshore from a
straight coast over 1000 km of deep water. The domain is periodic in x
(`isGlobal = .true.`), so the coast and the wave field are uniform along it,
without the lateral boundaries that would drain the waves traveling at an
angle to the wind in a narrow channel; the wind blows toward +y from the
land row at y = 0. Grid row n (1-based) is at a fetch of (n - 1) × 10 km.

After about 36 hours the wave field is steady. The significant wave height
grows from 2.7 m at 10 km to 5.4 m at 100 km and 10.1 m at 990 km, and the
dominant period from 6.3 s to 10 s and 16 s; at the far end, the sea is close
to the Pierson-Moskowitz full development for this wind, with a significant
wave height of about 9.9 m and a peak period of 15.7 s (Pierson and
Moskowitz, 1964). Spectra are written hourly at 50, 200, 500 and 990 km.

**`duration-limited`**: the same 20 m/s wind starts blowing over a sea
without fetch limitation. The wave field is then the same everywhere, and
the model reduces to the integration of the source functions at a point: the
domain is a single sea cell, so large (10<sup>6</sup> km) that propagation out of
it is negligible. The significant wave height grows to 6.0 m in 6 hours,
9.6 m in 24 hours and 11.9 m in 48 hours, with the dominant period growing
from 11 s to 18 s. After about 30 hours of wind the sea is as developed as at
1000 km fetch in the fetch-limited case. The spectrum is written hourly;
after 48 hours its peak, at 0.057 Hz, is three bins above the lowest
frequency, so lower `fmin` to run the case longer.

**`jonswap`**: fetch-limited growth over the fetches of the Joint North Sea
Wave Project (JONSWAP; Hasselmann et al., 1973), with a 10 m/s wind over 100
km of deep water at 1 km resolution, in the same periodic layout as
`fetch-limited`. After about 10 hours the wave field is steady:

| Fetch | Hs, UMWM | Hs, JONSWAP | Peak frequency, UMWM | Peak frequency, JONSWAP |
| ---: | ---: | ---: | ---: | ---: |
| 10 km | 0.71 m | 0.51 m | 0.28 Hz | 0.35 Hz |
| 20 km | 0.88 m | 0.72 m | 0.25 Hz | 0.28 Hz |
| 40 km | 1.08 m | 1.02 m | 0.20 Hz | 0.22 Hz |
| 80 km | 1.32 m | 1.45 m | 0.18 Hz | 0.18 Hz |

where the JONSWAP values are the fetch laws gH<sub>s</sub>/U<sup>2</sup> = 1.6 × 10<sup>-3</sup>
X<sup>1/2</sup> and f<sub>p</sub>U/g = 3.5 X<sup>-0.33</sup>, with X = gx/U<sup>2</sup> for the fetch x.
The spectra written at these four fetches show the peaked spectrum shifting
to lower frequencies with fetch.

## Laboratory wave tank

**`wave-tank`**: a 10 m long tank with 0.5 m of fresh water, under a 10 m/s
wind, at 20 cm resolution. The wind speed is the wind at 10 m, to which
laboratory winds are customarily extrapolated with the logarithmic profile
for comparison with field measurements. The tank is periodic in x, of
unlimited width. The waves are steady along the tank within a minute: the
significant wave height grows from 1 cm at 1 m to 3.0 cm at 9.8 m, and the
peak frequency falls from 3.2 to 1.8 Hz. The frequencies span 1 to 20 Hz,
through the gravity-capillary range: the phase speed is smallest at 13.5 Hz,
for a wavenumber of (g ρ<sub>w</sub>/σ)<sup>1/2</sup> = 367 rad/m and a wavelength of
1.7 cm. The gridded output and the spectra along the tank are written every
10 s, every global time step (`outgrid = outspec = -1`).

In a tank the waves are far younger than a fully developed sea, and their
spectral peak lies far above the default cutoff of the prognostic range,
f<sub>c</sub> = 0.53 g/U ≈ 4 f<sub>PM</sub>, which is 0.52 Hz for this wind; with the
default, the whole spectrum of the tank would be in the diagnostic
equilibrium range, at nearly the same height, about 8 cm, all along the
tank. The case sets `fcut_fac = 10`, which raises the cutoff to 9.8 Hz.

## Waves on other worlds

The planetary cases use the parameters of Schneck et al. (2026), from their
PlanetWaves model (v1.0.0, [doi:10.5281/zenodo.18906357](https://doi.org/10.5281/zenodo.18906357)),
which give the air density from the ideal gas law and the properties of the
liquids from the DIPPR database and, for Titan, the TitanPool model:

| Parameter | Namelist | Earth, ocean cases | Titan, Ontario Lacus | Mars, Jezero (0.5 bar) |
| --- | --- | ---: | ---: | ---: |
| Gravity [m/s<sup>2</sup>] | `g` | 9.80665 | 1.352 | 3.71 |
| Atmosphere | | | N<sub>2</sub>, 1.5 atm, 92 K | CO<sub>2</sub>, 0.5 bar, 288 K |
| Air density [kg/m<sup>3</sup>] | `rhoa0` | 1.2 | 5.564 | 0.9188 |
| Air kinematic viscosity [m<sup>2</sup>/s] | `nu_air` | 1.56e-5 | 1.156e-6 | 1.576e-5 |
| Liquid | | sea water | CH<sub>4</sub>-C<sub>2</sub>H<sub>6</sub>-N<sub>2</sub>, 92 K | fresh water, 288 K |
| Liquid density [kg/m<sup>3</sup>] | `rhow0` | 1030 | 588.15 | 998.21 |
| Liquid kinematic viscosity [m<sup>2</sup>/s] | `nu_water` | 0.9e-6 | 8.084e-7 | 1.156e-6 |
| Surface tension [N/m] | `sfct` | 0.07 | 0.032766 | 0.073836 |
| Smallest phase speed [Hz, cm] | | 13.7, 1.65 | 3.3, 4.0 | 6.5, 2.8 |
| Frequencies [Hz] | `fmin`, `fmax` | 0.04 to 14 | 0.02 to 5 | 0.05 to 10 |
| Frequency bins | `om` | 52 | 49 | 47 |

The geometry of the lakes is idealized as rectangular basins of uniform depth.
As in the Earth cases, the frequencies extend through the gravity-capillary
range, past the smallest phase speed, at the wavenumber (g ρ<sub>w</sub>/σ)<sup>1/2</sup>;
with the lower gravity of Titan and Mars, it lies at lower frequencies and
longer wavelengths than on Earth. The number of bins keeps the frequency
spacing of the default 37 bins from 0.0313 to 2 Hz.

**`titan-ontario-lacus`**: Ontario Lacus, the largest lake of Titan's south
polar region, as a 230 × 75 km basin, 20 m deep, under a 1.5 m/s wind along
its length, the average storm wind over Titan's lakes (Detelich et al.,
2026). The low gravity and the dense atmosphere make the waves slow, long,
and high for their wind. After a day, the significant wave height grows from
0.23 m at 5 km to 0.46 m at 100 km, with dominant periods from 5.6 to 10 s,
in the range of 15 to 48 cm and 6 to 10.5 s found by Detelich et al. (2026)
for such winds; farther down the lake the waves are limited by the duration
of the wind rather than by the fetch, and are still growing slowly. Spectra
are written hourly at 25, 100 and 225 km.

**`mars-jezero-lake`**: the lake that filled Jezero crater on early Mars, as a
40 × 40 km basin, 100 m deep, under a 5 m/s wind, with the low-pressure (0.5
bar) atmosphere of Schneck et al. (2026). After about 8 hours the waves are
steady, with the significant wave height growing from 0.21 m at 1 km to
0.57 m at 40 km and dominant periods from 3.2 to 5.6 s. Spectra are written
hourly at 10, 20 and 40 km. For the high-pressure (2 bar) atmosphere of
Schneck et al. (2026):

```
umwm run --case mars runs/mars-2bar -s forcing_constant.rhoa0=3.6752 -s physics.nu_air=3.9405e-6
```

Schneck et al. (2026) modeled Jezero with winds of 2 to 3 m/s, close to the
threshold of wave growth there; try `-s forcing_constant.wspd0=3`.

## Bathymetry and currents

These two cases, like the other ocean cases, model frequencies from 0.04 to
14 Hz in 52 bins, and have a steady, uniform 20 m/s wind blowing offshore
from a straight coast, as in `fetch-limited`: the domain is periodic in x,
and the wind blows toward +y from the land row at y = 0.

**`seamount`**: a Gaussian seamount rises from a 1000 m deep floor to 20 m
below the surface at the center of a 100 × 160 km domain, at 80 km fetch,
with a width (the standard deviation of the Gaussian) of 20 km, at 2.5 km
resolution. The bathymetry is read from `input/umwm.gridtopo`
(`topoFromFile`). After about 8 hours the waves over the seamount are
steady. Over its top, where the 10 s waves feel the bottom, the significant
wave height is 4.7 m, 7% lower than the 5.1 m away from it, and downwind the
seamount casts a weak shadow, 2% lower 20 km beyond its top. With the broad
directional spread of a wind sea, refraction over the seamount does not
focus the waves into a higher sea in its lee. Spectra are written upwind of,
on, and downwind of the seamount, and at the same fetch away from it.

**`eddy-field`**: beyond 60 km from the coast of a 200 × 300 km domain, at 4
km resolution, the waves cross 14 mesoscale eddies of Gaussian shape, 6
rotating clockwise and 8 counterclockwise, with radii of 13 to 28 km and
speeds up to 0.7 m/s (rms 0.3 m/s), which do not evolve in time. The
currents are read once, from `input/umwmin.nc`, and held constant
(`currents` and `steady` in `&FORCING`). The currents refract and advect
the waves: across the eddies, the significant wave height varies
by -8% to +11% around its mean at each fetch (standard deviation 3%), in
streaks downwind of the eddies, where it would be uniform without currents.
After 15 hours the pattern is steady, while the waves at the far end are
still growing slowly. Spectra are written upwind of the eddies and across
the domain beyond them.

The `make_input.py` script of each case writes its input file on the grid of
its `main.nml`, from parameters at its top: the depths and the width of the
seamount, and the number, sizes, speeds, and random seed of the eddies.
After changing them, `python make_input.py` in the case directory writes the
input file again.

## ERA5 cases

These cases start from, and are forced by, the ERA5 reanalysis, on a
longitude-latitude grid (`lon0`, `lat0`, `dlon` and `dlat` in `&GRID`), with
the bathymetry of the ERA5 wave model (`topoFromFile`). With `source =
'era5'` in `&FORCING`, the input files are prepared from the analysis-ready,
cloud-optimized ERA5 store of Google Cloud before the model runs, by the
`umwm` executable, by `umwm run`, or beforehand by `umwm forcing era5`: the
bathymetry, the hourly wind, air density and sea ice, and the wave
partitions of ERA5, the wind sea and three swells, from which the model
rebuilds its initial spectrum (`waves`), with the significant wave height,
mean period, mean direction and directional spread of each.

```
mkdir -p runs/helene/output && cd runs/helene
../../umwm --case era5-helene                        # prepares input/, then runs
umwm run --case era5-global                          # in ./era5-global
umwm init --case era5-helene runs/milton -s domain.starttimestr="2024-10-09 00:00:00" \
  -s domain.stoptimestr="2024-10-10 00:00:00"        # another hurricane
```

**`era5-helene`**: hurricane Helene crossing the Gulf of Mexico on 26
September 2024, at the 0.25 degrees of ERA5, for 24 hours. The domain is
regional, and its open boundaries, the Yucatan Channel, the Straits of
Florida and the western Caribbean, take the spectrum of the ERA5 partitions
throughout the run. The preparation reads about 1 GB of ERA5 fields, as the
ARCO store has a chunk of the whole globe per hour and variable. ERA5
resolves the core of a hurricane poorly: its wind peaks at 32 m/s, far below
the 60 m/s of Helene, and the significant wave height at 10.7 m in ERA5 and
11.0 m in the model. Against ERA5 interpolated to the grid of the model, the
significant wave height of the model differs by:

| Time | From ERA5 | From rest |
| ---: | ---: | ---: |
| 0 h | 0.007 m (0.4%) | 1.78 m |
| 6 h | 0.26 m (13%) | 0.72 m |
| 12 h | 0.40 m (17%) | 0.62 m |
| 24 h | 0.68 m (28%) | 0.72 m |

(root mean square difference over the sea cells, and relative to the mean
of ERA5), most under the strongest winds, where the model grows the waves by
its own physics; within three cells of the open boundaries, it differs by
0.15 m or less throughout.

**`era5-global`**: the global ocean at 1 degree, from 78S to 78N, for 24
hours from the same time, with 38500 sea cells. The domain is periodic in
longitude and has no open boundaries: its edges, the rows at 78S and 78N,
are land. The preparation reads 0.24 GB of ERA5 fields, the wave partitions
at the start time only. It runs in about 15 minutes in serial, and in 5
minutes on 4 processes (`mpiexec -n 4 ./umwm --case era5-global`), with the
same results. Against ERA5, the significant wave height of the model
differs by 0.009 m (0.4%) at the start, 0.37 m (17%) after 6 hours and
0.51 m (23%) after 24 hours, and is lower by 0.2 m on average from 6 hours
on: the total wave energy decreases by 7% in the first hour and 12% in six
hours, and then holds, as the model damps the swell that travels across or
against the local wind (the damping of `sin_diss1`; Donelan et al., 2012)
much faster than the ERA5 wave model, most in the trade winds of the
southern hemisphere, which the long swell of the Southern Ocean crosses.

## Adding a case

A case is a directory with a complete `main.nml`, which sets every
namelist variable, and optionally a `spectrum.nml`, listed in `cases.txt`
with a one-line description. Its input files, if any, go in its `input/`
directory; forcing read from file should be steady (`forcing.steady`), so
that one file of the start time is enough, unless the case prepares its
input files from a data source (`forcing.source`), as the ERA5 cases do. The
test suite of the Python package checks that every case listed is complete
and valid, and that its input files match its grid; it runs the ERA5 cases
only with `UMWM_TEST_ERA5=1`.

## References

- Donelan, M. A., M. Curcic, S. S. Chen, and A. K. Magnusson, 2012: Modeling
  waves and wind stress, *J. Geophys. Res. Oceans*, **117**, C00J23,
  doi:10.1029/2011JC007787.
- Detelich, C. E., et al., 2026: Modeling the seasonality of wind-driven
  hydrocarbon waves in Titan's polar lakes, *J. Geophys. Res. Planets*,
  **131**(5), e2026JE009693, doi:10.1029/2026JE009693.
- Hasselmann, K., et al., 1973: Measurements of wind-wave growth and swell
  decay during the Joint North Sea Wave Project (JONSWAP), *Ergänzungsheft
  zur Deutschen Hydrographischen Zeitschrift*, Reihe A, **8**(12), 95 pp.
- Pierson, W. J., and L. Moskowitz, 1964: A proposed spectral form for fully
  developed wind seas based on the similarity theory of S. A. Kitaigorodskii,
  *J. Geophys. Res.*, **69**(24), 5181-5190, doi:10.1029/JZ069i024p05181.
- Schneck, U. G., et al., 2026: Modeling wind-driven waves on other planets:
  Applications to Mars, Titan, and exoplanets, *J. Geophys. Res. Planets*,
  **131**(4), e2025JE009490, doi:10.1029/2025JE009490.
