# University of Miami Wave Model (UMWM)

A tiny, fast, parallel spectral ocean wave model.

This is the reference implementation of UMWM, 
described by [Donelan et al. (2012)](https://doi.org/10.1029/2011JC007787),
with later improvements and bug fixes.
UMWM solves the wave energy balance equation on a curvilinear grid.
It has been used to simulate:

* Global swell and windsea
* Waves in coastal and hurricane conditions
* Wave-induced material transport (Stokes drift)
* Ancient Martian seas and methane lakes on Titan
* Waves in laboratory settings such as wave tanks

UMWM was initially designed with a primary goal of accurate and conservative
momentum coupling with atmosphere and ocean circulation models.
That design goal remains a priority.
Further, UMWM takes a highly simplified approach to nonlinear downshifting
of wave energy; this allows it to run signinficantly faster than other
spectral wave models.

## Getting started

### Getting the code

```
git clone https://github.com/umwm/umwm
```

### System dependencies

* `make`
* A recent Fortran compiler (known to work with GNU, Intel, Cray, and IBM)
* NetCDF for I/O
* MPI for parallel processing (optional)

### Building UMWM

Type `make` to build the serial model and auxiliary tools:

```
make
```

The default serial compiler is `gfortran`. Select another compiler by setting
`FC`:

```
make FC=ifx
make FC=flang
```

The default `FCFLAGS` are aggressive GNU-style optimization flags:
`-Ofast -march=native -ffast-math -funroll-loops -Wall`. Override `FCFLAGS`
when using a compiler that needs different optimization options:

```
make FC=ifx FCFLAGS="-O3 -xHost -ipo -fp-model fast=2"
```

Build with MPI by setting `MPI=yes`. This defines `-DMPI` and defaults to the
`mpif90` wrapper unless `FC` is set explicitly:

```
make MPI=yes
make MPI=yes FC=mpiifx
```

NetCDF Fortran flags are discovered with `nf-config` by default. If you have
multiple NetCDF installs, point to the matching `nf-config` for your compiler:

```
make NF_CONFIG=/path/to/nf-config
```

As a fallback, set `NETCDF` to an install prefix with `include` and `lib`
subdirectories, or provide the flags manually:

```
make NETCDF=/path/to/netcdf
make NETCDF_FFLAGS="-I/path/include" NETCDF_FLIBS="-L/path/lib -lnetcdff -lnetcdf"
```

Use `make print-config` to see the resolved compiler and NetCDF settings.
Executable `umwm` will be built in the top-level directory and auxiliary tool
executables will be built in `tools/`.
PDF documentation can be built separately with `make docs`.

### Running UMWM

Running in serial mode:

```
./umwm
```

Running in parallel, for example on 16 cores:

```
mpiexec -n 16 ./umwm
```

You can read the full technical reference doc [here](DOCS.md).

## Papers

A list of papers that use UMWM is [here](PAPERS.md).

## Thanks

UMWM development has been supported by [NSF Award 1745384](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1745384),
the [Gulf of Mexico Research Initiative](http://gulfresearchinitiative.org/),
and the [National Oceanographic Partnership Program](https://www.nopp.org/).

UMWM has also been improved by a number of [open source contributors](CONTRIBUTORS.md).
