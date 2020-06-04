# University of Miami Wave Model (UMWM)

[![Build Status](https://travis-ci.org/umwm/umwm.svg?branch=master)](https://travis-ci.org/umwm/umwm)
[![GitHub issues](https://img.shields.io/github/issues/umwm/umwm.svg)](https://github.com/umwm/umwm/issues)

A third-generation spectral ocean wave model.

This is the reference implementation of UMWM, 
described by [Donelan et al. (2012)](https://github.com/milancurcic/publications/blob/master/Donelan_etal_JGR2012.pdf).
UMWM solves the wave energy balance equation on a curvilinear grid.
It has been used to simulate:

* Global swell and windsea
* Waves in coastal and hurricane conditions
* Wave-induced material transport (Stokes drift)
* Ancient Martian seas and methane lakes on Titan
* Waves in laboratory settings such as wave tanks

## Getting started

### Getting the code

```
git clone https://github.com/umwm/umwm
```

### System dependencies

* `make`
* GNU, Intel, or Cray Fortran compiler
* NetCDF for I/O
* MPI for parallel processing (optional)

### Building UMWM

Edit the following variables in the top-level Makefile:

* `FC`: Fortran compiler (e.g. `mpif90` for parallel builds)
* `FCFLAGS`: Flags to pass to the Fortran compiler
* `CPPFLAGS`: Pre-processor flags -- set to `-DMPI` if building for parallel execution, 
and leave blank for serial builds.

Path to the NetCDF library must be set as `NETCDF` environment variable.
If your library is installed in non-standard directories (something 
other than `$NETCDF/lib` for library files and `$NETCDF/include` for modules)
edit the `NETCDFLIB` and `NETCDFINC` variables in `src/Makefile`.

Type `make`. Executable `umwm` will be built in the top-level directory. 
Auxilliary tools executables will be built in `tools/`. 
Documentation will be built in `docs/`.

### Running UMWM

Running in serial mode:

```
./umwm
```

Running in parallel, for example on 16 cores:

```
mpiexec -n 16 ./umwm
```

Read the [docs](docs) for more information.

## Publications

See [publications](PUBLICATIONS.md) for a full list of publications.

## Thanks

UMWM development is currently supported by [NSF Award 1745384](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1745384).
It has previously been supported by the [Gulf of Mexico Research Initiative](http://gulfresearchinitiative.org/)
and the [National Oceanographic Partnership Program](https://www.nopp.org/).

UMWM has also been improved by a number of [open source contributors](CONTRIBUTORS.md).
