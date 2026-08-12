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

### Get the code

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

### Running tests

Run the test suite from the top-level directory:

```
make test
```

### Building with fpm

[Fortran Package Manager (fpm)](https://fpm.fortran-lang.org/) 0.13 or newer
can build the MPI-enabled model and run its test suite. The active compiler
must have matching MPI and netCDF Fortran installations discoverable through
their compiler wrappers and `pkg-config` files.

With GNU Fortran, OpenMPI, and a system netCDF installation:

```
fpm build --compiler gfortran
fpm test --compiler gfortran --runner ./tests/fpm-mpi-test-runner.sh
```

With Intel Fortran, Intel MPI, and a netCDF-Fortran installation built with
`ifx`, first initialize oneAPI and expose that netCDF installation:

```
source ~/intel/oneapi/setvars.sh
export PKG_CONFIG_PATH=/path/to/intel-netcdf/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/path/to/intel-netcdf/lib:$LD_LIBRARY_PATH
fpm build --compiler ifx --build-dir build-ifx
fpm test --compiler ifx --build-dir build-ifx \
  --runner ./tests/fpm-mpi-test-runner.sh
```

The MPI test runner launches one process by default; set `FPM_TEST_PROCESSES`
to override that count. It prepares the generated Rankine-vortex forcing
fixture and runs tests from the `tests` directory, preserving the paths used
by the Make suite.

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

Papers about or using UMWM:

* Detelich, C. E., U. G. Schneck, A. G. Hayes, M. Curcic, R. V. Palermo, A. D. Ashton, J. T. Perron, J. M. Lora, and J. Steckloff, 2026: Modeling the seasonality of wind-driven hydrocarbon waves in Titan's polar lakes, *J. Geophys. Res. Planets*, **131**(5), e2026JE009693, doi:10.1029/2026JE009693. [Link](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2026JE009693)

* Schneck, U. G., C. E. Detelich, M. Curcic, A. D. Ashton, A. G. Hayes, and J. T. Perron, 2026: Modeling wind-driven waves on other planets: Applications to Mars, Titan, and exoplanets, *J. Geophys. Res. Planets*, **131**(4), e2025JE009490, doi:10.1029/2025JE009490. [Link](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2025JE009490)

* Barr, B. W. and S. S. Chen, 2025: Impacts of seastate-dependent sea spray heat fluxes on tropical cyclone structure and intensity in fully coupled atmosphere-wave-ocean model simulations, *J. Adv. Model. Earth Syst.*, **17**(7), e2024MS004550, doi:10.1029/2024MS004550. [Link](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2024MS004550)

* Barr, B. W., S. S. Chen, and C. W. Fairall, 2023: Sea-state-dependent sea spray and air-sea heat fluxes in tropical cyclones: A new parameterization for fully coupled atmosphere-wave-ocean models, *J. Atmos. Sci.*, **80**(4), 933-960, doi:10.1175/JAS-D-22-0126.1. [Link](https://journals.ametsoc.org/view/journals/atsc/80/4/JAS-D-22-0126.1.xml)

* Haza A. C., N. Paldor, T. M. Özgökmen, M. Curcic, S. S. Chen, and G. Jacobs, 2019: Wind-based estimations of ocean surface currents from massive clusters of drifters in the Gulf of Mexico, *J. Geophys. Res. Oceans*, **124**, doi:10.1029/2018JC014813. [PDF](https://github.com/milancurcic/publications/blob/master/Haza_etal_JGR2019.pdf).

* Li, G., M. Curcic, M. Iskandarani, S. S. Chen, and O. M. Knio, 2019: Uncertainty propagation in coupled atmosphere-wave-ocean system: A study of Hurricane Earl (2010), *Mon. Wea. Rev.*, **147**, 221-245, doi:10.1175/MWR-D-17-0371.1. [PDF](https://github.com/milancurcic/publications/blob/master/Li_etal_MWR2019.pdf)

* Haza, A., E. D'Asaro, H. Cheng, S. S. Chen, M. Curcic, C. Guigand, H. S. Huntley, G. Jacobs, G. Novelli, T. M. Özgökmen, A. C. Poje, E. Ryan, and A. Scherbina, 2018: Drogue-loss detection for surface drifters during the Lagrangian Submesoscale Experiment (LASER), *J. Atmos. Oceanic Technol.*, **35(4)**, 705-725, doi:10.1175/JTECH-D-17-0143.1. [PDF](https://github.com/milancurcic/publications/blob/master/Haza_etal_JTECH2018.pdf)

* Dietrich, J. C., A. Muhammad, M. Curcic, A. Fathi, C. N. Dawson, S. S. Chen, and R. A. Luettich Jr., 2018: Sensitivity of storm surge predictions to atmospheric forcing during Hurricane Isaac, *J. Waterway, Port, Coastal, Ocean Eng.*, **144**(1): 04017035. [PDF](https://github.com/milancurcic/publications/blob/master/Dietrich_etal_WWENG2018.pdf)

* Kim, E., M. Lance, M. Curcic, S. S. Chen, C. Phillips, and P. Veers, 2016: On the use of coupled wind, wave, and current fields in the simulation of loads on bottom-supported offshore wind turbines during hurricanes, *Technical Report NREL/TP--5000-65283*, National Renewable Energy Lab. (NREL), Golden, CO, United States, doi:10.2172/1266702. [Link](http://www.osti.gov/scitech/biblio/1266702)

* Judt, F., S. S. Chen, and M. Curcic, 2016: Atmospheric forcing of the upper ocean transport in the Gulf of Mexico: From seasonal to diurnal scales, *J. Geophys. Res. Oceans*, **121**, 4416-4433, doi:10.1002/2015JC011555. [PDF](https://github.com/milancurcic/publications/blob/master/Judt_etal_JGR2016.pdf)

* Curcic, M., S. S. Chen, and T. M. Ozgökmen, 2016: Hurricane-induced ocean waves and Stokes drift and their impacts on surface transport and dispersion in the Gulf of Mexico, *Geophys. Res. Lett.*, **43**, 2773–2781, doi:10.1002/2015GL067619. [PDF](https://github.com/milancurcic/publications/blob/master/Curcic_etal_GRL2016.pdf)

* Zhu, P., Y. Wang, S. S. Chen, M. Curcic, and C. Gao, 2016: Impact of storm-induced cooling of sea surface temperature on large turbulent eddies and vertical turbulent transport in the atmospheric boundary layer of Hurricane Isaac, *J. Geophys. Res. Oceans*, **121**, 861–876, doi:10.1002/2015JC011320. [PDF](https://github.com/milancurcic/publications/blob/master/Zhu_etal_JGR2016.pdf)

* Chen, S. S. and M. Curcic, 2016: Ocean surface waves in Hurricane Ike (2008) and Superstorm Sandy (2012): Coupled modeling and observations, *Oce. Mod.*, **103**, 161-176. doi:10.1016/j.ocemod.2015.08.005. [PDF](https://github.com/milancurcic/publications/blob/master/Chen_and_Curcic_OM2016.pdf)

* Curcic, M., 2015: Explicit air-sea momentum exchange in coupled atmosphere-wave-ocean modeling of tropical cyclones, *Ph.D. Thesis*, University of Miami. [Link](http://scholarlyrepository.miami.edu/oa_dissertations/1512)

* Banfield, D., M. A. Donelan, and L. Cavaleri, 2015: Winds, waves and shorelines from ancient martian seas, *Icarus*, **250**, 368-383, doi:10.1016/j.icarus.2014.12.001. [Link](http://www.sciencedirect.com/science/article/pii/S0019103514006794)

* Reichl, B. G., T. Hara, and I. Ginis, 2014: Sea state dependence of the wind stress over the ocean under hurricane winds, *J. Geophys. Res. Oceans*, **119**, 30-51, doi:10.1002/2013JC009289. [Link](http://onlinelibrary.wiley.com/doi/10.1002/2013JC009289/full)

* Curcic M., E. Kim, L. Manuel, S. S. Chen, M. A. Donelan, J. Michalakes, 2013: Coupled atmosphere-wave-ocean modeling to characterize hurricane load cases for offshore wind turbines, *51st AIAA Aerospace Sciences Meeting*, January 2013, Grapevine TX, doi:10.2514/6.2013-198. [PDF](https://github.com/milancurcic/publications/blob/master/Curcic_etal_AIAA2013.pdf)

* Donelan, M. A., M. Curcic, S. S. Chen, and A. K. Magnusson, 2012: Modeling waves and wind stress, *J. Geophys. Res. Oceans*, **117**, C00J23, doi:10.1029/2011JC007787. [PDF](https://github.com/milancurcic/publications/blob/master/Donelan_etal_JGR2012.pdf)

## Thanks

UMWM development is currently supported by the [NSF Award 2543464](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2543464).

Previously, UMWM was supported by the [NSF Award 1745384](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1745384),
the [Gulf of Mexico Research Initiative](http://gulfresearchinitiative.org/),
and the [National Oceanographic Partnership Program](https://www.nopp.org/).

Many thanks as well to the open source contributors:

* Tim Campbell (Naval Research Laboratory)
* Milan Curcic (University of Miami)
* Anton Darmenov (NASA Goddard Space Flight Center)
* Mark Donelan (University of Miami)
* Michael Hirsch (Scivision Inc.)
* Edoardo Mazza (University of Washington)
* Aishwarya Raman (NASA Goddard Space Flight Center)
* Dalton Kei Sasaki (University of Sao Paolo)
* Martin Schmidt (Leibniz Institute for Baltic Sea Research)
* Ashwanth Srinivasan (Tendral LLC)
