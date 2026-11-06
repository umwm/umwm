# University of Miami Wave Model (UMWM)

A tiny, ultrafast, parallel spectral ocean wave model for CPUs and GPUs.

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
`-Ofast -march=native -ffast-math -funroll-loops -Wall`, or `-O3 -fast` when
`FC` is `nvfortran` or an MPI wrapper of it. Override `FCFLAGS` when using
another compiler that needs different optimization options:

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
`make shared` builds the shared library `src/shared/libumwm.so` with the C
API used by the Python package.
Executable `umwm` will be built in the top-level directory and auxiliary tool
executables will be built in `tools/`.
PDF documentation can be built separately with `make docs`.

### Building for NVIDIA GPUs

Build for NVIDIA GPUs by setting `GPU=yes`. This needs the
[NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk), for `nvfortran` with
OpenACC and its CUDA-aware HPC-X Open MPI, and a netCDF-Fortran built with
`nvfortran`:

```
module load nvhpc-hpcx-cuda13
unset FC
make GPU=yes MPI=yes NF_CONFIG=/path/to/nvfortran/netcdf-fortran/bin/nf-config
```

Load the HPC-X module of the SDK rather than `nvhpc`, so that MPI uses the
UCX of the SDK, which transfers GPU memory, and unset the `FC` that the module
exports, which would take precedence over `mpif90`. `GPU=yes` selects
`nvfortran` (or `mpif90` with `MPI=yes`) with the default `FCFLAGS`
`-O3 -fast`, and adds the OpenACC flags `-acc=gpu -gpu=$(GPU_ARCH),mem:unified`
also when `FCFLAGS` is set. `GPU_ARCH` is the compute capability of the GPU,
`cc120` (Blackwell RTX PRO) by default. The kernels of the time step, in
`src/umwm_domain_acc.f90`, are compiled with `ACC_KERNEL_FLAGS`, which inline
the functions of the spectral bins and limit the registers of a thread.

In the GPU build, the source time step, the forcing interpolation, the
diagnostics, the ocean stress, the Stokes drift, and most of the
initialization run on the GPU, and the input and output on the CPU. The
fields are in unified memory, which the CPU and the GPU share, as the code has
no OpenACC data directives; this needs a system where the GPU can address the
memory of the CPU, as with heterogeneous memory management (HMM), for which
`nvidia-smi -q` shows the addressing mode `HMM`. It supports at most 128
frequencies and 128 directions. Its results differ from those of the CPU
build at the level of single-precision rounding, because the GPU raises to
the dissipation power in single precision, where the CPU uses double
precision, and sums over the spectrum in another order. Run it on one GPU as
in serial, selecting the GPU with `CUDA_VISIBLE_DEVICES`:

```
CUDA_VISIBLE_DEVICES=0 ./umwm CONFIG -o OUTPUT
```

The GPU build is fast for large grids, of the order of 100,000 sea cells and
more; small grids, such as most example cases, use a small part of the GPU.

With MPI, the GPU build runs on several GPUs, with the same results as on
one. The processes of a node share its visible GPUs, one process per GPU, or
blocks of consecutive ranks if there are more processes than GPUs, and the
halo is exchanged while the cells inside the tiles are integrated. Map the
processes to the NUMA domains in rank order, so that each runs next to its
GPU on nodes where the GPUs are numbered by NUMA domain; for example, on a
node with four GPUs, two on each of two NUMA domains:

```
mpirun -np 4 --map-by ppr:2:numa --bind-to numa ./umwm CONFIG -o OUTPUT
```

The MPI library should be CUDA-aware, as the HPC-X of the SDK, which copies
the halo from and to the GPU; through other libraries, the halo pages
migrate to the CPU and back at each time step. On two GPUs, a time step of a
grid of 360,000 sea cells or more takes half the time of one GPU.

### Running tests

Run the test suite from the top-level directory:

```
make test
```

The tests generate their forcing with Python (`requirements.txt`); set
`PYTHON` to use a particular interpreter, for example `make test
PYTHON="uv run python"` with the environment of the Python package.

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

The model runs a configuration: a directory with the `main.nml` namelist,
the `spectrum.nml` points of the spectrum output, and the input files in
`input/`, such as the [example cases](examples/). Running in serial mode:

```
./umwm examples/jonswap -o runs/jonswap
```

Running in parallel, for example on 16 cores:

```
mpiexec -n 16 ./umwm examples/jonswap -o runs/jonswap
```

The model reads the input files from `CONFIG/input`, or the directory given
with `-i INPUT`, and writes the output and restart files, with copies of the
namelists, `umwmout_main.nml` and `umwmout_spectrum.nml`, to the directory
given with `-o OUTPUT`, created if needed, or else the current directory. A
restart reads its restart file from the output directory. The output
directory may be the configuration or the input directory, as the files the
model writes there are all named `umwmout_*`, `umwmspc_*` or `umwmrst_*`.
The configuration may also be given as its main namelist file, whose
directory holds the others. `./umwm --help` lists the options.

You can read the full technical reference doc [here](DOCS.md).

### Example cases

The [examples](examples/README.md) directory holds a gallery of idealized
cases with constant wind: fetch-limited and duration-limited growth, the
fetches of JONSWAP, a laboratory wave tank, lakes on Titan and Mars with
their planetary parameters, and waves over a seamount and through a field
of ocean eddies, whose bathymetry and currents come with the cases. Each
case is a configuration directory:

```
./umwm --list-cases
./umwm examples/jonswap -o runs/jonswap
./umwm examples/titan-ontario-lacus -o runs/titan
```

or with the Python command-line tool, which also lists the cases with their
directories, and copies one into a directory to edit:

```
umwm cases
umwm run examples/mars-jezero-lake -o runs/mars
umwm init --case mars-jezero-lake runs/mars-windy -s forcing_constant.wspd0=8
```

### Running from ERA5

UMWM starts from, and is forced by, the [ERA5](https://doi.org/10.1002/qj.3803)
reanalysis on a regional or global longitude-latitude grid, with no other
program to run: the namelist names the grid and the data source, and the
model prepares its missing input files from ERA5 before it runs, and again
those it prepared for another grid, such as for another case in the same
directory.

```
&GRID
  gridFromFile = .F.
  topoFromFile = .T.               ! the bathymetry of the ERA5 wave model
  lon0 = -98.5, lat0 = 17.5        ! cell (1,1) of a regular lon-lat grid
  dlon = 0.25,  dlat = 0.25        ! its spacing in degrees
  ...
/
&FORCING
  winds = .T., air_density = .T., seaice = .T.
  waves = .T.                      ! wave spectrum from the ERA5 wave partitions
  source = 'era5'                  ! prepare the missing input files from ERA5
  ...
/
```

The input files are prepared from the
[analysis-ready, cloud-optimized ERA5 store](https://cloud.google.com/storage/docs/public-datasets/era5)
of Google, read anonymously over HTTPS, on the grid of the namelist: the
bathymetry, the hourly wind, air density and sea ice, and the wave
partitions, the wind sea and three swells. The model rebuilds its initial
spectrum from the partitions, with the significant wave height, mean
period, mean direction and directional spread of each, and in a regional
domain it propagates the spectrum of the ERA5 partitions in from its open
boundaries throughout the run.

The `umwm` executable prepares the files with the Python package (below)
installed with the dependencies of ERA5, `pip install '.[era5]'`; set
`UMWM_PYTHON` to the Python interpreter that has it if that is not
`python3`. The Python command-line tool prepares them as well, when it runs
the model or beforehand, for example before a parallel run on nodes without
network access:

```
./umwm examples/era5-helene -o runs/helene -i runs/helene/input         # hurricane Helene
umwm run examples/era5-global-1deg -o runs/global -i runs/global/input  # the global ocean at 1 degree
umwm init --case era5-helene gulf                                       # a configuration to edit
umwm forcing era5 gulf                                                  # its input files, in gulf/input
mpiexec -n 8 ./umwm gulf -o runs/gulf                                   # and run it in parallel
```

The input files are prepared in the input directory: the directory given
with `-i`, or else the `input/` directory of the configuration, such as
`gulf/input`; with `-i`, the example cases are left as they are.

The example cases `era5-helene`, `era5-global-1deg` and `era5-global-0.25deg`
are templates for other regions and times. [DOCS.md](DOCS.md#wave-partitions-initial-state-and-open-boundaries)
describes the spectra of the partitions, the open boundaries, and how to add
another data source, such as NCEP GFS with its WAVEWATCH III wave model.

### Python interface

The `umwm` Python package configures and runs the model from Python, with
the model state and the output available as NumPy arrays and xarray
datasets, and provides the `umwm` command-line tool. The time integration
runs in the compiled Fortran library, so there is no performance cost.
Install it from the repository (it builds the library with the same
Makefile configuration, for example `UMWM_MPI=yes pip install .` for MPI or
`UMWM_GPU=yes pip install .` for NVIDIA GPUs),
or set up an environment with [uv](https://docs.astral.sh/uv/):

```
pip install .
uv sync && uv run umwm info
```

```python
import umwm

config = umwm.Config()                        # the model defaults
config.domain.stoptimestr = "2026-01-01 12:00:00"
config.forcing_constant.wspd0 = 15.0

with umwm.Domain(config, output_dir="runs/case1") as domain:
    domain.run()
    swh = domain.grid_field("swh")            # NumPy array (ny, nx)
    output = domain.output()                  # xarray.Dataset of the output files
```

```
umwm init case1 -s domain.stoptimestr="2026-01-01 12:00:00"
umwm run case1 -o runs/case1
umwm output runs/case1 --stats
umwm run examples/jonswap -o runs/jonswap     # an example case
```

See [python/README.md](python/README.md) for the documentation of the
package and the command-line tool.

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
