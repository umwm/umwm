# UMWM test Makefile

include ../mk/config.mk

.SUFFIXES:

SRC_DIR := ../src
TESTS := test_dispersion test_spectrum test_config test_grid test_forcing \
         test_source_functions test_modulation test_advection test_domain
MODEL_OBJS := datetime_module.o umwm_advection.o umwm_config.o umwm_constants.o \
              umwm_dispersion.o umwm_domain.o umwm_domain_acc.o umwm_domain_output.o umwm_env.o \
              umwm_forcing.o umwm_grid.o umwm_io.o umwm_modulation.o umwm_mpi.o umwm_physics.o \
              umwm_restart.o umwm_sheltering.o umwm_source_functions.o \
              umwm_spectrum.o umwm_stokes.o umwm_stress.o umwm_util.o
PYTHON ?= python
RANKINE_FORCING_STAMP := input/.rankine_forcing.stamp

.PHONY: test clean clean_all
.DEFAULT_GOAL := test

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(NETCDF_FFLAGS) -I. -I$(SRC_DIR) -c $<

%.o: $(SRC_DIR)/%.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(NETCDF_FFLAGS) -I. -I$(SRC_DIR) -c $<

$(SRC_DIR)/datetime_module.f90:
	$(MAKE) --directory=$(SRC_DIR) datetime_module.f90

test: $(TESTS)
	mkdir -p output
	./test_dispersion
	./test_spectrum
	./test_config
	./test_grid
	./test_forcing
	./test_source_functions
	./test_modulation
	./test_advection
	./test_domain

test_dispersion: test_dispersion.f90 tuff.o umwm_constants.o umwm_dispersion.o umwm_env.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(LDLIBS)

test_spectrum: test_spectrum.f90 tuff.o umwm_constants.o umwm_spectrum.o umwm_env.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(LDLIBS)

test_config: test_config.f90 tuff.o umwm_constants.o umwm_source_functions.o umwm_config.o umwm_env.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(LDLIBS)

test_grid: test_grid.f90 tuff.o umwm_constants.o umwm_source_functions.o umwm_config.o umwm_spectrum.o umwm_env.o umwm_util.o umwm_grid.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(NETCDF_FLIBS) $(LDLIBS)

test_forcing: test_forcing.f90 tuff.o umwm_constants.o umwm_source_functions.o umwm_config.o umwm_spectrum.o umwm_env.o umwm_util.o umwm_grid.o umwm_forcing.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(NETCDF_FLIBS) $(LDLIBS)

test_source_functions: test_source_functions.f90 tuff.o umwm_constants.o umwm_source_functions.o umwm_env.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(LDLIBS)

test_modulation: test_modulation.f90 tuff.o umwm_constants.o umwm_modulation.o umwm_source_functions.o umwm_env.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(LDLIBS)

test_advection: test_advection.f90 tuff.o umwm_constants.o umwm_dispersion.o umwm_advection.o umwm_env.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(LDLIBS)

test_domain: test_domain.f90 tuff.o $(MODEL_OBJS)
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(NETCDF_FLIBS) $(LDLIBS)

test_forcing test_domain: | $(RANKINE_FORCING_STAMP)

$(RANKINE_FORCING_STAMP): forcing/idealized_rankine_forcing.py Makefile
	mkdir -p input
	$(PYTHON) forcing/idealized_rankine_forcing.py \
	  --output-dir input \
	  --start-time 2026-01-01T00:00:00 \
	  --last-hour 1 \
	  --nx 21 \
	  --ny 11 \
	  --dx-m 10000 \
	  --dy-m 10000 \
	  --vmax 30 \
	  --rmw-m 50000 \
	  --initial-center-x-m 0 \
	  --initial-center-y-m 0 \
	  --translation-x-m-per-hour 10000 \
	  --translation-y-m-per-hour 0 \
	  --uc 0.25 \
	  --vc -0.5 \
	  --rhoa 1.3 \
	  --rhow 1025
	touch $@

umwm_config.o: umwm_constants.o umwm_source_functions.o
umwm_dispersion.o: umwm_constants.o
umwm_spectrum.o: umwm_constants.o
umwm_grid.o: umwm_config.o umwm_constants.o umwm_env.o umwm_spectrum.o umwm_util.o
umwm_forcing.o: umwm_config.o umwm_grid.o
umwm_mpi.o: umwm_env.o umwm_grid.o
umwm_io.o: datetime_module.o umwm_env.o umwm_grid.o umwm_mpi.o
umwm_restart.o: umwm_env.o umwm_grid.o umwm_io.o umwm_spectrum.o umwm_util.o
umwm_sheltering.o: umwm_constants.o
umwm_source_functions.o: umwm_constants.o
umwm_modulation.o:
umwm_physics.o: umwm_source_functions.o
umwm_stokes.o: umwm_constants.o
umwm_stress.o: umwm_constants.o
umwm_domain.o: datetime_module.o umwm_advection.o umwm_config.o umwm_constants.o \
               umwm_dispersion.o umwm_env.o umwm_forcing.o umwm_grid.o umwm_io.o \
               umwm_modulation.o umwm_mpi.o umwm_physics.o umwm_restart.o umwm_sheltering.o \
               umwm_source_functions.o umwm_spectrum.o umwm_stokes.o umwm_stress.o \
               umwm_util.o
umwm_domain_output.o: umwm_domain.o umwm_io.o
umwm_domain_acc.o: umwm_advection.o umwm_domain.o umwm_physics.o umwm_source_functions.o

# the inlining of the GPU build, as in ../src/Makefile
INLINE_SRC = umwm_advection.f90 umwm_physics.f90 umwm_sheltering.f90 \
             umwm_source_functions.f90 umwm_stress.f90
ACC_KERNEL_FLAGS ?= -Minline=lib:umwm_inline.lib -gpu=maxregcount:64
ifneq ($(UMWM_GPU_ENABLED),)
umwm_inline.lib: $(addprefix $(SRC_DIR)/,$(INLINE_SRC)) | $(INLINE_SRC:.f90=.o)
	$(RM) -r $@
	$(FC) $(CPPFLAGS) $(FCFLAGS) -Mextract=lib:$@ $^
umwm_domain_acc.o: $(SRC_DIR)/umwm_domain_acc.f90 umwm_inline.lib
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(ACC_KERNEL_FLAGS) $(NETCDF_FFLAGS) -I. -I$(SRC_DIR) -c $<
endif

clean:
	$(RM) -r umwm_inline.lib
	$(RM) $(TESTS) *.o *.mod *.smod $(RANKINE_FORCING_STAMP) \
	  input/umwmin_2026-01-01_00:00:00.nc \
	  input/umwmin_2026-01-01_01:00:00.nc output/umwmout.grid

clean_all: clean
