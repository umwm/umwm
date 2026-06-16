# UMWM test Makefile

include ../mk/config.mk

.SUFFIXES:

SRC_DIR := ../src
TESTS := test_dispersion test_spectrum test_config test_grid test_forcing
PYTHON ?= python
RANKINE_FORCING_STAMP := input/.rankine_forcing.stamp

.PHONY: test clean clean_all
.DEFAULT_GOAL := test

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(NETCDF_FFLAGS) -I. -I$(SRC_DIR) -c $<

%.o: $(SRC_DIR)/%.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(NETCDF_FFLAGS) -I. -I$(SRC_DIR) -c $<

test: $(TESTS)
	./test_dispersion
	./test_spectrum
	./test_config
	./test_grid
	./test_forcing

test_dispersion: test_dispersion.f90 tuff.o umwm_constants.o umwm_dispersion.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(LDLIBS)

test_spectrum: test_spectrum.f90 tuff.o umwm_constants.o umwm_spectrum.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(LDLIBS)

test_config: test_config.f90 tuff.o umwm_constants.o umwm_config.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(LDLIBS)

test_grid: test_grid.f90 tuff.o datetime_module.o umwm_constants.o umwm_config.o umwm_spectrum.o umwm_module.o umwm_env.o umwm_util.o umwm_grid.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(NETCDF_FLIBS) $(LDLIBS)

test_forcing: test_forcing.f90 tuff.o datetime_module.o umwm_constants.o umwm_config.o umwm_spectrum.o umwm_module.o umwm_util.o umwm_grid.o umwm_forcing.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(NETCDF_FLIBS) $(LDLIBS)

test_forcing: | $(RANKINE_FORCING_STAMP)

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

umwm_module.o: datetime_module.o umwm_constants.o
umwm_env.o: umwm_module.o
umwm_config.o: umwm_constants.o
umwm_dispersion.o: umwm_constants.o
umwm_grid.o: umwm_config.o umwm_constants.o umwm_module.o umwm_spectrum.o umwm_util.o
umwm_forcing.o: umwm_config.o umwm_grid.o umwm_module.o
umwm_spectrum.o: umwm_constants.o
umwm_util.o: umwm_module.o umwm_spectrum.o

clean:
	$(RM) $(TESTS) *.o *.mod $(RANKINE_FORCING_STAMP) \
	  input/umwmin_2026-01-01_00:00:00.nc \
	  input/umwmin_2026-01-01_01:00:00.nc

clean_all: clean
