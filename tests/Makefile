# UMWM test Makefile

include ../mk/config.mk

.SUFFIXES:

SRC_DIR := ../src
TESTS := test_dispersion test_spectrum test_config test_grid test_forcing

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

umwm_module.o: datetime_module.o umwm_constants.o
umwm_env.o: umwm_module.o
umwm_config.o: umwm_constants.o
umwm_dispersion.o: umwm_constants.o
umwm_grid.o: umwm_config.o umwm_constants.o umwm_module.o umwm_spectrum.o umwm_util.o
umwm_forcing.o: umwm_config.o umwm_grid.o umwm_module.o
umwm_spectrum.o: umwm_constants.o
umwm_util.o: umwm_module.o umwm_spectrum.o

clean:
	$(RM) $(TESTS) *.o *.mod

clean_all: clean
