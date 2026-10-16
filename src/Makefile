# University of Miami Wave Model Makefile

include ../mk/config.mk

.SUFFIXES:

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(NETCDF_FFLAGS) -c $<

# Object list
OBJS = datetime_module.o \
       umwm_advection.o \
       umwm_config.o \
       umwm_constants.o \
       umwm_diagnostics.o \
       umwm_dispersion.o \
       umwm_domain.o \
       umwm_domain_output.o \
       umwm_env.o \
       umwm_forcing.o \
       umwm_grid.o \
       umwm_io.o \
       umwm_mpi.o \
       umwm_physics.o \
       umwm_util.o \
       umwm_restart.o \
       umwm_spectrum.o \
       umwm_sheltering.o \
       umwm_source_functions.o \
       umwm_stokes.o \
       umwm_stress.o

.PHONY: all test clean clean_all print-config
.DEFAULT_GOAL := all

all: umwm

# Main target
umwm: umwm.o $(OBJS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ umwm.o $(OBJS) $(NETCDF_FLIBS) $(LDLIBS)
	cp $@ ../.
	$(AR) rcs libumwm.a $(OBJS)

# Object            Source              Dependencies
umwm.o:             umwm.f90            umwm_domain.o umwm_env.o
datetime_module.o:  datetime_module.f90
umwm_domain.o:      umwm_domain.f90     datetime_module.o umwm_advection.o umwm_config.o umwm_constants.o umwm_dispersion.o umwm_env.o umwm_forcing.o umwm_grid.o umwm_io.o umwm_mpi.o umwm_physics.o umwm_restart.o umwm_sheltering.o umwm_source_functions.o umwm_spectrum.o umwm_stokes.o umwm_stress.o umwm_util.o
umwm_domain_output.o: umwm_domain_output.f90 umwm_domain.o umwm_io.o
umwm_constants.o:   umwm_constants.f90
umwm_config.o:      umwm_config.f90     umwm_constants.o
umwm_env.o:         umwm_env.f90
umwm_util.o:        umwm_util.f90
umwm_grid.o:        umwm_grid.f90       umwm_config.o umwm_constants.o umwm_env.o umwm_spectrum.o umwm_util.o
umwm_diagnostics.o: umwm_diagnostics.f90 umwm_constants.o
umwm_dispersion.o:  umwm_dispersion.f90  umwm_constants.o
umwm_spectrum.o:    umwm_spectrum.f90    umwm_constants.o
umwm_mpi.o:         umwm_mpi.f90        umwm_env.o umwm_grid.o
umwm_physics.o:     umwm_physics.f90    umwm_source_functions.o
umwm_advection.o:   umwm_advection.f90
umwm_io.o:          umwm_io.f90         datetime_module.o umwm_env.o umwm_grid.o umwm_mpi.o
umwm_forcing.o:     umwm_forcing.f90    umwm_config.o umwm_grid.o
umwm_restart.o:     umwm_restart.f90    umwm_env.o umwm_grid.o umwm_io.o umwm_spectrum.o umwm_util.o
umwm_sheltering.o:  umwm_sheltering.f90 umwm_constants.o
umwm_stokes.o:      umwm_stokes.f90     umwm_constants.o
umwm_stress.o:      umwm_stress.f90     umwm_constants.o
umwm_source_functions.o: umwm_source_functions.f90

test:
	$(MAKE) --directory=../tests test

datetime_module.f90:
	curl -O https://raw.githubusercontent.com/wavebitscientific/datetime-fortran/master/src/datetime_module.f90

clean:
	$(RM) umwm *.o *.mod *.smod libumwm.a

clean_all:
	$(MAKE) clean
