# University of Miami Wave Model Makefile

include ../mk/config.mk

.SUFFIXES:

%.o: %.F90
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(NETCDF_FFLAGS) -c $<

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(NETCDF_FFLAGS) -c $<

# Object list
OBJS = datetime_module.o \
       umwm_advection.o \
       umwm_constants.o \
       umwm_diagnostics.o \
       umwm_env.o \
       umwm_forcing.o \
       umwm_init.o \
       umwm_io.o \
       umwm_module.o \
       umwm_mpi.o \
       umwm_physics.o \
       umwm_util.o \
       umwm_restart.o \
       umwm_sheltering.o \
       umwm_source_functions.o \
       umwm_stokes.o \
       umwm_stress.o \
       umwm_top.o

.PHONY: all clean clean_all print-config
.DEFAULT_GOAL := all

all: umwm

# Main target
umwm: umwm.o $(OBJS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ umwm.o $(OBJS) $(NETCDF_FLIBS) $(LDLIBS)
	cp $@ ../.
	$(AR) rcs libumwm.a $(OBJS)

# Object            Source              Dependencies
umwm.o:            umwm.F90         umwm_module.o umwm_top.o
datetime_module.o:  datetime_module.f90
umwm_top.o:         umwm_top.F90 datetime_module.o umwm_constants.o umwm_env.o umwm_module.o umwm_mpi.o umwm_forcing.o umwm_physics.o umwm_advection.o umwm_io.o umwm_init.o umwm_restart.o umwm_source_functions.o umwm_stress.o
umwm_module.o:      umwm_module.F90     datetime_module.o umwm_constants.o
umwm_constants.o:   umwm_constants.F90
umwm_diagnostics.o: umwm_diagnostics.F90 umwm_constants.o
umwm_mpi.o:         umwm_mpi.F90
umwm_init.o:        umwm_init.F90       umwm_module.o umwm_mpi.o umwm_util.o umwm_io.o umwm_constants.o
umwm_physics.o:     umwm_physics.F90    umwm_module.o umwm_io.o  umwm_advection.o
umwm_advection.o:   umwm_advection.F90  umwm_module.o umwm_io.o  umwm_constants.o
umwm_io.o:          umwm_io.F90         umwm_module.o umwm_mpi.o umwm_util.o umwm_stokes.o
umwm_forcing.o:     umwm_forcing.F90    umwm_module.o umwm_io.o  umwm_util.o
umwm_util.o:        umwm_util.F90       umwm_module.o
umwm_restart.o:     umwm_restart.F90    umwm_module.o umwm_io.o
umwm_sheltering.o:  umwm_sheltering.F90 umwm_constants.o
umwm_stokes.o:      umwm_stokes.F90     umwm_module.o umwm_constants.o
umwm_stress.o:      umwm_stress.F90     umwm_advection.o umwm_module.o umwm_constants.o umwm_stokes.o
umwm_source_functions.o: umwm_source_functions.F90 umwm_constants.o umwm_io.o umwm_sheltering.o

datetime_module.f90:
	curl -O https://raw.githubusercontent.com/wavebitscientific/datetime-fortran/master/src/datetime_module.f90

clean:
	$(RM) umwm *.o *.mod libumwm.a 

clean_all:
	$(MAKE) clean
