# University of Miami Wave Model Makefile
#
# FC, FCFLAGS, and CPPFLAGS come from the top-level
# Makefile. NETCDF comes from the environment.
#
# If your NetCDF libraries are in non-standard directory layout
# (something other than $NETCDF/include and $NETCDF/lib, edit
# NETCDFINC and NETCDFLIB below.

NETCDFINC ?= -I$(NETCDF)/include
NETCDFLIB ?= -L$(NETCDF)/lib -lnetcdf -lnetcdff

# Rules, object list and dependencies start here.

.SUFFIXES: .F90 .f90 .o

# Rule
.F90.o:
	$(FC) -c $(CPPFLAGS) $(FCFLAGS) $(NETCDFINC) $<

.f90.o:
	$(FC) -c $(CPPFLAGS) $(FCFLAGS) $(NETCDFINC) $<

# Object list
OBJS = datetime_module.o \
       umwm_advection.o \
       umwm_constants.o \
       umwm_diagnostics.o \
       umwm_domain.o \
       umwm_forcing.o \
       umwm_grid.o \
       umwm_init.o \
       umwm_io.o \
       umwm_module.o \
       umwm_mpi.o \
       umwm_physics.o \
       umwm_util.o \
       umwm_restart.o \
       umwm_sheltering.o \
       umwm_source_functions.o \
       umwm_spectrum.o \
       umwm_stokes.o \
       umwm_stress.o \
       umwm_top.o

all: umwm

# Main target
umwm: umwm.F90 $(OBJS)
	$(FC) -c $(CPPFLAGS) $(FCFLAGS) $(NETCDFINC) $<
	$(FC) $(FCFLAGS) -o $@ umwm.o $(OBJS) $(NETCDFLIB)
	cp $@ ../.
	ar r libumwm.a $(OBJS)

%.o: %.mod

# Object           Source              Dependencies
datetime_module.o: datetime_module.f90
umwm_domain.o: umwm_domain.F90 datetime_module.o umwm_constants.o umwm_grid.o umwm_spectrum.o
umwm_top.o:        umwm_top.F90 datetime_module.o umwm_constants.o umwm_module.o umwm_mpi.o umwm_forcing.o umwm_physics.o umwm_advection.o umwm_io.o umwm_init.o umwm_restart.o umwm_source_functions.o umwm_stress.o
umwm_module.o: umwm_module.F90    
umwm_constants.o: umwm_constants.F90
umwm_diagnostics.o: umwm_diagnostics.F90 umwm_constants.o
umwm_grid.o: umwm_grid.F90 umwm_constants.o
umwm_mpi.o:        umwm_mpi.F90
umwm_init.o:       umwm_init.F90       umwm_module.o umwm_mpi.o umwm_util.o umwm_io.o umwm_constants.o
umwm_physics.o:    umwm_physics.F90    umwm_module.o umwm_io.o  umwm_advection.o
umwm_advection.o:  umwm_advection.F90  umwm_module.o umwm_io.o  umwm_constants.o
umwm_io.o:         umwm_io.F90         umwm_module.o umwm_mpi.o umwm_util.o umwm_stokes.o
umwm_forcing.o:    umwm_forcing.F90    umwm_module.o umwm_io.o  umwm_util.o
umwm_util.o:       umwm_util.F90       umwm_module.o
umwm_restart.o:    umwm_restart.F90    umwm_module.o umwm_io.o
umwm_sheltering.o: umwm_sheltering.F90 umwm_constants.o
umwm_spectrum.o: umwm_spectrum.F90 umwm_constants.o
umwm_stokes.o:     umwm_stokes.F90     umwm_module.o umwm_constants.o
umwm_stress.o:     umwm_stress.F90     umwm_advection.o umwm_module.o umwm_constants.o umwm_stokes.o
umwm_source_functions.o: umwm_source_functions.F90 umwm_constants.o umwm_io.o umwm_sheltering.o

datetime_module.f90:
	curl -O https://raw.githubusercontent.com/wavebitscientific/datetime-fortran/master/src/datetime_module.f90

.PHONY:
clean:
	$(RM) umwm *.o *.mod libumwm.a 

clean_all:
	$(MAKE) clean
