# University of Miami Wave Model Makefile

# This Makefile may be run from another directory with -f, as the shared
# target does, so the build configuration is located relative to it.
UMWM_SRC_MAKEFILE := $(abspath $(lastword $(MAKEFILE_LIST)))
UMWM_SRC_DIR := $(dir $(UMWM_SRC_MAKEFILE))

include $(UMWM_SRC_DIR)../mk/config.mk

# Directory of the sources when building out of tree, as the shared target
# does; only the sources are searched there, not the objects.
SRC_VPATH ?=
ifneq ($(strip $(SRC_VPATH)),)
vpath %.f90 $(SRC_VPATH)
endif

.SUFFIXES:

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(NETCDF_FFLAGS) -c $<

# Object list
OBJS = datetime_module.o \
       umwm_advection.o \
       umwm_capi.o \
       umwm_config.o \
       umwm_constants.o \
       umwm_diagnostics.o \
       umwm_dispersion.o \
       umwm_domain.o \
       umwm_domain_acc.o \
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

.PHONY: all shared test clean clean_all print-config
.DEFAULT_GOAL := all

all: umwm

# Main target
umwm: umwm.o umwm_cli.o $(OBJS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ umwm.o umwm_cli.o $(OBJS) $(NETCDF_FLIBS) $(LDLIBS)
	cp $@ ../.
	$(AR) rcs libumwm.a $(OBJS)

# Shared library with the C API, used by the Python package. The objects
# are compiled with -fPIC in SHARED_DIR, apart from those of the executable,
# and the library is linked with the run path of the netCDF libraries.
SHARED_DIR ?= shared
SHARED_LIB ?= libumwm.so
comma := ,
SHARED_RPATH := $(patsubst -L%,-Wl$(comma)-rpath$(comma)%,$(filter -L%,$(NETCDF_FLIBS)))

shared:
	mkdir -p $(SHARED_DIR)
	$(MAKE) --directory=$(SHARED_DIR) -f $(UMWM_SRC_MAKEFILE) SRC_VPATH=$(UMWM_SRC_DIR) \
	  FCFLAGS="$(FCFLAGS) -fPIC" $(SHARED_LIB)

$(SHARED_LIB): $(OBJS)
	$(FC) -shared $(FCFLAGS) $(LDFLAGS) -o $@ $^ $(NETCDF_FLIBS) $(SHARED_RPATH) $(LDLIBS)

# Object            Source              Dependencies
umwm.o:             umwm.f90            umwm_cli.o umwm_domain.o umwm_env.o
umwm_cli.o:         umwm_cli.f90        umwm_env.o
datetime_module.o:  datetime_module.f90
umwm_domain.o:      umwm_domain.f90     datetime_module.o umwm_advection.o umwm_config.o umwm_constants.o umwm_dispersion.o umwm_env.o umwm_forcing.o umwm_grid.o umwm_io.o umwm_mpi.o umwm_physics.o umwm_restart.o umwm_sheltering.o umwm_source_functions.o umwm_spectrum.o umwm_stokes.o umwm_stress.o umwm_util.o
umwm_domain_output.o: umwm_domain_output.f90 umwm_domain.o umwm_io.o
umwm_domain_acc.o:  umwm_domain_acc.f90 umwm_advection.o umwm_domain.o umwm_physics.o umwm_source_functions.o
umwm_capi.o:        umwm_capi.f90       umwm_domain.o umwm_env.o umwm_io.o
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

# In the GPU build, the elemental functions of the spectral bins are inlined
# into the kernels of umwm_domain_acc from a library extracted from their
# modules; nvfortran does not inline them across files otherwise. The
# kernels run faster with fewer registers for more threads at a time.
INLINE_SRC = umwm_advection.f90 umwm_physics.f90 umwm_sheltering.f90 \
             umwm_source_functions.f90 umwm_stress.f90
ACC_KERNEL_FLAGS ?= -Minline=lib:umwm_inline.lib -gpu=maxregcount:64
ifneq ($(UMWM_GPU_ENABLED),)
umwm_inline.lib: $(INLINE_SRC) | $(INLINE_SRC:.f90=.o)
	$(RM) -r $@
	$(FC) $(CPPFLAGS) $(FCFLAGS) -Mextract=lib:$@ $^
umwm_domain_acc.o: umwm_domain_acc.f90 umwm_inline.lib
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(ACC_KERNEL_FLAGS) $(NETCDF_FFLAGS) -c $<
endif

test:
	$(MAKE) --directory=../tests test

datetime_module.f90:
	curl -O https://raw.githubusercontent.com/wavebitscientific/datetime-fortran/master/src/datetime_module.f90

clean:
	$(RM) umwm *.o *.mod *.smod libumwm.a $(SHARED_LIB)
	$(RM) -r umwm_inline.lib
	$(RM) -r $(SHARED_DIR)

clean_all:
	$(MAKE) clean
