# Shared UMWM build configuration.
ifndef UMWM_CONFIG_INCLUDED
UMWM_CONFIG_INCLUDED := 1

MPI ?= no
UMWM_MPI_ENABLED := $(filter yes YES true TRUE 1 on ON,$(strip $(MPI)))

ifeq ($(UMWM_MPI_ENABLED),)
  UMWM_DEFAULT_FC := gfortran
else
  UMWM_DEFAULT_FC := mpif90
  override CPPFLAGS += -DMPI
endif

ifeq ($(origin FC),default)
  FC := $(UMWM_DEFAULT_FC)
else ifeq ($(origin FC),undefined)
  FC := $(UMWM_DEFAULT_FC)
endif

ifeq ($(origin FCFLAGS),default)
  FCFLAGS := -Ofast -march=native -ffast-math -funroll-loops -Wall
else ifeq ($(origin FCFLAGS),undefined)
  FCFLAGS := -Ofast -march=native -ffast-math -funroll-loops -Wall
endif

LDFLAGS ?=
LDLIBS ?=
NF_CONFIG ?= nf-config

UMWM_NO_NETCDF_GOALS := clean clean_all docs help
ifeq ($(strip $(MAKECMDGOALS)),)
  UMWM_NEED_NETCDF := yes
else ifneq ($(filter-out $(UMWM_NO_NETCDF_GOALS),$(MAKECMDGOALS)),)
  UMWM_NEED_NETCDF := yes
else
  UMWM_NEED_NETCDF :=
endif

NF_CONFIG_PATH := $(shell command -v $(NF_CONFIG) 2>/dev/null)

ifeq ($(origin NETCDF_FFLAGS),undefined)
  ifneq ($(strip $(NF_CONFIG_PATH)),)
    NETCDF_FFLAGS := $(shell $(NF_CONFIG) --fflags 2>/dev/null)
    NETCDF_FFLAGS_SOURCE := nf-config
  else ifneq ($(strip $(NETCDF)),)
    NETCDF_FFLAGS := -I$(NETCDF)/include
    NETCDF_FFLAGS_SOURCE := NETCDF
  endif
else
  NETCDF_FFLAGS_SOURCE := manual
endif

ifeq ($(origin NETCDF_FLIBS),undefined)
  ifneq ($(strip $(NF_CONFIG_PATH)),)
    NETCDF_FLIBS := $(shell $(NF_CONFIG) --flibs 2>/dev/null)
    NETCDF_FLIBS_SOURCE := nf-config
  else ifneq ($(strip $(NETCDF)),)
    NETCDF_FLIBS := -L$(NETCDF)/lib -lnetcdff -lnetcdf
    NETCDF_FLIBS_SOURCE := NETCDF
  endif
else
  NETCDF_FLIBS_SOURCE := manual
endif

ifneq ($(UMWM_NEED_NETCDF),)
  ifeq ($(strip $(NETCDF_FFLAGS_SOURCE)),)
    $(error NetCDF Fortran compile flags not found. Set NF_CONFIG=/path/to/nf-config, NETCDF=/prefix, or NETCDF_FFLAGS/NETCDF_FLIBS manually)
  endif
  ifeq ($(strip $(NETCDF_FLIBS)),)
    $(error NetCDF Fortran link flags not found. Set NF_CONFIG=/path/to/nf-config, NETCDF=/prefix, or NETCDF_FLIBS manually)
  endif
endif

.PHONY: print-config
print-config:
	@printf 'MPI=%s\n' '$(MPI)'
	@printf 'FC=%s\n' '$(FC)'
	@printf 'FCFLAGS=%s\n' '$(FCFLAGS)'
	@printf 'CPPFLAGS=%s\n' '$(CPPFLAGS)'
	@printf 'LDFLAGS=%s\n' '$(LDFLAGS)'
	@printf 'LDLIBS=%s\n' '$(LDLIBS)'
	@printf 'NETCDF_FFLAGS=%s\n' '$(NETCDF_FFLAGS)'
	@printf 'NETCDF_FLIBS=%s\n' '$(NETCDF_FLIBS)'

endif
