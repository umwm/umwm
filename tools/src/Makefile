# University of Miami Wave Model tools Makefile
#
# FC, FCFLAGS, and CPPFLAGS come from the top-level
# Makefile. NETCDF comes from the environment.
#
# If your NetCDF libraries are in non-standard directory layout
# (something other than $NETCDF/include and $NETCDF/lib, edit
# NETCDFINC and NETCDFLIB below.

NETCDFINC = -I$(NETCDF)/include
NETCDFLIB = -L$(NETCDF)/lib -lnetcdf -lnetcdff

# Rules:

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) $(NETCDFINC) $<

all: umwm_gridgen umwm_topogen wrf2umwmgrid wrf2umwmin 

umwm_gridgen: umwm_gridgen.f90 umwm_gridgen.o umwm_utils.o
	$(FC) $(FCFLAGS) -o $@ umwm_gridgen.o umwm_utils.o $(NETCDFLIB)
	cp $@ ../.

umwm_topogen: umwm_topogen.f90 umwm_topogen.o umwm_utils.o
	$(FC) $(FCFLAGS) -o $@ umwm_topogen.o umwm_utils.o $(NETCDFLIB)
	cp $@ ../.

wrf2umwmgrid: wrf2umwmgrid.f90 wrf2umwmgrid.o umwm_utils.o grid_utils.o
	$(FC) $(FCFLAGS) -o $@ wrf2umwmgrid.o umwm_utils.o grid_utils.o $(NETCDFLIB)
	cp $@ ../.

wrf2umwmin: wrf2umwmin.f90 wrf2umwmin.o umwm_utils.o
	$(FC) $(FCFLAGS) -o $@ wrf2umwmin.o umwm_utils.o $(NETCDFLIB)
	cp $@ ../.

# Dependencies
umwm_gridgen.o: umwm_gridgen.f90  umwm_utils.o
umwm_topogen.o: umwm_topogen.f90  umwm_utils.o
wrf2umwmgrid.o: wrf2umwmgrid.f90  umwm_utils.o grid_utils.o
wrf2umwmin.o:   wrf2umwmin.f90 umwm_utils.o
umwm_utils.o:   umwm_utils.f90
grid_utils.o:   grid_utils.f90

.PHONY:
clean:
	rm -f wrf2umwmin wrf2umwmgrid umwm_gridgen umwm_topogen *.o *.mod 
