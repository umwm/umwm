# University of Miami Wave Model tools Makefile

include ../../mk/config.mk

.SUFFIXES:

.PHONY: all clean print-config
.DEFAULT_GOAL := all

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) $(NETCDF_FFLAGS) -c $<

all: umwm_gridgen umwm_topogen wrf2umwmgrid wrf2umwmin 

umwm_gridgen: umwm_gridgen.f90 umwm_gridgen.o umwm_utils.o
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ umwm_gridgen.o umwm_utils.o $(NETCDF_FLIBS) $(LDLIBS)
	cp $@ ../.

umwm_topogen: umwm_topogen.f90 umwm_topogen.o umwm_utils.o
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ umwm_topogen.o umwm_utils.o $(NETCDF_FLIBS) $(LDLIBS)
	cp $@ ../.

wrf2umwmgrid: wrf2umwmgrid.f90 wrf2umwmgrid.o umwm_utils.o grid_utils.o
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ wrf2umwmgrid.o umwm_utils.o grid_utils.o $(NETCDF_FLIBS) $(LDLIBS)
	cp $@ ../.

wrf2umwmin: wrf2umwmin.f90 wrf2umwmin.o umwm_utils.o
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ wrf2umwmin.o umwm_utils.o $(NETCDF_FLIBS) $(LDLIBS)
	cp $@ ../.

# Dependencies
umwm_gridgen.o: umwm_gridgen.f90  umwm_utils.o
umwm_topogen.o: umwm_topogen.f90  umwm_utils.o
wrf2umwmgrid.o: wrf2umwmgrid.f90  umwm_utils.o grid_utils.o
wrf2umwmin.o:   wrf2umwmin.f90 umwm_utils.o
umwm_utils.o:   umwm_utils.f90
grid_utils.o:   grid_utils.f90

clean:
	$(RM) wrf2umwmin wrf2umwmgrid umwm_gridgen umwm_topogen *.o *.mod
