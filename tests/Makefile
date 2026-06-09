# UMWM test Makefile

include ../mk/config.mk

.SUFFIXES:

SRC_DIR := ../src
TESTS := test_dispersion test_spectrum

.PHONY: test clean clean_all
.DEFAULT_GOAL := test

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -c $<

%.o: $(SRC_DIR)/%.F90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -c $<

%.o: $(SRC_DIR)/%.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -c $<

test: $(TESTS)
	./test_dispersion
	./test_spectrum

test_dispersion: test_dispersion.f90 tuff.o umwm_constants.o umwm_dispersion.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(LDLIBS)

test_spectrum: test_spectrum.f90 tuff.o umwm_constants.o umwm_spectrum.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -I. -I$(SRC_DIR) -o $@ $^ $(LDLIBS)

umwm_dispersion.o: umwm_constants.o
umwm_spectrum.o: umwm_constants.o

clean:
	$(RM) $(TESTS) *.o *.mod

clean_all: clean
