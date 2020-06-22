# University of Miami Wave Model

# Top-level Makefile

FC ?= mpif90
FCFLAGS ?= -O3 -Wall -std=f2018
CPPFLAGS ?= -DMPI

# don't use make's implicit rule for FC
ifeq ($(FC), f77)
  FC = mpif90
endif

.PHONY: all umwm docs tools clean clean_all

all: umwm tools docs

umwm:
	FC=$(FC) CPPFLAGS='$(CPPFLAGS)' FCFLAGS='$(FCFLAGS)' $(MAKE) --directory=src

docs:
	cd docs && pdflatex umwm_manual_v2.tex
	cd docs && pdflatex umwm_manual_v2.tex
	cd docs && pdflatex umwm_manual_v2.tex
	$(RM) docs/*.{aux,log,toc}

tools:
	FC=$(FC) FCFLAGS='$(FCFLAGS)' $(MAKE) --directory=tools/src

clean:
	$(RM) umwm
	$(RM) tools/umwm_gridgen
	$(RM) tools/umwm_topogen
	$(RM) tools/wrf2umwmgrid
	$(RM) tools/wrf2umwmin
	$(MAKE) --directory=src clean
	$(MAKE) --directory=tools/src clean

clean_all:
	$(MAKE) clean
	$(MAKE) --directory=src clean_all
