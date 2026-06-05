# University of Miami Wave Model

# Top-level Makefile

include mk/config.mk

.PHONY: all umwm docs tools clean clean_all print-config
.DEFAULT_GOAL := all

all: umwm tools

umwm:
	$(MAKE) --directory=src

docs:
	cd docs && pdflatex umwm_manual_v2.tex
	cd docs && pdflatex umwm_manual_v2.tex
	cd docs && pdflatex umwm_manual_v2.tex
	$(RM) docs/*.aux docs/*.log docs/*.toc

tools:
	$(MAKE) --directory=tools/src

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
