RSCRIPT = Rscript
RUN_WALD_SIM ?= 1

all: visualization

clean:
	find figures -type f -not -name ".gitkeep" -delete

visualization:
	RUN_WALD_SIM=$(RUN_WALD_SIM) $(RSCRIPT) R/main.R

requirements:
	$(RSCRIPT) -e "renv::restore()"


.PHONY: requirements visualization all clean test
