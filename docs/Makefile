
# Need to set PYTHONPATH so that we pick up the local pyslim
PYPATH=${PWD}/..
PYSLIM_VERSION:=$(shell PYTHONPATH=${PYPATH} \
   python3 -c 'import pyslim; print(pyslim.__version__.split("+")[0])')

dev: prereqs
	if (cd _static && md5sum -c pedigrees.ink.svg.md5 | grep -l "FAILED"); then echo "remaking figures with inkscape"; $(MAKE) -C _static; else echo "not remaking figures"; fi
	PYTHONPATH=${PYPATH} ./build.sh

dist: prereqs
	@echo Building distribution for pyslim version ${PYSLIM_VERSION}
	sed -i s/__PYSLIM_VERSION__/${PYSLIM_VERSION}/g _config.yml
	$(MAKE) -C _static
	PYTHONPATH=${PYPATH} ./build.sh

clean:
	rm -fR _build
	rm -f *.trees

example_sim.trees : example_sim.slim
	slim -s 23 $<

prereqs: example_sim.trees
