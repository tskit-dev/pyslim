
# Need to set PYTHONPATH so that we pick up the local pyslim
PYPATH=${PWD}/..
PYSLIM_VERSION:=$(shell PYTHONPATH=${PYPATH} \
   python3 -c 'import pyslim; print(pyslim.__version__.split("+")[0])')

dev: prereqs
	$(MAKE) -C _static
	PYTHONPATH=${PYPATH} ./build.sh

dist: prereqs
	@echo Building distribution for pyslim version ${PYSLIM_VERSION}
	sed -i s/__PYSLIM_VERSION__/${PYSLIM_VERSION}/g _config.yml
	$(MAKE) -C _static
	PYTHONPATH=${PYPATH} ./build.sh

clean:
	rm -fR _build

prereqs:
	slim -s 23 example_sim.slim
