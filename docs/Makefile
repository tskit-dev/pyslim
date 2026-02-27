
all: prereqs
	if (cd _static && md5sum -c pedigrees.ink.svg.md5 | grep -l "FAILED"); then echo "remaking figures with inkscape"; $(MAKE) -C _static; else echo "not remaking figures"; fi
	./build.sh

clean:
	rm -fR _build
	rm -f *.trees

example_sim.trees : example_sim.slim
	slim -s 23 $<

prereqs: example_sim.trees
