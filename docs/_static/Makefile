.PHONY : all

# NOTE: all layers must be visible in pedigrees.ink.svg for this to work.

all: pedigree01.png pedigree2.png pedigree3.png pedigree_recapitate.png pedigree_simplify.png pedigree_mutate.png

pedigree01.png : pedigree0.png pedigree1.png
	convert +append $^ $@

pedigree0.pdf : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		ts_indivs simp_indivs pedigree simp_pedigree chromosomes simplify_removes \
		extra_chromosomes extra_indivs label0 \
		> $@

pedigree1.pdf : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		time_label ts_indivs simp_indivs pedigree simp_pedigree chromosomes simplify_removes \
		label1 \
		> $@

pedigree2.pdf : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		time_label ts_indivs simp_indivs pedigree simp_pedigree chromosomes simplify_removes \
		indiv_types \
		> $@

pedigree3.pdf : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		time_label ts_indivs simp_indivs pedigree simp_pedigree chromosomes simplify_removes \
		labels \
		> $@

pedigree_recapitate.pdf : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		ts_indivs simp_indivs pedigree simp_pedigree chromosomes simplify_removes \
		recap_pedigree extra_chromosomes tskit_time recap_label recap \
		> $@

pedigree_simplify.pdf : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		ts_indivs pedigree chromosomes \
		recap_pedigree tskit_time recap \
		> $@

pedigree_mutate.pdf : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		ts_indivs pedigree chromosomes \
		recap_pedigree tskit_time recap \
		mutations \
		> $@

## 
# Graphics whatnot

# save inkscape svg files as .ink.svg and this'll do the right thing
%.svg : %.ink.svg
	inkscape $< --export-plain-svg=$@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.png : %.svg
	convert -density 300 $< -flatten $@
