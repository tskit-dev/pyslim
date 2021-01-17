.PHONY : all

# NOTE: all layers must be visible in pedigrees.ink.svg for this to work.

all: pedigree01.png pedigree2.png pedigree3.png pedigree_recapitate.png pedigree_simplify.png pedigree_mutate.png

pedigree01.png : pedigree0.png pedigree1.png
	convert +append $^ $@

pedigree0.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		ts_indivs simp_indivs pedigree simp_pedigree chromosomes simplify_removes \
		extra_chromosomes extra_indivs label0 \
		> $@

pedigree1.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		time_label ts_indivs simp_indivs pedigree simp_pedigree chromosomes simplify_removes \
		label1 \
		> $@

pedigree2.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		time_label ts_indivs simp_indivs pedigree simp_pedigree chromosomes simplify_removes \
		indiv_types \
		> $@

pedigree3.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		time_label ts_indivs simp_indivs pedigree simp_pedigree chromosomes simplify_removes \
		labels \
		> $@

pedigree_recapitate.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		ts_indivs simp_indivs pedigree simp_pedigree chromosomes simplify_removes \
		recap_pedigree extra_chromosomes tskit_time recap_label recap \
		> $@

pedigree_simplify.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		ts_indivs pedigree chromosomes \
		recap_pedigree tskit_time recap \
		> $@

pedigree_mutate.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		ts_indivs pedigree chromosomes \
		recap_pedigree tskit_time recap \
		mutations \
		> $@

## 
# Graphics whatnot

%.png : %.svg
	inkscape $<  --export-type="png" --export-dpi=300

clean :
	-rm pedigree*.png