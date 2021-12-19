export SHELL=/bin/bash
.PRECIOUS : recipe_17.4.trees

check : recipe_17.1_overlaid.log recipe_17.3.log recipe_17.4.png recipe_17.5.png recipe_17.7.log recipe_17.8_III.log recipe_17.9.log recipe_nucleotides_II.log recipe_nucleotides_III.log
	echo "-------------"
	echo "17.1 + 17.2"
	cat recipe_17.1_overlaid.log
	echo "-------------"
	echo "17.3"
	cat recipe_17.3.log
	echo "-------------"
	echo "17.4"
	eog recipe_17.4.png
	echo "-------------"
	echo "17.5"
	eog recipe_17.5.png
	echo "-------------"
	echo "17.7"
	cat recipe_17.7.log
	echo "-------------"
	echo "17.8"
	head  recipe_17.8_III.log
	echo "(snipped $$(wc -l recipe_17.8_III.log) lines)"
	echo "-------------"
	echo "17.9"
	cat recipe_17.9.log 
	echo "-------------"
	echo "18.13 II"
	cat recipe_nucleotides_II.log
	echo "-------------"
	echo "18.13 III"
	head recipe_nucleotides_III.log
	echo "(snipped $$(wc -l recipe_nucleotides_III.log) lines)"
	echo "-------------"

clean :
	-rm -f *.trees *.log *.png

sanitize : 
	for x in *; do mv -v "$x" "${x/:/-}"; done
	for x in *; do mv -v "$x" "${x/(/-}"; done
	for x in *; do mv -v "$x" "${x/)/-}"; done

recipe_17.1.trees : Recipe\ 17.1\ -\ A\ minimal\ tree-seq\ model.txt
	slim -s 5 "$<" &> $@.log

recipe_17.1_overlaid.trees : Recipe\ 17.2\ -\ Overlaying\ neutral\ mutations.py recipe_17.1.trees
	python3 "$<"

recipe_17.1_overlaid.log : recipe_17.1_overlaid.trees
	./check_trees.py $< > $@

slim_2345_FIXED.trees : Recipe\ 17.3\ -\ Simulation\ conditional\ upon\ fixation\ of\ a\ sweep,\ preserving\ ancestry\ II.txt
	slim -s 2345 "$<" &> $@.log

recipe_17.3.log : slim_2345_FIXED.trees
	./check_trees.py $< > $@

recipe_17.4.trees : Recipe\ 17.4\ -\ Detecting\ the\ dip\ in\ diversity\ (analyzing\ tree\ heights\ in\ Python)\ I.txt
	slim -s 2345 '$<' &> $@.log

recipe_17.4.png : Recipe\ 17.4\ -\ Detecting\ the\ dip\ in\ diversity\ (analyzing\ tree\ heights\ in\ Python)\ II.py recipe_17.4.trees
	python3 '$<'

recipe_17.5.trees: Recipe\ 17.5\ -\ Mapping\ admixture\ (analyzing\ ancestry\ in\ Python)\ I.txt
	slim -s 2345 "$<" &> $@.log

recipe_17.5.png : Recipe\ 17.5\ -\ Mapping\ admixture\ (analyzing\ ancestry\ in\ Python)\ II.py recipe_17.5.trees
	python3 "$<"

recipe_17.7.trees : Recipe\ 17.7\ -\ Analyzing\ selection\ coefficients\ in\ Python\ with\ pyslim\ I.txt
	slim -s 2345 "$<" &> $@.log

recipe_17.7.log : Recipe\ 17.7\ -\ Analyzing\ selection\ coefficients\ in\ Python\ with\ pyslim\ II.py recipe_17.7.trees
	python3 "$<" &> $@

recipe_17.8.trees : Recipe\ 17.8\ -\ Starting\ a\ hermaphroditic\ WF\ model\ with\ a\ coalescent\ history\ I.py
	python3 "$<"

recipe_17.8_II.trees : Recipe\ 17.8\ -\ Starting\ a\ hermaphroditic\ WF\ model\ with\ a\ coalescent\ history\ II.txt recipe_17.8.trees
	slim -s 2345 "$<" &> $@.log

recipe_17.8_III.log : Recipe\ 17.8\ -\ Starting\ a\ hermaphroditic\ WF\ model\ with\ a\ coalescent\ history\ III.py recipe_17.8_II.trees
	python3 "$<" &> $@

recipe_17.9.trees : Recipe\ 17.9\ -\ Starting\ a\ sexual\ nonWF\ model\ with\ a\ coalescent\ history\ I.py
	python3 "$<"

recipe_17.9.log : Recipe\ 17.9\ -\ Starting\ a\ sexual\ nonWF\ model\ with\ a\ coalescent\ history\ II.txt recipe_17.9.trees
	slim -s 2345 "$<" &> $@

recipe_nucleotides.trees : Recipe\ 18.13\ -\ Tree-sequence\ recording\ and\ nucleotide-based\ models\ I.txt
	slim -s 2345 "$<" &> $@.log

recipe_nucleotides_II.log : Recipe\ 18.13\ -\ Tree-sequence\ recording\ and\ nucleotide-based\ models\ II.py recipe_nucleotides.trees
	python3 "$<" &> $@

recipe_nucleotides_III.log : Recipe\ 18.13\ -\ Tree-sequence\ recording\ and\ nucleotide-based\ models\ III.py recipe_nucleotides.trees
	python3 "$<" &> $@
