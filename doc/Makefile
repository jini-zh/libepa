.PHONY: clean

notes.pdf: notes.tex
	pdflatex $<
	@if grep -q 'Warning.* Rerun ' notes.log; then pdflatex $<; fi

clean:
	rm -vf notes.aux notes.log notes.out notes.pdf
