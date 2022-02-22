.PHONY: clean

notes.pdf: notes.tex
	pdflatex $<
	@if grep -q 'Warning.* Rerun ' notes.log; then pdflatex $<; fi

clean:
	rm notes.{aux,log,out,pdf} 2> /dev/null; true
