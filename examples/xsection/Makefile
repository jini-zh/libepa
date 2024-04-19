CXXFLAGS ?= -O2 -pipe -march=native -fno-stack-protector

cxx := $(CXX) $(CXXFLAGS) -std=gnu++17

.NOTPARALLEL: xsection.pdf

xsection.pdf: xsection.plot ../epslatex2pdf ratio.txt
	gnuplot $<
	../epslatex2pdf xsection-plot.tex -o $@
	@rm xsection-plot.tex xsection-plot.pdf

xsection.txt: xsection
	LD_LIBRARY_PATH=../.. ./$< -v -f 90 -t 250 -n 32 > $@ || rm $@

xsection_b.txt: xsection
	LD_LIBRARY_PATH=../.. ./$< -v -f 90 -t 250 -n 32 -S --rs-max=3e3 > $@ || rm $@

ratio.txt: ../ratio xsection.txt .WAIT xsection_b.txt
	$< xsection_b.txt xsection.txt > $@ || rm $@

xsection: xsection.cpp ../common.cpp ../../libepa.so
	$(cxx) -I ../../include $< -o $@ -L ../.. -lepa `pkg-config --libs gsl` -lpthread

clean:
	rm -f xsection xsection.pdf
	rm -f xsection-plot.tex xsection-plot.pdf
	rm -f xsection.txt xsection_b.txt ratio.txt
