CXXFLAGS ?= -O2 -pipe -march=native -fno-stack-protector

cxx := $(CXX) $(CXXFLAGS) -std=gnu++17

.NOTPARALLEL:

atlas.pdf: atlas.plot ../epslatex2pdf differential.txt integrated.txt ratio.txt differential_b.txt integrated_b.txt ratio_b.txt
	gnuplot $<
	../epslatex2pdf atlas-plot.tex -o $@
	@rm atlas-plot.tex atlas-plot.pdf

differential.txt: atlas
	LD_LIBRARY_PATH=../.. ./$< -n 11 -v > $@

differential_b.txt: atlas
	LD_LIBRARY_PATH=../.. ./$< -n 11 -v -S > $@

integrated.txt: integrate differential.txt
	./$^ > $@

integrated_b.txt: integrate differential_b.txt
	./$^ > $@

ratio.txt: ratio integrated.txt
	./$^ > $@

ratio_b.txt: ratio integrated_b.txt
	./$^ > $@

atlas: atlas.cpp ../common.cpp ../../libepa.so
	$(cxx) -I ../../include $< -o $@ -L ../.. -lepa `pkg-config --libs gsl` -lpthread

clean:
	rm -f atlas atlas.pdf
	rm -f atlas-plot.tex atlas-plot.pdf
	rm -f differential.txt differential_b.txt
	rm -f integrated.txt integrated_b.txt
	rm -f ratio.txt ratio_b.txt
