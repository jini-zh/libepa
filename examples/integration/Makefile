CXXFLAGS ?= -O2 -pipe -march=native -fno-stack-protector

cxx := $(CXX) $(CXXFLAGS) -std=gnu++17

integration.pdf: integration.plot ../epslatex2pdf integration.txt
	gnuplot $<
	../epslatex2pdf integration-plot.tex -o $@
	@rm integration-plot.tex integration-plot.pdf

integration.txt: integration
	LD_LIBRARY_PATH=../.. ./$< -v -f 1 -t 100 -n 100 > $@

integration: integration.cpp ../common.cpp ../../libepa.so
	$(cxx) -I ../../include $< -o $@ -L ../.. -lepa `pkg-config --libs gsl` -lpthread

clean:
	rm -f integration integration.pdf integration.txt integration-plot.tex integration-plot.pdf
