CXXFLAGS ?= -O2 -pipe -march=native -fno-stack-protector

cxx := $(CXX) $(CXXFLAGS) -std=gnu++17

luminosity := LD_LIBRARY_PATH=../.. ./luminosity -v -f 1 -t 3e3 -n 1000

.NOTPARALLEL: survival.pdf

survival.pdf: survival.plot ../epslatex2pdf survival.txt .WAIT parallel.txt .WAIT perpendicular.txt
	gnuplot survival.plot
	../epslatex2pdf survival-plot.tex -o $@
	@rm survival-plot.tex survival-plot.pdf

luminosity.txt: luminosity
	$(luminosity) > $@

luminosity_b.txt: luminosity
	$(luminosity) -S > $@

parallel.txt: luminosity
	$(luminosity) -p parallel > $@

perpendicular.txt: luminosity
	$(luminosity) -p perpendicular > $@

survival.txt: ../ratio luminosity.txt .WAIT luminosity_b.txt
	$< luminosity_b.txt luminosity.txt > $@

luminosity: luminosity.cpp ../common.cpp ../../libepa.so
	$(cxx) -I ../../include $< -o $@ -L ../.. -lepa `pkg-config --libs gsl` -lpthread

clean:
	rm -f luminosity survival.pdf
	rm -f survival-plot.tex survival-plot.pdf
	rm -f survival.txt luminosity.txt luminosity_b.txt parallel.txt perpendicular.txt
