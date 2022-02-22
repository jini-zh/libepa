CXXFLAGS ?= -O2 -pipe -march=native -fno-stack-protector

cxx := $(CXX) $(CXXFLAGS) -std=gnu++17

luminosity := LD_LIBRARY_PATH=../.. ./luminosity -v -f 1 -t 3e3 -n 1000

.NOTPARALLEL: %.txt

survival.pdf: survival.plot ../epslatex2pdf ../reformat-logscale survival.txt luminosity.txt luminosity_b.txt parallel.txt perpendicular.txt
	gnuplot survival.plot
	../reformat-logscale survival-plot.tex
	../epslatex2pdf survival-plot.tex -o $@
	rm survival-plot.{tex,pdf}

luminosity.txt: luminosity
	$(luminosity) > $@

luminosity_b.txt: luminosity
	$(luminosity) -S > $@

parallel.txt: luminosity
	$(luminosity) -p parallel > $@

perpendicular.txt: luminosity
	$(luminosity) -p perpendicular > $@

survival.txt: ../ratio luminosity.txt luminosity_b.txt
	$< luminosity_b.txt luminosity.txt > $@

luminosity: luminosity.cpp ../common.cpp ../../libepa.so
	$(cxx) -I ../../include $< -o $@ -L ../.. -lepa `pkg-config --libs gsl` -lpthread

clean:
	rm luminosity survival.pdf survival-plot.{tex,pdf} {survival,luminosity{,_b},parallel,perpendicular}.txt 2>/dev/null; true
