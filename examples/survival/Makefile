CXXFLAGS ?= -O2 -pipe -march=native -fno-stack-protector

cxx := $(CXX) $(CXXFLAGS) -std=gnu++17

.NOTPARALLEL: %.txt

survival.pdf: survival.plot ../epslatex2pdf ../reformat-logscale survival.txt luminosity.txt luminosity_b.txt parallel.txt perpendicular.txt
	gnuplot survival.plot
	../reformat-logscale survival-plot.tex
	../epslatex2pdf survival-plot.tex -o $@
	rm survival-plot.{tex,pdf}

luminosity.txt luminosity_b.txt parallel.txt perpendicular.txt: survival
	LD_LIBRARY_PATH=../.. ./survival -v -f 1 -t 3e3 -n 1000 --$(basename $@) $@

survival.txt: luminosity.txt luminosity_b.txt
	paste $^ | perl -ae 'die "Grid misaligned: $$F[0] != $$F[2] at line $$.\n" if $$F[0] != $$F[2]; printf "%s %19.12e\n", $$F[0], $$F[3] / $$F[1]' > $@

survival: survival.cpp ../common.cpp ../../libepa.so
	$(cxx) -iquote ../../src $< -o $@ -L ../.. -lepa `pkg-config --libs gsl` -lpthread

clean:
	rm survival{,.pdf} survival-plot.{tex,pdf} {survival,luminosity{,_b},parallel,perpendicular}.txt 2>/dev/null; true
