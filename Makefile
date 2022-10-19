CXXFLAGS ?= -O2 -pipe -march=native -fno-stack-protector

cxx = $(CXX) $(CXXFLAGS) -std=gnu++20 -I include

.PHONY: clean test all test_all doc ffi

objects := gsl algorithms epa proton
objects := $(foreach object,$(objects),src/$(object).o)

ffi := ffi epa proton
ffi := $(foreach object,$(ffi),ffi/c/$(object).o)

all: libepa.so doc ffi

libepa.so: $(objects) $(ffi)
	$(cxx) -shared $^ -o $@

src/gsl.o: include/epa/gsl.hpp
src/epa.o: include/epa/epa.hpp include/epa/gsl.hpp include/epa/algorithms.hpp
src/proton.o: include/epa/proton.hpp include/epa/epa.hpp include/epa/gsl.hpp \
	include/epa/algorithms.hpp
src/algorithms.o: include/epa/algorithms.hpp

ffi: $(ffi) ffi/python/epa/_epa_cffi.so

ffi/c/epa.o: ffi/c/epa.cpp ffi/c/epa.h ffi/c/epa_vars.h ffi/c/ffi.hpp ffi/c/ffi.h
ffi/c/proton.o: ffi/c/proton.cpp ffi/c/ffi.hpp ffi/c/ffi.h
ffi/c/ffi.o: ffi/c/ffi.cpp ffi/c/ffi.hpp ffi/c/ffi.h

ffi/python/epa/_epa_cffi.so: ffi/python/build.py ffi/c/ffi.h ffi/c/epa.h \
	ffi/c/epa_vars.h ffi/c/proton.h | libepa.so
	cd ffi/python && LD_LIBRARY_PATH=../.. ./build.py
	mv -v ffi/python/epa/_epa_cffi.cpython-*.so ffi/python/epa/_epa_cffi.so

%.o: %.cpp
	$(cxx) -fPIC -c $< -o $@

test_all: test/test
	LD_LIBRARY_PATH=. $< -l success -t '*'

test: test/test
	LD_LIBRARY_PATH=. $< -l success

test/test: test/test.o libepa.so
	$(cxx) $< -o $@ -L . -lepa `pkg-config --libs gsl` -lboost_unit_test_framework

test/test.o: test/test.cpp test/a1.cpp include/epa/proton.hpp \
	include/epa/epa.hpp include/epa/gsl.hpp include/epa/algorithms.hpp
	$(cxx) -iquote test -c $< -o $@

test/a1.cpp: test/make-a1-form-factor test/a1.dat
	$< test/a1.dat > $@

doc: doc/notes.pdf

doc/notes.pdf: doc/notes.tex
	make -C doc

clean:
	rm $(objects) libepa.so $(ffi) ffi/python/epa/{_epa_cffi.*,_epa_functions.py,_epa_vars.py} 2> /dev/null; true
	rm -r ffi/python/epa/__pycache__ 2> /dev/null; true
	make -C doc clean
