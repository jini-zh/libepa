CXXFLAGS ?= -O2 -pipe -march=native -fno-stack-protector

cxx = $(CXX) $(CXXFLAGS) -std=gnu++17 -I include

.PHONY: clean test all test_all

objects := gsl algorithms epa proton
objects := $(foreach object,$(objects),src/$(object).o)

all: libepa.so

libepa.so: $(objects)
	$(cxx) -shared $^ -o $@

src/gsl.o: include/epa/gsl.hpp
src/epa.o: include/epa/epa.hpp include/epa/gsl.hpp include/epa/algorithms.hpp
src/proton.o: include/epa/proton.hpp include/epa/epa.hpp include/epa/gsl.hpp \
	include/epa/algorithms.hpp
src/algorithms.o: include/epa/algorithms.hpp

%.o: %.cpp
	$(cxx) -fPIC -c $< -o $@

test_all: test/test
	LD_LIBRARY_PATH=. $< -l success -t '*'

test: test/test
	LD_LIBRARY_PATH=. $< -l success

test/test: test/test.o libepa.so
	$(cxx) $< -o $@ -L . -lepa `pkg-config --libs gsl` -lboost_unit_test_framework-mt

test/test.o: test/test.cpp test/a1.cpp include/epa/proton.hpp \
	include/epa/epa.hpp include/epa/gsl.hpp include/epa/algorithms.hpp
	$(cxx) -iquote test -c $< -o $@

test/a1.cpp: test/make-a1-form-factor test/a1.dat
	$< test/a1.dat > $@

clean:
	rm $(objects) libepa.so 2> /dev/null; true
