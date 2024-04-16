This example calculates the cross section for fermion pair production in
ultraperipheral proton-proton collisions with respect to the fermion mass. See
section 4.2 of [2311.01353] for the discussion.


Two versions of the program are provided, one is written in C++, the other is
in Python.

To compile the C++ version, build libepa.so and install it or set the
environment variable LD_LIBRARY_PATH as

export LD_LIBRARY_PATH="$PWD/../.."

Then execute `make xsection`.

To use the Python version, build libepa.so and the python FFI and install the
or set the environment variables LD_LIBRARY_PATH --- as above --- and
PYTHONPATH --- as

export PYTHONPATH="$PWD/../../ffi/python"

To run the example, execute `./xsection -f 90 -t 250 -n 32`, where `-f` (from)
and `-t` (to) define the range of invariant mass in GeV, and `-n` defines the
number of points. The program will print two columns of numbers: the fermion
mass in GeV and the cross section in barns. By default the cross section is
calculated neglecting non-electromagnetic interactions. To compute the cross
section taking into account non-electromagnetic interactions, add the `-S`
key. This will take a long time. With the `-v` key, the program will print
intermediate results to the standard error, allowing you to monitor the
process. See `./xsection --help` for more options.

Run `make` to produce `xsection.pdf` with the cross sections and their ratio
(texlive and gnuplot are require; note that calculation of the cross section
taking into account non-electromagnetic interactions will take a long time).
