This program is an example of multiple integral calculation with the help of libepa.
It calculates the following integral in a given range of parameter a:
        a  sqrt(1 - (x/a)^2)  sqrt(1 - (x/a)^2 - y^2)
        /         /                      /          x * y * z         x * (x + 1/2)
 I(a) = | dx      | dy                   | dz --------------------- = -------------
        /         /                      /    sqrt(x^2 + y^2 + z^2)     (x + 1)^2
        0         0                      0

This is the example provided in section 4.1 of [2311.01353].

Two versions of the program are provided, one is written in C++, the other is
in Python.

To compile the C++ version, build libepa.so and install it or set the
environment variable LD_LIBRARY_PATH as

export LD_LIBRARY_PATH="$PWD/../.."

Then execute `make integration`.

To use the Python version, build libepa.so and the python FFI and install the
or set the environment variables LD_LIBRARY_PATH --- as above --- and
PYTHONPATH --- as

export PYTHONPATH="$PWD/../../ffi/python"

To run the example, execute `./integration -f 1 -t 100 -n 100`
(`./integration.py -f 1 -t 100 -n 100` for the Python version), where `-f`
(from) and `-t` (to) define the integration domain, and `-n` defines the
number of points to compute. See `./integration --help` for a short
description and the list of options.

The output of the program is two columns: a and I(a).

Run `make` to produce `integration.pdf` with a plot of the integral computed
as above and compared to the analytical calculation (texlive and gnuplot are
required).
