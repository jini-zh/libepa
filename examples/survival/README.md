This example calculates the survival factor in proton-proton collisions. The
survival factor is a ratio of two luminosities (neglecting the photon fusion
cross section dependence on the photons polarizations): the luminosity taking
into account non-electromagnetic interactions divided by the luminosity
neglecting them. The program provided here calculates the luminosities, and
`../ratio` calculates the ratio. See section 4.3 of [2311.01353] for the
discussion of this example; see [2106.14842] for the discussion of the
survival factor in proton-proton collisions.

Two versions of the program are provided, one is written in C++, the other is
in Python.

To compile the C++ version, build `libepa.so`, and install it or set the
environment variable `LD_LIBRARY_PATH` as

    export LD_LIBRARY_PATH="$PWD/../.."

Then execute `make luminosity`.

To use the Python version, build `libepa.so`, and the python FFI and install
them or set the environment variables `LD_LIBRARY_PATH` &mdash; as above
&mdash; and `PYTHONPATH` &mdash; as

    export PYTHONPATH="$PWD/../../ffi/python"

To compute a luminosity neglecting non-electromagnetic interactions, execute
`./luminosity -f 1 -t 3e3 -n 1000` (`./luminosity.py -f 1 -t 3e3 -n 1000` for
the Python version), where `-f` (from) and `-t` (to) define the range of the
photon-photon invariant masses in GeV and `-n` defines the number of points to
compute. The program will print two columns of output: the invariant mass of
the photons in GeV and the luminosity in GeV^{-1}. To compute the luminosity
taking into account non-electromagnetic interactions, add the `-S` key. This
will take a long time. With the `-v` key, the program will print intermediate
results to the standard error, allowing you to monitor the process. See
`./luminosity --help` for other options. 

Run `make` to produce `survival.pdf` with the luminosities and their ratio
(texlive and gnuplot are required; note that calculation of the luminosity
taking into account non-electromagnetic interactions will take a long time).

[2106.14842]: https://arxiv.org/abs/2106.14842
[2311.01353]: https://arxiv.org/abs/2311.01353
