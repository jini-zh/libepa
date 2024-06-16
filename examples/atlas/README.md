This example calculates the fiducial differential cross section for muon pair
production in ultraperipheral proton-proton collisions.  The fiducial region is
defined by the cuts on transverse momentum pT and pseudorapidity eta of each
muon:

     pT >  6 GeV for 12 GeV < sqrt(s) < 30 GeV,
     pT > 10 GeV for 30 GeV < sqrt(s) < 70 GeV;
     abs(eta) < 2.4.

The goal of the calculation is to reproduce the experimental data in
[1708.04053]. See section 4.4 in [2311.01353] for the discussion and
[2106.14842], [1806.07238] for previous calculations.

Two versions of the program are provided, one is written in C++, the other is
in Python.

To compile the C++ version, build libepa.so, and install it or set the
environment variable `LD_LIBRARY_PATH` as

    export LD_LIBRARY_PATH="$PWD/../.."

Then execute `make atlas`.

To use the Python version, build `libepa.so` and the python FFI, and install
them or set the environment variables `LD_LIBRARY_PATH` --- as above
--- and `PYTHONPATH` --- as

    export PYTHONPATH="$PWD/../../ffi/python"

To run the example, execute `./atlas`. The program will print two columns of
numbers and a commented line to the standard output. The first column is the
invariant mass of the muon pair, in GeV. The second column is the differential
cross section, in barn/GeV. The commented line provides the integrated cross
section. The calculation by default does not take into account
non-electromagnetic interactions. To take into account non-electromagnetic
interactions, execute `./atlas -S`. This will take a long time. With the `-v`
key, the program will print intermediate results to the standard error, allowing
you to monitor the process.  The keys `-n` and `-s` allow changing the number of
points in the output. For example, `./atlas -n 10` will split the range of
invariant masses [12:70] GeV into 10 segments. The mass 30 GeV will always
appear twice, for the two different cuts imposed on the transverse momentum pT
(not implemented in the Python version).

Run `./atlas --help` for a short description of the program and its arguments.

Run `make` to produce `atlas.pdf` with the plot of the calculation (texlive and
gnuplot are required; also it will take a long time to compute the cross
sections with non-electromagnetic interactions taken into account).

[1708.04053]: https://arxiv.org/abs/1708.04053
[1806.07238]: https://arxiv.org/abs/1806.07238
[2106.14842]: https://arxiv.org/abs/2106.14842
[2311.01353]: https://arxiv.org/abs/2311.01353
