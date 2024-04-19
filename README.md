A charged particle at rest is surrounded by electric field which can be
represented as a bunch of virtual photons with zero time component of the
photons momenta. If the particle is moving with the speed close to the speed of
light, momenta of the photons of its field acquire longitudinal component, so
that the photons virtuality q<sup>2</sup> &#8810; ω<sup>2</sup>, where q is the
photon 4-momentum, ω is the photon energy. Under the equivalent photon
approximation (EPA) these photons are threated as real and are assumed to be
distributed according to a known spectrum.

In many experiments in high energy physics, collisions of two ultrarelativistic
charged particles are studied. Particles often do not collide head-on, but
rather pass at some distance from each other and collide with their
electromagnetic fields. If both particles survive in the collision, such
collision is called ultraperipheral. Ultraperipheral collisions (UPC) can be
considered as photon-photon collisions with the photons meeting the requirements
of the equivalent photon approximation.

libepa is a library for calculations of cross sections of ultraperipheral
collisions under the equivalent photon approximation. The cross sections are
calculated as integrals over the convolution of the photon EPA spectra with the
photon-photon cross section provided by the user. By default, the integration
is performed with the help of GNU Scientific Library [(GSL)][GSL]. The user can
supply their own integration routines.

The documentation is available as [doc/libepa.html][doc]. An overview of the
theory of EPA and UPC can be found in papers[[1-3]](#references). Paper [3] in
particular is devoted to the libepa approach and discusses the examples that can
be found in the `examples` directory. Cite paper [3] if you want to make a
reference to libepa.

[doc]: https://jini-zh.org/libepa/libepa.html

A Python API is available.

# Installation

Requirements:
* A C++ compiler supporting the C++17 standard.
* GNU Make.
* GNU Scientific Library ([GSL][GSL]).
* Boost [test][boost.test] for tests.
* Python module [cffi][python-cffi] for the Python interface.
* TeX Live for the documentation (should be possible to compile manually with
  any other TeX distribution).
* TeX Live, Perl and Gnuplot for the plots in examples.

[GSL]: https://www.gnu.org/software/gsl/
[boost.test]:  https://www.boost.org/doc/libs/1_84_0/libs/test/doc/html/index.html
[python-cffi]: https://pypi.org/project/cffi/

To compile, execute `make` in the project directory. `Makefile` supports the
following targets:
* `default`: compiles `libepa.so` and the Python interface.
* `all`: also compiles `doc/notes.pdf`.
* `libepa.so`: the libepa shared object file.
* `ffi`: foreign function interface (FFI) for Python.
* `test`: compile and run inexpensive tests; see `test/README.md`.
* `test_all`: compile and run all tests; see `test/README.md`.
* `install`: compile and install `libepa.so`, the headers and the Python
  modules.
* `uninstall`: do the reverse of `install`: delete `libepa.so`, the headers and
  the Python modules from the system.

Also `Makefile` supports the standard GNU Make conventions:
* [`prefix`][make-dirs]: installation root. Used to construct the variables
  below. Default is `/usr/local`.
* [`exec_prefix`][make-dirs]: prefix for machine-specific files. Used for
  `libdir`. Default is `$(prefix)`.
* [`includedir`][make-dirs]: the directory for installing header files. Default
  is `$(prefix)/include`.
* [`libdir`][make-dirs]: the directory for the library shared object file.
  Default is `$(exec_prefix)/lib`.
* [`DESTDIR`][make-destdir]: prefix prepended to each file name during
  installation.

[make-dirs]: https://www.gnu.org/prep/standards/html_node/Directory-Variables.html
[make-destdir]: https://www.gnu.org/prep/standards/html_node/DESTDIR.html

So, to compile and install the library and the FFI to the default locations,
execute

    make install

To compile and install the library and the FFI to `~/usr`, execute

    make prefix=$HOME/usr install

To compile and install the library to `/usr/lib64` and the header files to
`/usr/include`, execute

    make prefix=/usr libdir=/usr/lib64 install

When installing system-wide, do remember to execute `ldconfig` as root; the
target directory (`/usr/local/lib` by default) should be included in
`/etc/ld.so.conf`.

Note that Python modules are installed via simple file copying to
[`platlib`][python-sysconfig] with `$(prefix)` as `platbase` (this usually
resolves to something like `/usr/lib/python3.11/site-packages`). If you want to
have them installed into a different location, modify `ffi/python/install.py` or
do it manually (see the script for the list of files). If you want to have them
installed by `pip`, you can use `ffi/python/setup.py`.

[python-sysconfig]: https://docs.python.org/3/library/sysconfig.html

To uninstall, execute `make uninstall`.

# Directory structure

* `doc`: library documentation.
* `examples`: examples of how to use the library.
* `ffi`: foreign function interface implementations.
* `include`: headers.
* `src`: sources.
* `test`: tests.

See `README.md` in each of the directories for more.

# Acknowledgements

libepa development was supported by the Russian Science Foundation grant No 19-12-00123.

# References

1. M. I. Vysotsky, E. V. Zhemchugov.
   Equivalent photons in proton-proton and ion-ion collisions at the LHC.
   [Physics Uspekhi 62, 910 (2019)](http://dx.doi.org/10.3367/UFNe.2018.07.038389);
   [arXiv:1806.07238](https://arxiv.org/abs/1806.07238).
2. S. I. Godunov, V. A. Novikov, A. N. Rozanov, M. I. Vysotsky, E. V. Zhemchugov.
   Production of heavy charged particles in proton-proton ultraperipheral collisions at the Large Hadron Collider: survival factor.
   [Journal of High Energy Physics 2021, 234 (2021)](https://doi.org/10.1007/JHEP10%282021%29234);
   [arXiv:2106.14842](https://arxiv.org/abs/2106.14842).
3. E. V. Zhemchugov, S. I. Godunov, E. K. Karkaryan, V. A. Novikov, A. N. Rozanov, M. I. Vysotsky.
   libepa &mdash; a C++/Python library for calculations of cross sections of ultraperipheral collisions.
   [arxiv:2311.01353](https://arxiv.org/abs/2311.01353).

[GSL]: https://www.gnu.org/software/gsl
[notes]: https://jini-zh.org/libepa/notes.pdf
