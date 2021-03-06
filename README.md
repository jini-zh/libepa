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

The documentation is under development. Currently the main formulas are
collected in [doc/notes.pdf][notes]. An overview of the theory of EPA and UPC
can be found in papers [[1,2]](#references). libepa is based on the code
developed to perform calculations in these papers. Sample code that calculates
some of the plots published in [2] can be found in the `examples/` directory
(actual values can be slightly different because libepa uses more recent
values of physical constants, in particular the proton charge radius).

A Python API is planned but not implemented yet.

# Installation

Requirements:
* A C++ compiler from the GNU Compiler Collection (GCC) supporting the C++17 standard.
* GNU Make.
* GNU Scientific Library.
* TeX Live for the documentation (should be possible to compile manually with
  any other TeX distribution).
* TeX Live, Perl and Gnuplot for the plots in examples.

To compile, simply execute `make` in the project directory. To build and run an
example, execute `make` in its directory. Note that calculation of most examples
takes time and will use the number of threads equal to the number of logical
CPU cores at the system.

The `make install` target is not implemented yet.

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

[GSL]: https://www.gnu.org/software/gsl
[notes]: https://jini-zh.org/libepa/notes.pdf
