#pragma once

#include <memory>

#include "algorithms.hpp"
#include "gsl.hpp"

namespace epa {

const double pi              = M_PI;
const double alpha           = 7.2973525693e-3;        // [1] (fine structure constant)
const double planck          = 0.1973269804;           // GeV fm [1]
const double barn            = sqr(planck) * 1e-2;     // (GeV * barn)^{-2}
const double fm              = 1. / (10 * sqrt(barn)); // (GeV fm)^{-1}
const double electron_charge = 1.602176634e-19;        // C [1]
const double light_speed     = 299792458;              // m/s [1]
const double amu
  = (1.66053906660e-27 * sqr(light_speed)) / electron_charge * 1e-9; // GeV, [1]
const double infinity        = gsl::infinity;

/*
 [1] E. Tiesinga, P. J. Mohr, D. B. Newell, B. N. Taylor.
     CODATA recommended values of the fundamental physical constants: 2018.
     Reviews of Modern Physics 93, 025010 (2020).
*/

extern bool print_backtrace;

#define EPA_TRY try {
#define EPA_BACKTRACE(message, ...) \
  } catch (std::exception&) { \
    if (print_backtrace) fprintf(stderr, message "\n", __VA_ARGS__); \
    throw; \
  }

// default parameters for GSL routines for quadrature integration
extern double default_absolute_error;
extern double default_relative_error;
extern double default_error_step;
extern size_t default_integration_limit;
extern size_t default_cquad_integration_limit;
extern gsl::integration::QAGMethod default_integration_method;

// Function: f, a, b -> integral of f from a to b
typedef std::function<
          double (const std::function<double (double)>&, double, double)
        > Integrator;

// GSL quadratic adaptive integrator with default initialization
Integrator gsl_integrator(
    double absolute_error = default_absolute_error,
    double relative_error = default_relative_error,
    gsl::integration::QAGMethod = gsl::integration::GAUSS41,
    std::shared_ptr<gsl::integration::QAGWorkspace> = nullptr
);

struct gsl_integrator_keys {
  double absolute_error = default_absolute_error;
  double relative_error = default_relative_error;
  gsl::integration::QAGMethod method = gsl::integration::GAUSS41;
  std::shared_ptr<gsl::integration::QAGWorkspace> workspace;
};

// Helper function --- with keyword parameters
Integrator gsl_integrator(const gsl_integrator_keys&);

// GSL integrator with relative_error = default_relative_error *
// default_error_step ** level
Integrator gsl_integrator(unsigned level = 0);

// Default integrator
const std::function<Integrator (unsigned)> default_integrator
  = static_cast<Integrator (*)(unsigned)>(gsl_integrator);

// GSL CQUAD integrator with default initialization
Integrator gsl_cquad_integrator(
    double absolute_error = default_absolute_error,
    double relative_error = default_relative_error,
    std::shared_ptr<gsl::integration::CQuadWorkspace> = nullptr
);

struct gsl_cquad_integrator_keys {
  double absolute_error = default_absolute_error;
  double relative_error = default_relative_error;
  std::shared_ptr<gsl::integration::CQuadWorkspace> workspace;
};

Integrator gsl_cquad_integrator(const gsl_cquad_integrator_keys&);

Integrator gsl_cquad_integrator(unsigned level = 0);

// Integrator that takes one extra parameter: the value of the integral
// calculated so far. It can be used to avoid extra work in calculations that
// don't need that much accuracy. This kind of integrator is currently used in
// spectrum_b_function1d* functions to evaluate the integral of the form factor
// in the region after the one described by Function1d. 
typedef std::function<
          double (const std::function<double (double)>&, double, double, double)
        > Integrator_I;

Integrator_I gsl_integrator_i(
    double relative_error = default_relative_error,
    gsl::integration::QAGMethod = gsl::integration::GAUSS41,
    std::shared_ptr<gsl::integration::QAGWorkspace> = nullptr
);

struct gsl_integrator_i_keys {
  double relative_error = default_relative_error;
  gsl::integration::QAGMethod method = gsl::integration::GAUSS41;
  std::shared_ptr<gsl::integration::QAGWorkspace> workspace;
};

Integrator_I gsl_integrator_i(const gsl_integrator_keys&);
Integrator_I gsl_integrator_i(unsigned level = 0);

const std::function<Integrator_I (unsigned)> default_integrator_i
  = static_cast<Integrator_I (*)(unsigned)>(gsl_integrator_i);

// Electromagnetic form factor of a particle. Q2 is the photon 3-momentum
// squared
typedef std::function<double (double /* Q2 */)> FormFactor;

FormFactor form_factor_monopole(double lambda2);
FormFactor form_factor_dipole(double lambda2);

// Equivalent photon spectrum integrated across transversal plane. w is the
// photon energy
typedef std::function<double (double /* w */)> Spectrum;

// Equivalent photon spectrum of a particle with charge Ze moving with the
// Lorentz factor gamma.  This is the generic function. See below for optimized
// versions for specific form factors
Spectrum spectrum(
    unsigned Z,
    double gamma,
    FormFactor,
    Integrator = default_integrator(0)
);

// EPA spectrum for monopole form factor with the parameter lambda^2
Spectrum spectrum_monopole(unsigned Z, double gamma, double lambda2);

// EPA spectrum for dipole form factor with the parameter lambda^2
Spectrum spectrum_dipole(unsigned Z, double gamma, double lambda2);

// Equivalent photon spectrum at distance b from the source particle in the
// transversal plane. w is the photon energy
typedef std::function<double (double /* b */, double /* w */)> Spectrum_b;

// Equivalent photon spectrum of a particle with charge Ze moving with the
// Lorentz factor gamma. This is the generic function. See below for optimized
// versions for specific form factors.
//
// Note that the form factor is weighted with an oscillating function (the
// Bessel function J1) and integrated over semi-infinite interval. Unless the
// form factor falls very rapidly, this integration is very difficult to
// perform numerically. Take care with the default quadrature integrator.
Spectrum_b spectrum_b(
    unsigned Z,
    double gamma,
    FormFactor,
    Integrator = default_integrator(0)
);

// EPA spectrum for point-like particle
Spectrum_b spectrum_b_point(unsigned Z, double gamma);

// EPA spectrum for monopole form factor
Spectrum_b spectrum_b_monopole(unsigned Z, double gamma, double lambda2);

// EPA spectrum for dipole form factor
Spectrum_b spectrum_b_dipole(unsigned Z, double gamma, double lambda2);

// EPA spectrum for form factor given by Function1d as a set of points (q2, ff)
// for 0 <= q2 <= q2_max. "g" stands for "global": the form factor is integrated
// from 0 to q2_max as a regular function of one variable,
// with special consideration given only to the endpoint q2_max (cf.
// spectrum_b_function1d_qs).
//
// rest_form_factor: form factor to use for q2 > q2_max (optional). Multiplied
// by (*form_factor)(q2_max) / rest_form_factor(q2_max) so that the form factor
// is continuous.
//
// rest_spectrum: EPA spectrum for rest_form_factor (optional). Can be used to
// improve convergence. If provided, rest_form_factor is integrated from 0 to
// q2_max, otherwise --- from q2_max to infinity.
//
// b_max: for b > b_max > 0, assume that the source particle is point-like and
// use spectrum_b_pointlike.
//
// This function uses slightly different integration interface:
// integrate parameters are:
//   the function to integrate
//   integration lower bound
//   integration upper bound
//   absolute value of already calculated part of the integral. It can be used
//     to speed up convergence by requiring that absolute error of numerical
//     integration is something like 1e-3 times the value of this parameter.
Spectrum_b
spectrum_b_function1d_g(
    unsigned z,
    double gamma,
    std::shared_ptr<Function1d> form_factor,
    FormFactor rest_form_factor = FormFactor(),
    Spectrum_b rest_spectrum    = Spectrum_b(),
    double b_max = 0,
    Integrator_I = default_integrator_i(0)
);

// EPA spectrum for form factor given by Function1d as a set of points (q2, ff)
// for 0 <= q2 <= q2_max. "s" stands for "segmented": the form factor is
// integrated in between each known point. This calculation takes a long time
// compared to spectrum_b_function1d_g but has better convergence.
//
// See spectrum_b_function1d_g for the description of arguments.
Spectrum_b
spectrum_b_function1d_s(
    unsigned z,
    double gamma,
    std::shared_ptr<Function1d> form_factor,
    FormFactor rest_form_factor = FormFactor(),
    Spectrum_b rest_spectrum    = Spectrum_b(),
    double b_max = 0,
    Integrator_I = default_integrator_i(0)
);

// EPA spectrum for form factor given by Function1d as a set of points (q2, ff)
// for 0 <= q2 <= q2_max. Uses spectrum_b_function1d_g and falls back to
// spectrum_b_function1d_s in the case of convergency failure.
Spectrum_b
spectrum_b_function1d(
    unsigned z,
    double gamma,
    std::shared_ptr<Function1d> form_factor,
    FormFactor rest_form_factor = FormFactor(),
    Spectrum_b rest_spectrum    = Spectrum_b(),
    double b_max = 0,
    Integrator_I = default_integrator_i(0)
);

// Photon-photon luminosity with non-electromagnetic interactions neglected
typedef std::function<double (double /* s */)> Luminosity;
// differentiated with respect to rapidity of the system
typedef std::function<double (double /* s */, double /* y */)> Luminosity_y;
// for calculating fiducial cross section
typedef std::function<double (
    double /* s */, double /* ymin */, double /* ymax */
)> Luminosity_fid;

// Photon-photon luminosity in ultraperipheral collisions of particles with EPA
// spectra nA and nB with non-electromagnetic interactions neglected
Luminosity luminosity(
    Spectrum nA, Spectrum nB, Integrator = default_integrator(0)
);
// Same for identical particles.
Luminosity luminosity(Spectrum, Integrator = default_integrator(0));

Luminosity_y luminosity_y(Spectrum nA, Spectrum nB);
Luminosity_y luminosity_y(Spectrum);

Luminosity_fid luminosity_fid(
    Spectrum nA, Spectrum nB, Integrator = default_integrator(0)
);
Luminosity_fid luminosity_fid(Spectrum, Integrator = default_integrator(0));

struct Polarization {
  double parallel;
  double perpendicular;
};

// Photon-photon luminosity with non-electromagnetic interactions respected.
// `Polarization' are the weights to multiply components with parallel or
// perpendicular photons polarizations. Pass parallel = 1, perpendicular = 0 to
// compute only for parallel photon polarization, parallel = 1, perpendicular =
// 1 to compute for the sum of polarizations.
typedef std::function<double (double /* s */, Polarization)> Luminosity_b;
// differentiated with respect to rapidity of the system
typedef std::function<
  double (double /* s */, double /* y */, Polarization)
> Luminosity_b_y;
// for calculating fiducial cross section
typedef std::function<
  double (double /* s */, Polarization, double /* ymin */, double /* ymax */)
> Luminosity_b_fid;

// Photon-photon luminosity in ultraperipheral collisions of particles with EPA
// spectra nA and nB and the probability to survive upc_probability(b) where b
// is the impact parameter of the collision.
Luminosity_b
luminosity_b(
    Spectrum_b nA,
    Spectrum_b nB,
    std::function<double (double)> upc_probability,
    const std::function<Integrator (unsigned)>& = default_integrator
);

Luminosity_b
luminosity_b(
    Spectrum_b,
    std::function<double (double)> upc_probability,
    const std::function<Integrator (unsigned)>& = default_integrator
);

// differentiated with respect to rapidity of the system
Luminosity_b_y
luminosity_b_y(
    Spectrum_b nA,
    Spectrum_b nB,
    std::function<double (double)> upc_probability,
    const std::function<Integrator (unsigned)>& = default_integrator
);

Luminosity_b_y
luminosity_b_y(
    Spectrum_b,
    std::function<double (double)> upc_probability,
    const std::function<Integrator (unsigned)>& = default_integrator
);

// for calculating fiducial cross section
Luminosity_b_fid
luminosity_b_fid(
    Spectrum_b nA,
    Spectrum_b nB,
    std::function<double (double)> upc_probability,
    const std::function<Integrator (unsigned)>& = default_integrator
);
Luminosity_b_fid
luminosity_b_fid(
    Spectrum_b,
    std::function<double (double)> upc_probability,
    const std::function<Integrator (unsigned)>& = default_integrator
);

// Generic differential cross section. sqrt{s} is the invariant mass of the
// system.
typedef std::function<double (double /* s */)> XSection;

// Photon fusion cross section taking into account photons polarization.
// Returns two values: cross sections for the parallel and perpendicular
// photons polarizations respectively.
typedef std::function<Polarization (double /* s */)> XSection_b;

// Differential cross section for a process occuring in ultraperipheral
// collisions with non-electromagnetic interactions neglected
XSection
xsection(
    XSection,  // photon fusion cross section
    Luminosity // photon-photon luminosity
);
// same with non-electromagnetic interactions respected
XSection
xsection_b(
    XSection_b   xsection,
    Luminosity_b luminosity
);

// Fiducial cross section for the production of a pair of charged particles in
// ultraperipheral collision with non-electromagnetic interactions neglected
// and with the constraints on the phase space pT > pT_max, abs(eta) < eta_max,
// where pT and eta are the transverse momentum and pseudorapidity of each
// particle.
XSection
xsection_fid(
    // differential with respect to pT cross section for the pair production in
    // fusion of real photons
    std::function<double (double /* s */, double /* pT */)> xsection_pT,
    Luminosity_fid,
    double mass, // the mass of the charged particle
    double pT_min,
    double eta_max,
    Integrator = default_integrator(0)
);

// Fiducial cross section for the production of a pair of charged particles in
// ultraperipheral collision with non-electromagnetic interactions respected
// and with the constraints on the phase space pT > pT_max, abs(eta) < eta_max,
// where pT and eta are the transverse momentum and pseudorapidity of each
// particle.
//
// Note that the photon-photon differential cross section is usually divergent
// at pT = 0. GSL CQUAD integration method converges better in this case than
// GSL QAG.
XSection
xsection_b_fid(
    // differential with respect to pT cross section for the pair production in
    // fusion of real photons
    std::function<Polarization (double /* s */, double /* pT */)> xsection_pT,
    Luminosity_b_fid,
    double mass, // the mass of the charged particle
    double pT_min,
    double eta_max,
    Integrator = gsl_cquad_integrator(0u)
);

// Cross section for the production of a pair of fermions in photon-photon
// collisions (the Breit-Wheeler cross section) differentiated with respect to
// pT --- fermion transverse momentum. `mass' and `charge' are the fermion mass
// and charge.
std::function<double (double /* s */, double /* pT */)>
photons_to_fermions_pT(double mass, unsigned charge = 1);
// same for polarized photons
std::function<Polarization (double /* s */, double /* pT */)>
photons_to_fermions_pT_b(double mass, unsigned charge = 1);

}; // namespace epa
