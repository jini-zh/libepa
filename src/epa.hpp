#pragma once

#include <memory>

#include "algorithms.hpp"
#include "gsl.hpp"

namespace epa {

const double pi              = M_PI;
const double alpha           = 7.2973525693e-3; // [1] (fine structure constant)
const double planck          = 0.1973269804;           // GeV fm [1]
const double barn            = sqr(planck) * 1e-2;     // (GeV * barn)^{-2}
const double fm              = 1. / (10 * sqrt(barn)); // (GeV fm)^{-1}
const double electron_charge = 1.602176634e-19;        // C [1]
const double light_speed     = 299792458;              // m/s [1]
const double amu
  = (1.66053906660e-27 * sqr(light_speed)) / electron_charge * 1e-9; // GeV, [1]
/*
 [1] E. Tiesinga, P. J. Mohr, D. B. Newell, B. N. Taylor.
     CODATA recommended values of the fundamental physical constants: 2018.
     Reviews of Modern Physics 93, 025010 (2020).
*/

// default parameters for GSL routines for quadrature integration
extern double default_absolute_error;
extern double default_relative_error;
extern double default_error_step;
extern size_t default_integration_limit;
extern gsl::integration::QAGMethod default_integration_method;

// Function: f, a, b -> integral of f from a to b
typedef std::function<
          double (const std::function<double (double)>&, double, double)
        > IntegratorAB;
// Function: f, a -> integral of f from a to infinity
typedef std::function<double (const std::function<double (double)>&, double)>
        IntegratorAInf;
// Function: f -> integral of f from a to infinity
typedef std::function<double (const std::function<double (double)>&)>
        IntegratorInf;

// GSL integrators  with default initialization
IntegratorAB gsl_integrator_ab(
    double absolute_error = default_absolute_error,
    double relative_error = default_relative_error,
    size_t limit          = default_integration_limit,
    gsl::integration::QAGMethod = gsl::integration::GAUSS41,
    std::shared_ptr<gsl::integration::Workspace> = nullptr
);
IntegratorAInf gsl_integrator_ainf(
    double absolute_error = default_absolute_error,
    double relative_error = default_relative_error,
    size_t limit          = default_integration_limit,
    std::shared_ptr<gsl::integration::Workspace> = nullptr
);
IntegratorInf gsl_integrator_inf(
    double absolute_error = default_absolute_error,
    double relative_error = default_relative_error,
    size_t limit          = default_integration_limit,
    std::shared_ptr<gsl::integration::Workspace> = nullptr
);

struct gsl_integrator_ab_keys {
  double absolute_error = default_absolute_error;
  double relative_error = default_relative_error;
  size_t limit          = default_integration_limit;
  gsl::integration::QAGMethod method = gsl::integration::GAUSS41;
  std::shared_ptr<gsl::integration::Workspace> workspace;
};

struct gsl_integrator_inf_keys {
  double absolute_error = default_absolute_error;
  double relative_error = default_relative_error;
  size_t limit          = default_integration_limit;
  std::shared_ptr<gsl::integration::Workspace> workspace;
};

// Helper functions --- with keyword parameters
IntegratorAB   gsl_integrator_ab  (const gsl_integrator_ab_keys&);
IntegratorAInf gsl_integrator_ainf(const gsl_integrator_inf_keys&);
IntegratorInf  gsl_integrator_inf (const gsl_integrator_inf_keys&);

// GSL integrator with relative_error = default_relative_error *
// default_error_step ** level
IntegratorAB   gsl_integrator_ab  (unsigned level = 0);
IntegratorAInf gsl_integrator_ainf(unsigned level = 0);
IntegratorInf  gsl_integrator_inf (unsigned level = 0);

// Default integrators
const std::function<IntegratorAB (unsigned)> default_integrator_ab
  = static_cast<IntegratorAB (*)(unsigned)>(gsl_integrator_ab);
const std::function<IntegratorAInf (unsigned)> default_integrator_ainf
  = static_cast<IntegratorAInf (*)(unsigned)>(gsl_integrator_ainf);
const std::function<IntegratorInf (unsigned)> default_integrator_inf
  = static_cast<IntegratorInf (*)(unsigned)>(gsl_integrator_inf);

// Integrators that take one extra parameter: the value of the integral
// calculated so far. It can be used to avoid extra work in calculations that
// don't need that much accuracy. These integrators are currently used in
// spectrum_b_function1d* functions to evaluate the integral of the form factor
// in the region after the one described by Function1d. 
typedef std::function<
          double (const std::function<double (double)>&, double, double, double)
        > IntegratorAB_I;
typedef std::function<
          double (const std::function<double (double)>&, double, double)
        > IntegratorAInf_I;
IntegratorAB_I gsl_integrator_ab_i(
    double relative_error = default_relative_error,
    size_t limit          = default_integration_limit,
    gsl::integration::QAGMethod = gsl::integration::GAUSS41,
    std::shared_ptr<gsl::integration::Workspace> = nullptr
);
IntegratorAInf_I gsl_integrator_ainf_i(
    double relative_error = default_relative_error,
    size_t limit          = default_integration_limit,
    std::shared_ptr<gsl::integration::Workspace> = nullptr
);
struct gsl_integrator_ab_i_keys {
  double relative_error = default_relative_error;
  size_t limit          = default_integration_limit;
  gsl::integration::QAGMethod method = gsl::integration::GAUSS41;
  std::shared_ptr<gsl::integration::Workspace> workspace;
};
struct gsl_integrator_ainf_i_keys {
  double relative_error = default_relative_error;
  size_t limit          = default_integration_limit;
  std::shared_ptr<gsl::integration::Workspace> workspace;
};
IntegratorAB_I   gsl_integrator_ab_i  (const gsl_integrator_ab_keys&);
IntegratorAInf_I gsl_integrator_ainf_i(const gsl_integrator_inf_keys&);
IntegratorAB_I   gsl_integrator_ab_i  (unsigned level = 0);
IntegratorAInf_I gsl_integrator_ainf_i(unsigned level = 0);
const std::function<IntegratorAB_I (unsigned)> default_integrator_ab_i
  = static_cast<IntegratorAB_I (*)(unsigned)>(gsl_integrator_ab_i);
const std::function<IntegratorAInf_I (unsigned)> default_integrator_ainf_i
  = static_cast<IntegratorAInf_I (*)(unsigned)>(gsl_integrator_ainf_i);

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
    IntegratorAInf = default_integrator_ainf(0)
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
    IntegratorAInf = default_integrator_ainf(0)
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
// integrate_ab parameters are:
//   the function to integrate
//   integration lower bound
//   integration upper bound
//   absolute value of already calculated part of the integral. It can be used
//     to speed up convergence by requiring that absolute error of numerical
//     integration is something like 1e-3 times the value of this parameter.
// integrate_ainf parameters are similar:
//   the function to integrate
//   integration lower bound; upper bound is infinity
//   absolute value of already calculated part of the integral
Spectrum_b
spectrum_b_function1d_g(
    unsigned z,
    double gamma,
    std::shared_ptr<Function1d> form_factor,
    FormFactor rest_form_factor = FormFactor(),
    Spectrum_b rest_spectrum    = Spectrum_b(),
    double b_max = 0,
    IntegratorAB_I   = default_integrator_ab_i(0),
    IntegratorAInf_I = default_integrator_ainf_i(0)
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
    IntegratorAB_I   = default_integrator_ab_i(0),
    IntegratorAInf_I = default_integrator_ainf_i(0)
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
    IntegratorAB_I   = default_integrator_ab_i(0),
    IntegratorAInf_I = default_integrator_ainf_i(0)
);

// Photon-photon luminosity with non-electromagnetic interactions neglected
typedef std::function<double (double /* s */)> Luminosity;
// same differentiated with respect to rapidity of the system
typedef std::function<double (double /* s */, double /* y */)> Luminosity_y;

// Photon-photon luminosity in ultraperipheral collisions of particles with EPA
// spectra nA and nB differentiated with respect to rapidity of
// the system.
Luminosity_y luminosity_y(Spectrum nA, Spectrum nB);
// Same for identical particles.
Luminosity_y luminosity_y(Spectrum);

// Photon-photon luminosity
Luminosity luminosity(
    Spectrum nA, Spectrum nB, IntegratorAInf = default_integrator_ainf(0)
);
Luminosity luminosity(Spectrum, IntegratorAInf);
// This function takes advantage of the symmetry of the system
Luminosity luminosity(Spectrum, IntegratorAB = default_integrator_ab(0));

// Photon-photon luminosity with non-electromagnetic interactions taken into account


}; // namespace epa
