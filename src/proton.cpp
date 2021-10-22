#include "proton.hpp"

namespace epa {

FormFactor proton_dipole_form_factor(double lambda2) {
  return [=](double q2) -> double {
    double tau = q2 / sqr(2 * proton_mass);
    return (1 + (proton_magnetic_moment - 1) * tau / (tau + 1))
         / sqr(1 + q2 / lambda2);
  };
};

Spectrum
proton_dipole_spectrum(double energy, double lambda2) {
  double c = alpha / pi;
  double b = sqr(2 * proton_mass) / lambda2;
  const double mu = proton_magnetic_moment;
  double l01 = 4 - 2 * (mu - 1) / b;
  double l10;
  double l11;
  {
    double c1 = (mu - 1) / pow(b - 1, 4);
    double c2 = (mu - 1) / (b - 1);
    l10 = c1 * (c2 * (1 + 3 * b) - 2);
    l11 = c1 * (c2 * 4 - 2 / b);
  };
  double a0 = -17. / 6
            + (mu - 1) * (2 * b * b - 7 * b + 11) / (3 * pow(b - 1, 3))
            + sqr(mu - 1) * (b * b - 8 * b - 17) / (6 * pow(b - 1, 4));
  double a1 = -7
            + (mu - 1) * (3 * b * b - 9 * b + 10) / pow(b - 1, 3)
            - sqr(mu - 1) * (b + 7) / pow(b - 1, 4);
  double a2 = -4
            + (mu - 1) * 2 * (b * b - 3 * b + 3) / pow(b - 1, 3)
            - 4 * sqr(mu - 1) / pow(b - 1, 4);
  double lg2 = lambda2 * sqr(energy / proton_mass);
  return [=](double w) -> double {
    double a = sqr(w) / lg2;
    double n1 = (1 + l01 * a) * log(1 + 1. / a)
              + (l10 + l11 * a) * log((a + b) / (a + 1));
    double n2 = (a0 + a1 * a + a2 * a * a) / sqr(a + 1);
    double n = n1 + n2;
//    if (abs(n / n1) < 1e3 * a * std::numeric_limits<double>().epsilon()) return 0;
    return c / w * n;
  };
};

Spectrum_b
proton_dipole_spectrum_b(double energy, double lambda2) {
  double c = alpha / sqr(pi);
  double m2 = sqr(2 * proton_mass);
  double gamma = energy / proton_mass;
  const double mu = proton_magnetic_moment;
  double x = lambda2 / m2;
  double k12 = (mu - 1) * sqr(x / (1 - x));
  double k11 = 1 + k12;
  double k00 = (1 - mu * x) / (1 - x) * lambda2 / 2;
  return [=](double b, double w) -> double {
    EPA_TRY
      double wg = w / gamma;
      double wg2 = sqr(wg);
      double rl = sqrt(lambda2 + wg2);
      double rm = sqrt(m2      + wg2);
      return c / w * sqr(
            wg * gsl::bessel_K1(b * wg)
          - k11 * rl * gsl::bessel_K1(b * rl)
          + k12 * rm * gsl::bessel_K1(b * rm)
          - k00 * b  * gsl::bessel_K0(b * rl)
      );
    EPA_BACKTRACE(
        "lambda (b, w) %e, %e\n  defined in proton_dipole_spectrum_b(%e, %e)",
        b, w, energy, lambda2
    );
  };
};

double pp_elastic_slope(double collision_energy) {
  const double B0 = 12;    // GeV^{-2}
  const double B1 = -0.22; // +/- 0.17 GeV^{-2}
  const double B2 = 0.037; // +/- 0.006 GeV^{-2}
  double l = log(collision_energy /* / 1 GeV */);
  return B0 + 2 * (B1 + 2 * B2 * l) * l;
};

std::function<double (double)>
pp_upc_probability(double collision_energy) {
  // see hep-ph/0608271
  double B = 2 * pp_elastic_slope(collision_energy);
  return [B](double b) -> double {
    return sqr(1 - exp(-sqr(b) / B));
  };
};

}; // namespace epa
