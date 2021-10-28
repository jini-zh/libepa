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

Luminosity_y
pp_luminosity_y(double collision_energy) {
  return luminosity_y(proton_dipole_spectrum(collision_energy / 2));
};

Luminosity
pp_luminosity(double collision_energy, Integrator integrate) {
  return luminosity(proton_dipole_spectrum(collision_energy / 2), integrate);
};

Luminosity_fid
pp_luminosity_fid(double collision_energy, Integrator integrate) {
  return luminosity_fid(
      proton_dipole_spectrum(collision_energy / 2), integrate
  );
};

static double ppx_luminosity_internal(
    double b1, double b2, double B, double psum, double pdifference
) {
  double e = exp(-0.5 * sqr(b1 - b2) / B);
  double z = b1 * b2 / B;
  double i01 = psum * gsl::bessel_I0_scaled(z);
  double i02 = psum * gsl::bessel_I0_scaled(2 * z);
  double i21, i22;
  if (pdifference != 0) {
    i21 = pdifference * gsl::bessel_In_scaled(2, z);
    i22 = pdifference * gsl::bessel_In_scaled(2, 2 * z);
  } else {
    i21 = 0;
    i22 = 0;
  };
  return -2 * e * (i01 + i21) + e * e * (i02 + i22);
};

Luminosity_y_b
ppx_luminosity_y_b(
    Spectrum_b n,
    double B,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  struct Env {
    double w1;
    double w2;
    double b1;
    double psum;
    double pdifference;
  };

  auto env = std::make_shared<Env>();

  auto fb2 = [env, B, n = std::move(n)](double b2) -> double {
    EPA_TRY
      return b2 * n(env->b1, env->w1) * n(b2, env->w2)
           * (env->psum
              + ppx_luminosity_internal(
                  env->b1, b2, B, env->psum, env->pdifference
                )
             );
    EPA_BACKTRACE("lambda (b2) %e", b2);
  };

  auto fb1 = [env, fb2 = std::move(fb2), integrate = integrator(level + 1)](
      double b1
  ) -> double {
    EPA_TRY
      env->b1 = b1;
      return b1 * integrate(fb2, 0, infinity);
    EPA_BACKTRACE("lambda (b1) %e", b1);
  };

  return [env, fb1 = std::move(fb1), integrate = integrator(level)](
      double s, double y, Polarization polarization
  ) -> double {
    EPA_TRY
      double rs = sqrt(s);
      double rx = exp(y);
      env->w1          = rs * rx;
      env->w2          = rs / rx;
      env->psum        = polarization.parallel + polarization.perpendicular;
      env->pdifference = polarization.parallel - polarization.perpendicular;
      return 0.5 * sqr(pi) * integrate(fb1, 0, infinity);
    EPA_BACKTRACE(
        "lambda (s, y, polarization) %e, %e, {%e, %e}\n"
        "  defined in ppx_luminosity_y",
        s, y, polarization.parallel, polarization.perpendicular
    );
  };
};

Luminosity_fid_b
ppx_luminosity_fid_b(
    Spectrum n,
    Spectrum_b n_b,
    double B,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  struct Env {
    double rs;
    double xmin;
    double xmax;
    double b1;
    double b2;
    double psum;
    double pdifference;
  };

  auto env = std::make_shared<Env>();

  auto fx = [env, n_b = std::move(n_b)](double x) -> double {
    EPA_TRY
      double rx = sqrt(x);
      return n_b(env->b1, env->rs * rx) * n_b(env->b2, env->rs / rx) / x;
    EPA_BACKTRACE("lambda (x) %e; y = %e", x, 0.5 * log(x));
  };

  auto fb2 = [
    env,
    B,
    fx        = std::move(fx),
    integrate = integrator(level + 2),
    one       = n ? 1 : 0
  ](double b2) -> double {
    EPA_TRY
      env->b2 = b2;
      return b2 * integrate(fx, env->xmin, env->xmax)
           * ppx_luminosity_internal(
               env->b1, b2, B, env->psum, env->pdifference
             );
    EPA_BACKTRACE("lambda (b2) %e", b2);
  };

  auto fb1 = [env, fb2 = std::move(fb2), integrate = integrator(level + 1)](
      double b1
  ) -> double {
    EPA_TRY
      env->b1 = b1;
      return b1 * integrate(fb2, 0, infinity);
    EPA_BACKTRACE("lambda (b1) %e", b1);
  };

  return [
    env,
    l   = n ? luminosity_fid(std::move(n), integrator(0)) : Luminosity_fid(),
    fb1 = std::move(fb1),
    integrate = integrator(level)
  ](
      double s, Polarization polarization, double ymin, double ymax
  ) -> double {
    EPA_TRY
      env->rs          = 0.5 * sqrt(s);
      env->xmin        = exp(2 * ymin);
      env->xmax        = exp(2 * ymax);
      env->psum        = polarization.parallel + polarization.perpendicular;
      env->pdifference = polarization.parallel - polarization.perpendicular;
      return (l ? 0.5 * env->psum * l(s, ymin, ymax) : 0)
             + 0.25 * sqr(pi) * integrate(fb1, 0, infinity);
    EPA_BACKTRACE(
        "lambda (s, polarization, ymin, ymax) %e, {%e, %e}, %e, %e\n"
        "  defined in ppx_luminosity_fid",
        s, polarization.parallel, polarization.perpendicular, ymin, ymax
    );
  };
};

Luminosity_b ppx_luminosity_b(
    Spectrum n,
    Spectrum_b n_b,
    double B,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  return [
    l = ppx_luminosity_fid_b(std::move(n), std::move(n_b), B, integrator, level)
  ](double s, Polarization polarization) -> double {
    return 2 * l(s, polarization, -infinity, 0);
  };
};

Luminosity_y_b pp_luminosity_y_b(
    double collision_energy,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  return ppx_luminosity_y_b(
      proton_dipole_spectrum_b(0.5 * collision_energy),
      pp_elastic_slope(collision_energy),
      integrator,
      level
  );
};

Luminosity_fid_b pp_luminosity_fid_b(
    double collision_energy,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  double energy = 0.5 * collision_energy;
  return ppx_luminosity_fid_b(
      proton_dipole_spectrum(energy),
      proton_dipole_spectrum_b(energy),
      pp_elastic_slope(collision_energy),
      integrator,
      level
  );
};

Luminosity_b pp_luminosity_b(
    double collision_energy,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  double energy = 0.5 * collision_energy;
  return ppx_luminosity_b(
      proton_dipole_spectrum(energy),
      proton_dipole_spectrum_b(energy),
      pp_elastic_slope(collision_energy),
      integrator,
      level
  );
};

}; // namespace epa
