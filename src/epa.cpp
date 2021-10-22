#include <cmath>

#include "epa.hpp"

namespace epa {

bool print_backtrace = true;

double default_absolute_error     = 0;
double default_relative_error     = 1e-3;
double default_error_step         = 1e-1;
size_t default_integration_limit  = 1000;
gsl::integration::QAGMethod default_integration_method = gsl::integration::GAUSS41;

Integrator gsl_integrator(
    double absolute_error,
    double relative_error,
    size_t limit,
    gsl::integration::QAGMethod method,
    std::shared_ptr<gsl::integration::Workspace> workspace
) {
  if (!workspace)
    workspace = std::make_shared<gsl::integration::Workspace>(limit);
  return [=](const std::function<double (double)>& f, double a, double b)
         -> double {
    return gsl::integrate(
        f, a, b, absolute_error, relative_error, limit, method, *workspace
    ).first;
  };
};

Integrator gsl_integrator(const gsl_integrator_keys& keys) {
  return gsl_integrator(
      keys.absolute_error,
      keys.relative_error,
      keys.limit,
      keys.method,
      keys.workspace
  );
};

Integrator gsl_integrator(unsigned level) {
  return gsl_integrator(
      default_absolute_error,
      default_relative_error * pow(default_error_step, level)
  );
};

Integrator_I gsl_integrator_i(
    double relative_error,
    size_t limit,
    gsl::integration::QAGMethod method,
    std::shared_ptr<gsl::integration::Workspace> workspace
) {
  if (!workspace)
    workspace = std::make_shared<gsl::integration::Workspace>(limit);
  return [=](
      const std::function<double (double)>& f, double a, double b, double I
  ) -> double {
    return gsl::integrate(
        f,
        a,
        b,
        relative_error * abs(I),
        relative_error,
        limit,
        method,
        *workspace
    ).first;
  };
};

Integrator_I gsl_integrator_i(const gsl_integrator_keys& keys) {
  return gsl_integrator_i(
      keys.relative_error,
      keys.limit,
      keys.method,
      keys.workspace
  );
};

Integrator_I gsl_integrator_i(unsigned level) {
  return gsl_integrator_i(default_relative_error * pow(default_error_step, level));
};

FormFactor form_factor_monopole(double lambda2) {
  return [=](double q2) -> double {
    return 1. / (1 + q2 / lambda2);
  };
};

FormFactor form_factor_dipole(double lambda2) {
  return [=](double q2) -> double {
    return 1. / sqr(1 + q2 / lambda2);
  };
};

Spectrum
spectrum(
    unsigned Z, double gamma, FormFactor form_factor, Integrator integrate
) {
  auto wg2 = std::make_shared<double>();

  auto iqt = [=, F = std::move(form_factor)](double qt) -> double {
    EPA_TRY
      double q2 = sqr(qt) + *wg2;
      return qt * sqr(qt / q2 * F(q2));
    EPA_BACKTRACE("lambda (qt) %e", qt);
  };

  double c = 2 * sqr(Z) * alpha / pi;

  return [=, integrate = std::move(integrate)](double w) -> double {
    EPA_TRY
      *wg2 = sqr(w / gamma);
      return c / w * integrate(iqt, 0, infinity);
    EPA_BACKTRACE("lambda (w) %e\n  defined in epa::spectrum(%u, %e)", w, Z, gamma);
  };
};

Spectrum
spectrum_monopole(unsigned Z, double gamma, double lambda2) {
  double c = sqr(Z) * alpha / pi;
  return [=](double w) -> double {
    double x = sqr(w / gamma) / lambda2;
//    if (x > 1.4e4) return 0; // prevent roundoff errors
    return c / w * ((1 + 2 * x) * log(1 + 1. / x) - 2);
  };
};

Spectrum
spectrum_dipole(unsigned Z, double gamma, double lambda2) {
  double c = sqr(Z) * alpha / pi;
  return [=](double w) -> double {
    double x = sqr(w / gamma) / lambda2;
//    if (x > 200) return 0; // prevent roundoff errors
    return c / w * (
       (1 + 4 * x) * log(1 + 1. / x)
       - ((24 * x + 42) * x + 17) / (6 * sqr(x + 1))
    );
  };
};

Spectrum_b
spectrum_b(
    unsigned Z, double gamma, FormFactor form_factor, Integrator integrate
) {
  struct Env {
    double b;
    double wg2;
  };

  auto env = std::make_shared<Env>();

  auto iqt = [=, F = std::move(form_factor)](double qt) -> double {
    EPA_TRY
      double qt2 = sqr(qt);
      double q2  = qt2 + env->wg2;
      return qt2 / q2 * F(q2) * gsl::bessel_J1(env->b * qt);
    EPA_BACKTRACE("lambda (qt) %e", qt);
  };

  double c = alpha * sqr(Z / pi);

  return [=, integrate = std::move(integrate)](double b, double w) -> double {
    EPA_TRY
      env->b   = b;
      env->wg2 = sqr(w / gamma);
      return c / w * sqr(integrate(iqt, 0, infinity));
    EPA_BACKTRACE(
        "lambda (b, w) %e, %e\n  defined in epa::spectrum_b(%u, %e)",
        b, w, Z, gamma
    );
  };
};

Spectrum_b
spectrum_b_point(unsigned Z, double gamma) {
  double c = alpha * sqr(Z / pi / gamma);
  return [=](double b, double w) -> double {
    EPA_TRY
      return c * w * sqr(gsl::bessel_K1(b * w / gamma));
    EPA_BACKTRACE(
        "lambda (b, w) %e, %e\n  defined in epa::spectrum_b_point(%u, %e)",
        b, w, Z, gamma
    );
  };
};

// x * K1(x) - 1 for small x
// x K1(x) - 1 = (x/2)^2 (2 ln(x/2) + 2γ - 1) + (x/2)^4 (ln(x/2) + γ - 5/4)
// where γ is the Euler constant
static double xk1_1(double x) {
  const double euler = 0.5772156649015329;
  double le = log(0.5 * x) + euler;
  double x2 = x * x;
  return x2 * (0.5 * (le - 0.5) + (0.0625 * x2 * (le - 1.25)));
};

Spectrum_b
spectrum_b_monopole(unsigned Z, double gamma, double lambda2) {
  double c = alpha * sqr(Z / pi);
  return [=](double b, double w) -> double {
    EPA_TRY
      double wg = w / gamma;
      double u = b * wg;
      double a = lambda2 / sqr(wg);
      double d;
      if (a < 1e-6)
        d = 0.5 * b * lambda2 * gsl::bessel_K0(u);
      else {
        double v = sqrt(b * b * lambda2 + sqr(u));
        if (v < 1e-2)
          d = xk1_1(u) - xk1_1(v);
        else
          d = u * gsl::bessel_K1(u) - v * gsl::bessel_K1(v);
        d /= b;
      };
      return c / w * sqr(d);
    EPA_BACKTRACE(
        "lambda (b, w) %e, %e\n  defined in epa::spectrum_b_monopole(%u, %e, %e)",
        b, w, Z, gamma, lambda2
    );
  };
};

Spectrum_b
spectrum_b_dipole(unsigned Z, double gamma, double lambda2) {
  double c = alpha * sqr(Z / pi);
  return [=](double b, double w) -> double {
    EPA_TRY
      double wg = w / gamma;
      double r = sqrt(lambda2 + sqr(wg));
      return c / w * sqr(
            wg * gsl::bessel_K1(b*wg)
          - r * gsl::bessel_K1(b*r)
          - 0.5 * b * lambda2 * gsl::bessel_K0(b*r)
      );
    EPA_BACKTRACE(
        "lambda (b, w) %e, %e\n  defined in epa::spectrum_b_dipole(%u, %e, %e)",
        b, w, Z, gamma, lambda2
    );
  };
};

static
Spectrum_b
spectrum_b_function1d_x(
    unsigned Z,
    double gamma,
    std::shared_ptr<Function1d> form_factor,
    FormFactor rest_form_factor,
    Spectrum_b rest_spectrum,
    double b_max,
    std::function<double (double, double, double)>&& integral_qt_max,
    Integrator_I integrate
) {
  if (!form_factor || form_factor->points.size() < 2)
    if (rest_spectrum)
      return rest_spectrum;
    else if (rest_form_factor)
      return spectrum_b(Z, gamma, rest_form_factor);
    else if (!form_factor)
      throw std::invalid_argument(
          "epa::spectrum_b_function1d_x: null form factor"
      );
    else
      throw std::invalid_argument(
          "epa::spectrum_b_function1d_x: too few points in the form factor"
      );

  struct Env {
    double b;
    double wg2;
  };

  auto env = std::make_shared<Env>();

  double c = alpha * sqr(Z / pi);

  double norm;
  if (rest_form_factor) {
    auto& last = form_factor->points.back();
    norm = last.second / rest_form_factor(last.first);
  };

  Spectrum_b n0;
  if (b_max > 0) n0 = spectrum_b_point(Z, gamma);

  auto rqt = [env, ff = std::move(rest_form_factor)](double qt) -> double {
    EPA_TRY
      double qt2 = sqr(qt);
      double q2  = qt2 + env->wg2;
      return qt2 / q2 * ff(q2) * gsl::bessel_J1(env->b * qt);
    EPA_BACKTRACE("lambda (qt) %e", qt);
  };

  return [
    =,
    integral_qt_max = std::move(integral_qt_max),
    integrate       = std::move(integrate),
    n0              = std::move(n0)
  ](double b, double w) -> double {
    EPA_TRY
      if (b_max > 0 && b > b_max) return n0(b, w);
      double wg2 = sqr(w / gamma);
      env->b   = b;
      env->wg2 = wg2;

      double qt_max = form_factor->points.back().first - wg2;
      qt_max = qt_max > 0 ? sqrt(qt_max) : 0;

      double I = qt_max == 0 ? 0 : integral_qt_max(b, wg2, qt_max);

      double C = c / w;
      if (rest_form_factor)
        I += norm * (
            rest_spectrum
            ? sqrt(rest_spectrum(b, w) / C) - integrate(rqt, 0, qt_max, I)
            : integrate(rqt, qt_max, infinity, I)
        );
      return C * sqr(I);
    EPA_BACKTRACE(
        "lambda (b, w) %e, %e\n"
        "  defined in epa::spectrum_b_function1d_x(%u, %e, %p, %s, %s, %e)",
        b, w,
        Z, gamma, form_factor.get(),
        rest_form_factor ? "<rest_form_factor>" : "nullptr",
        rest_spectrum    ? "<rest_spectrum>"    : "nullptr",
        b_max
    );
  };
};

Spectrum_b
spectrum_b_function1d_g(
    unsigned Z,
    double gamma,
    std::shared_ptr<Function1d> form_factor,
    FormFactor rest_form_factor,
    Spectrum_b rest_spectrum,
    double b_max,
    Integrator_I integrate
) {
  struct Env {
    double b;
    double wg2;
  };

  auto env = std::make_shared<Env>();

  auto fqt = [env, form_factor](double qt) -> double {
    double qt2 = sqr(qt);
    double q2  = qt2 + env->wg2;
    return qt2 / q2 * (*form_factor)(q2) * gsl::bessel_J1(env->b * qt);
  };

  return spectrum_b_function1d_x(
      Z,
      gamma,
      form_factor,
      rest_form_factor,
      rest_spectrum,
      b_max,
      [env, fqt = std::move(fqt), integrate]
      (double b, double wg2, double qt_max) -> double {
        EPA_TRY
          env->b = b;
          env->wg2 = wg2;
          return integrate(fqt, 0, qt_max, 0);
        EPA_BACKTRACE(
            "lambda (b, sqr(w/gamma), qt_max) %e, %e, %e\n"
            "  defined in epa::spectrum_b_function1d_g",
            b, wg2, qt_max
        );
      },
      integrate
  );
};

Spectrum_b
spectrum_b_function1d_s(
    unsigned Z,
    double gamma,
    std::shared_ptr<Function1d> form_factor,
    FormFactor rest_form_factor,
    Spectrum_b rest_spectrum,
    double b_max,
    Integrator_I integrate
) {
  struct Env {
    double b;
    double wg2;
  };

  auto env = std::make_shared<Env>();

  auto fj1 = [env](double qt) -> double {
    EPA_TRY
      return gsl::bessel_J1(env->b * qt) / (sqr(qt) + env->wg2);
    EPA_BACKTRACE("lambda (qt) %e", qt);
  };

  return spectrum_b_function1d_x(
      Z,
      gamma,
      form_factor,
      rest_form_factor,
      rest_spectrum,
      b_max,
      [form_factor, env, fj1 = std::move(fj1), integrate]
      (double b, double wg2, double qt_max) -> double {
        EPA_TRY
          env->b   = b;
          env->wg2 = wg2;
          std::vector<std::pair<double, double>>::const_iterator left
            = form_factor->locate(wg2);
          double start = left->first < wg2 ? 0 : sqrt(left->first - wg2);
          double I = 0;
          while (true) {
            auto right = left + 1;
            if (right == form_factor->points.end()) break;
            double end = sqrt(right->first - wg2);
            double A = (right->second - left->second)
                     / sqrt(right->first - left->first);
            double B = left->second - A * left->first;
            I += A / b
                 * ( sqr(end)   * gsl::bessel_Jn(2, b * end)
                   - sqr(start) * gsl::bessel_Jn(2, b * start)
                   )
               + B / b
                 * (gsl::bessel_J0(b * start) - gsl::bessel_J0(b * end))
               - B * wg2 * integrate(fj1, start, end, 0);
            left = right;
            start = end;
          };
          return I;
        EPA_BACKTRACE(
            "lambda (b, sqr(w/gamma), qt_max) %e, %e, %e\n"
            "  defined in epa::spectrum_b_function1d_s",
            b, wg2, qt_max
        );
      },
      integrate
  );
};

Spectrum_b
spectrum_b_function1d(
    unsigned Z,
    double gamma,
    std::shared_ptr<Function1d> form_factor,
    FormFactor rest_form_factor,
    Spectrum_b rest_spectrum,
    double b_max,
    Integrator_I integrate
) {
  auto global = spectrum_b_function1d_g(
      Z,
      gamma,
      form_factor,
      rest_form_factor,
      rest_spectrum,
      b_max,
      integrate
  );
  auto segmented = spectrum_b_function1d_s(
      Z,
      gamma,
      form_factor,
      rest_form_factor,
      rest_spectrum,
      b_max,
      integrate
  );

  return [global = std::move(global), segmented = std::move(segmented)]
         (double b, double w) -> double {
    // stub
    return global(b, w);
  };
};

Luminosity_y luminosity_y(Spectrum nA, Spectrum nB) {
  return [nA = std::move(nA), nB = std::move(nB)]
         (double s, double y) -> double {
    EPA_TRY
      s = 0.5 * sqrt(s);
      y = exp(y);
      return 0.25 * nA(s * y) * nB(s / y);
    EPA_BACKTRACE(
        "lambda (s, y) %e, %e\n  defined in epa::luminosity_y", s, y
    );
  };
};

Luminosity_y luminosity_y(Spectrum n) {
  return luminosity_y(n, n);
};

Luminosity_fid luminosity_fid(Spectrum nA, Spectrum nB, Integrator integrate) {
  auto rs = std::make_shared<double>();
  auto fx = [rs, nA = std::move(nA), nB = std::move(nB)](double x) -> double {
    EPA_TRY
      double rx = sqrt(x);
      return nA(*rs * rx) * nB(*rs / rx) / x;
    EPA_BACKTRACE("lambda (x) %e; y = %e", x, 0.5 * log(x));
  };

  return [rs, fx = std::move(fx), integrate = std::move(integrate)](
      double s, double ymin, double ymax
  ) -> double {
    EPA_TRY
      *rs = sqrt(s) / 2;
      return 0.125 * integrate(fx, exp(2 * ymin), exp(2 * ymax));
    EPA_BACKTRACE(
        "lambda (s, ymin, ymax) %e, %e, %e\n  defined in epa::luminosity_fid",
        s, ymin, ymax
    );
  };
};

Luminosity_fid luminosity_fid(Spectrum n, Integrator integrate) {
  return [l = luminosity_fid(n, n, std::move(integrate))](
      double s, double ymin, double ymax
  ) -> double {
    return ymin == -ymax ? 2 * l(s, ymin, 0) : l(s, ymin, ymax);
  };
};

Luminosity luminosity(Spectrum nA, Spectrum nB, Integrator integrate) {
  return [l = luminosity_fid(nA, nB, std::move(integrate))](
      double s
  ) -> double {
    return l(s, -infinity, infinity);
  };
};

Luminosity luminosity(Spectrum n, Integrator integrate) {
  return [l = luminosity_fid(n, n, std::move(integrate))](double s) -> double {
    return 2 * l(s, -infinity, 0);
  };
};

Luminosity_b_y
luminosity_b_y(
    Spectrum_b nA,
    Spectrum_b nB,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator
) {
  struct Env {
    double rs;
    double rx;
    double b1;
    double b2;
    Polarization polarization;
  };

  auto env = std::make_shared<Env>();

  auto fphi = [env, upc = std::move(upc)](double phi) -> double {
    EPA_TRY
      double c = cos(phi);
      double s = sin(phi);
      return upc(sqrt(sqr(env->b1) + sqr(env->b2) - 2 * env->b1 * env->b2 * c))
           * (
               env->polarization.parallel * sqr(c)
             + env->polarization.perpendicular * sqr(s)
             );
    EPA_BACKTRACE("lambda (phi) %e", phi);
  };

  auto fb2 = [
    env,
    nA   = std::move(nA),
    nB   = std::move(nB),
    fphi = std::move(fphi),
    integrate = integrator(2)
  ](double b2) -> double {
    EPA_TRY
      env->b2 = b2;
      return b2
             * nA(env->b1, env->rs * env->rx)
             * nB(b2, env->rs / env->rx)
             * integrate(fphi, 0, 2 * pi);
    EPA_BACKTRACE("lambda (b2) %e", b2);
  };

  auto fb1 = [env, fb2 = std::move(fb2), integrate = integrator(1)](
      double b1
  ) -> double {
    EPA_TRY
      env->b1 = b1;
      return b1 * integrate(fb2, 0, infinity);
    EPA_BACKTRACE("lambda (b1) %e", b1);
  };

  return [env, fb1 = std::move(fb1), integrate = integrator(0)](
      double s, double y, Polarization polarization
  ) -> double {
    EPA_TRY
      env->rs = 0.5 * sqrt(s);
      env->rx = exp(y);
      env->polarization = polarization;
      return 0.25 * pi / sqr(env->rx) * integrate(fb1, 0, infinity);
    EPA_BACKTRACE(
        "lambda (s, y, polarization) %e, %e, {%e, %e}\n"
        "  defined in luminosity_b_y",
        s, y, polarization.parallel, polarization.perpendicular
    );
  };
};

Luminosity_b_y
luminosity_b_y(
    Spectrum_b n,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator
) {
  return luminosity_b_y(n, n, std::move(upc), integrator);
};

Luminosity_b_fid
luminosity_b_fid(
    Spectrum_b nA,
    Spectrum_b nB,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator
) {
  // luminosity_b_y is not used here because it would result in an unoptimal
  // order of integration
  struct Env {
    double rs;
    double xmin;
    double xmax;
    double b1;
    double b2;
    Polarization polarization;
  };

  auto env = std::make_shared<Env>();

  auto fx = [env, nA = std::move(nA), nB = std::move(nB)](double x) -> double {
    EPA_TRY
      double rx = sqrt(x);
      return nA(env->b1, env->rs * rx) * nB(env->b2, env->rs / rx) / x;
    EPA_BACKTRACE("lambda (x) %e; y = %e", x, 0.5 * log(x));
  };

  auto fphi = [env, upc = std::move(upc)](double phi) -> double {
    EPA_TRY
      double c = cos(phi);
      double s = sin(phi);
      return upc(sqrt(sqr(env->b1) + sqr(env->b2) - 2 * env->b1 * env->b2 * c))
           * (
               env->polarization.parallel * sqr(c)
             + env->polarization.perpendicular * sqr(s)
             );
    EPA_BACKTRACE("lambda (phi) %e", phi);
  };

  auto fb2 = [
    env, fx = std::move(fx), fphi = std::move(fphi), integrate = integrator(2)
  ](double b2) -> double {
    EPA_TRY
      env->b2 = b2;
      return b2 * integrate(fx, env->xmin, env->xmax) * integrate(fphi, 0, 2 * pi);
    EPA_BACKTRACE("lambda (b2) %e", b2);
  };

  auto fb1 = [env, fb2 = std::move(fb2), integrate = integrator(1)](
      double b1
  ) -> double {
    EPA_TRY
      env->b1 = b1;
      return b1 * integrate(fb2, 0, infinity);
    EPA_BACKTRACE("lambda (b1) %e", b1);
  };

  return [env, fb1 = std::move(fb1), integrate = integrator(0)](
        double s,
        Polarization polarization,
        double ymin,
        double ymax
  ) -> double {
    EPA_TRY
      env->rs            = 0.5 * sqrt(s);
      env->xmin          = exp(2 * ymin);
      env->xmax          = exp(2 * ymax);
      env->polarization  = polarization;
      return 0.25 * pi * integrate(fb1, 0, infinity);
    EPA_BACKTRACE(
        "lambda (s, polarization, ymin, ymax) %e, {%e, %e}, %e, %e\n"
        "  defined in luminosity_b_fid",
        s, polarization.parallel, polarization.perpendicular, ymin, ymax
    );
  };
};

Luminosity_b_fid
luminosity_b_fid(
    Spectrum_b n,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator
) {
  return [l = luminosity_b_fid(n, n, std::move(upc), integrator)](
      double s,
      Polarization polarization,
      double ymin,
      double ymax
  ) -> double {
    return ymin == -ymax
         ? 2 * l(s, polarization, ymin, 0)
         : l(s, polarization, ymin, ymax);
  };
};

Luminosity_b
luminosity_b(
    Spectrum_b nA,
    Spectrum_b nB,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator
) {
  return [
    l = luminosity_b_fid(
            std::move(nA),
            std::move(nB),
            std::move(upc),
            integrator
        )
  ](double s, Polarization polarization) -> double {
    return l(s, polarization, -infinity, infinity);
  };
};

Luminosity_b
luminosity_b(
    Spectrum_b n,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator
) {
  return [l = luminosity_b_fid(n, n, std::move(upc), integrator)](
      double s, Polarization polarization
  ) -> double {
    return 2 * l(s, polarization, -infinity, 0);
  };
};

}; // namespace epa
