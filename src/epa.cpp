#include <cmath>

#include <epa/epa.hpp>

namespace epa {

bool print_backtrace = true;

double default_absolute_error          = 0;
double default_relative_error          = 1e-3;
double default_error_step              = 1e-1;
size_t default_integration_limit       = 1000;
size_t default_cquad_integration_limit = 100;
gsl::integration::QAGMethod default_integration_method
  = gsl::integration::GAUSS41;

std::function<Integrator (unsigned)> default_integrator
  = static_cast<Integrator (*)(unsigned)>(qag_integrator);

std::function<Integrator_I (unsigned)> default_integrator_i
  = static_cast<Integrator_I (*)(unsigned)>(qag_integrator_i);

Integrator qag_integrator(
    double absolute_error,
    double relative_error,
    gsl::integration::QAGMethod method,
    std::shared_ptr<gsl::integration::QAGWorkspace> workspace
) {
  if (!workspace)
    workspace = std::make_shared<gsl::integration::QAGWorkspace>(
        default_integration_limit
    );
  return [=](const std::function<double (double)>& f, double a, double b)
         -> double {
    return gsl::integrate(
        f,
        a,
        b,
        absolute_error,
        relative_error,
        workspace->limit(),
        method,
        *workspace
    ).result;
  };
};

Integrator qag_integrator(const qag_integrator_keys& keys) {
  return qag_integrator(
      keys.absolute_error,
      keys.relative_error,
      keys.method,
      keys.workspace
  );
};

Integrator qag_integrator(unsigned level) {
  return qag_integrator(
      default_absolute_error,
      default_relative_error * pow(default_error_step, level)
  );
};

std::function<Integrator (unsigned)>
qag_integrator_generator(double relative_error, double error_step) {
  return [=](unsigned level) -> Integrator {
    return qag_integrator(0, relative_error * pow(error_step, level));
  };
};

Integrator cquad_integrator(
    double absolute_error,
    double relative_error,
    std::shared_ptr<gsl::integration::CQuadWorkspace> workspace
) {
  if (!workspace)
    workspace = std::make_shared<gsl::integration::CQuadWorkspace>(
        default_cquad_integration_limit
    );
  return [=](
      const std::function<double (double)>& f, double a, double b
  ) -> double {
    return gsl::integration::cquad(
        f,
        a,
        b,
        absolute_error,
        relative_error,
        *workspace
    ).result;
  };
};

Integrator cquad_integrator(const cquad_integrator_keys& keys) {
  return cquad_integrator(
      keys.absolute_error,
      keys.relative_error,
      keys.workspace
  );
};

Integrator cquad_integrator(unsigned level) {
  return cquad_integrator(
      default_absolute_error,
      default_relative_error * pow(default_error_step, level)
  );
};

Integrator_I qag_integrator_i(
    double relative_error,
    gsl::integration::QAGMethod method,
    std::shared_ptr<gsl::integration::QAGWorkspace> workspace
) {
  if (!workspace)
    workspace = std::make_shared<gsl::integration::QAGWorkspace>(
        default_integration_limit
    );
  return [=](
      const std::function<double (double)>& f, double a, double b, double I
  ) -> double {
    return gsl::integrate(
        f,
        a,
        b,
        relative_error * abs(I),
        relative_error,
        workspace->limit(),
        method,
        *workspace
    ).result;
  };
};

Integrator_I qag_integrator_i(const qag_integrator_keys& keys) {
  return qag_integrator_i(
      keys.relative_error,
      keys.method,
      keys.workspace
  );
};

Integrator_I qag_integrator_i(unsigned level) {
  return qag_integrator_i(default_relative_error * pow(default_error_step, level));
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
  if (!form_factor || form_factor->points->size() < 2)
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
    auto& last = form_factor->points->back();
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

      double qt_max = form_factor->points->back().first - wg2;
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
            if (right == form_factor->points->end()) break;
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
         (double rs, double y) -> double {
    EPA_TRY
      double E = 0.5 * rs;
      double x = exp(y);
      return E * nA(E * x) * nB(E / x);
    EPA_BACKTRACE(
        "lambda (rs, y) %e, %e\n  defined in epa::luminosity_y", rs, y
    );
  };
};

Luminosity_y luminosity_y(Spectrum n) {
  return luminosity_y(n, n);
};

Luminosity_fid luminosity_fid(Spectrum nA, Spectrum nB, Integrator integrate) {
  auto E = std::make_shared<double>();

  auto fx = [E, nA = std::move(nA), nB = std::move(nB)](double x) -> double {
    EPA_TRY
      double rx = sqrt(x);
      return nA(*E * rx) * nB(*E / rx) / x;
    EPA_BACKTRACE("lambda (x) %e; y = %e", x, 0.5 * log(x));
  };

  return [E, fx = std::move(fx), integrate = std::move(integrate)](
      double rs, double y_min, double y_max
  ) -> double {
    EPA_TRY
      *E  = 0.5 * rs;
      return 0.25 * rs * integrate(fx, exp(2 * y_min), exp(2 * y_max));
    EPA_BACKTRACE(
        "lambda (rs, y_min, y_max) %e, %e, %e\n  defined in epa::luminosity_fid",
        rs, y_min, y_max
    );
  };
};

Luminosity_fid luminosity_fid(Spectrum n, Integrator integrate) {
  return [l = luminosity_fid(n, n, std::move(integrate))](
      double rs, double y_min, double y_max
  ) -> double {
    return y_min == -y_max ? 2 * l(rs, y_min, 0) : l(rs, y_min, y_max);
  };
};

Luminosity luminosity(Spectrum nA, Spectrum nB, Integrator integrate) {
  return [l = luminosity_fid(nA, nB, std::move(integrate))](
      double rs
  ) -> double {
    return l(rs, -infinity, infinity);
  };
};

Luminosity luminosity(Spectrum n, Integrator integrate) {
  return [l = luminosity_fid(n, n, std::move(integrate))](double rs) -> double {
    return 2 * l(rs, -infinity, 0);
  };
};

Luminosity_y_b
luminosity_y_b(
    Spectrum_b nA,
    Spectrum_b nB,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  struct Env {
    double E;
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
    integrate = integrator(level + 2)
  ](double b2) -> double {
    EPA_TRY
      env->b2 = b2;
      return b2
             * nA(env->b1, env->E * env->rx)
             * nB(b2, env->E / env->rx)
             * integrate(fphi, 0, 2 * pi);
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
      double rs, double y, Polarization polarization
  ) -> double {
    EPA_TRY
      env->E = 0.5 * rs;
      env->rx = exp(y);
      env->polarization = polarization;
      return env->E * pi / sqr(env->rx) * integrate(fb1, 0, infinity);
    EPA_BACKTRACE(
        "lambda (rs, y, polarization) %e, %e, {%e, %e}\n"
        "  defined in luminosity_y_b",
        rs, y, polarization.parallel, polarization.perpendicular
    );
  };
};

Luminosity_y_b
luminosity_y_b(
    Spectrum_b n,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  return luminosity_y_b(n, n, std::move(upc), integrator);
};

Luminosity_fid_b
luminosity_fid_b(
    Spectrum_b nA,
    Spectrum_b nB,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  // luminosity_y_b is not used here because it would result in an unoptimal
  // order of integration
  struct Env {
    double E;
    double x_min;
    double x_max;
    double b1;
    double b2;
    Polarization polarization;
  };

  auto env = std::make_shared<Env>();

  auto fx = [env, nA = std::move(nA), nB = std::move(nB)](double x) -> double {
    EPA_TRY
      double rx = sqrt(x);
      return nA(env->b1, env->E * rx) * nB(env->b2, env->E / rx) / x;
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
    env,
    fx        = std::move(fx),
    fphi      = std::move(fphi),
    integrate = integrator(level + 2)
  ](double b2) -> double {
    EPA_TRY
      env->b2 = b2;
      return b2 * integrate(fx, env->x_min, env->x_max) * integrate(fphi, 0, 2 * pi);
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
        double rs,
        Polarization polarization,
        double y_min,
        double y_max
  ) -> double {
    EPA_TRY
      env->E             = 0.5 * rs;
      env->x_min         = exp(2 * y_min);
      env->x_max         = exp(2 * y_max);
      env->polarization  = polarization;
      return env->E * pi * integrate(fb1, 0, infinity);
    EPA_BACKTRACE(
        "lambda (rs, polarization, y_min, y_max) %e, {%e, %e}, %e, %e\n"
        "  defined in luminosity_fid_b",
        rs, polarization.parallel, polarization.perpendicular, y_min, y_max
    );
  };
};

Luminosity_fid_b
luminosity_fid_b(
    Spectrum_b n,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  return [l = luminosity_fid_b(n, n, std::move(upc), integrator, level)](
      double rs,
      Polarization polarization,
      double y_min,
      double y_max
  ) -> double {
    return y_min == -y_max
         ? 2 * l(rs, polarization, y_min, 0)
         : l(rs, polarization, y_min, y_max);
  };
};

Luminosity_b
luminosity_b(
    Spectrum_b nA,
    Spectrum_b nB,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  return [
    l = luminosity_fid_b(
            std::move(nA), std::move(nB), std::move(upc), integrator, level
        )
  ](double rs, Polarization polarization) -> double {
    return l(rs, polarization, -infinity, infinity);
  };
};

Luminosity_b
luminosity_b(
    Spectrum_b n,
    std::function<double (double)> upc,
    const std::function<Integrator (unsigned)>& integrator,
    unsigned level
) {
  return [l = luminosity_fid_b(n, n, std::move(upc), integrator, level)](
      double rs, Polarization polarization
  ) -> double {
    return 2 * l(rs, polarization, -infinity, 0);
  };
};

XSection
xsection(XSection xsection, Luminosity luminosity) {
  return [
    xsection   = std::move(xsection),
    luminosity = std::move(luminosity)
  ](double rs) -> double {
    return xsection(rs) * luminosity(rs);
  };
};

XSection
xsection_b(
    XSection_b   xsection,
    Luminosity_b luminosity
) {
  return [
    xsection   = std::move(xsection),
    luminosity = std::move(luminosity)
  ](double rs) -> double {
    return luminosity(rs, xsection(rs));
  };
};

static
XSection
xsection_fid_x(
    std::function<
      double (
        double /* sqrt(s) */,
        double /* pT */,
        double /* y_min */,
        double /* y_max */
      )
    > xl, // cross section * luminosity
    double           mass,
    double           pT_min,
    double           eta_max,
    double           w1_min,
    double           w1_max,
    double           w2_min,
    double           w2_max,
    Integrator       integrate
) {
  struct Env {
    double energy;
    double y_min;
    double y_max;
  };

  auto env = std::make_shared<Env>();

  double m2 = sqr(mass);
  double sinh_eta = sinh(eta_max);
  double cosh_eta = cosh(eta_max);
  double cosh2_eta = sqr(cosh_eta);
  double E_min = sqrt(w1_min * w2_min);
  double E_max = sqrt(w1_max * w2_max);

  auto fpT = [=, xl = std::move(xl)](double pT) -> double {
    EPA_TRY
      double pT2 = sqr(pT);
      double r = 1 - (pT2 + m2) / sqr(env->energy);
      if (r <= 0) return infinity;
      double y = log(
            pT / env->energy
          * (sinh_eta + sqrt(cosh2_eta + m2 / pT2)) / (1 + sqrt(r))
      );
      double y_min = std::max(-y, env->y_min);
      double y_max = std::min( y, env->y_max);
      if (y_min >= y_max) return 0;
      return xl(2 * env->energy, pT, y_min, y_max);
    EPA_BACKTRACE("lambda (pT) %e", pT);
   };

  return [=, fpT = std::move(fpT)](double rs) -> double {
    EPA_TRY
      double E = 0.5 * rs;
      if (E <= E_min || E >= E_max) return 0;

      // Transverse momentum integration limits (u <= pT <= v)
      double v = sqrt(sqr(E) - m2);
      // when pT < v / cosh_eta, the integration limit y computed in fpT above
      // is negative, and the integration domain is empty
      double u = std::max(pT_min, v / cosh_eta);
      if (u >= v) return 0;

      env->energy = E;
      env->y_min = log(std::max(w1_min / E, E / w2_max));
      env->y_max = log(std::min(w1_max / E, E / w2_min));

      double y = std::max(env->y_min, -env->y_max);
      if (y > 0) {
        // when y > eta_max, the integration limit y computed in fpT above is
        // negative, and the integration domain is empty
        if (y >= eta_max) return 0;

        double r = 1 - m2 / sqr(E) * (1 + sqr(sinh(y) / cosh_eta));
        if (r > 0) {
          r = sqrt(r);
          double u1 = 0.5 * E * (
              (1 + r) / cosh(y - eta_max) - (1 - r) / cosh(y + eta_max)
          );
          if (u1 > u) {
            u = u1;
            if (u >= v) return 0;
          };
        };
      };
      // when y < 0, u1 < E / cosh(eta_max) <= u

      return integrate(fpT, u, v);
    EPA_BACKTRACE("lambda (rs) %e\n  defined in xsection_fid_x", rs);
  };
};

XSection
xsection_fid(
    std::function<double (double, double)> xsection,
    Luminosity_fid luminosity,
    double mass,
    double pT_min,
    double eta_max,
    double w1_min,
    double w1_max,
    double w2_min,
    double w2_max,
    Integrator integrate
) {
  return xsection_fid_x(
      [xsection = std::move(xsection), luminosity = std::move(luminosity)]
      (double rs, double pT, double y_min, double y_max) -> double {
        return xsection(rs, pT)
             * (
                 y_min == -y_max
                 ? 2 * luminosity(rs, y_min, 0)
                 : luminosity(rs, y_min, y_max)
               );
      },
      mass,
      pT_min, eta_max,
      w1_min, w1_max, w2_min, w2_max,
      std::move(integrate)
  );
};

XSection
xsection_fid_b(
    std::function<Polarization (double, double)> xsection,
    Luminosity_fid_b luminosity,
    double mass,
    double pT_min,
    double eta_max,
    double w1_min,
    double w1_max,
    double w2_min,
    double w2_max,
    Integrator integrate
) {
  return xsection_fid_x(
      [xsection = std::move(xsection), luminosity = std::move(luminosity)]
      (double rs, double pT, double y_min, double y_max) -> double {
        return y_min == -y_max
             ? 2 * luminosity(rs, xsection(rs, pT), y_min, 0)
             : luminosity(rs, xsection(rs, pT), y_min, y_max);
      },
      mass,
      pT_min, eta_max,
      w1_min, w1_max, w2_min, w2_max,
      std::move(integrate)
  );
};

std::function<double (double)>
photons_to_fermions(double mass, double charge) {
  double c = 4 * pi * sqr(sqr(charge) * alpha) * barn;
  double m2 = sqr(mass);
  return [=](double rs) -> double {
    double s = sqr(rs);
    double x  = m2 / s;
    double x2 = sqr(x);
    double r  = sqrt(1 - 4 * x);
    return c / s * (
        (1 + 4 * x - 8 * x2) * log((1 + r) / (1 - r)) - (1 + 4 * x) * r
    );
  };
};

std::function<Polarization (double)>
photons_to_fermions_b(double mass, double charge) {
  double c = 4 * pi * sqr(sqr(charge) * alpha) * barn;
  double m2 = sqr(mass);
  return [=](double rs) -> Polarization {
    double s = sqr(rs);
    double x  = m2 / s;
    double x2 = sqr(x);
    double r  = sqrt(1 - 4 * x);
    double l  = log((1 + r) / (1 - r));
    double C  = c / s;
    return {
      C * ((1 + 4 * x - 12 * x2) * l - (1 + 6 * x) * r),
      C * ((1 + 4 * x -  4 * x2) * l - (1 + 2 * x) * r)
    };
  };
};

static
std::function<double (double, double)>
photons_to_fermions_pT_x(
    double mass, double charge, double polarization_factor
) {
  double c = 8 * pi * sqr(sqr(charge) * alpha) * barn;
  double m2 = sqr(mass);
  double m4 = sqr(m2);
  return [=](double rs, double pT) -> double {
    double s = sqr(rs);
    double pT2 = sqr(pT);
    double z = pT2 + m2;
    double sz = s * z;
    return c * pT / sz
         * (1 - 2 * (sqr(pT2) + polarization_factor * m4) / sz)
         / sqrt(1 - 4 * z / s);
  };
};

std::function<double (double, double)>
photons_to_fermions_pT(double mass, double charge) {
  return photons_to_fermions_pT_x(mass, charge, 1);
};

std::function<Polarization (double, double)>
photons_to_fermions_pT_b(double mass, double charge) {
  return [
    parallel      = photons_to_fermions_pT_x(mass, charge, 2),
    perpendicular = photons_to_fermions_pT_x(mass, charge, 0)
  ](double rs, double pT) -> Polarization {
    return { parallel(rs, pT), perpendicular(rs, pT) };
  };
};

}; // namespace epa
