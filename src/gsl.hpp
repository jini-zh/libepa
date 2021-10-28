#pragma once

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include <functional>
#include <utility>
#include <stdexcept>
#include <limits>

namespace gsl {

class Error: public std::exception {
  public:
    Error(int err): err_(err) {};

    const char* what() const throw();
    int err() const { return err_; };
  private:
    int err_;
};

const double infinity = std::numeric_limits<double>().infinity();

const auto bessel_J0 = gsl_sf_bessel_J0;
const auto bessel_J1 = gsl_sf_bessel_J1;
const auto bessel_Jn = gsl_sf_bessel_Jn;

const auto bessel_I0_scaled = gsl_sf_bessel_I0_scaled;
const auto bessel_In_scaled = gsl_sf_bessel_In_scaled;

// avoid underflow
inline double bessel_K0(double x) {
  if (x > 700) return 0;
  return gsl_sf_bessel_K0(x);
};

inline double bessel_K1(double x) {
  if (x > 700) return 0;
  return gsl_sf_bessel_K1(x);
};

void init();

namespace integration {

class QAGWorkspace {
  public:
    QAGWorkspace(size_t limit = 1000);
    QAGWorkspace(const QAGWorkspace&) = delete;
    QAGWorkspace(QAGWorkspace&&);
    ~QAGWorkspace();

    QAGWorkspace& operator=(const QAGWorkspace&) = delete;

    gsl_integration_workspace* get() const {
      return workspace;
    }

    size_t limit() const {
      return limit_;
    };

  private:
    gsl_integration_workspace* workspace;
    size_t limit_;
};

enum QAGMethod {
  GAUSS15 = GSL_INTEG_GAUSS15,
  GAUSS21 = GSL_INTEG_GAUSS21,
  GAUSS31 = GSL_INTEG_GAUSS31,
  GAUSS41 = GSL_INTEG_GAUSS41,
  GAUSS51 = GSL_INTEG_GAUSS51,
  GAUSS61 = GSL_INTEG_GAUSS61
};

struct QAGResult {
  double result;
  double abserr;
};

// Unified interface for quadrature integration (qag, qagiu, qagil, qagi). Pass
// the infinity constant to integrate to infinity.
QAGResult qag(
    const std::function<double (double)>& f,
    double a,
    double b,
    double epsabs,
    double epsrel,
    size_t limit,
    QAGMethod,
    const QAGWorkspace&
);

}; // namespace integration

const auto integrate = integration::qag;

}; // namespace gsl
