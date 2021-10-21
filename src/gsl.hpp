#pragma once

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include <functional>
#include <utility>
#include <stdexcept>

namespace gsl {

class Error: public std::exception {
  public:
    Error(int err): err_(err) {};

    const char* what() const throw();
    int err() const { return err_; };
  private:
    int err_;
};

const auto bessel_J0 = gsl_sf_bessel_J0;
const auto bessel_J1 = gsl_sf_bessel_J1;
const auto bessel_Jn = gsl_sf_bessel_Jn;
const auto bessel_K0 = gsl_sf_bessel_K0;
const auto bessel_K1 = gsl_sf_bessel_K1;

void init();

namespace integration {

class Workspace {
  public:
    Workspace(size_t limit = 1000);
    Workspace(const Workspace&) = delete;
    Workspace(Workspace&&);
    ~Workspace();

    Workspace& operator=(const Workspace&) = delete;

    gsl_integration_workspace*
    get() const {
      return workspace;
    }

    size_t
    limit() const {
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

std::pair<double, double> qag(
    const std::function<double (double)>& f,
    double a,
    double b,
    double epsabs,
    double epsrel,
    size_t limit,
    QAGMethod,
    const Workspace&
);

std::pair<double, double> qagiu(
    const std::function<double (double)>& f,
    double a,
    double epsabs,
    double epsrel,
    size_t limit,
    const Workspace&
);

std::pair<double, double> qagi(
    const std::function<double (double)>& f,
    double epsabs,
    double epsrel,
    size_t limit,
    const Workspace&
);

}; // namespace integration

}; // namespace gsl
