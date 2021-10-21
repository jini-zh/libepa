#include <gsl/gsl_errno.h>

#include "gsl.hpp"

namespace gsl {

const char* Error::what() const throw() {
  return gsl_strerror(err_);
};

void init() {
  gsl_set_error_handler(
      [](const char* reason, const char* file, int line, int gsl_errno) {
         throw Error(gsl_errno);
      }
  );
};

namespace integration {

static double closure_trampoline(double x, void* data) {
  auto f = reinterpret_cast<std::function<double (double)>*>(data);
  return (*f)(x);
};

Workspace::Workspace(size_t limit):
  workspace(gsl_integration_workspace_alloc(limit))
{};

Workspace::Workspace(Workspace&& w) {
  workspace = w.workspace;
  limit_ = w.limit_;
  w.workspace = nullptr;
};

Workspace::~Workspace() {
  if (workspace) gsl_integration_workspace_free(workspace);
};

std::pair<double, double> qag(
    const std::function<double (double)>& f,
    double a,
    double b,
    double epsabs,
    double epsrel,
    size_t limit,
    QAGMethod method,
    const Workspace& workspace
) {
  gsl_function F;
  F.function = closure_trampoline;
  F.params = const_cast<std::function<double (double)>*>(&f);

  std::pair<double, double> result;
  gsl_integration_qag(
      &F,
      a,
      b,
      epsabs,
      epsrel,
      limit,
      method,
      workspace.get(),
      &result.first,
      &result.second
  );

  return result;
};

std::pair<double, double> qagiu(
    const std::function<double (double)>& f,
    double a,
    double epsabs,
    double epsrel,
    size_t limit,
    const Workspace& workspace
) {
  gsl_function F;
  F.function = closure_trampoline;
  F.params = const_cast<std::function<double (double)>*>(&f);

  std::pair<double, double> result;
  gsl_integration_qagiu(
      &F,
      a,
      epsabs,
      epsrel,
      limit,
      workspace.get(),
      &result.first,
      &result.second
  );

  return result;
};

std::pair<double, double> qagi(
    const std::function<double (double)>& f,
    double epsabs,
    double epsrel,
    size_t limit,
    const Workspace& workspace
) {
  gsl_function F;
  F.function = closure_trampoline;
  F.params = const_cast<std::function<double (double)>*>(&f);

  std::pair<double, double> result;
  gsl_integration_qagi(
      &F,
      epsabs,
      epsrel,
      limit,
      workspace.get(),
      &result.first,
      &result.second
  );

  return result;
};

}; // namespace integration

}; // namespace gsl
