#include <epa/epa.hpp>

#include "ffi.hpp"

using namespace epa;
using namespace ffi;

extern "C" void epa_init() {
  gsl::init();
};

#define defconst(type, name) extern "C" type epa_ ## name() { return name; }
#define defvar(type, name) \
  extern "C" type epa_get_ ## name() { return name; }; \
  extern "C" void epa_set_ ## name(type value_) { name = value_; }

#include "epa_vars.h"

#undef defvar
#undef defconst

template <typename Workspace>
static inline
std::shared_ptr<Workspace>*
epa_make_integration_workspace(size_t limit) {
  Workspace* w = nullptr;
  try {
    w = new Workspace(limit);
    return new std::shared_ptr<Workspace>(w);
  } catch (std::exception&) {
    if (w) delete w;
    set_error();
    return nullptr;
  };
};

extern "C"
std::shared_ptr<gsl::integration::QAGWorkspace>*
epa_make_qag_integration_workspace(size_t limit) {
  return epa_make_integration_workspace<gsl::integration::QAGWorkspace>(limit);
};

extern "C"
void
epa_destroy_qag_integration_workspace(
    std::shared_ptr<gsl::integration::QAGWorkspace>* workspace
) {
  delete workspace;
};

extern "C"
Function*
epa_qag_integrator(
    double absolute_error,
    double relative_error,
    gsl::integration::QAGMethod method,
    std::shared_ptr<gsl::integration::QAGWorkspace>* workspace
) {
  try {
    return lift(
        qag_integrator(
          absolute_error, relative_error, method, workspace ? *workspace : nullptr
        )
    );
  } FFI_CATCH;
};

extern "C"
std::shared_ptr<gsl::integration::CQuadWorkspace>*
epa_make_cquad_integration_workspace(size_t limit) {
  return epa_make_integration_workspace<gsl::integration::CQuadWorkspace>(
      limit
  );
};

extern "C"
void
epa_destroy_cquad_integration_workspace(
    std::shared_ptr<gsl::integration::CQuadWorkspace>* workspace
) {
  delete workspace;
};

extern "C"
Function*
epa_cquad_integrator(
    double absolute_error,
    double relative_error,
    std::shared_ptr<gsl::integration::CQuadWorkspace>* workspace
) {
  try {
    return lift(
        cquad_integrator(
          absolute_error, relative_error, workspace ? *workspace : nullptr
        )
    );
  } FFI_CATCH;
};


extern "C" int epa_get_default_integration_method() {
  return default_integration_method;
};

extern "C" void epa_set_default_integration_method(int method) {
  default_integration_method = static_cast<gsl::integration::QAGMethod>(method);
};

extern "C" Function* epa_form_factor_monopole(double lambda2) {
  try {
    return lift(form_factor_monopole(lambda2));
  } FFI_CATCH;
};

extern "C" Function* epa_form_factor_dipole(double lambda2) {
  try {
    return lift(form_factor_dipole(lambda2));
  } FFI_CATCH;
};

extern "C"
Function*
epa_spectrum(
    unsigned Z, double gamma, Function* form_factor, Function* integrator
) {
  try {
    return lift(
        spectrum(
          Z,
          gamma,
          lower<FormFactor>(form_factor),
          integrator ? lower<Integrator>(integrator) : default_integrator(0)
        )
    );
  } FFI_CATCH;
};

extern "C"
Function*
epa_spectrum_monopole(unsigned Z, double gamma, double lambda2) {
  try {
    return lift(spectrum_monopole(Z, gamma, lambda2));
  } FFI_CATCH;
};

extern "C"
Function*
epa_spectrum_dipole(unsigned Z, double gamma, double lambda2) {
  try {
    return lift(spectrum_dipole(Z, gamma, lambda2));
  } FFI_CATCH;
};

extern "C"
Function* epa_spectrum_b(
    unsigned Z, double gamma, Function* form_factor, Function* integrator
) {
  try {
    return lift(
        spectrum_b(
          Z,
          gamma,
          lower<FormFactor>(form_factor),
          integrator ? lower<Integrator>(integrator) : default_integrator(0)
        )
    );
  } FFI_CATCH;
};

extern "C"
Function*
epa_spectrum_b_point(unsigned Z, double gamma) {
  try {
    return lift(spectrum_b_point(Z, gamma));
  } FFI_CATCH;
};

extern "C"
Function*
epa_spectrum_b_monopole(unsigned Z, double gamma, double lambda2) {
  try {
    return lift(spectrum_b_monopole(Z, gamma, lambda2));
  } FFI_CATCH;
};

extern "C"
Function*
epa_spectrum_b_dipole(unsigned Z, double gamma, double lambda2) {
  try {
    return lift(spectrum_b_dipole(Z, gamma, lambda2));
  } FFI_CATCH;
};

template <typename Luminosity>
static inline Function* epa_luminosity_(
    Luminosity (*luminosity1)(Spectrum, Integrator),
    Luminosity (*luminosity2)(Spectrum, Spectrum, Integrator),
    Function* spectrum1,
    Function* spectrum2,
    Function* integrator
) {
  try {
    auto i = integrator ? lower<Integrator>(integrator) : default_integrator(0);
    auto n1 = lower<Spectrum>(spectrum1);
    if (spectrum1 == spectrum2) return lift(luminosity1(std::move(n1), i));
    auto n2 = lower<Spectrum>(spectrum2);
    return lift(luminosity2(std::move(n1), std::move(n2), i));
  } FFI_CATCH;
};

extern "C"
Function*
epa_luminosity(Function* spectrum1, Function* spectrum2, Function* integrator) {
  return epa_luminosity_(
      luminosity, luminosity, spectrum1, spectrum2, integrator
  );
};

extern "C"
Function*
epa_luminosity_y(Function* spectrum1, Function* spectrum2) {
  try {
    auto n1 = lower<Spectrum>(spectrum1);
    if (spectrum1 == spectrum2) return lift(luminosity_y(std::move(n1)));
    auto n2 = lower<Spectrum>(spectrum1);
    return lift(luminosity_y(std::move(n1), std::move(n2)));
  } FFI_CATCH;
};

extern "C"
Function*
epa_luminosity_fid(
    Function* spectrum1, Function* spectrum2, Function* integrator
) {
  return epa_luminosity_(
      luminosity_fid, luminosity_fid, spectrum1, spectrum2, integrator
  );
};

template <typename Luminosity>
static inline Function* epa_luminosity_b_(
    Luminosity (*luminosity1)(
      Spectrum_b,
      std::function<double (double)>,
      const std::function<Integrator (unsigned)>&,
      unsigned
    ),
    Luminosity (*luminosity2)(
      Spectrum_b,
      Spectrum_b,
      std::function<double (double)>,
      const std::function<Integrator (unsigned)>&,
      unsigned
    ),
    Function* spectrum1,
    Function* spectrum2,
    Function* upc_probability,
    Function* integrator_generator,
    unsigned  integration_level
) {
  try {
    auto n1  = lower<Spectrum_b>(spectrum1);
    auto upc = lower<double (double)>(upc_probability);

    auto ig  = integrator_generator
             ? lower<Integrator (unsigned)>(integrator_generator)
             : default_integrator;

    if (spectrum1 == spectrum2)
      return lift(
          luminosity1(
            std::move(n1),
            std::move(upc),
            std::move(ig),
            integration_level
          )
      );

    auto n2 = lower<Spectrum_b>(spectrum2);
    return lift(
        luminosity2(
          std::move(n1),
          std::move(n2),
          std::move(upc),
          std::move(ig),
          integration_level
        )
    );
  } catch (ForeignError&) {
    return nullptr;
  } catch (std::exception&) {
    set_error();
    return nullptr;
  };
};

extern "C"
Function*
epa_luminosity_b(
    Function* spectrum1,
    Function* spectrum2,
    Function* upc_probability,
    Function* integrator_generator,
    unsigned  integration_level
) {
  return epa_luminosity_b_(
      luminosity_b,
      luminosity_b,
      spectrum1,
      spectrum2,
      upc_probability,
      integrator_generator,
      integration_level
  );
};

extern "C"
Function*
epa_luminosity_y_b(
    Function* spectrum1,
    Function* spectrum2,
    Function* upc_probability,
    Function* integrator_generator,
    unsigned integration_level
) {
  return epa_luminosity_b_(
      luminosity_y_b,
      luminosity_y_b,
      spectrum1,
      spectrum2,
      upc_probability,
      integrator_generator,
      integration_level
  );
};

extern "C"
Function*
epa_luminosity_fid_b(
    Function* spectrum1,
    Function* spectrum2,
    Function* upc_probability,
    Function* integrator_generator,
    unsigned  integration_level
) {
  return epa_luminosity_b_(
      luminosity_fid_b,
      luminosity_fid_b,
      spectrum1,
      spectrum2,
      upc_probability,
      integrator_generator,
      integration_level
  );
};

extern "C"
Function*
epa_xsection(Function* photons_xsection, Function* luminosity) {
  try {
    return lift(
        xsection(
          lower<XSection>(photons_xsection),
          lower<Luminosity>(luminosity)
        )
    );
  } FFI_CATCH;
};

extern "C"
Function*
epa_xsection_b(Function* photons_xsection, Function* luminosity) {
  try {
    return lift(
        xsection_b(
          lower<XSection_b>(photons_xsection),
          lower<Luminosity_b>(luminosity)
        )
    );
  } FFI_CATCH;
};

extern "C"
Function*
epa_xsection_fid(
    Function* photons_xsection_pT,
    Function* luminosity,
    double mass,
    double pT_min,
    double eta_max,
    double w1_min,
    double w1_max,
    double w2_min,
    double w2_max,
    Function* integrator
) {
  try {
    return lift(
        xsection_fid(
          lower<double (double, double)>(photons_xsection_pT),
          lower<Luminosity_fid>(luminosity),
          mass,
          pT_min,
          eta_max,
          w1_min,
          w1_max,
          w2_min,
          w2_max,
          integrator ? lower<Integrator>(integrator) : default_integrator(0)
        )
    );
  } FFI_CATCH;
};

extern "C"
Function*
epa_xsection_fid_b(
    Function* photons_xsection_pT,
    Function* luminosity,
    double mass,
    double pT_min,
    double eta_max,
    double w1_min,
    double w1_max,
    double w2_min,
    double w2_max,
    Function* integrator
) {
  try {
    return lift(
        xsection_fid_b(
          lower<Polarization (double, double)>(photons_xsection_pT),
          lower<Luminosity_fid_b>(luminosity),
          mass,
          pT_min,
          eta_max,
          w1_min,
          w1_max,
          w2_min,
          w2_max,
          integrator ? lower<Integrator>(integrator) : cquad_integrator(0u)
        )
     );
  } FFI_CATCH;
};

extern "C" Function* epa_photons_to_fermions(double mass, double charge) {
  try {
    return lift(photons_to_fermions(mass, charge));
  } FFI_CATCH;
};

extern "C" Function* epa_photons_to_fermions_pT(double mass, double charge) {
  try {
    return lift(photons_to_fermions_pT(mass, charge));
  } FFI_CATCH;
};

extern "C" Function* epa_photons_to_fermions_b(double mass, double charge) {
  try {
    return lift(photons_to_fermions_b(mass, charge));
  } FFI_CATCH;
};

extern "C" Function* epa_photons_to_fermions_pT_b(double mass, double charge) {
  try {
    return lift(photons_to_fermions_pT_b(mass, charge));
  } FFI_CATCH;
};
