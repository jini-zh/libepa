#include <epa/proton.hpp>

#include "ffi.hpp"

using namespace epa;
using namespace ffi;

extern "C" double epa_proton_mass() {
  return proton_mass;
};

extern "C" double epa_proton_magnetic_moment() {
  return proton_magnetic_moment;
};

extern "C" double epa_proton_dipole_form_factor_lambda2() {
  return proton_dipole_form_factor_lambda2;
};

extern "C" Function* epa_proton_dipole_form_factor(double lambda2) {
  try {
    return lift(proton_dipole_form_factor(lambda2));
  } FFI_CATCH;
};

extern "C" Function* epa_proton_dipole_spectrum(double energy, double lambda2) {
  try {
    return lift(proton_dipole_spectrum(energy, lambda2));
  } FFI_CATCH;
};

extern "C"
Function*
epa_proton_dipole_spectrum_Dirac(double energy, double lambda2) {
  try {
    return lift(proton_dipole_spectrum_Dirac(energy, lambda2));
  } FFI_CATCH;
};

extern "C"
Function*
epa_proton_dipole_spectrum_b_Dirac(double energy, double lambda2) {
  try {
    return lift(proton_dipole_spectrum_b_Dirac(energy, lambda2));
  } FFI_CATCH;
};

extern "C" double epa_pp_elastic_slope(double collision_energy) {
  try {
    return pp_elastic_slope(collision_energy);
  } FFI_CATCH_R(0)
};

extern "C" Function* epa_pp_upc_probability(double collision_energy) {
  try {
    return lift(pp_upc_probability(collision_energy));
  } FFI_CATCH;
};

extern "C"
Function*
epa_pp_luminosity(double collision_energy, Function* integrator) {
  try {
    return lift(
        pp_luminosity(
          collision_energy,
          integrator ? lower<Integrator>(integrator) : default_integrator(0)
        )
    );
  } FFI_CATCH;
};

extern "C" Function* epa_pp_luminosity_y(double collision_energy) {
  try {
    return lift(pp_luminosity_y(collision_energy));
  } FFI_CATCH;
};

extern "C"
Function*
epa_pp_luminosity_fid(double collision_energy, Function* integrator) {
  try {
    return lift(
        pp_luminosity_fid(
          collision_energy,
          integrator ? lower<Integrator>(integrator) : default_integrator(0)
        )
    );
  } FFI_CATCH;
};

template <typename Luminosity>
static inline Function* epa_ppx_luminosity_b_(
    Luminosity (*luminosity)(
      Spectrum,
      Spectrum_b,
      double,
      const std::function<Integrator (unsigned)>&,
      unsigned
    ),
    Function* spectrum,
    Function* spectrum_b,
    double    B,
    Function* integrator_generator,
    unsigned  integration_level
) {
  try {
    return lift(
        luminosity(
          spectrum ? lower<Spectrum>(spectrum) : Spectrum(),
          lower<Spectrum_b>(spectrum_b),
          B,
          integrator_generator
          ? lower<Integrator (unsigned)>(integrator_generator)
          : default_integrator,
          integration_level
        )
    );
  } FFI_CATCH;
};

extern "C"
Function*
epa_ppx_luminosity_b(
    Function* spectrum,
    Function* spectrum_b,
    double    B,
    Function* integrator_generator,
    unsigned  integration_level
) {
  return epa_ppx_luminosity_b_(
      ppx_luminosity_b,
      spectrum,
      spectrum_b,
      B,
      integrator_generator,
      integration_level
  );
};

extern "C"
Function*
epa_ppx_luminosity_y_b(
    Function* spectrum_b,
    double    B,
    Function* integrator_generator,
    unsigned  integration_level
) {
  try {
    return lift(
        ppx_luminosity_y_b(
          lower<Spectrum_b>(spectrum_b),
          B,
          integrator_generator
          ? lower<Integrator (unsigned)>(integrator_generator)
          : default_integrator,
          integration_level
        )
    );
  } FFI_CATCH;
};

extern "C"
Function*
epa_ppx_luminosity_fid_b(
    Function* spectrum,
    Function* spectrum_b,
    double    B,
    Function* integrator_generator,
    unsigned  integration_level
) {
    return epa_ppx_luminosity_b_(
      ppx_luminosity_fid_b,
      spectrum,
      spectrum_b,
      B,
      integrator_generator,
      integration_level
    );
};

template <typename Result>
static inline Function* epa_pp_luminosity_b_(
    Result (*luminosity)(
      double, const std::function<Integrator (unsigned)>&, unsigned
    ),
    double    collision_energy,
    Function* integrator_generator,
    unsigned  integration_level
) {
  try {
    return lift(
        luminosity(
          collision_energy,
          integrator_generator
          ? lower<Integrator (unsigned)>(integrator_generator)
          : default_integrator,
          integration_level
        )
    );
  } FFI_CATCH;
};

extern "C"
Function*
epa_pp_luminosity_b(
    double    collision_energy,
    Function* integrator_generator,
    unsigned  integration_level
) {
  return epa_pp_luminosity_b_(
      pp_luminosity_b,
      collision_energy,
      integrator_generator,
      integration_level
  );
};

extern "C"
Function*
epa_pp_luminosity_y_b(
    double    collision_energy,
    Function* integrator_generator,
    unsigned  integration_level
) {
  return epa_pp_luminosity_b_(
      pp_luminosity_y_b,
      collision_energy,
      integrator_generator,
      integration_level
  );
};

extern "C"
Function*
epa_pp_luminosity_fid_b(
    double    collision_energy,
    Function* integrator_generator,
    unsigned  integration_level
) {
  return epa_pp_luminosity_b_(
      pp_luminosity_fid_b,
      collision_energy,
      integrator_generator,
      integration_level
  );
};

extern "C"
Function*
epa_pp_to_ppll(
    double collision_energy,
    double mass,
    double charge,
    double pT_min,
    double eta_max,
    Function* integration_generator,
    unsigned  integration_level
) {
  try {
    return lift(
        pp_to_ppll(
          collision_energy,
          mass,
          charge,
          pT_min,
          eta_max,
          integration_generator
          ? lower<Integrator (unsigned)>(integration_generator)
          : default_integrator,
          integration_level
        )
    );
  } FFI_CATCH;
};

extern "C"
Function*
epa_pp_to_ppll_b_integrator(unsigned level) {
  try {
    return lift(pp_to_ppll_b_integrator(level));
  } FFI_CATCH;
};

extern "C"
Function*
epa_pp_to_ppll_b(
    double    collision_energy,
    double    mass,
    double    charge,
    double    pT_min,
    double    eta_max,
    Function* integrator_generator,
    unsigned  integration_level
) {
  try {
    return lift(
        pp_to_ppll_b(
          collision_energy,
          mass,
          charge,
          pT_min,
          eta_max,
          integrator_generator
          ? lower<Integrator (unsigned)>(integrator_generator)
          : pp_to_ppll_b_integrator(0),
          integration_level
        )
    );
  } FFI_CATCH;
};
