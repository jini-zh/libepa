#ifndef EPA_EPA_H
#define EPA_EPA_H

#include "ffi.h"

typedef struct epa_qag_workspace   epa_qag_workspace;
typedef struct epa_cquad_workspace epa_cquad_workspace;

typedef struct epa_polarization {
  double parallel;
  double perpendicular;
} epa_polarization;

#define defun(name, result, ...) \
  typedef struct name { \
    result (*function)(__VA_ARGS__ __VA_OPT__(,) void*); \
    void* data; \
    void (*destructor)(void*); \
  } name

defun(epa_function1d, double, double);
defun(epa_function2d, double, double, double);
defun(epa_function3d, double, double, double, double);

defun(epa_integrator, double, epa_function1d*, double, double);
defun(epa_integrator_generator, epa_integrator*, unsigned);

defun(epa_luminosity_b_f,     double, double, epa_polarization);
defun(epa_luminosity_y_b_f,   double, double, double, epa_polarization);
defun(epa_luminosity_fid_b_f, double, double, epa_polarization, double, double);

defun(epa_xsection_b_f,  epa_polarization, double);
defun(epa_xsection_pT_b, epa_polarization, double, double);

#undef defun

#define defconst(type, name) type epa_ ## name()
#define defvar(type, name) \
  type epa_get_##name(); \
  void epa_set_##name(type)

#include "epa_vars.h"

#undef defvar
#undef defconst

void epa_init(void);
void epa_version(int* major, int* minor, int* patch);

int epa_get_default_integration_method(void);
void epa_set_default_integration_method(int);

epa_qag_workspace* epa_make_qag_integration_workspace(size_t limit);
void epa_destroy_qag_integration_workspace(epa_qag_workspace*);

epa_integrator*
epa_qag_integrator(
    double absolute_error,
    double relative_error,
    int method,
    epa_qag_workspace*
);

epa_cquad_workspace* epa_make_cquad_integration_workspace(size_t limit);
void epa_destroy_cquad_integration_workspace(epa_cquad_workspace*);

epa_integrator*
epa_cquad_integrator(
    double absolute_error,
    double relative_error,
    epa_cquad_workspace*
);

epa_function1d* epa_form_factor_monopole(double lambda2);
epa_function1d* epa_form_factor_dipole  (double lambda2);

epa_function1d*
epa_spectrum(
    unsigned Z, double gamma, epa_function1d* form_factor, epa_integrator*
);

epa_function1d* epa_spectrum_monopole(unsigned Z, double gamma, double lambda2);
epa_function1d* epa_spectrum_dipole  (unsigned Z, double gamma, double lambda2);

epa_function2d*
epa_spectrum_b(
    unsigned Z, double gamma, epa_function1d* form_factor, epa_integrator*
);

epa_function2d* epa_spectrum_b_point(unsigned Z, double gamma);

epa_function2d*
epa_spectrum_b_monopole(unsigned Z, double gamma, double lambda2);

epa_function2d*
epa_spectrum_b_dipole(unsigned Z, double gamma, double lambda2);

epa_function1d*
epa_luminosity(
    epa_function1d* spectrum1,
    epa_function1d* spectrum2,
    epa_integrator*
);

epa_function2d*
epa_luminosity_y(
    epa_function1d* spectrum1,
    epa_function1d* spectrum2
);

epa_function3d*
epa_luminosity_fid(
    epa_function1d* spectrum1,
    epa_function1d* spectrum2,
    epa_integrator*
);

epa_luminosity_b_f*
epa_luminosity_b(
    epa_function2d* spectrum1,
    epa_function2d* spectrum2,
    epa_function1d* upc_probability,
    epa_integrator_generator*,
    unsigned        level
);

epa_luminosity_y_b_f*
epa_luminosity_y_b(
    epa_function2d* spectrum1,
    epa_function2d* spectrum2,
    epa_function1d* upc_probability,
    epa_integrator_generator*,
    unsigned        level
);

epa_luminosity_fid_b_f*
epa_luminosity_fid_b(
    epa_function2d* spectrum1,
    epa_function2d* spectrum2,
    epa_function1d* upc_probability,
    epa_integrator_generator*,
    unsigned        level
);

epa_function1d*
epa_xsection(
    epa_function1d* photons_xsection,
    epa_function1d* luminosity
);

epa_function1d*
epa_xsection_b(
    epa_xsection_b_f*   photons_xsection,
    epa_luminosity_b_f* luminosity
);

epa_function1d*
epa_xsection_fid(
    epa_function2d* photons_xsection_pT,
    epa_function3d* luminosity_fid,
    double          mass,
    double          pT_min,
    double          eta_max,
    double          w1_min,
    double          w1_max,
    double          w2_min,
    double          w2_max,
    epa_integrator*
);

epa_function1d*
epa_xsection_fid_b(
    epa_xsection_pT_b*      photons_xsection_pT,
    epa_luminosity_fid_b_f* luminosity,
    double                  mass,
    double                  pT_min,
    double                  eta_max,
    double                  w1_min,
    double                  w1_max,
    double                  w2_min,
    double                  w2_max,
    epa_integrator*
);

epa_function1d*    epa_photons_to_fermions     (double mass, double charge);
epa_function2d*    epa_photons_to_fermions_pT  (double mass, double charge);
epa_xsection_b_f*  epa_photons_to_fermions_b   (double mass, double charge);
epa_xsection_pT_b* epa_photons_to_fermions_pT_b(double mass, double charge);


#endif // EPA_EPA_H
