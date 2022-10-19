#ifndef EPA_PROTON_H
#define EPA_PROTON_H

#include "epa.h"

double epa_proton_mass(void);
double epa_proton_magnetic_moment(void);
double epa_proton_dipole_form_factor_lambda2(void);

epa_function1d* epa_proton_dipole_form_factor(double lambda2);

epa_function1d* epa_proton_dipole_spectrum(double energy, double lambda2);
epa_function1d* epa_proton_dipole_spectrum_Dirac(double energy, double lambda2);

epa_function2d*
epa_proton_dipole_spectrum_b_Dirac(double energy, double lambda2);

double epa_pp_elastic_slope(double collision_energy);

epa_function1d* epa_pp_upc_probability(double collision_energy);

epa_function1d* epa_pp_luminosity(double collision_energy, epa_integrator*);
epa_function2d* epa_pp_luminosity_y(double collision_energy);
epa_function3d* epa_pp_luminosity_fid(double collision_energy, epa_integrator*);

epa_luminosity_b_f*
epa_ppx_luminosity_b(
    epa_function1d* spectrum,
    epa_function2d* spectrum_b,
    double          B,
    epa_integrator_generator*,
    unsigned        integration_level
);

epa_luminosity_y_b_f*
epa_ppx_luminosity_y_b(
    epa_function2d* spectrum_b,
    double          B,
    epa_integrator_generator*,
    unsigned        integration_level
);

epa_luminosity_fid_b_f*
epa_ppx_luminosity_fid_b(
    epa_function1d* spectrum,
    epa_function2d* spectrum_b,
    double          B,
    epa_integrator_generator*,
    unsigned        integration_level
);

epa_luminosity_b_f*
epa_pp_luminosity_b(
    double collision_energy,
    epa_integrator_generator*,
    unsigned integration_level
);

epa_luminosity_y_b_f*
epa_pp_luminosity_y_b(
    double collision_energy,
    epa_integrator_generator*,
    unsigned integration_level
);

epa_luminosity_fid_b_f*
epa_pp_luminosity_fid_b(
    double collision_energy,
    epa_integrator_generator*,
    unsigned integration_level
);

epa_function1d*
epa_pp_to_ppll(
    double collision_energy,
    double mass,
    double charge,
    double pT_min,
    double eta_max,
    epa_integrator_generator*,
    unsigned integration_level
);

epa_integrator_generator* epa_pp_to_ppll_b_integrator(unsigned level);

epa_function1d*
epa_pp_to_ppll_b(
    double collision_energy,
    double mass,
    double charge,
    double pT_min,
    double eta_max,
    epa_integrator_generator*,
    unsigned integration_level
);

#endif // EPA_PROTON_H
