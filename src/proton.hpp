#pragma once

#include "epa.hpp"

namespace epa {

const double proton_mass            = 0.93827208816; // GeV [1]
const double proton_radius          = 0.8414 * fm;   // GeV^{-1} [1]
const double proton_magnetic_moment = 2.79284734463; // [1]
/*
 [1] E. Tiesinga, P. J. Mohr, D. B. Newell, B. N. Taylor.
     CODATA recommended values of the fundamental physical constants: 2018.
     Reviews of Modern Physics 93, 025010 (2020).
*/

const double proton_dipole_form_factor_lambda2
  = 12.0 / sqr(proton_radius);

// Proton form factor in dipole approximation.
FormFactor
proton_dipole_form_factor(double lambda2 = proton_dipole_form_factor_lambda2);

// Equivalent photon spectra for a proton of given eneregy (dipole approximation).
Spectrum
proton_dipole_spectrum(
    double energy,
    double lambda = proton_dipole_form_factor_lambda2
);
Spectrum_b
proton_dipole_spectrum_b(
    double energy,
    double lambda = proton_dipole_form_factor_lambda2
);

// The slope of the cross section for elastic scattering of two protons
// (parameter B in [1112.3243])
double pp_elastic_slope(double collision_energy);

// The probability to avoid non-electromagnetic interactions in an
// ultraperipheral collision of two protons with the impact parameter b
std::function<double (double /* b */)>
pp_upc_probability(double collision_energy);

// Photon-photon luminosity in ultraperipheral proton-proton collisions with
// non-electromagnetic interactions neglected
Luminosity     pp_luminosity(
    double collision_energy, Integrator = default_integrator(0)
);
Luminosity_y   pp_luminosity_y(double collision_energy);
Luminosity_fid pp_luminosity_fid(
    double collision_energy, Integrator = default_integrator(0)
);

// Photon-photon luminosity in ultraperipheral collisions of two identical
// particles with the same probability to avoid non-electromagnetic
// interactions as proton but with arbitrary spectrum. Use this function if you
// want to try non-default proton spectrum.
Luminosity_b ppx_luminosity_b(
    Spectrum, // optional --- when provided might speed up calculations
    Spectrum_b,
    double B, // = pp_elastic_slope(collision_energy)
    const std::function<Integrator (unsigned)>& = default_integrator
);
Luminosity_b_y ppx_luminosity_b_y(
    Spectrum_b,
    double B,
    const std::function<Integrator (unsigned)>& = default_integrator
);
Luminosity_b_fid ppx_luminosity_b_fid(
    Spectrum,
    Spectrum_b,
    double B,
    const std::function<Integrator (unsigned)>& = default_integrator
);

// Photon-photon luminosity in proton-proton ultraperipheral collisions
Luminosity_b pp_luminosity_b(
    double collision_energy,
    const std::function<Integrator (unsigned)>& = default_integrator
);
Luminosity_b_y pp_luminosity_b_y(
    double collision_energy,
    const std::function<Integrator (unsigned)>& = default_integrator
);
Luminosity_b_fid pp_luminosity_b_fid(
    double collision_energy,
    const std::function<Integrator (unsigned)>& = default_integrator
);

}; // namespace epa
