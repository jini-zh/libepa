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

}; // namespace epa
