#pragma once

#include <epa/epa.hpp>

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

// Equivalent photon spectrum for a proton of given energy (dipole approximation).
// Both Dirac and Pauli current terms are taken into account (F_1(Q^2) \psi
// \gamma^\mu \psi + F_2(Q^2) / 2 m_p * \psi \sigma^{\mu \nu} q_\nu \psi).
// This is the correct spectrum, but it has no _b counterpart.
Spectrum
proton_dipole_spectrum(
    double energy,
    double lambda2 = proton_dipole_form_factor_lambda2
);

// Equivalent photon spectrum for a proton of given energy (dipole
// approximation, only the Dirac current term (F_1(Q^2) \psi \gamma^\mu \psi)
// is taken into account).
Spectrum
proton_dipole_spectrum_Dirac(
    double energy,
    double lambda2 = proton_dipole_form_factor_lambda2
);
Spectrum_b
proton_dipole_spectrum_b_Dirac(
    double energy,
    double lambda2 = proton_dipole_form_factor_lambda2
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
    const std::function<Integrator (unsigned)>& = default_integrator,
    unsigned integration_level = 0
);
Luminosity_y_b ppx_luminosity_y_b(
    Spectrum_b,
    double B,
    const std::function<Integrator (unsigned)>& = default_integrator,
    unsigned integration_level = 0
);
Luminosity_fid_b ppx_luminosity_fid_b(
    Spectrum,
    Spectrum_b,
    double B,
    const std::function<Integrator (unsigned)>& = default_integrator,
    unsigned integration_level = 0
);

// Photon-photon luminosity in proton-proton ultraperipheral collisions
Luminosity_b pp_luminosity_b(
    double collision_energy,
    const std::function<Integrator (unsigned)>& = default_integrator,
    unsigned integration_level = 0
);
Luminosity_y_b pp_luminosity_y_b(
    double collision_energy,
    const std::function<Integrator (unsigned)>& = default_integrator,
    unsigned integration_level = 0
);
Luminosity_fid_b pp_luminosity_fid_b(
    double collision_energy,
    const std::function<Integrator (unsigned)>& = default_integrator,
    unsigned integration_level = 0
);

// Convenience functions to calculate cross sections for the production of a
// fermion pair in ultraperipheral collisions of protons. It is assumed that
// the fermions only interact electromagnetically with the protons
// (approximately like leptons).

// Total cross section with non-electromagnetic interactions between protons
// neglected
XSection pp_to_ppll(
    double collision_energy,
    double mass,   // of the produced particle
    double charge, // of the produced particle
    Integrator = default_integrator(0)
);
// Same for charge = 1
XSection pp_to_ppll(
    double collision_energy,
    double mass,
    Integrator = default_integrator(0)
);
// Fiducial cross section with non-electromagnetic interactions between protons
// neglected
XSection pp_to_ppll(
    double collision_energy,
    double mass,    // of the produced particle
    double charge,  // of the produced particle
    double pT_min,  // minimal transverse momentum of the produced particle
    double eta_max, // maximal pseudorapidity of the produced particle
    const std::function<Integrator (unsigned)>& = default_integrator,
    unsigned integration_level = 0
);
// Same for charge = 1
XSection pp_to_ppll(
    double collision_energy,
    double mass,
    double pT_min,
    double eta_max,
    const std::function<Integrator (unsigned)>& = default_integrator,
    unsigned integration_level = 0
);
// Total cross section with non-electromagnetic interactions between protons
// respected
XSection pp_to_ppll_b(
    double collision_energy,
    double mass,   // of the produced particle
    double charge, // of the produced particle
    const std::function<Integrator (unsigned)>& = default_integrator,
    unsigned integration_level = 0
);
// Same for charge = 1
XSection pp_to_ppll_b(
    double collision_energy,
    double mass,
    const std::function<Integrator (unsigned)>& = default_integrator,
    unsigned integration_level = 0
);
// Default integrator for the following function.
// Returns gsl_cquad_integrator for level `level'
std::function<Integrator (unsigned)>
pp_to_ppll_b_integrator(unsigned level = 0);
// Fiducial cross section with non-electromagnetic interactions between protons
// respected
XSection pp_to_ppll_b(
    double collision_energy,
    double mass,    // of the produced particle
    double charge,  // of the produced particle
    double pT_min,  // minimal transverse momentum of the produced particle
    double eta_max, // maximal pseudorapidity of the produced particle
    const std::function<Integrator (unsigned)>& = pp_to_ppll_b_integrator(0),
    unsigned integration_level = 0
);
// Same for charge = 1
XSection pp_to_ppll_b(
    double collision_energy,
    double mass,
    double pT_min,
    double eta_max,
    const std::function<Integrator (unsigned)>& = pp_to_ppll_b_integrator(0),
    unsigned integration_level = 0
);

}; // namespace epa
