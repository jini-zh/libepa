#include <filesystem>
#include <fstream>

#define BOOST_TEST_MODULE epa
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <epa/proton.hpp>
#include "a1.cpp"

namespace epa {

namespace test {

struct InitFixture {
  void setup() {
    gsl::init();
  };
};

BOOST_TEST_GLOBAL_FIXTURE(InitFixture);

BOOST_AUTO_TEST_SUITE(epa_test);

BOOST_AUTO_TEST_CASE(epa_form_factors, *boost::unit_test::tolerance(1e-5)) {
  BOOST_TEST(
      form_factor_monopole(sqr(80e-3))(1e1) == 6.3959066197633518e-4
  );
  BOOST_TEST(
      form_factor_dipole(proton_dipole_form_factor_lambda2)(1e1)
      == 3.833392459597915e-3
  );
};

BOOST_AUTO_TEST_CASE(epa_spectra, *boost::unit_test::tolerance(1e-5)) {
  BOOST_TEST(
      spectrum(
        82,
        5.02e3 / 2 / amu,
        [](double) -> double { return 1; },
        qag_integrator({ .relative_error = 1e-5 })
      )(1e2)
      == 111.73111054940296
  );

  BOOST_TEST(
      spectrum_monopole(82, 5.02e3 / 2 / amu, sqr(80e-3))(1e2)
      == 0.074371161172122113
  );

  BOOST_TEST(
      spectrum_dipole(
        1, 13e3 / 2 / proton_mass, proton_dipole_form_factor_lambda2
      )(1e2)
      == 0.00012165634320545708
  );

  BOOST_CHECK_CLOSE_FRACTION(
      spectrum_b(
        82,
        5.02e3 / 2 / amu,
        form_factor_monopole(sqr(80e-3)),
        qag_integrator({ .relative_error = 1e-3 })
      )(1e1 * fm, 1e2),
      1.4675442690292509e-06,
      1e-3
  );

  BOOST_TEST(
      spectrum_b_point(82, 5.02e3 / 2 / amu)(1e1 * fm, 1e2)
      == 1.8377284965685942e-06
  );

  {
    auto n1 = spectrum_b_monopole(82, 5.02e3 / 2 / amu, sqr(80e-3));
    BOOST_TEST(n1(10 * fm, 1e2) == 1.4677997339675521e-6);
    BOOST_CHECK_CLOSE_FRACTION(n1(0.1 * fm, 2.2e5), 2.5853814193024425e-48, 1e-3);
    BOOST_TEST(n1(0.02 * fm, 1e2) == 1.3855209874282745e-7);
  };

  BOOST_TEST(
      spectrum_b_dipole(
        1, 13e3 / 2 / proton_mass, proton_dipole_form_factor_lambda2
      )(1.5 * fm, 1e2)
      == 1.1706570931517264e-07
  );
};

struct A1_fixture {
  std::shared_ptr<Function1d> form_factor;
  A1_fixture() {
    size_t npoints = sizeof(A1_FORM_FACTOR) / sizeof(double) / 3;
    form_factor = std::make_shared<Function1d>();
    form_factor->points->reserve(npoints);
    double m2 = sqr(2 * proton_mass);
    const double* a1 = A1_FORM_FACTOR;
    for (size_t i = 0; i < npoints; ++i) {
      double q2       = *a1++;
      double electric = *a1++;
      double magnetic = *a1++;
      double tau = q2 / m2;
      form_factor->points->push_back({
          q2,
          (electric + tau * proton_magnetic_moment * magnetic) / (1 + tau)
      });
    };
  };
};

BOOST_FIXTURE_TEST_SUITE(a1_form_factor, A1_fixture);

BOOST_AUTO_TEST_CASE(epa_spectra) {
  // Note the difference in the values in the following two tests. These
  // methods have poor accuracy due to oscillating nature of the integrand.  It
  // is strongly advised to use analytical calculation for the EPA spectrum.
  BOOST_TEST(
      spectrum_b_function1d_g(
        1,
        13e3 / 2 / proton_mass,
        form_factor
      )(fm, 1e2) == 2.5985895117901073e-07
  );

  BOOST_TEST(
      spectrum_b_function1d_s(
        1,
        13e3 / 2 / proton_mass,
        form_factor
      )(fm, 1e2) == 2.6014110821783035e-07
  );

  // after ~0.8 GeV^2 the error becomes too large
  size_t i = std::upper_bound(
                 form_factor->points->begin(),
                 form_factor->points->end(),
                 0.8,
                 [](double x, const std::pair<double, double>& point) {
                   return x < point.first;
                 }
             )
            - form_factor->points->begin();
  form_factor->points->resize(i);

  // This calculation has poor accuracy and depends on where the form factor is
  // cut off (the value of i above)
  BOOST_TEST(
      spectrum_b_function1d_g(
        1,
        13e3 / 2 / proton_mass,
        form_factor,
        proton_dipole_form_factor(0.67),
        proton_dipole_spectrum_b_Dirac(13e3 / 2, 0.67),
        5 * proton_radius
      )(fm, 1e2) == 2.7996649347083567e-07
  );
  
  BOOST_TEST(
      spectrum_b_function1d_s(
        1,
        13e3 / 2 / proton_mass,
        form_factor,
        proton_dipole_form_factor(0.67),
        proton_dipole_spectrum_b_Dirac(13e3 / 2, 0.67),
        5 * proton_radius
      )(fm, 1e2) == 2.8030979763022275e-07
  );

  auto n = spectrum_b_function1d(
      1,
      13e3 / 2 / proton_mass,
      form_factor
  );
  BOOST_TEST(n(fm, 1e2) == 2.7996649347083567e-07);
  // TODO: find such b, w that n fails to calculate
};

BOOST_AUTO_TEST_SUITE_END(); // fixture

BOOST_AUTO_TEST_CASE(epa_luminosity) {
  BOOST_TEST(
      luminosity_y(proton_dipole_spectrum(13e3 / 2))(100, 1)
      == 4.8783899969641711e-06
  );

  BOOST_TEST(
      luminosity_y(proton_dipole_spectrum_Dirac(13e3 / 2))(100, 1)
      == 4.6213647977030971e-06
  );

  BOOST_TEST(
      luminosity(proton_dipole_spectrum(13e3 / 2))(100)
      == 2.6904531638939847e-05
  );

  BOOST_TEST(
      luminosity(proton_dipole_spectrum_Dirac(13e3 / 2))(100)
      == 2.5099166887786259e-05
  );
};

BOOST_AUTO_TEST_CASE(epa_xsection) {
  BOOST_TEST(photons_to_fermions_pT(100)(250, 15) == 1.4291814382449728e-14);

  BOOST_TEST(
      xsection_fid(
        photons_to_fermions_pT(100),
        pp_luminosity_fid(13e3),
        100,
        10,
        2.5
      )(250) == 1.2127955941113192e-17
  );

  {
    auto x = photons_to_fermions_pT_b(100)(250, 15);
    BOOST_TEST(x.parallel      == 7.7801624393910958e-15);
    BOOST_TEST(x.perpendicular == 2.0803466325508358e-14);
  };
};

BOOST_AUTO_TEST_SUITE_END(); // epa_test


BOOST_AUTO_TEST_CASE(proton_test, *boost::unit_test::tolerance(1e-5)) {
  BOOST_TEST(
      proton_dipole_form_factor()(1e1) == 0.008916207850484446
  );
  BOOST_TEST(
      proton_dipole_spectrum(13e3 / 2)(1e2) == 0.00013072158743362693
  );
  BOOST_TEST(
      proton_dipole_spectrum_Dirac(13e3 / 2)(1e2) == 0.00012672479265553522
  );
  BOOST_TEST(
      proton_dipole_spectrum_b_Dirac(13e3 / 2)(fm, 1e2) == 2.2984963313647681e-07
  );

  BOOST_TEST(pp_luminosity(13e3)(100) == 2.6904531638939847e-05);

  BOOST_TEST(pp_to_ppll(13e3, 100, 10, 2.5)(250) == 1.2127955941113193e-17);
};

BOOST_AUTO_TEST_SUITE(expensive, *boost::unit_test::disabled())

BOOST_AUTO_TEST_CASE(test_luminosity_b, *boost::unit_test::tolerance(1e-5)) {
  BOOST_TEST(
      luminosity_b(proton_dipole_spectrum_b_Dirac(13e3 / 2), pp_upc_probability(13e3))(
        100, { 1, 1 }
      ) == 2.2902623308910459e-05
  );
};

BOOST_AUTO_TEST_CASE(test_pp_luminosity_b, *boost::unit_test::tolerance(1e-5)) {
  BOOST_TEST(pp_luminosity_b(13e3)(100, { 1, 1 }) == 2.290215747968901e-05);
};

BOOST_AUTO_TEST_CASE(test_xsection_fid_b, *boost::unit_test::tolerance(1e-5)) {
  BOOST_TEST(
      xsection_fid_b(
        photons_to_fermions_pT_b(100),
        pp_luminosity_fid_b(13e3),
        100,
        10,
        2.5,
        cquad_integrator(0, 1e-2)
      )(250) == 9.7800068779352812e-18
  );
};

BOOST_AUTO_TEST_CASE(test_ppll_b, *boost::unit_test::tolerance(1e-5)) {
  BOOST_TEST(
      pp_to_ppll_b(
        13e3, 100, 10, 2.5,
        [](unsigned level) -> Integrator {
           return level == 0
                  ? cquad_integrator(0, 1e-2)
                  : qag_integrator(level - 1);
        }
      )(250) == 9.7800054406137899e-18
  );
};

BOOST_AUTO_TEST_SUITE_END(); // expensive

}; // namespace test

}; // namespace epa
