#include "../common.cpp"

#include <getopt.h>

// TODO
void usage(const char* argv0) {
  std::cout
    << "This program calculates the fiducial differential cross section for "
       "muon pair production in ultraperipheral proton-proton collisions. The "
       "fiducial region is defined by the cuts on transverse momentum pT and "
       "pseudorapidity eta of each muon: "
       "pT > 6 GeV for 12 GeV < sqrt(s) < 30 GeV, "
       "pT > 10 GeV for 30 GeV < sqrt(s) < 70 GeV; "
       "abs(eta) < 2.4. "
       "Parameters of this calculation are the same as were set in the ATLAS "
       "experiment for the paper 1708.04053.\n"
       "Usage: " << argv0 << "options...\n"
       "Allowed options:\n"
       "  -h or --help:     print\n"
       "  -n or --npoints:  number of points in the invariant mass range\n"
       "  -s or --step:     points step in the invariant mass range (GeV)\n"
       "  -v or --verbose:  be verbose while calculating\n"
       "  -S or --survival: account for non-electromagnetic interactions\n"
  ;
};

int main(int argc, char** argv) {
  option options[] = {
    { "help",     0, nullptr, 'h' },
    { "npoints",  1, nullptr, 'n' },
    { "step",     1, nullptr, 's' },
    { "verbose",  0, nullptr, 'v' },
    { "survival", 0, nullptr, 'S' },
    { nullptr,    0, nullptr, 0   }
  };

  std::string optstring;
  {
    std::stringstream ss;
    for (auto& o: options)
      if (o.val > 0) {
        ss << static_cast<char>(o.val);
        if (o.has_arg) ss << ':';
      };
    optstring = ss.str();
  };

  unsigned npoints = 0;
  double step = 0;
  bool verbose = false;
  bool survival = false;

  while (true) {
    int c = getopt_long(argc, argv, optstring.c_str(), options, nullptr);
    if (c == -1) break;
    switch (c) {
      case 'h':
        usage(argv[0]);
        return 0;

      case 'n':
        npoints = parse_unsigned(optarg, "number of points");
        break;

      case 's':
        step = parse_energy(optarg, "invariant mass step");
        break;

      case 'v':
        verbose = true;
        break;

      case 'S':
        survival = true;
        break;

      case '?':
        return 1;
    };
  };

  if (npoints != 0 && step != 0) {
    std::cerr << "--npoints and --step parameters are incompatible\n";
    return 0;
  };

  if (npoints == 0) npoints = step == 0 ? 100 : (70. - 30.) / step;

  auto grid = make_grid(12, 70, npoints);
  {
    auto i = std::lower_bound(
        grid.begin(),
        grid.end(),
        30.,
        [](const std::pair<double, double>& point, double x) -> bool {
          return point.first < x;
        }
    );
    if (i->first != 30) grid.insert(i, { 30, 0 });
  };

  const double muon_mass = 105.6583745e-3;

  gsl::init();
  epa::default_error_step = 1./3;

  auto xsection = survival
                ? make_function1d_async(
                    [=]() -> std::function<double (double)> {
                      return [
                        x0 = epa::pp_to_ppll_b(13e3, muon_mass, 6,  2.4),
                        x1 = epa::pp_to_ppll_b(13e3, muon_mass, 10, 2.4)
                      ](double rs) -> double {
                        return rs < 30 ? x0(rs) : x1(rs);
                      };
                    },
                    grid,
                    verbose ? "xsection_b " : nullptr
                  )
                : make_function1d(
                    [
                      x0 = epa::pp_to_ppll(13e3, muon_mass, 6,  2.4),
                      x1 = epa::pp_to_ppll(13e3, muon_mass, 10, 2.4)
                    ](double rs) -> double {
                      return rs < 30 ? x0(rs) : x1(rs);
                    },
                    grid
                  );
  xsection.dump(std::cout);

  if (verbose) {
    auto integrate = epa::gsl_integrator(0, 1e-2);
    std::cerr
      << "Integrated cross section: "
      << (  integrate(xsection, grid.front().first, 30)
          + integrate(xsection, 30, grid.back().first))
      << " b\n";
  };

  return 0;
};
