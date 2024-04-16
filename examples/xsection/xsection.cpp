#include "../common.cpp"

#include <getopt.h>

void usage(const char* argv0) {
  std::cout
    << "This program plots the cross section for the production of a pair of "
       "fermions in ultraperipheral proton-proton collisions with respect to "
       "the fermion mass\n"
       "Usage: " << argv0 << " options...\n"
       "Allowed options:\n"
       "  -h or --help:             print this message and exit (no argument)\n"
       "  -E or --collision-energy: specify proton-proton collision energy (in GeV). Default is 13e3\n"
       "  -f or --from:             beginning of the mass range (GeV)\n"
       "  -t or --to:               ending of the mass range (GeV)\n"
       "  -n or --npoints:          number of points in the mass range\n"
       "  -s or --step:             points step in the mass range (GeV)\n"
       "  -l or --log:              assume logarithmic scale in the mass range (no argument)\n"
       "  -v or --verbose:          be verbose while calculating\n"
       "  -S or --survival:         account for non-electromagnetic interactions\n"
       "  --rs-max:                 cut off integration over invariant mass (sqrt{s}) at this value (GeV)\n"
  ;
};

int main(int argc, char** argv) {
  gsl::init();
  epa::default_relative_error = 1e-2;

  option options[] = {
    { "help",             0, nullptr, 'h' },
    { "collision-energy", 1, nullptr, 'E' },
    { "from",             1, nullptr, 'f' },
    { "to",               1, nullptr, 't' },
    { "npoints",          1, nullptr, 'n' },
    { "step",             1, nullptr, 's' },
    { "log",              0, nullptr, 'l' },
    { "verbose",          0, nullptr, 'v' },
    { "survival",         0, nullptr, 'S' },
    { "rs-max",           1, nullptr,  -2 },
    { nullptr,            0,       0,   0 }
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

  double collision_energy = 13e3;
  double from = -1;
  double to = -1;
  double step = 0;
  bool log = false;
  unsigned npoints = 0;
  bool verbose = false;
  bool survival = false;
  double rs_max = 0;

  while (true) {
    int c = getopt_long(argc, argv, optstring.c_str(), options, nullptr);
    if (c == -1) break;
    switch (c) {
      case 'h':
        usage(argv[0]);
        return 0;

      case 'E':
        collision_energy = parse_energy(optarg, "collision energy");
        break;

      case 'f':
        from = parse_energy(optarg, "beginning mass (from)");
        break;

      case 't':
        to = parse_energy(optarg, "ending mass (to)");
        break;

      case 's':
        step = parse_energy(optarg, "mass step");
        break;

      case 'n':
        npoints = parse_unsigned(optarg, "number of points");
        break;

      case 'l':
        log = true;
        break;

      case 'v':
        verbose = true;
        break;

      case 'S':
        survival = true;
        break;

      case -2:
        rs_max = parse_energy(optarg, "invariant mass cutoff");
        break;

      case '?':
        return 1;
    };
  };

  if ((from > 0) + (to > 0) + (step > 0) + (npoints > 0) != 3) {
    std::cerr << "Exactly 3 of --from, --to, --step, --npoints must be specified\n";
    return 1;
  };

  if (from < 0) {
    from = to - step * npoints;
    if (from < 0) {
      std::cerr << "Beginning of the invariant mass range cannot be negative or zero\n";
      return 1;
    };
  };

  if (!npoints) npoints = (to - from) / step + 1;

  if (rs_max == 0) rs_max = collision_energy / 2;

  auto grid = make_grid(from, to, npoints, log);

  auto xsection = survival
                ? make_function1d_async(
                    [=]() -> std::function<double (double)> {
                      return [
                        =,
                        integrate  = epa::default_integrator(0),
                        luminosity = epa::pp_luminosity_b(
                            collision_energy,
                            epa::default_integrator,
                            1
                        )
                      ](double mass) -> double {
                        return integrate(
                            epa::xsection_b(epa::photons_to_fermions_b(mass), luminosity),
                            2 * mass,
                            rs_max
                        );
                      };
                    },
                    grid,
                    verbose ? std::string("cross section ") : std::string()
                  )
                : make_function1d(
                    [
                      =,
                      integrate = epa::default_integrator(0),
                      luminosity = epa::pp_luminosity(
                        collision_energy, epa::default_integrator(1)
                      )
                    ](double mass) -> double {
                      return integrate(
                          epa::xsection(epa::photons_to_fermions(mass), luminosity),
                          2 * mass,
                          rs_max
                      );
                    },
                    grid
                  );

  xsection.dump(std::cout);

  return 0;
};
