#include "../common.cpp"

#include <getopt.h>

void usage(const char* argv0) {
  std::cout
    << "This program calculates photon-photon luminosity in ultraperipheral proton-proton collisions\n"
       "Usage: " << argv0 << " options... > output\n"
       "Allowed options:\n"
       "  -h or --help:             print this message and exit (no argument)\n"
       "  -E or --collision-energy: specify proton-proton collision energy (in GeV). Default is 13e3\n"
       "  -f or --from:             beginning of the invariant mass range (GeV)\n"
       "  -t or --to:               ending of the invariant mass range (GeV)\n"
       "  -n or --npoints:          number of points in the invariant mass range\n"
       "  -s or --step:             points step in the invariant mass range (GeV)\n"
       "  -l or --log:              assume logarithmic scale in the invariant mass range (no argument)\n"
       "  -T or --nthreads:         use this many threads for calculation. Default is " << get_nprocs() << '\n'
    << "  -v or --verbose:          be verbose while calculating.\n"
       "  -S or --survival:         account for non-electromagnetic interactions.\n"
       "  -p or --polarization parallel|perpendicular: assume polarized photons (no argument). Implies -S.\n";
};

epa::Function1d pp_luminosity_b(
    double collision_energy,
    epa::Polarization polarization,
    std::vector<std::pair<double, double>> grid,
    const std::string& verbose = std::string()
) {
  return make_function1d_async(
      [&]() -> std::function<double (double)> {
        return [=, l = epa::pp_luminosity_b(collision_energy)](
            double rs
        ) -> double {
          return l(rs, polarization);
        };
      },
      std::move(grid),
      verbose
  );
};

int main(int argc, char** argv) {
  option options[] = {
    { "help",             0, nullptr, 'h' },
    { "collision-energy", 1, nullptr, 'E' },
    { "from",             1, nullptr, 'f' },
    { "to",               1, nullptr, 't' },
    { "npoints",          1, nullptr, 'n' },
    { "step",             1, nullptr, 's' },
    { "log",              0, nullptr, 'l' },
    { "nthreads",         1, nullptr, 'T' },
    { "verbose",          0, nullptr, 'v' },
    { "survival",         0, nullptr, 'S' },
    { "polarization",     1, nullptr, 'p' },
    { 0,                  0, 0,        0  }
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
  unsigned nthreads = 0;
  bool verbose = false;
  bool survival = false;
  epa::Polarization polarization = { 1, 1 };

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
        from = parse_energy(optarg, "beginning invariant mass (from)");
        break;

      case 't':
        to = parse_energy(optarg, "ending invariant mass (to)");
        break;

      case 's':
        step = parse_energy(optarg, "invariant mass step");
        break;

      case 'n':
        npoints = parse_unsigned(optarg, "number of points");
        break;

      case 'l':
        log = true;
        break;

      case 'T':
        nthreads = parse_unsigned(optarg, "number of threads");
        break;

      case 'v':
        verbose = true;
        break;

      case 'S':
        survival = true;
        break;

      case 'p':
        if (strcmp(optarg, "parallel") == 0)
          polarization = { 1, 0 };
        else if (strcmp(optarg, "perpendicular") == 0)
          polarization = { 0, 1 };
        else {
          std::cerr
            << "Invalid polarization `" << optarg << "', expected `parallel' or `perpendicular'\n";
          return 1;
        };
        survival = true;
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

  if (nthreads == 0) nthreads = get_nprocs();

  auto grid = make_grid(from, to, npoints, log);

  auto luminosity = survival
                  ? pp_luminosity_b(
                      collision_energy,
                      polarization,
                      grid,
                      verbose ? "luminosity_b " : nullptr
                    )
                  : make_function1d(epa::pp_luminosity(collision_energy), grid);
  luminosity.dump(std::cout);

  return 0;
};
