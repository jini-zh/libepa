#include "../common.cpp"

#include <getopt.h>

void usage(const char* argv0) {
  std::cout
    << "This program calculates photon-photon luminosities and the survival "
       "factor in proton-proton collisions\n"
       "Usage: " << argv0 << " options...\n"
       "Allowed options:\n"
       "  -h or --help:             print this message and exit (no argument)\n"
       "  -E or --collision-energy: specify proton-proton collision energy (in GeV). Default is 13e3\n"
       "  -f or --from:             beginning of the invariant mass range (GeV)\n"
       "  -t or --to:               ending of the invariant mass range (GeV)\n"
       "  -n or --npoints:          number of points in the invariant mass range\n"
       "  -s or --step:             points step in the invariant mass range (GeV)\n"
       "  -l or --log:              assume logarithmic scale in the invariant mass range (no argument)\n"
    << "  -v or --verbose:          be verbose while calculating\n"
       "  --survival:               output file for the survival factor\n"
       "  --luminosity:             output file for the luminosity with non-electromagnetic interactions between the protons neglected\n"
       "  --luminosity_b:           output file for the luminosity with non-electromagnetic interactions between the protons respected\n"
       "  --parallel:               output file for the luminosity of photons with parallel polarization and non-electromagnetic interactions between the protons respected\n"
       "  --perpendicular:          output file for the luminosity of photons with perpendicular polarization and non-electromagnetic interactions between the protons respected\n";
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
  const char* survival      = nullptr;
  const char* luminosity    = nullptr;
  const char* luminosity_b  = nullptr;
  const char* parallel      = nullptr;
  const char* perpendicular = nullptr;

  option options[] = {
    { "help",             0, nullptr, 'h' },
    { "collision-energy", 1, nullptr, 'E' },
    { "from",             1, nullptr, 'f' },
    { "to",               1, nullptr, 't' },
    { "npoints",          1, nullptr, 'n' },
    { "step",             1, nullptr, 's' },
    { "log",              0, nullptr, 'l' },
    { "verbose",          0, nullptr, 'v' },
    { "luminosity",       1, nullptr, -2 },
    { "luminosity_b",     1, nullptr, -3 },
    { "parallel",         1, nullptr, -4 },
    { "perpendicular",    1, nullptr, -5 },
    { "survival",         1, nullptr, -6 },
    { 0,                  0, 0,        0 }
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

      case 'v':
        verbose = true;
        break;

      case -2:
        luminosity = optarg;
        break;

      case -3:
        luminosity_b = optarg;
        break;

      case -4:
        parallel = optarg;
        break;

      case -5:
        perpendicular = optarg;
        break;

      case -6:
        survival = optarg;
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

  auto grid = make_grid(from, to, npoints, log);

  if (verbose)
    std::cerr
      << "Calculating photon-photon luminosities in ultraperipheral collisions "
         "of protons with the energy " << collision_energy << " GeV for "
         "invariant masses from " << from << " GeV to " << to << " GeV ("
      << npoints << " points)\n";

  epa::Function1d f_luminosity;
  if (luminosity || survival) {
    if (verbose)
      std::cerr
        << "  with non-electromagnetic interactions between protons "
           "neglected... ";
    f_luminosity = make_function1d(epa::pp_luminosity(collision_energy), grid);
    if (verbose) std::cerr << "done.";
    if (luminosity) {
      f_luminosity.save(luminosity);
      if (verbose)
        std::cerr << " Result is written to `" << luminosity << "'\n";
    };
  };

  epa::Function1d f_luminosity_b;
  if (luminosity_b || survival) {
    if (verbose)
      std::cerr
        << "  with non-electromagnetic interactions between protons "
           "respected:\n";
    f_luminosity_b = pp_luminosity_b(
        collision_energy, {1, 1}, grid, verbose ? "luminosity_b " : nullptr
    );
    if (verbose) std::cerr << "Done.";
    f_luminosity_b.save(luminosity_b);
    if (verbose)
      std::cerr << " Result is written to `" << luminosity_b << "'\n";
  };

  if (survival) {
    if (verbose) std::cerr << "Calculating ratio... ";
    epa::Function1d f_survival;
    f_survival.points.resize(grid.size());
    for (size_t i = 0; i < grid.size(); ++i) {
      f_survival.points[i].first = grid[i].first;
      f_survival.points[i].second = f_luminosity_b.points[i].second
                                  / f_luminosity.points[i].second;
    };
    if (verbose) std::cerr << "done.";
    f_survival.save(survival);
    if (verbose) std::cerr << " Result is written to `" << survival << "'\n";
  };

  if (parallel) {
    if (verbose) std::cerr << "  for parallel polarization:\n";
    pp_luminosity_b(
        collision_energy, {1, 0}, grid, verbose ? "parallel " : nullptr
    ).save(parallel);
    if (verbose)
      std::cerr << "Done. Result is written to `" << parallel << "'\n";
  };

  if (perpendicular) {
    if (verbose) std::cerr << "  for perpendicular polarization:\n";
    pp_luminosity_b(
        collision_energy, {0, 1}, grid, verbose ? "perpendicular " : nullptr
    ).save(perpendicular);
    if (verbose)
      std::cerr << "Done. Result is written to `" << perpendicular << "'\n";
  };

  return 0;
};
