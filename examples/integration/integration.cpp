#include "../common.cpp"

#include <getopt.h>

using epa::sqr;

epa::Integrator integrator(unsigned level) {
  std::cout << "integrator " << static_cast<int>(level) << '\n';
  return epa::default_integrator(level);
};

void usage(const char* argv0) {
  std::cout
    << "This program is an example of multiple integral calculation with the help of libepa.\n"
       "It calculates the following integral in a given range of parameter a:\n"
       "        a  sqrt(1 - (x/a)^2)  sqrt(1 - (x/a)^2 - y^2)\n"
       "        /         /                      /          x * y * z         x * (x + 1/2)\n"
       " I(a) = | dx      | dy                   | dz --------------------- = -------------\n"
       "        /         /                      /    sqrt(x^2 + y^2 + z^2)     (x + 1)^2\n"
       "        0         0                      0\n"
       "Usage: " << argv0 << " options...\n"
       "Allowed options:\n"
       "  -h or --help:             print this message and exit (no argument)\n"
       "  -f or --from:             beginning of the parameter range\n"
       "  -t or --to:               ending of the parameter range\n"
       "  -n or --npoints:          number of points in the parameter range\n"
       "  -s or --step:             points step in the parameter range\n"
       "  -l or --log:              assume logarithmic scale in the parameter range (no argument)\n"
       "  -v or --verbose:          be verbose while calculating\n"
  ;
};

int main(int argc, char** argv) {

  gsl::init();
  epa::default_relative_error = 1e-2;

  // Parse parameters
  
  option options[] = {
    { "help",             0, nullptr, 'h' },
    { "from",             1, nullptr, 'f' },
    { "to",               1, nullptr, 't' },
    { "npoints",          1, nullptr, 'n' },
    { "step",             1, nullptr, 's' },
    { "log",              0, nullptr, 'l' },
    { "verbose",          0, nullptr, 'v' },
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

  // Define integral

  auto integral = []() -> std::function<double (double)> {
    return [
      integrate_x = epa::default_integrator(0),
      integrate_y = epa::default_integrator(1),
      integrate_z = epa::default_integrator(2)
    ](double a) -> double {
      return (15 / a) * integrate_x([=](double x) -> double {
        return x * integrate_y([=](double y) -> double {
          return y * integrate_z([=](double z) -> double {
            return z / sqrt(x*x + y*y + z*z);},
            0, sqrt(1 - pow(x/a, 2) - y*y));},
          0, sqrt(1 - pow(x/a, 2)));},
        0, a);
    };
  };

  // Calculate for the requested set of parameters
  
  auto grid = make_grid(from, to, npoints, log);
  auto values = make_function1d_async(integral, grid, verbose ? std::string("I(a): ") : std::string());
  values.dump(std::cout);

  return 0;
};
