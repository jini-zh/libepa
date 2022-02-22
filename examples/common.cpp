#include <epa/proton.hpp>

#include <iostream>
#include <list>
#include <mutex>
#include <thread>

#include <string.h>

#include <sys/sysinfo.h>

std::vector<std::pair<double, double>> make_grid(
    double   from,
    double   to,
    unsigned n,
    bool     log = false
) {
  double x = from;
  double step = log ? pow(to / from, 1. / n) : (to - from) / n;
  std::vector<std::pair<double, double>> grid(n);
  --n;
  for (size_t i = 0; i < n; ++i) {
    grid[i].first = x;
    x = log ? x * step : x + step;
  };
  grid[n].first = to;
  return grid;
};

epa::Function1d make_function1d(
    const std::function<double (double)>& f,
    std::vector<std::pair<double, double>> grid
) {
  epa::Function1d F(std::move(grid));
  auto& points = *F.points;
  for (size_t i = 0; i < points.size(); ++i)
    points[i].second = f(points[i].first);
  return F;
};

epa::Function1d make_function1d_async(
    const std::function<std::function<double (double)> ()>& generator,
    std::vector<std::pair<double, double>> grid,
    const std::string& verbose = std::string()
) {
  epa::Function1d F(std::move(grid));
  std::mutex m_queue;
  std::mutex m_output;
  std::list<std::thread> threads;
  size_t i = 0;
  int nprocs = get_nprocs();
  if (!verbose.empty()) std::cerr << "Using " << nprocs << " threads.\n";
  auto& points = *F.points;
  for (int t = 0; t < nprocs; ++t)
    threads.emplace_back(
        [&]() {
          auto f = generator();
          while (true) {
            size_t j;
            {
              std::unique_lock lock(m_queue);
              j = i++;
            };
            if (j >= points.size()) break;
            points[j].second = f(points[j].first);
            if (!verbose.empty()) {
              std::unique_lock lock(m_output);
              std::cerr
                << verbose
                << points[j].first << " => " << points[j].second
                << '\n';
            };
          };
        }
    );

  for (auto& t: threads) t.join();
  return F;
};

double parse_energy(const char* arg, const char* value) {
  char* end;
  errno = 0;
  double result = strtod(arg, &end);
  if (!*end && !errno && result > 0) return result;
  std::cerr << "Invalid " << value << ": " << arg;
  if (errno) std::cerr << " (" << strerror(errno) << ')';
  else if (result <= 0) std::cerr << " (not positive)";
  std::cerr << '\n';
  std::exit(1);
};

unsigned parse_unsigned(const char* arg, const char* value) {
  size_t end;
  int result = std::stoi(arg, &end);
  if (!arg[end] && result > 0) return result;
  std::cerr << "Invalid " << value << ": " << arg;
  if (result < 0) std::cerr << " (not positive)";
  return result;
};
