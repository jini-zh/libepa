#pragma once

#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epa {

static inline double sqr(double x) {
  return x * x;
};

struct Function1d {
  class OutOfBounds: public std::exception {
    public:
      std::string message;

      OutOfBounds(const Function1d&, double x);

      const char* what() const throw ();
  };

  std::vector<std::pair<double, double>> points;

  double operator()(double x) const;
  std::vector<std::pair<double, double>>::const_iterator locate(double x) const;

  static Function1d load(const std::filesystem::path&);
  void save(const std::filesystem::path&) const;
};

}; // namespace epa
