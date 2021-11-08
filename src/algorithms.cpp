#include <algorithm>
#include <fstream>
#include <sstream>

#include "algorithms.hpp"

namespace epa {

Function1d::OutOfBounds::OutOfBounds(const Function1d& f, double x) {
  std::stringstream ss;
  ss
    << x
    << " is out of range ["
    << f.points.front().first
    << ", "
    << f.points.back().first
    << ']';
  message = ss.str();
};

const char* Function1d::OutOfBounds::what() const throw () {
  return message.c_str();
};

Function1d::Function1d(Function1d&& f): points(std::move(f.points)) {};

Function1d::Function1d(std::vector<std::pair<double, double>>&& points):
  points(points)
{};

Function1d& Function1d::operator=(Function1d&& f) {
  points = std::move(f.points);
  return *this;
};

double Function1d::operator()(double x) const {
  auto i = locate(x);
  auto j = i+1;
  if (x < i->first || j == points.end()) throw OutOfBounds(*this, x);
  return i->second
       + (x - i->first) * (j->second - i->second) / (j->first - i->first);
};

std::vector<std::pair<double, double>>::const_iterator
Function1d::locate(double x) const {
  auto result = std::upper_bound(
      points.begin(),
      points.end(),
      x,
      [](double x, const std::pair<double, double>& point) {
        return x < point.first;
      }
  );
  if (result != points.begin() && result != points.end()) --result;
  return result;
};

Function1d Function1d::load(const std::filesystem::path& file) {
  std::ifstream f(file);
  Function1d result;
  while (true) {
    std::pair<double, double> point;
    f >> point.first >> point.second;
    if (!f) return result;
    result.points.push_back(point);
    f.ignore(1, '\n');
  };
};

void Function1d::dump(std::ostream& output) const {
  auto precision = output.precision(12);
  auto flags     = output.setf(std::ios_base::scientific);
  for (auto& point: points)
    output << point.first << ' ' << point.second << '\n';
};

void Function1d::save(const std::filesystem::path& file) const {
  std::ofstream f(file);
  dump(f);
};

}; // namespace epa
