#pragma once

#include "optimization/export.hpp"

#include <atomic>
#include <cstddef>
#include <tuple>
#include <vector>

class OPTIMIZATION_API Minimizer
{
public:
  Minimizer(long long maxiter,
        double rho,
        double tol,
        double dt,
        double dmin_people,
        double dmin_walls,
        bool nesterov);
  ~Minimizer();

  std::tuple<std::vector<double>, std::vector<double>,
             std::vector<double>, long long>
  run(std::vector<double>& xyrv,
      std::vector<double>& contacts,
      std::vector<double>& Vd,
      std::vector<double>& Mass);

  // Zero selects std::thread::hardware_concurrency().
  void set_thread_count(std::size_t thread_count) noexcept;
  std::size_t thread_count() const noexcept;

  void print();

private:
  std::size_t _maxiter;
  double _rho;
  double _tol;
  double _dt;
  double _dmin_people;
  double _dmin_walls;
  bool _nesterov;
  std::atomic<std::size_t> _thread_count{0};
};
