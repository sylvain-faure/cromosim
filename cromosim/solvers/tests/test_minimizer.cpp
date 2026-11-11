#include "optimization/minimizer.hpp"

#include <atomic>
#include <cmath>
#include <future>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
int failures = 0;

#define CHECK(condition) check((condition), #condition, __FILE__, __LINE__)

void check(bool condition, const char* expression, const char* file, int line)
{
  if (!condition) {
    std::cerr << file << ':' << line << ": CHECK failed: " << expression << '\n';
    ++failures;
  }
}

void check_near(double actual, double expected, double tolerance,
                const char* file, int line)
{
  if (std::abs(actual - expected) > tolerance) {
    std::cerr << file << ':' << line << ": expected " << expected
              << ", got " << actual << " (tolerance " << tolerance << ")\n";
    ++failures;
  }
}

#define CHECK_NEAR(actual, expected, tolerance) \
  check_near((actual), (expected), (tolerance), __FILE__, __LINE__)

template<class Function>
void expect_invalid_argument(Function&& function, const char* file, int line)
{
  try {
    function();
    std::cerr << file << ':' << line << ": expected std::invalid_argument\n";
    ++failures;
  } catch (const std::invalid_argument&) {
  } catch (...) {
    std::cerr << file << ':' << line << ": expected std::invalid_argument, got another exception\n";
    ++failures;
  }
}

#define EXPECT_INVALID(expression) \
  expect_invalid_argument([&] { expression; }, __FILE__, __LINE__)

using Result = std::tuple<std::vector<double>, std::vector<double>,
                          std::vector<double>, long long>;

Result solve(Minimizer& solver,
             std::vector<double> xyrv,
             std::vector<double> contacts,
             std::vector<double> velocity,
             std::vector<double> mass)
{
  return solver.run(xyrv, contacts, velocity, mass);
}

void check_vectors(const std::vector<double>& left,
                   const std::vector<double>& right,
                   double tolerance)
{
  CHECK(left.size() == right.size());
  if (left.size() != right.size()) {
    return;
  }
  for (std::size_t index = 0; index < left.size(); ++index) {
    CHECK_NEAR(left[index], right[index], tolerance);
  }
}

void test_no_contacts()
{
  Minimizer solver(100, 0.5, 1e-10, 0.1, 1.0, 0.5, false);
  auto result = solve(solver, {0.0, 0.0, 0.5, 1.0}, {}, {2.0, -3.0}, {4.0});
  check_vectors(std::get<0>(result), {2.0, -3.0}, 0.0);
  check_vectors(std::get<1>(result), {0.0}, 0.0);
  CHECK(std::get<2>(result).empty());
  CHECK(std::get<3>(result) == 0);
}

void test_wall_contact(bool nesterov)
{
  Minimizer solver(1000, 0.5, 1e-9, 0.1, 1.0, 0.5, nesterov);
  auto result = solve(solver,
                      {0.0, 0.0, 0.5, 1.0},
                      {0.0, -1.0, 0.4, 1.0, 0.0},
                      {0.0, 0.0}, {1.0});
  CHECK(std::get<3>(result) > 0);
  CHECK_NEAR(std::get<0>(result)[0], 1.0, 2e-7);
  CHECK_NEAR(std::get<0>(result)[1], 0.0, 1e-12);
  const double expected_lambda = nesterov ? 10.0 : 1.0;
  CHECK_NEAR(std::get<2>(result)[0], expected_lambda, 2e-6);
  CHECK_NEAR(std::get<1>(result)[0], 3.0 * expected_lambda / 3.14159265358979323846, 3e-6);
}

void test_people_contact()
{
  Minimizer solver(1000, 0.25, 1e-9, 0.1, 1.0, 0.5, false);
  auto result = solve(solver,
                      {0.0, 0.0, 0.5, 1.0, 0.8, 0.0, 0.5, 1.0},
                      {0.0, 1.0, 0.8, 1.0, 0.0},
                      {0.0, 0.0, 0.0, 0.0}, {1.0, 1.0});
  CHECK(std::get<3>(result) > 0);
  check_vectors(std::get<0>(result), {-1.0, 0.0, 1.0, 0.0}, 2e-7);
  CHECK_NEAR(std::get<2>(result)[0], 1.0, 2e-6);
  CHECK_NEAR(std::get<1>(result)[0], 3.0 / 3.14159265358979323846, 3e-6);
  CHECK_NEAR(std::get<1>(result)[1], 3.0 / 3.14159265358979323846, 3e-6);
}

void test_inactive_and_iteration_limit()
{
  Minimizer inactive(10, 0.5, 1e-9, 0.1, 1.0, 0.5, false);
  auto inactive_result = solve(inactive,
                               {0.0, 0.0, 0.5, 1.0},
                               {0.0, -1.0, 0.6, 1.0, 0.0},
                               {0.0, 0.0}, {1.0});
  check_vectors(std::get<0>(inactive_result), {0.0, 0.0}, 0.0);
  check_vectors(std::get<1>(inactive_result), {0.0}, 0.0);
  check_vectors(std::get<2>(inactive_result), {0.0}, 0.0);
  CHECK(std::get<3>(inactive_result) > 0);

  Minimizer limited(1, 1e-6, 0.0, 0.1, 1.0, 0.5, false);
  auto limited_result = solve(limited,
                              {0.0, 0.0, 0.5, 1.0},
                              {0.0, -1.0, 0.4, 1.0, 0.0},
                              {0.0, 0.0}, {1.0});
  CHECK(std::get<3>(limited_result) == -1);
}

void test_validation()
{
  Minimizer solver(10, 0.5, 1e-6, 0.1, 1.0, 0.5, false);
  EXPECT_INVALID(solve(solver, {0.0, 0.0}, {}, {}, {}));
  EXPECT_INVALID(solve(solver, {0.0, 0.0, 0.5, 1.0}, {0.0}, {0.0, 0.0}, {1.0}));
  EXPECT_INVALID(solve(solver, {0.0, 0.0, 0.5, 1.0}, {}, {0.0}, {1.0}));
  EXPECT_INVALID(solve(solver, {0.0, 0.0, 0.5, 1.0}, {}, {0.0, 0.0}, {0.0}));
  EXPECT_INVALID(solve(solver, {0.0, 0.0, 0.5, 1.0},
                       {0.5, -1.0, 0.4, 1.0, 0.0}, {0.0, 0.0}, {1.0}));
  EXPECT_INVALID(solve(solver, {0.0, 0.0, 0.5, 1.0},
                       {0.0, -2.0, 0.4, 1.0, 0.0}, {0.0, 0.0}, {1.0}));
  auto nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_INVALID(solve(solver, {0.0, 0.0, 0.5, 1.0}, {}, {nan, 0.0}, {1.0}));
  EXPECT_INVALID(Minimizer(0, 0.5, 1e-6, 0.1, 1.0, 0.5, false));
}

void make_star(std::size_t particle_count,
               std::vector<double>& xyrv,
               std::vector<double>& contacts,
               std::vector<double>& velocity,
               std::vector<double>& mass)
{
  xyrv.reserve(4 * particle_count);
  velocity.assign(2 * particle_count, 0.0);
  mass.assign(particle_count, 1.0);
  for (std::size_t particle = 0; particle < particle_count; ++particle) {
    xyrv.insert(xyrv.end(), {static_cast<double>(particle), 0.0,
                             particle == particle_count - 1 ? 0.4 : 0.5, 1.0});
  }
  for (std::size_t particle = 1; particle < particle_count; ++particle) {
    contacts.insert(contacts.end(), {0.0, static_cast<double>(particle),
                                     0.9, 1.0, 0.0});
  }
}

void test_parallel_and_reentrant()
{
  std::vector<double> xyrv, contacts, velocity, mass;
  make_star(128, xyrv, contacts, velocity, mass);
  Minimizer solver(80, 0.002, 1e-7, 0.1, 1.0, 0.5, false);
  solver.set_thread_count(1);
  auto serial = solve(solver, xyrv, contacts, velocity, mass);
  solver.set_thread_count(4);
  auto parallel = solve(solver, xyrv, contacts, velocity, mass);
  check_vectors(std::get<0>(serial), std::get<0>(parallel), 0.0);
  check_vectors(std::get<1>(serial), std::get<1>(parallel), 0.0);
  check_vectors(std::get<2>(serial), std::get<2>(parallel), 0.0);
  CHECK(std::get<3>(serial) == std::get<3>(parallel));

  auto first = std::async(std::launch::async, [&] { return solve(solver, xyrv, contacts, velocity, mass); });
  auto second = std::async(std::launch::async, [&] { return solve(solver, xyrv, contacts, velocity, mass); });
  check_vectors(std::get<0>(first.get()), std::get<0>(second.get()), 0.0);
}
} // namespace

int main()
{
  test_no_contacts();
  test_wall_contact(false);
  test_wall_contact(true);
  test_people_contact();
  test_inactive_and_iteration_limit();
  test_validation();
  test_parallel_and_reentrant();
  if (failures != 0) {
    std::cerr << failures << " test assertion(s) failed\n";
    return 1;
  }
  std::cout << "All optimization unit tests passed\n";
  return 0;
}
