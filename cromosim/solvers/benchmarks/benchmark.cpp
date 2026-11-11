#include "optimization/minimizer.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

namespace
{
struct Options
{
  std::size_t particles = 20000;
  long long iterations = 50;
  std::size_t warmup = 2;
  std::size_t repetitions = 9;
  std::size_t threads = 0;
};

std::size_t parse_size(const char* value, const char* option)
{
  const std::string text(value);
  std::size_t consumed = 0;
  const auto parsed = std::stoull(text, &consumed);
  if (consumed != text.size() || parsed == 0) {
    throw std::invalid_argument(std::string(option) + " expects a positive integer");
  }
  return static_cast<std::size_t>(parsed);
}

Options parse_options(int argc, char** argv)
{
  Options options;
  for (int index = 1; index < argc; index += 2) {
    if (index + 1 == argc) {
      throw std::invalid_argument(std::string(argv[index]) + " requires a value");
    }
    const std::string option(argv[index]);
    if (option == "--particles") options.particles = parse_size(argv[index + 1], argv[index]);
    else if (option == "--iterations") options.iterations = static_cast<long long>(parse_size(argv[index + 1], argv[index]));
    else if (option == "--warmup") options.warmup = parse_size(argv[index + 1], argv[index]);
    else if (option == "--repetitions") options.repetitions = parse_size(argv[index + 1], argv[index]);
    else if (option == "--threads") options.threads = parse_size(argv[index + 1], argv[index]);
    else throw std::invalid_argument("unknown option: " + option);
  }
  if (options.particles < 2) {
    throw std::invalid_argument("--particles must be at least 2");
  }
  return options;
}

struct Data
{
  std::vector<double> xyrv;
  std::vector<double> contacts;
  std::vector<double> velocity;
  std::vector<double> mass;
};

Data make_data(std::size_t particle_count)
{
  Data data;
  data.xyrv.reserve(4 * particle_count);
  data.contacts.reserve(5 * particle_count);
  data.velocity.resize(2 * particle_count);
  data.mass.resize(particle_count);
  for (std::size_t particle = 0; particle < particle_count; ++particle) {
    data.xyrv.insert(data.xyrv.end(), {static_cast<double>(particle), 0.0,
                                      0.45 + 0.01 * static_cast<double>(particle % 3), 1.0});
    data.velocity[2 * particle] = static_cast<double>(particle % 5) * 0.001;
    data.velocity[2 * particle + 1] = 0.0;
    data.mass[particle] = 1.0 + 0.1 * static_cast<double>(particle % 4);
    data.contacts.insert(data.contacts.end(), {static_cast<double>(particle), -1.0,
                                               0.4, 1.0, 0.0});
  }
  return data;
}

using Result = std::tuple<std::vector<double>, std::vector<double>,
                          std::vector<double>, long long>;

struct Measurement
{
  std::vector<double> milliseconds;
  Result result;
};

Measurement measure(const Options& options, const Data& original, std::size_t threads)
{
  Minimizer solver(options.iterations, 0.001, 0.0, 0.1, 1.0, 0.5, false);
  solver.set_thread_count(threads);
  Measurement measurement;
  for (std::size_t repetition = 0; repetition < options.warmup + options.repetitions; ++repetition) {
    auto data = original;
    const auto start = std::chrono::steady_clock::now();
    measurement.result = solver.run(data.xyrv, data.contacts, data.velocity, data.mass);
    const auto stop = std::chrono::steady_clock::now();
    if (repetition >= options.warmup) {
      measurement.milliseconds.push_back(
        std::chrono::duration<double, std::milli>(stop - start).count());
    }
  }
  std::sort(measurement.milliseconds.begin(), measurement.milliseconds.end());
  return measurement;
}

double percentile(const std::vector<double>& sorted, double fraction)
{
  const double position = fraction * static_cast<double>(sorted.size() - 1);
  const auto lower = static_cast<std::size_t>(std::floor(position));
  const auto upper = static_cast<std::size_t>(std::ceil(position));
  const double weight = position - static_cast<double>(lower);
  return sorted[lower] * (1.0 - weight) + sorted[upper] * weight;
}

double checksum(const Result& result)
{
  double sum = 0.0;
  for (double value : std::get<0>(result)) sum += value;
  for (double value : std::get<1>(result)) sum += value;
  for (double value : std::get<2>(result)) sum += value;
  return sum;
}

bool equivalent(const Result& left, const Result& right)
{
  for (std::size_t part = 0; part < 3; ++part) {
    const auto& a = part == 0 ? std::get<0>(left) : part == 1 ? std::get<1>(left) : std::get<2>(left);
    const auto& b = part == 0 ? std::get<0>(right) : part == 1 ? std::get<1>(right) : std::get<2>(right);
    if (a.size() != b.size()) return false;
    for (std::size_t index = 0; index < a.size(); ++index) {
      const double tolerance = 1e-12 * (1.0 + std::abs(a[index]));
      if (std::abs(a[index] - b[index]) > tolerance) return false;
    }
  }
  return std::get<3>(left) == std::get<3>(right);
}

void print_measurement(const char* label, std::size_t threads, const Measurement& measurement)
{
  std::cout << label << ",threads=" << threads
            << ",median_ms=" << percentile(measurement.milliseconds, 0.5)
            << ",p25_ms=" << percentile(measurement.milliseconds, 0.25)
            << ",p75_ms=" << percentile(measurement.milliseconds, 0.75)
            << ",checksum=" << checksum(measurement.result) << '\n';
}
} // namespace

int main(int argc, char** argv)
{
  try {
    const auto options = parse_options(argc, argv);
    const auto hardware = std::max(1u, std::thread::hardware_concurrency());
    const auto parallel_threads = options.threads == 0
      ? std::min<std::size_t>(4, hardware) : options.threads;
    const auto data = make_data(options.particles);
    std::cout << std::fixed << std::setprecision(3)
              << "workload,particles=" << options.particles
              << ",contacts=" << options.particles
              << ",max_iterations=" << options.iterations
              << ",warmup=" << options.warmup
              << ",repetitions=" << options.repetitions << '\n';
    const auto serial = measure(options, data, 1);
    const auto parallel = measure(options, data, parallel_threads);
    if (!equivalent(serial.result, parallel.result)) {
      std::cerr << "serial and parallel outputs differ\n";
      return 2;
    }
    print_measurement("serial", 1, serial);
    print_measurement("parallel", parallel_threads, parallel);
    const double serial_median = percentile(serial.milliseconds, 0.5);
    const double parallel_median = percentile(parallel.milliseconds, 0.5);
    std::cout << "speedup=" << serial_median / parallel_median << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "benchmark error: " << error.what() << '\n';
    return 1;
  }
}
