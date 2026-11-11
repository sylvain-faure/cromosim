#include "optimization/minimizer.hpp"

#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <functional>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <iostream>
#include <thread>
#include <utility>

namespace
{
constexpr std::size_t particle_stride = 4;
constexpr std::size_t contact_stride = 5;
constexpr std::size_t dimension = 2;
constexpr double pi = 3.141592653589793238462643383279502884;

class ParallelExecutor
{
public:
  ParallelExecutor(std::size_t requested_threads, std::size_t maximum_work)
  {
    const auto hardware = static_cast<std::size_t>(std::thread::hardware_concurrency());
    const auto desired = requested_threads == 0 ? std::max<std::size_t>(hardware, 1) : requested_threads;
    thread_count_ = std::max<std::size_t>(1, std::min(desired, std::max<std::size_t>(maximum_work, 1)));
    workers_.reserve(thread_count_ - 1);
    for (std::size_t id = 1; id < thread_count_; ++id) {
      workers_.emplace_back([this, id] { worker_loop(id); });
    }
  }

  ~ParallelExecutor()
  {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stopping_ = true;
      ++generation_;
    }
    work_ready_.notify_all();
    for (auto& worker : workers_) {
      worker.join();
    }
  }

  ParallelExecutor(const ParallelExecutor&) = delete;
  ParallelExecutor& operator=(const ParallelExecutor&) = delete;

  template<class Function>
  void for_each_block(std::size_t count, Function&& function)
  {
    if (count == 0) {
      return;
    }
    if (workers_.empty()) {
      function(0, count);
      return;
    }

    {
      std::lock_guard<std::mutex> lock(mutex_);
      count_ = count;
      task_ = std::forward<Function>(function);
      remaining_ = workers_.size();
      error_ = nullptr;
      ++generation_;
    }
    work_ready_.notify_all();

    try {
      function(0, block_end(0, count));
    } catch (...) {
      std::lock_guard<std::mutex> lock(mutex_);
      error_ = std::current_exception();
    }

    std::unique_lock<std::mutex> lock(mutex_);
    work_done_.wait(lock, [this] { return remaining_ == 0; });
    if (error_) {
      std::rethrow_exception(error_);
    }
  }

private:
  std::size_t block_begin(std::size_t id, std::size_t count) const noexcept
  {
    return count * id / thread_count_;
  }

  std::size_t block_end(std::size_t id, std::size_t count) const noexcept
  {
    return count * (id + 1) / thread_count_;
  }

  void worker_loop(std::size_t id)
  {
    std::size_t observed_generation = 0;
    for (;;) {
      std::function<void(std::size_t, std::size_t)> task;
      std::size_t count = 0;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        work_ready_.wait(lock, [this, observed_generation] {
          return stopping_ || generation_ != observed_generation;
        });
        if (stopping_) {
          return;
        }
        observed_generation = generation_;
        task = task_;
        count = count_;
      }

      try {
        task(block_begin(id, count), block_end(id, count));
      } catch (...) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!error_) {
          error_ = std::current_exception();
        }
      }

      {
        std::lock_guard<std::mutex> lock(mutex_);
        --remaining_;
      }
      work_done_.notify_one();
    }
  }

  std::size_t thread_count_ = 1;
  std::vector<std::thread> workers_;
  std::mutex mutex_;
  std::condition_variable work_ready_;
  std::condition_variable work_done_;
  std::function<void(std::size_t, std::size_t)> task_;
  std::exception_ptr error_;
  std::size_t generation_ = 0;
  std::size_t remaining_ = 0;
  std::size_t count_ = 0;
  bool stopping_ = false;
};

struct CsrMatrix
{
  std::size_t rows = 0;
  std::size_t columns = 0;
  std::vector<std::size_t> row_offsets;
  std::vector<std::size_t> column_indices;
  std::vector<double> values;
};

CsrMatrix transpose(const CsrMatrix& source)
{
  CsrMatrix result;
  result.rows = source.columns;
  result.columns = source.rows;
  result.row_offsets.assign(result.rows + 1, 0);
  result.column_indices.resize(source.column_indices.size());
  result.values.resize(source.values.size());

  for (const auto column : source.column_indices) {
    ++result.row_offsets[column + 1];
  }
  std::partial_sum(result.row_offsets.begin(), result.row_offsets.end(), result.row_offsets.begin());
  auto next = result.row_offsets;
  for (std::size_t row = 0; row < source.rows; ++row) {
    for (std::size_t offset = source.row_offsets[row]; offset < source.row_offsets[row + 1]; ++offset) {
      const auto target = next[source.column_indices[offset]]++;
      result.column_indices[target] = row;
      result.values[target] = source.values[offset];
    }
  }
  return result;
}

void multiply(const CsrMatrix& matrix,
              const std::vector<double>& input,
              std::vector<double>& output,
              ParallelExecutor& executor)
{
  executor.for_each_block(matrix.rows, [&](std::size_t begin, std::size_t end) {
    for (std::size_t row = begin; row < end; ++row) {
      double sum = 0.0;
      for (std::size_t offset = matrix.row_offsets[row]; offset < matrix.row_offsets[row + 1]; ++offset) {
        sum += matrix.values[offset] * input[matrix.column_indices[offset]];
      }
      output[row] = sum;
    }
  });
}

bool is_finite(double value) noexcept
{
  return std::isfinite(value) != 0;
}

void require_finite(const std::vector<double>& values, const char* name)
{
  if (!std::all_of(values.begin(), values.end(), is_finite)) {
    throw std::invalid_argument(std::string(name) + " must contain only finite values");
  }
}

std::size_t contact_index(double value, std::size_t particle_count, bool wall_allowed)
{
  if (!is_finite(value) || std::trunc(value) != value) {
    throw std::invalid_argument("contact indices must be finite integers");
  }
  if (wall_allowed && value == -1.0) {
    return particle_count;
  }
  if (value < 0.0 || value >= static_cast<double>(particle_count)) {
    throw std::invalid_argument("contact index is outside the particle array");
  }
  return static_cast<std::size_t>(value);
}

double dot(const std::vector<double>& left, const std::vector<double>& right)
{
  return std::inner_product(left.begin(), left.end(), right.begin(), 0.0);
}

double maximum(const std::vector<double>& values)
{
  return *std::max_element(values.begin(), values.end());
}
} // namespace

Minimizer::Minimizer(long long maxiter,
             double rho,
             double tol,
             double dt,
             double dmin_people,
             double dmin_walls,
             bool nesterov)
  : _maxiter(0),
    _rho(rho),
    _tol(tol),
    _dt(dt),
    _dmin_people(dmin_people),
    _dmin_walls(dmin_walls),
    _nesterov(nesterov)
{
  if (maxiter <= 0) {
    throw std::invalid_argument("maxiter must be positive");
  }
  if (!is_finite(rho) || rho <= 0.0 || !is_finite(tol) || tol < 0.0 ||
      !is_finite(dt) || dt <= 0.0 || !is_finite(dmin_people) || dmin_people < 0.0 ||
      !is_finite(dmin_walls) || dmin_walls < 0.0) {
    throw std::invalid_argument("solver parameters are outside their valid range");
  }
  _maxiter = static_cast<std::size_t>(maxiter);
}

Minimizer::~Minimizer() = default;

void Minimizer::set_thread_count(std::size_t thread_count) noexcept
{
  _thread_count.store(thread_count, std::memory_order_relaxed);
}

std::size_t Minimizer::thread_count() const noexcept
{
  return _thread_count.load(std::memory_order_relaxed);
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, long long>
Minimizer::run(std::vector<double>& xyrv,
           std::vector<double>& contacts,
           std::vector<double>& Vd,
           std::vector<double>& Mass)
{
  // if (xyrv.size() % particle_stride != 0) {
  //   throw std::invalid_argument("xyrv size must be a multiple of 4");
  // }
  // if (contacts.size() % contact_stride != 0) {
  //   throw std::invalid_argument("contacts size must be a multiple of 5");
  // }
  const auto particle_count = xyrv.size() / particle_stride;
  const auto contact_count = contacts.size() / contact_stride;
  // if (Vd.size() != dimension * particle_count) {
  //   throw std::invalid_argument("Vd must contain two values per particle");
  // }
  // if (Mass.size() != particle_count) {
  //   throw std::invalid_argument("Mass must contain one value per particle");
  // }
  // require_finite(xyrv, "xyrv");
  // require_finite(contacts, "contacts");
  // require_finite(Vd, "Vd");
  // require_finite(Mass, "Mass");
  // for (std::size_t particle = 0; particle < particle_count; ++particle) {
  //   if (xyrv[particle * particle_stride + 2] <= 0.0) {
  //     throw std::invalid_argument("particle radii must be positive");
  //   }
  //   if (Mass[particle] <= 0.0) {
  //     throw std::invalid_argument("particle masses must be positive");
  //   }
  // }

  std::vector<std::size_t> first_particles(contact_count);
  std::vector<std::size_t> second_particles(contact_count);
  for (std::size_t contact = 0; contact < contact_count; ++contact) {
    const auto base = contact * contact_stride;
    first_particles[contact] = contact_index(contacts[base], particle_count, false);
    second_particles[contact] = contact_index(contacts[base + 1], particle_count, true);
    const double normal_squared = contacts[base + 3] * contacts[base + 3] +
                                  contacts[base + 4] * contacts[base + 4];
    if (normal_squared <= std::numeric_limits<double>::epsilon()) {
      throw std::invalid_argument("contact normals must be non-zero");
    }
  }

  std::vector<double> pressure(particle_count, 0.0);
  std::vector<double> lambda(contact_count, 0.0);
  if (contact_count == 0) {
    return {Vd, std::move(pressure), std::move(lambda), 0};
  }

  CsrMatrix B;
  B.rows = contact_count;
  B.columns = dimension * particle_count;
  B.row_offsets.reserve(contact_count + 1);
  B.column_indices.reserve(4 * contact_count);
  B.values.reserve(4 * contact_count);
  B.row_offsets.push_back(0);
  for (std::size_t contact = 0; contact < contact_count; ++contact) {
    const auto base = contact * contact_stride;
    const auto first = first_particles[contact];
    const auto second = second_particles[contact];
    const auto wall = second == particle_count;
    const double scale = _nesterov ? _dt : 1.0;
    const double first_sign = _nesterov ? -1.0 : 1.0;
    B.column_indices.push_back(dimension * first);
    B.values.push_back(first_sign * scale * contacts[base + 3]);
    B.column_indices.push_back(dimension * first + 1);
    B.values.push_back(first_sign * scale * contacts[base + 4]);
    if (!wall) {
      B.column_indices.push_back(dimension * second);
      B.values.push_back(-first_sign * scale * contacts[base + 3]);
      B.column_indices.push_back(dimension * second + 1);
      B.values.push_back(-first_sign * scale * contacts[base + 4]);
    } else if (!_nesterov) {
      B.values[B.values.size() - 2] *= -1.0;
      B.values[B.values.size() - 1] *= -1.0;
    } else {
      B.values[B.values.size() - 2] *= -1.0;
      B.values[B.values.size() - 1] *= -1.0;
    }
    B.row_offsets.push_back(B.values.size());
  }
  const auto Bt = transpose(B);
  ParallelExecutor executor(_thread_count.load(std::memory_order_relaxed),
                            std::max(contact_count, dimension * particle_count));

  std::vector<double> distances(contact_count);
  std::vector<double> minimum_distances(contact_count);
  double minimum_radius = std::numeric_limits<double>::infinity();
  for (std::size_t particle = 0; particle < particle_count; ++particle) {
    minimum_radius = std::min(minimum_radius, xyrv[particle * particle_stride + 2]);
  }
  for (std::size_t contact = 0; contact < contact_count; ++contact) {
    distances[contact] = contacts[contact * contact_stride + 2];
    minimum_distances[contact] = second_particles[contact] == particle_count
                                   ? _dmin_walls : _dmin_people;
  }

  std::vector<double> velocity(Vd.size(), 0.0);
  std::vector<double> residual(contact_count, 0.0);
  long long info = -1;

  if (!_nesterov) {
    std::vector<double> transpose_product(Vd.size(), 0.0);
    double residual_max = std::numeric_limits<double>::infinity();
    std::size_t iteration = 0;
    while (_dt * residual_max > _tol * 2.0 * minimum_radius && iteration < _maxiter) {
      multiply(Bt, lambda, transpose_product, executor);
      executor.for_each_block(velocity.size(), [&](std::size_t begin, std::size_t end) {
        for (std::size_t index = begin; index < end; ++index) {
          velocity[index] = Vd[index] - transpose_product[index] / Mass[index / dimension];
        }
      });
      multiply(B, velocity, residual, executor);
      executor.for_each_block(contact_count, [&](std::size_t begin, std::size_t end) {
        for (std::size_t contact = begin; contact < end; ++contact) {
          residual[contact] -= (distances[contact] - minimum_distances[contact]) / _dt;
          lambda[contact] = std::max(lambda[contact] + _rho * residual[contact], 0.0);
        }
      });
      residual_max = maximum(residual);
      if (!is_finite(residual_max)) {
        throw std::runtime_error("optimization iteration produced a non-finite residual");
      }
      ++iteration;
    }
    if (_dt * residual_max <= _tol * 2.0 * minimum_radius) {
      std::cout << "Number of iterations = " << iteration << std::endl;             
      info = static_cast<long long>(iteration);
    }
  } else {
    std::vector<double> temporary_velocity(Vd.size(), 0.0);
    std::vector<double> A_result(contact_count, 0.0);
    const auto apply_A = [&](const std::vector<double>& input, std::vector<double>& output) {
      multiply(Bt, input, temporary_velocity, executor);
      executor.for_each_block(temporary_velocity.size(), [&](std::size_t begin, std::size_t end) {
        for (std::size_t index = begin; index < end; ++index) {
          temporary_velocity[index] /= Mass[index / dimension];
        }
      });
      multiply(B, temporary_velocity, output, executor);
    };

    std::vector<double> C(contact_count, 0.0);
    multiply(B, Vd, C, executor);
    for (std::size_t contact = 0; contact < contact_count; ++contact) {
      C[contact] += distances[contact] - minimum_distances[contact];
    }

    std::vector<double> probe(contact_count, -1.0);
    apply_A(probe, A_result);
    const double probe_norm = std::sqrt(dot(probe, probe));
    double lipschitz = std::sqrt(dot(A_result, A_result)) / probe_norm;
    if (!is_finite(lipschitz) || lipschitz <= std::numeric_limits<double>::epsilon()) {
      throw std::runtime_error("cannot estimate a positive Nesterov Lipschitz constant");
    }
    double step = 1.0 / lipschitz;
    double theta_previous = 1.0;
    std::vector<double> previous_lambda(lambda);
    std::vector<double> y(lambda);
    std::vector<double> gradient(contact_count, 0.0);
    std::vector<double> Ay(contact_count, 0.0);
    std::vector<double> difference(contact_count, 0.0);

    std::size_t iteration = 0;
    bool converged = false;
    while (iteration < _maxiter) {
      apply_A(y, gradient);
      for (std::size_t contact = 0; contact < contact_count; ++contact) {
        gradient[contact] += C[contact];
        lambda[contact] = std::max(y[contact] - step * gradient[contact], 0.0);
      }

      apply_A(y, Ay);
      const double gy = 0.5 * dot(y, Ay) + dot(y, C);
      for (std::size_t line_search = 0; line_search < 10; ++line_search) {
        apply_A(lambda, A_result);
        for (std::size_t contact = 0; contact < contact_count; ++contact) {
          difference[contact] = lambda[contact] - y[contact];
        }
        const double objective = 0.5 * dot(lambda, A_result) + dot(lambda, C);
        const double bound = gy + dot(gradient, difference) +
                             0.5 * lipschitz * dot(difference, difference);
        if (objective <= bound) {
          break;
        }
        lipschitz *= 2.0;
        step = 1.0 / lipschitz;
        for (std::size_t contact = 0; contact < contact_count; ++contact) {
          lambda[contact] = std::max(y[contact] - step * gradient[contact], 0.0);
        }
      }

      const double theta = 0.5 * (-theta_previous * theta_previous +
                           std::sqrt(theta_previous * theta_previous *
                                     (theta_previous * theta_previous + 4.0)));
      const double beta = theta_previous * (1.0 - theta_previous) /
                          (theta_previous * theta_previous + theta);
      for (std::size_t contact = 0; contact < contact_count; ++contact) {
        y[contact] = lambda[contact] + beta * (lambda[contact] - previous_lambda[contact]);
      }
      double next_theta = theta;
      for (std::size_t contact = 0; contact < contact_count; ++contact) {
        difference[contact] = lambda[contact] - previous_lambda[contact];
      }
      if (dot(gradient, difference) > 0.0) {
        y = lambda;
        next_theta = 1.0;
      }
      previous_lambda = lambda;
      theta_previous = next_theta;
      lipschitz *= 0.9;
      step = 1.0 / lipschitz;
      ++iteration;

      apply_A(lambda, A_result);
      for (std::size_t contact = 0; contact < contact_count; ++contact) {
        residual[contact] = -A_result[contact] - C[contact];
      }
      const double residual_max = maximum(residual);
      if (!is_finite(residual_max)) {
        throw std::runtime_error("Nesterov iteration produced a non-finite residual");
      }
      if (residual_max <= _tol * 2.0 * minimum_radius) {
        std::cout << "Number of iterations = " << iteration << std::endl;
        converged = true;
        info = static_cast<long long>(iteration);
        break;
      }
    }

    multiply(Bt, lambda, temporary_velocity, executor);
    executor.for_each_block(velocity.size(), [&](std::size_t begin, std::size_t end) {
      for (std::size_t index = begin; index < end; ++index) {
        velocity[index] = Vd[index] + temporary_velocity[index] / Mass[index / dimension];
      }
    });
    if (!converged) {
      info = -1;
    }
  }

  for (std::size_t contact = 0; contact < contact_count; ++contact) {
    const auto first = first_particles[contact];
    const auto second = second_particles[contact];
    const double first_radius = xyrv[first * particle_stride + 2];
    pressure[first] += 3.0 * lambda[contact] / (4.0 * pi * first_radius * first_radius);
    if (second != particle_count) {
      const double second_radius = xyrv[second * particle_stride + 2];
      pressure[second] += 3.0 * lambda[contact] / (4.0 * pi * second_radius * second_radius);
    }
  }
  return {std::move(velocity), std::move(pressure), std::move(lambda), info};
}

void Minimizer::print()
{
  // Retained for source compatibility. Solver workspaces are now local to run().
}
