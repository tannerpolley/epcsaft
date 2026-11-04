#include "regression/fit.hpp"

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"
#include "equilibrium/stoichiometry.hpp"

#include <ceres/ceres.h>
#include <Eigen/Dense>
#ifndef _WIN32
#include <dlfcn.h>
#include <fcntl.h>
#include <poll.h>
#include <signal.h>
#include <sys/resource.h>
#ifndef __APPLE__
#include <sys/syscall.h>
#include <sys/prctl.h>
#endif
#include <sys/socket.h>
#include <sys/wait.h>
#include <unistd.h>
#ifdef __APPLE__
#include <crt_externs.h>
#include <spawn.h>
#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0  // SO_NOSIGPIPE on each worker socket
#endif
#endif
#else
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0A00
#endif
#define NOMINMAX
#include <windows.h>
#include <sddl.h>
#include <atomic>
#endif

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <chrono>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>

#if !defined(SYS_close_range) && defined(__linux__) && defined(__x86_64__)
#define SYS_close_range 436  // The package targets linux-64, including the glibc 2.17 sysroot.
#endif

namespace epcsaft::regression {
namespace {

const char* compute_cpu_directory = std::getenv("EPCSAFT_FIT_CPU_PROFILE");

void require(bool condition, const char* what) {
    if (!condition) throw std::invalid_argument(std::string("fit ") + what);
}

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); });
}

// The same point-local D maps raw observation columns and covariance.
void chain(const FitRequest& request, const std::vector<double>& map, double slope, const double* gradient, double* out) {
    const auto columns = request.coordinates.coordinates.size();
    std::fill(out, out + columns, 0.0);
    for (std::size_t r = 0; r < map.size() / columns; ++r)
        for (std::size_t c = 0; c < columns; ++c) out[c] += slope * gradient[r] * map[r * columns + c];
}

// Each input row names the first row with the same state; this is the sole membership map.
std::vector<std::size_t> state_groups(const std::vector<Observation>& observations, const std::vector<std::size_t>* membership = nullptr) {
    std::vector<std::size_t> groups;
    // ponytail: quadratic declaration comparison; index exact state records if fit datasets grow large.
    for (std::size_t i = 0; i < observations.size(); ++i) {
        const auto prior = std::find_if(observations.begin(), observations.begin() + static_cast<std::ptrdiff_t>(i),
            [&](const auto& row) { return (!membership || (*membership)[i] == (*membership)[&row - observations.data()])
                && same_state(observations[i], row); });
        const auto j = static_cast<std::size_t>(std::distance(observations.begin(), prior));
        groups.push_back(j == i ? i : groups[j]);
    }
    return groups;
}

#include "fit_measurement.hpp"

struct Trial final {
    std::unique_ptr<const HelmholtzModel> model;  // declared first: compiled centers are destroyed first
    std::vector<Observation> observations;
    std::vector<Prediction> centers;
    std::vector<double> optimizer;
    ObservationEvaluation evaluation;
    bool jacobian = false;
    EvaluatedCoordinates chart;
    ReferenceCache references;
};

// The serial sweep's first failure: `row` orders failures across worker processes; `observation` is
// the reported index (the observation count for a refused trial).
struct Failure final {
    std::size_t row = 0, observation = 0;
    ObservationStatus status = ObservationStatus::Nonfinite;
    bool jacobian = false;
    std::string message;
    std::vector<ResidualEnergyAction> parameter_actions{};
    std::optional<core::Convergence> first_numerical_failure;
};

FitHost fit_host;

bool interrupted() { return fit_host.interrupted && fit_host.interrupted(); }
#include "fit_worker.hpp"

}  // namespace

void set_fit_host(FitHost host, const void* extension) {
#ifdef _WIN32
    auto image = open_image(extension);
    if (!host.worker.empty() && (!image || !native_image)) throw std::runtime_error("fit worker images unavailable");
    extension_image = std::move(image);
#else
    const int image = open_image(extension);
    if (!host.worker.empty() && (image < 0 || native_image < 0)) {
        if (image >= 0) close(image);
        throw std::runtime_error("fit worker images unavailable");
    }
    if (extension_image >= 0) close(extension_image);
    extension_image = image;

#endif
    fit_host = std::move(host);
}

// One successful point owns its immutable model and compiled central states through later promotion.
// Worker processes each evaluate a fixed share of the rows and keep that share's starts and centers;
// every process publishes only when the whole evaluation succeeds.
struct ObservationEvaluator::Impl final {
    Impl(FitRequest request, FitResult& result, const std::vector<std::size_t>& groups,
              const std::vector<double>& sqrt_weights, std::optional<MeasurementChart> measurement = {})
        : request_([&] { request.controls.solver_options = request.controls.solver_options.resolved({});
            request.controls.solver_options.nested_solve = true;
            for (auto& o : request.observations) {
                o.problem.solver_options = o.problem.solver_options.resolved(request.controls.solver_options);
                o.final_problem.solver_options = o.final_problem.solver_options.resolved(request.controls.solver_options);
            } return request; }()), result_(result), measurement_(measurement ? std::move(*measurement) : MeasurementChart(request_)),
          groups_(state_groups(request_.observations, request_.uncertainty ? &measurement_.membership : nullptr)), sqrt_weights_(sqrt_weights),
          starts_(groups.size()), owner_(groups.size(), 0) {
        if (request.uncertainty) for (std::size_t i = 0; i < sqrt_weights_.size(); ++i)
            if (sqrt_weights_[i] > 0.0) sqrt_weights_[i] = request.observations[i].scale / request.uncertainty->output_sigma[i];
    }
    std::size_t width() const { return request_.coordinates.coordinates.size() + measurement_.groups.size(); }
    std::pair<std::vector<double>, std::vector<double>> bounds() const {
        auto box = request_.coordinates.optimizer_bounds();
        box.first.insert(box.first.end(), measurement_.lower.begin(), measurement_.lower.end());
        box.second.insert(box.second.end(), measurement_.upper.begin(), measurement_.upper.end());
        return box;
    }

    // Fixed partition: worker k evaluates every count-th state group in row order. Electrolyte rows share
    // one reference cache and stay here, as does every row under a trace capture or in a native host.
    void start_workers(int requested) {
        interruptible_ = true;
        const std::size_t count = core::current_trace() || fit_host.worker.empty() ? 1
            : requested > 0 ? static_cast<std::size_t>(requested) :
#ifdef _WIN32
                std::max(1ul, GetActiveProcessorCount(ALL_PROCESSOR_GROUPS));
#else
                std::max(1u, std::thread::hardware_concurrency());
#endif
        std::size_t rank = 0;
        for (std::size_t i = 0; i < groups_.size(); ++i)
            owner_[i] = groups_[i] != i ? owner_[groups_[i]]
                : electrolyte_observation(request_.observations[i].kind) ? 0 : rank++ % count;
#ifdef _WIN32
        for (std::size_t k = 1; k < std::min(count, rank); ++k) workers_.start();
#else
        std::vector<char*> command;
        for (auto& word : fit_host.worker) command.push_back(word.data());
        command.push_back(nullptr);
        char** arguments = command.data();
#ifdef __APPLE__
        for (std::size_t k = 1; k < std::min(count, rank); ++k) {
            int pair[2];
            if (socketpair(AF_UNIX, SOCK_STREAM, 0, pair) != 0) throw std::runtime_error("fit worker socket unavailable");
            workers_.sockets.push_back(pair[0]);
            const int on = 1;
            // As on Linux: socket end as fd 0, images as fd 3/4, the standard streams, its own process group.
            posix_spawn_file_actions_t actions;
            posix_spawnattr_t attributes;
            posix_spawn_file_actions_init(&actions);
            posix_spawnattr_init(&attributes);
            pid_t process = 0;
            const int error = fcntl(pair[0], F_SETFD, FD_CLOEXEC) || fcntl(pair[1], F_SETFD, FD_CLOEXEC)
                || setsockopt(pair[0], SOL_SOCKET, SO_NOSIGPIPE, &on, sizeof on) || setsockopt(pair[1], SOL_SOCKET, SO_NOSIGPIPE, &on, sizeof on)
                || posix_spawn_file_actions_adddup2(&actions, pair[1], 0) || posix_spawn_file_actions_adddup2(&actions, extension_image, 3)
                || posix_spawn_file_actions_adddup2(&actions, native_image, 4) || posix_spawn_file_actions_addinherit_np(&actions, 1)
                || posix_spawn_file_actions_addinherit_np(&actions, 2) || posix_spawnattr_setflags(&attributes, POSIX_SPAWN_CLOEXEC_DEFAULT | POSIX_SPAWN_SETPGROUP)
                || posix_spawnattr_setpgroup(&attributes, 0) || posix_spawn(&process, arguments[0], &actions, &attributes, arguments, *_NSGetEnviron());
            posix_spawn_file_actions_destroy(&actions);
            posix_spawnattr_destroy(&attributes);
            close(pair[1]);
            if (error != 0) throw std::runtime_error("fit worker process unavailable");
            workers_.processes.push_back(process);
        }
#else
        rlimit descriptors;
        if (getrlimit(RLIMIT_NOFILE, &descriptors) != 0) throw std::runtime_error("fit worker process unavailable");
        // linux-64: nr_open bounds descriptor allocation, even when the soft limit was lowered.
        rlim_t nr_open = std::numeric_limits<int>::max();
        std::ifstream("/proc/sys/fs/nr_open") >> nr_open;
        const rlim_t descriptor_limit = std::min(descriptors.rlim_max, nr_open);
        for (std::size_t k = 1; k < std::min(count, rank); ++k) {
            int pair[2];
            if (socketpair(AF_UNIX, SOCK_STREAM | SOCK_CLOEXEC, 0, pair) != 0) throw std::runtime_error("fit worker socket unavailable");
            workers_.sockets.push_back(pair[0]);
            // The worker keeps its socket end as fd 0, the extension/native images as fd 3/4 and the standard streams. Its own
            // process group keeps terminal signals away: the fit, not its workers, answers an interrupt.
            pid_t process = 0;
            const int error = (process = fork()) < 0;
            if (process == 0) {
                if (setpgid(0, 0) != 0) _exit(127);
                const int sources[] = {fcntl(pair[1], F_DUPFD_CLOEXEC, 5), fcntl(extension_image, F_DUPFD_CLOEXEC, 5), fcntl(native_image, F_DUPFD_CLOEXEC, 5)};
                const int targets[] = {0, 3, 4};
                for (int i = 0; i < 3; ++i)
                    if (sources[i] < 0 || dup2(sources[i], targets[i]) < 0 || fcntl(targets[i], F_SETFD, 0) < 0) _exit(127);
                if (syscall(SYS_close_range, 5, ~0U, 0) != 0) {
                    const int directory = static_cast<int>(syscall(SYS_open, "/proc/self/fd", O_RDONLY | O_DIRECTORY | O_CLOEXEC));
                    if (directory >= 0) {
                        char buffer[4096];
                        long bytes;
                        while ((bytes = syscall(SYS_getdents64, directory, buffer, sizeof buffer)) > 0) {
                            for (long at = 0; at < bytes;) {
                                // linux_dirent64 has d_reclen at byte 16 and d_name at byte 19 on linux-64.
                                unsigned short length;
                                std::memcpy(&length, buffer + at + 16, sizeof length);
                                int fd = 0;
                                for (const char* digit = buffer + at + 19; *digit >= '0' && *digit <= '9'; ++digit) fd = 10 * fd + *digit - '0';
                                if (fd >= 5 && fd != directory && close(fd) != 0) _exit(127);
                                at += length;
                            }
                        }
                        if (bytes < 0 || close(directory) != 0) _exit(127);
                    } else for (rlim_t fd = 5; fd < descriptor_limit; ++fd)
                        if (close(static_cast<int>(fd)) != 0 && errno != EBADF) _exit(127);
                }
                execve(arguments[0], arguments, environ);
                _exit(127);
            }
            close(pair[1]);
            if (error != 0) throw std::runtime_error("fit worker process unavailable");
            workers_.processes.push_back(process);
        }
#endif
#endif
        for (std::size_t k = 1; k <= workers_.sockets.size(); ++k) {
            put(workers_.sockets[k - 1], request_.parameters, request_.thermochemistry, request_.coordinates.coordinates,
                request_.observations, request_.holds, request_.uncertainty, request_.controls.solver_options.values, measurement_, sqrt_weights_);
            put(workers_.sockets[k - 1], owner_, k, deadline_);
        }
    }

    bool evaluate(const std::vector<double>& optimizer, bool jacobian, ObservationEvaluation& out) {
        core::Scope scope("fit_evaluation");
        scope.probe("optimizer", optimizer);
        scope.probe("jacobian", jacobian ? 1.0 : 0.0);
        const bool same_point = cached_ && cached_->optimizer == optimizer;
        const bool hit = same_point && (!jacobian || cached_->jacobian);
        scope.probe("reuse", hit ? (cached_->jacobian ? 2.0 : 1.0) : same_point ? 3.0 : 0.0);
        if (hit) {
            out = cached_->evaluation;
            return true;
        }
        out = ObservationEvaluation{};
        out.solver_options = request_.controls.solver_options;
        ++trials_;
        (jacobian ? result_.jacobian_evaluations : result_.residual_evaluations) += 1;
        for (const Channel fd : workers_.sockets) put(fd, optimizer, jacobian);
        const std::size_t n = request_.observations.size(), width = request_.coordinates.tie().size() + request_.holds.size();
        Trial fresh;
        Trial& trial = same_point ? *cached_ : fresh;
        auto failure = compute(optimizer, jacobian, out, trial, same_point, scope);
        for (std::size_t k = 1; k <= workers_.sockets.size(); ++k) {
            std::optional<Failure> theirs;
            ObservationEvaluation part;
            get(workers_.sockets[k - 1], theirs, part.predictions, part.residuals, part.statuses, part.gradient, part.input_gradient,
                part.state_densities, part.reference_values, part.saturation_predictions, part.first_numerical_failure, part.first_numerical_failure_observation);
            // A trial refused before its rows sends none; any other reply carries every row.
            const bool rows = part.statuses.size() == n && part.predictions.size() == n && part.residuals.size() == n
                && part.gradient.size() == (jacobian ? n * width : 0)
                && part.input_gradient.size() == (jacobian ? n * measurement_.groups.size() : 0);
            if (!rows && !(theirs && part.statuses.empty())) throw std::runtime_error("fit worker reply size");
            for (std::size_t i = 0; rows && out.statuses.size() == n && i < n; ++i) if (owner_[i] == k) {
                out.predictions[i] = part.predictions[i];
                out.residuals[i] = part.residuals[i];
                out.statuses[i] = part.statuses[i];
                if (!part.saturation_predictions.empty()) out.saturation_predictions[i] = std::move(part.saturation_predictions[i]);
                if (request_.uncertainty) { out.state_densities[i] = part.state_densities.at(i); out.reference_values[i] = part.reference_values.at(i); }
                if (jacobian) std::copy_n(part.gradient.begin() + i * width, width, out.gradient.begin() + i * width);
                if (jacobian && !measurement_.groups.empty()) std::copy_n(part.input_gradient.begin() + i * measurement_.groups.size(),
                    measurement_.groups.size(), out.input_gradient.begin() + i * measurement_.groups.size());
            }
            if (theirs && (!failure || theirs->row < failure->row)) failure = std::move(theirs);
            if (part.first_numerical_failure && (!out.first_numerical_failure
                || part.first_numerical_failure_observation < out.first_numerical_failure_observation)) {
                out.first_numerical_failure = part.first_numerical_failure;
                out.first_numerical_failure_observation = part.first_numerical_failure_observation;
            }
        }
        if (out.first_numerical_failure && failure && out.first_numerical_failure_observation > failure->row)
            out.first_numerical_failure.reset();
        if (failure) core::retain_first(out.first_numerical_failure, failure->first_numerical_failure);
        core::retain_first(result_.first_numerical_failure, out.first_numerical_failure);
        // The serial sweep stops at its first failure: later rows stay unevaluated.
        for (std::size_t i = failure ? failure->row + 1 : n; i < out.statuses.size(); ++i) {
            out.predictions[i] = out.residuals[i] = 0.0;
            out.statuses[i] = ObservationStatus::NotEvaluated;
            if (!out.saturation_predictions.empty()) out.saturation_predictions[i].reset();
            if (jacobian) std::fill_n(out.gradient.begin() + i * width, width, 0.0);
        }
        const auto reject = [&](std::size_t index, ObservationStatus status, std::string message, bool flagged, std::vector<ResidualEnergyAction> actions = {}) {
            for (const Channel fd : workers_.sockets) put(fd, false);
            if (!out.saturation_predictions.empty()) result_.saturation_predictions = out.saturation_predictions;
            return fail(index, optimizer, status, std::move(message), flagged, std::move(actions), out.first_numerical_failure);
        };
        if (failure) return reject(failure->observation, failure->status, std::move(failure->message), failure->jacobian, std::move(failure->parameter_actions));
        out.weighted.resize(n);
        for (std::size_t i = 0; i < n; ++i) out.weighted[i] = sqrt_weights_[i] * out.residuals[i];
        // One half of the row-ordered squared sum; reject its first nonfinite operand or sum.
        double sum = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double weighted = out.weighted[i];
            sum += weighted * weighted;
            if (!std::isfinite(weighted) || !std::isfinite(sum))
                return reject(i, ObservationStatus::Nonfinite,
                    std::isfinite(weighted) ? "accumulated weighted cost not finite" : "weighted residual not finite", jacobian);
        }
        out.cost = 0.5 * sum;
        const std::size_t parameters = request_.coordinates.coordinates.size(), latent = measurement_.groups.size(), columns = this->width();
        if (jacobian) out.optimizer_gradient.resize((n + latent) * columns);
        for (std::size_t i = 0; jacobian && i < n; ++i) {
            if (sqrt_weights_[i] == 0.0) continue;
            chain(request_, trial.chart.map, sqrt_weights_[i] * residual_slope(request_.observations[i], out.predictions[i]),
                  out.gradient.data() + i * width, out.optimizer_gradient.data() + i * columns);
            for (std::size_t c = 0; c < latent; ++c) out.optimizer_gradient[i * columns + parameters + c]
                = sqrt_weights_[i] * residual_slope(request_.observations[i], out.predictions[i])
                    * out.input_gradient[i * latent + c] * measurement_.sigma[c];
            if (!Eigen::Map<const Eigen::VectorXd>(out.optimizer_gradient.data() + i * columns,
                static_cast<Eigen::Index>(columns)).allFinite())
                return reject(i, ObservationStatus::Nonfinite, "optimizer Jacobian not finite", true);
        }
        if (request_.uncertainty) {
            auto& reconciled = out.reconciliation.emplace();
            reconciled.output_cost = out.cost;
            reconciled.row_groups = measurement_.membership;
            reconciled.column_groups = measurement_.groups; reconciled.column_inputs = measurement_.inputs;
            for (std::size_t c = 0; c < latent; ++c) {
                const double u = optimizer[parameters + c];
                reconciled.values.push_back(measurement_.measured[c] + measurement_.sigma[c] * u);
                reconciled.corrections.push_back(u); reconciled.penalties.push_back(u*u);
                reconciled.row_groups.push_back(measurement_.groups[c]);
                reconciled.input_cost += 0.5 * u*u; out.weighted.push_back(u);
                if (jacobian) out.optimizer_gradient[(n + c) * columns + parameters + c] = 1.0;
            }
            out.cost += reconciled.input_cost;
            if (!std::isfinite(out.cost)) return reject(n, ObservationStatus::Nonfinite, "input penalty not finite", jacobian);
            if (!measurement_.factor.empty() && !whiten(out, columns))
                return reject(n, ObservationStatus::Nonfinite, "whitened residual or Jacobian not finite", jacobian);
            if (jacobian) reconciled.jacobian = out.optimizer_gradient;
        }
        for (const Channel fd : workers_.sockets) put(fd, true);
        publish(optimizer, jacobian, out, trial, fresh, same_point);
        return true;
    }

    // z = L^-1 e on the compact rows of the residual and of every complete column; excluded rows stay zero.
    // Output-first order makes the output cost marginal and the input cost conditional on the outputs.
    bool whiten(ObservationEvaluation& out, std::size_t columns) const {
        const auto& rows = measurement_.rows;
        const auto m = static_cast<Eigen::Index>(rows.size()), inputs = static_cast<Eigen::Index>(measurement_.groups.size());
        const std::size_t width = out.optimizer_gradient.empty() ? 0 : columns;
        Eigen::MatrixXd e(m, static_cast<Eigen::Index>(1 + width));
        for (std::size_t k = 0; k < rows.size(); ++k) {
            e(static_cast<Eigen::Index>(k), 0) = out.weighted[rows[k]];
            for (std::size_t c = 0; c < width; ++c) e(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(1 + c)) = out.optimizer_gradient[rows[k] * columns + c];
        }
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(measurement_.factor.data(), m, m)
            .triangularView<Eigen::Lower>().solveInPlace(e);
        for (std::size_t k = 0; k < rows.size(); ++k) {
            out.weighted[rows[k]] = e(static_cast<Eigen::Index>(k), 0);
            for (std::size_t c = 0; c < width; ++c) out.optimizer_gradient[rows[k] * columns + c] = e(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(1 + c));
        }
        auto& reconciled = *out.reconciliation;
        reconciled.output_cost = 0.5 * e.col(0).head(m - inputs).squaredNorm();
        reconciled.input_cost = 0.5 * e.col(0).tail(inputs).squaredNorm();
        for (Eigen::Index c = 0; c < inputs; ++c) reconciled.penalties[static_cast<std::size_t>(c)] = e(m - inputs + c, 0) * e(m - inputs + c, 0);
        out.cost = reconciled.output_cost + reconciled.input_cost;
        const bool finite = e.allFinite() && std::isfinite(out.cost);
        if (!finite) {
            core::Convergence record; record.residual = out.cost; record.stop = core::Stop::nonfinite;
            record.cause = e.allFinite() ? "whitened cost not finite" : "whitened residual or Jacobian not finite";
            core::retain_first(out.first_numerical_failure, record);
            core::Scope scope("fit_whitening");
            scope.probe("solved_operands", e.data(), e.size());
            scope.probe("cost", out.cost); scope.finish(record);
        }
        return finite;
    }

    // A worker answers each evaluation with its rows and first failure, then publishes on the fit's verdict.
    [[noreturn]] void serve(Channel fd) {
        for (std::vector<double> optimizer;;) {
            bool jacobian = false, commit = false;
            get(fd, optimizer, jacobian);
            core::Scope scope("fit_evaluation");
            const bool same_point = cached_ && cached_->optimizer == optimizer;
            Trial fresh;
            Trial& trial = same_point ? *cached_ : fresh;
            ObservationEvaluation out;
            const auto failure = compute(optimizer, jacobian, out, trial, same_point, scope);
            put(fd, failure, out.predictions, out.residuals, out.statuses, out.gradient, out.input_gradient,
                out.state_densities, out.reference_values, out.saturation_predictions, out.first_numerical_failure, out.first_numerical_failure_observation);
            get(fd, commit);
            if (commit) publish(optimizer, jacobian, out, trial, fresh, same_point);
        }
    }

    bool gradient_at(const std::vector<double>& point, std::vector<double>& out) {
        ObservationEvaluation trial;
        if (!evaluate(point, true, trial)) return false;
        const std::size_t width = this->width();
        const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> j(
            trial.optimizer_gradient.data(), static_cast<Eigen::Index>(trial.weighted.size()), static_cast<Eigen::Index>(width));
        const Eigen::Map<const Eigen::VectorXd> r(trial.weighted.data(), static_cast<Eigen::Index>(trial.weighted.size()));
        const Eigen::VectorXd gradient = j.transpose() * r;
        out.assign(gradient.data(), gradient.data() + gradient.size());
        return finite(out) || fail(trial.weighted.size(), point, ObservationStatus::Nonfinite, "weighted objective gradient not finite", true);
    }

    bool fail(std::size_t observation, const std::vector<double>& optimizer, ObservationStatus status,
              std::string message, bool jacobian, std::vector<ResidualEnergyAction> actions = {}, std::optional<core::Convergence> first = {}) {
        core::retain_failure(first, message, "fit evaluation");
        core::retain_first(result_.first_numerical_failure, first);
        result_.failures.push_back({trials_, observation, optimizer, status, std::move(message), jacobian, std::move(actions), first});
        ++result_.trial_failures;
        return false;
    }

    // This process's rows of one trial; `i`, the sweep position, orders a refused trial among row failures.
    std::optional<Failure> compute(const std::vector<double>& optimizer, bool jacobian, ObservationEvaluation& out,
                                   Trial& trial, bool same_point, core::Scope& scope) {
        if (!compute_cpu_directory) return compute_rows(optimizer, jacobian, out, trial, same_point, scope);
        // #405/#498: process CPU inside compute, excluding output setup and writing.
        std::ofstream file(std::string(compute_cpu_directory) + "/" + std::to_string(
#ifdef _WIN32
            GetCurrentProcessId()
#else
            getpid()
#endif
        ) + ".csv", std::ios::app);
        if (!file) throw std::runtime_error("fit compute CPU output unavailable");
        const auto cpu = [] {
#ifdef _WIN32
            FILETIME created, exited, kernel, user;
            if (!GetProcessTimes(GetCurrentProcess(), &created, &exited, &kernel, &user)) throw std::runtime_error("fit compute CPU clock unavailable");
            const auto ticks = [](FILETIME t) { return (static_cast<unsigned long long>(t.dwHighDateTime) << 32) | t.dwLowDateTime; };
            return double(ticks(kernel) + ticks(user)) * 1e-7;
#else
            timespec value{};
            if (clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &value)) throw std::runtime_error("fit compute CPU clock unavailable");
            return double(value.tv_sec) + double(value.tv_nsec) * 1e-9;
#endif
        };
        const double start = cpu();
        auto result = compute_rows(optimizer, jacobian, out, trial, same_point, scope);
        const double elapsed = cpu() - start;
        file << ++profile_evaluations_ << ',' << process_ << ',' << jacobian << ',' << std::setprecision(17) << elapsed << '\n';
        file.flush();
        if (!file) throw std::runtime_error("fit compute CPU output failed");
        return result;
    }

    std::optional<Failure> compute_rows(const std::vector<double>& optimizer, bool jacobian, ObservationEvaluation& out,
                                       Trial& trial, bool same_point, core::Scope& scope) {
        const std::size_t n = request_.observations.size(), width = request_.coordinates.tie().size() + request_.holds.size();
        std::size_t i = 0;
        const auto retain = [&](const auto& first) {
            if (!out.first_numerical_failure && first) {
                out.first_numerical_failure = first; out.first_numerical_failure_observation = i;
            }
        };
        try {
            require(optimizer.size() == this->width(), "complete optimizer width");
            const std::vector<double> parameter_point(optimizer.begin(), optimizer.begin() + request_.coordinates.coordinates.size());
            auto observations = request_.observations;
            auto reaction_laws = request_.coordinates.apply(observations, parameter_point);
            measurement_.apply(request_, optimizer, observations);
            if (!same_point) {
                trial.chart = evaluate_coordinates(request_.parameters, request_.coordinates, parameter_point, request_.holds, request_.thermochemistry, request_.controls.solver_options);
                trial.chart.reaction_laws = std::move(reaction_laws);
                trial.model = std::make_unique<HelmholtzModel>(trial.chart.parameters, trial.chart.active, request_.thermochemistry, request_.controls.solver_options);
                trial.centers.assign(groups_.size(), Prediction{});
            }
            trial.observations = std::move(observations);
            const Trial* seed = nullptr;
            if (!same_point && cached_ && trial.chart.physical.size() == cached_->chart.physical.size()) {
                bool near = true;
                for (std::size_t k = 0; k < trial.chart.physical.size(); ++k) {
                    const double a = trial.chart.physical[k], b = cached_->chart.physical[k];
                    near &= std::isfinite(a) && std::isfinite(b)
                        && std::abs(a - b) / std::max({1.0, std::abs(a), std::abs(b)}) <= policy::FIT_DENSITY_PARAMETER_FRACTION.value;
                }
                if (near) seed = &*cached_;
            }
            scope.probe("density_seed_chart", seed ? 1.0 : 0.0);
            const auto& model = *trial.model;
            out.predictions.assign(n, 0.0);
            const bool capped = std::any_of(trial.observations.begin(), trial.observations.end(),
                [](const auto& o) { return o.kind == ObservationKind::SaturationCappedAdditionEnthalpy; });
            if (capped) out.saturation_predictions.resize(n);
            out.residuals.assign(n, 0.0);
            out.statuses.assign(n, ObservationStatus::NotEvaluated);
            if (request_.uncertainty) { out.state_densities.resize(n); out.reference_values.resize(n); }
            if (jacobian) out.gradient.assign(n * width, 0.0);
            if (jacobian) out.input_gradient.assign(n * measurement_.groups.size(), 0.0);
            std::vector<std::vector<PredictionGradient>> gradients(n);
            std::vector<std::size_t> position(n, 0);
            for (i = 0; i < n; ++i) {
                if (sqrt_weights_[i] == 0.0 || owner_[i] != process_) continue;
                if (interruptible_ && interrupted()) throw FitInterrupted{};
                const auto& observation = trial.observations[i];
                const std::size_t g = groups_[i];
                auto& prediction = trial.centers[g];
                if (!prediction.value) {
                    core::Scope scope("fit_group_prediction");
                    scope.probe("observation", static_cast<double>(i));
                    scope.probe("anchor", static_cast<double>(g));
                    scope.probe("targets", static_cast<double>(std::count(groups_.begin(), groups_.end(), g)));
                    const auto& prior = seed ? seed->observations[i] : observation;
                    const bool unchanged = !request_.uncertainty || (pressure_closed(observation.kind)
                        ? std::tie(observation.temperature, observation.pressure, observation.composition, observation.branch)
                            == std::tie(prior.temperature, prior.pressure, prior.composition, prior.branch)
                        : same_state(observation, prior));
                    prediction = predict(model, observation, starts_[g], std::move(prediction.compiled), &trial.references,
                        seed && unchanged ? &seed->centers[g] : nullptr, seed && unchanged ? seed->model.get() : nullptr,
                        seed && unchanged && trial.chart.physical == seed->chart.physical, deadline_);
                    if (prediction.status != ObservationStatus::Available && starts_[g] && observation.kind != ObservationKind::SaturationCappedAdditionEnthalpy) {
                        auto retry = predict(model, observation, std::nullopt, std::move(prediction.compiled), &trial.references);
                        core::retain_first(retry.first_numerical_failure, prediction.first_numerical_failure, true);
                        prediction = std::move(retry);
                    }
                } else if (!pressure_closed(observation.kind)) project_prediction(model, observation, prediction);
                retain(prediction.first_numerical_failure);
                if (prediction.saturation) out.saturation_predictions[i] = prediction;
                out.statuses[i] = prediction.status;
                if (prediction.status != ObservationStatus::Available)
                    return Failure{i, i, prediction.status, jacobian, prediction.message, {}, prediction.first_numerical_failure};
                out.predictions[i] = *prediction.value;
                out.residuals[i] = *prediction.residual;
                if (request_.uncertainty) {
                    if (prediction.activity) out.state_densities[i].push_back(prediction.activity->state.molar_density);
                    if (prediction.state) out.state_densities[i].push_back(prediction.state->molar_density);
                    if (prediction.solved) out.state_densities[i] = prediction.solved->molar_densities;
                    if (prediction.reference) {
                        out.reference_values[i].assign(prediction.reference->log_fugacity.begin(), prediction.reference->log_fugacity.end());
                        if (prediction.reference->pure_solvent_log_fugacity) out.reference_values[i].push_back(*prediction.reference->pure_solvent_log_fugacity);
                        if (prediction.reference->solvent_state) out.state_densities[i].push_back(prediction.reference->solvent_state->molar_density);
                    }
                }
                if (!jacobian) continue;
                if (gradients[g].empty()) {
                    core::Scope scope("fit_group_differentiation");
                    scope.probe("anchor", static_cast<double>(g));
                    std::vector<const Observation*> rows;
                    for (std::size_t row = 0; row < n; ++row)
                        if (groups_[row] == g && sqrt_weights_[row] > 0.0) rows.push_back(&trial.observations[row]);
                    const auto directions = measurement_.directions(request_, i).first;
                    gradients[g] = pressure_closed(observation.kind) && directions.empty()
                        ? std::vector<PredictionGradient>{differentiate(model, observation, prediction)}
                        : differentiate_group(model, rows, prediction, request_.coordinates, directions);
                }
                const auto& gradient = gradients[g][position[g]++];
                if (prediction.saturation) {
                    auto& detail = *out.saturation_predictions[i];
                    const auto map_loading = [&](const auto& physical, auto& mapped) {
                        if (physical.empty()) return;
                        mapped.resize(this->width()); chain(request_, trial.chart.map, 1, physical.data(), mapped.data());
                    };
                    map_loading(gradient.saturation_loading_gradient, detail.saturation_loading_gradient);
                    map_loading(gradient.endpoint_loading_gradient, detail.endpoint_loading_gradient);
                }
                retain(gradient.first_numerical_failure);
                if (gradient.status != ObservationStatus::Available)
                    return Failure{i, i, gradient.status, true, gradient.message, gradient.parameter_actions, gradient.first_numerical_failure};
                std::copy_n(gradient.value.begin(), std::min(width, gradient.value.size()),
                          out.gradient.begin() + static_cast<std::ptrdiff_t>(i * width));
                std::size_t local = width;
                for (std::size_t c = 0; c < measurement_.groups.size(); ++c) if (measurement_.membership[i] == measurement_.groups[c])
                    out.input_gradient[i * measurement_.groups.size() + c] = gradient.value.at(local++);
            }
            return std::nullopt;
        } catch (const std::bad_alloc&) {
            throw;
        } catch (const std::exception& error) {
            if (!out.first_numerical_failure) out.first_numerical_failure_observation = i;
            retain_density_failure(out.first_numerical_failure, error, "fit trial");
            return Failure{i, n, ObservationStatus::InadmissibleTrial, jacobian, error.what(), {}, out.first_numerical_failure};
        }
    }

    // Neither a failed row nor an overflowing weighted/coordinate chain publishes starts or a cache.
    void publish(const std::vector<double>& optimizer, bool jacobian, ObservationEvaluation& out, Trial& trial,
                 Trial& fresh, bool same_point) {
        for (std::size_t g = 0; g < groups_.size(); ++g)
            if (trial.centers[g].solved) starts_[g] = trial.centers[g].solved;
        out.chart = trial.chart;
        trial.optimizer = optimizer;
        trial.evaluation = out;
        trial.jacobian = jacobian;
        if (!same_point) {
            cached_.emplace(std::move(fresh));
        }
    }

    const FitRequest request_;
    std::chrono::steady_clock::time_point deadline_ = request_.controls.maximum_elapsed_time_seconds > 0
        ? std::chrono::steady_clock::now() + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
            std::chrono::duration<double>(request_.controls.maximum_elapsed_time_seconds))
        : std::chrono::steady_clock::time_point::max();
    FitResult& result_;
    const MeasurementChart measurement_;
    const std::vector<std::size_t> groups_;
    std::vector<double> sqrt_weights_;
    std::vector<std::optional<equilibrium::EquilibriumResult>> starts_;
    std::optional<Trial> cached_;
    std::size_t trials_ = 0;
    std::size_t profile_evaluations_ = 0;
    std::vector<std::size_t> owner_;  // the process that evaluates each row; 0 is the fit's own
    std::size_t process_ = 0;
    bool interruptible_ = false;  // a least-squares fit's own process polls cancellation between rows
    Workers workers_;
};

// A worker never outlives the thread that started it; it reads its fit's request once, then serves it.
// The socket names the fit's process (its creator), so a fit that died before the death signal was armed
// shows as a different parent, including a subreaper's adoption, and the worker leaves without reading.
#ifdef _WIN32
void serve_fit(std::uintptr_t channel) {
    const Channel fd{reinterpret_cast<HANDLE>(channel)};
#else
void serve_fit(int fd) {
#ifdef __APPLE__
    // No death signal on macOS: launchd (pid 1) adopts an orphan, so the worker leaves at entry or within 0.1 s.
    if (getppid() == 1) _exit(1);
    std::thread([] { while (getppid() != 1) std::this_thread::sleep_for(std::chrono::milliseconds(100)); _exit(1); }).detach();
#else
    ucred fit{};
    socklen_t size = sizeof fit;
    if (getsockopt(fd, SOL_SOCKET, SO_PEERCRED, &fit, &size) != 0 || prctl(PR_SET_PDEATHSIG, SIGKILL) != 0
        || getppid() != fit.pid) _exit(1);
#endif
#endif
    FitRequest request;
    MeasurementChart measurement;
    std::vector<double> sqrt_weights;
    get(fd, request.parameters, request.thermochemistry, request.coordinates.coordinates, request.observations,
        request.holds, request.uncertainty, request.controls.solver_options.values, measurement, sqrt_weights);
    const auto groups = state_groups(request.observations);
    FitResult result; result.solver_options = request.controls.solver_options.resolved({});
    ObservationEvaluator::Impl evaluator(request, result, groups, sqrt_weights, std::move(measurement));
    get(fd, evaluator.owner_, evaluator.process_, evaluator.deadline_);
    evaluator.serve(fd);
}

ObservationEvaluator::ObservationEvaluator(const FitRequest& request, FitResult& failures, const std::vector<std::size_t>& rows) {
    for (const auto& coordinate : request.coordinates.coordinates)
        if (coordinate.reaction && coordinate.reaction->coefficient == ReactionCoefficient::EnthalpyShift) {
            validate(request.coordinates, request.parameters, request.observations, request.thermochemistry);
            break;
        }
    std::vector<double> selection(request.observations.size(), 0.0);
    for (const auto row : rows) { require(row < selection.size(), "evaluator row"); selection[row] = 1.0; }
    auto canonical = request;
    if (canonical.uncertainty && canonical.uncertainty->groups.empty() && canonical.uncertainty->output_sigma.empty()
        && !canonical.uncertainty->covariance && !canonical.uncertainty->correlation) canonical.uncertainty.reset();
    impl_ = std::make_shared<Impl>(canonical, failures, state_groups(canonical.observations), selection);
}
bool ObservationEvaluator::evaluate(const std::vector<double>& optimizer, bool columns, ObservationEvaluation& out) {
    out.residual_optimizer_gradient.clear();
    if (!impl_->evaluate(optimizer, columns, out)) return false;
    if (columns) out.residual_optimizer_gradient = out.optimizer_gradient;
    const auto& measurement = impl_->measurement_;
    const std::size_t parameters = impl_->request_.coordinates.coordinates.size(), latent = measurement.groups.size();
    if (columns) for (std::size_t i = 0; i < out.predictions.size(); ++i) if (out.statuses[i] == ObservationStatus::Available) {
        double* row = out.optimizer_gradient.data() + i * optimizer.size();
        if (measurement.factor.empty()) {
            const double slope = residual_slope(impl_->request_.observations[i], out.predictions[i]) * impl_->sqrt_weights_[i];
            for (std::size_t c = 0; c < optimizer.size(); ++c) row[c] /= slope;
            continue;
        }
        // Whitening mixes rows, so raw dy/d(v,u) comes from the physical actions, never from objective rows.
        const std::size_t width = impl_->request_.coordinates.tie().size() + impl_->request_.holds.size();
        chain(impl_->request_, out.chart.map, 1.0, out.gradient.data() + i * width, row);
        for (std::size_t c = 0; c < latent; ++c) row[parameters + c] = out.input_gradient[i * latent + c] * measurement.sigma[c];
    }
    if (columns) out.optimizer_gradient.resize(out.predictions.size() * optimizer.size());
    return true;
}

namespace {
std::vector<int> select_held(const std::vector<double>& point, const std::vector<double>& gradient,
                            const std::vector<double>& lower, const std::vector<double>& upper) {
    std::vector<int> held;
    for (std::size_t c = 0; c < point.size(); ++c)
        if ((point[c] - lower[c] < policy::FIT_ACTIVE_BOUND_DISTANCE.value && gradient[c] > 0.0)
            || (upper[c] - point[c] < policy::FIT_ACTIVE_BOUND_DISTANCE.value && gradient[c] < 0.0))
            held.push_back(static_cast<int>(c));
    return held;
}

// The optimizer's chart of the box, in place of Ceres' bounds: held coordinates are fixed; a step leaving
// the box is shortened along its own direction to the first bound it meets (clamping one component keeps
// the rest of a step computed as if that bound could be crossed), clamping coordinates already on it. As
// the callback it ends a round on accepted new bound contact or a coordinate joining the held set, not on release.
class BoxChart final : public ceres::Manifold, public ceres::IterationCallback {
public:
    BoxChart(const std::vector<double>& u, const std::vector<double>& lower, const std::vector<double>& upper,
             const std::vector<int>& held,
             ObservationEvaluator::Impl& evaluator)
        : u_(u), lower_(lower), upper_(upper), held_(held), evaluator_(evaluator) {
        for (int c = 0; c < int(lower.size()); ++c)
            if (std::find(held.begin(), held.end(), c) == held.end()) free_.push_back(c);
        for (const int c : free_) at_bound_.push_back(side(u[c], c));
    }
    int side(double value, int c) const { return value == lower_[c] ? -1 : value == upper_[c] ? 1 : 0; }  // Plus writes bounds exactly
    ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override {
        core::Scope scope("fit_iteration");
        scope.probe("accepted", summary.step_is_successful ? 1.0 : 0.0, summary.iteration);
        scope.probe("cost", summary.cost, summary.iteration);
        scope.probe("ratio", summary.relative_decrease, summary.iteration);
        scope.probe("cost_change", summary.cost_change, summary.iteration);
        if (summary.relative_decrease != 0.0) scope.probe("predicted_cost_change", summary.cost_change / summary.relative_decrease, summary.iteration);
        scope.probe("trust_radius", summary.trust_region_radius, summary.iteration);
        if (summary.iteration > 0 && !summary.step_is_successful) return ceres::SOLVER_CONTINUE;
        std::vector<double> gradient;
        // Ceres evaluates the accepted Jacobian and updates u_ before user callbacks: this is a cache hit.
        if (!evaluator_.gradient_at(u_, gradient)) return ceres::SOLVER_ABORT;
        const auto next = select_held(u_, gradient, lower_, upper_);
        scope.observe([&] {
            std::vector<double> mask(u_.size(), 0.0), projected(u_.size());
            for (const int c : next) mask[c] = 1.0;
            for (std::size_t c = 0; c < u_.size(); ++c)
                projected[c] = u_[c] - std::clamp(u_[c] - gradient[c], lower_[c], upper_[c]);
            scope.probe("optimizer", u_, summary.iteration);
            scope.probe("held", mask, summary.iteration);
            scope.probe("weighted_gradient", gradient, summary.iteration);
            scope.probe("projected_gradient", projected, summary.iteration);
        });
        bool contact = false;
        for (std::size_t k = 0; summary.iteration > 0 && summary.step_is_successful && k < free_.size(); ++k) {
            contact = contact || (side(u_[free_[k]], free_[k]) != 0 && side(u_[free_[k]], free_[k]) != at_bound_[k]);
            at_bound_[k] = side(u_[free_[k]], free_[k]);
        }
        const bool join = !std::includes(held_.begin(), held_.end(), next.begin(), next.end());
        return summary.iteration > 0 && (contact || join)
            ? ceres::SOLVER_TERMINATE_SUCCESSFULLY : ceres::SOLVER_CONTINUE;
    }
    int AmbientSize() const override { return static_cast<int>(lower_.size()); }
    int TangentSize() const override { return static_cast<int>(free_.size()); }
    bool Plus(const double* x, const double* delta, double* out) const override {
        double alpha = 1.0;
        std::size_t limit = free_.size();
        for (std::size_t k = 0; k < free_.size(); ++k) {
            const int c = free_[k];
            const double bound = delta[k] < 0.0 ? lower_[c] : upper_[c];
            if (std::abs(bound - x[c]) >= policy::FIT_ACTIVE_BOUND_DISTANCE.value
                && (x[c] + delta[k] - bound) * delta[k] > 0.0 && (bound - x[c]) / delta[k] < alpha) {
                alpha = (bound - x[c]) / delta[k];
                limit = k;
            }
        }
        std::copy(x, x + lower_.size(), out);
        for (std::size_t k = 0; k < free_.size(); ++k)
            out[free_[k]] = std::clamp(x[free_[k]] + alpha * delta[k], lower_[free_[k]], upper_[free_[k]]);
        if (limit < free_.size()) out[free_[limit]] = delta[limit] < 0.0 ? lower_[free_[limit]] : upper_[free_[limit]];
        return true;
    }
    bool PlusJacobian(const double*, double* jacobian) const override {
        std::fill(jacobian, jacobian + lower_.size() * free_.size(), 0.0);
        for (std::size_t k = 0; k < free_.size(); ++k) jacobian[free_[k] * free_.size() + k] = 1.0;
        return true;
    }
    bool Minus(const double* y, const double* x, double* out) const override {
        for (std::size_t k = 0; k < free_.size(); ++k) out[k] = y[free_[k]] - x[free_[k]];
        return true;
    }
    bool MinusJacobian(const double*, double* jacobian) const override {
        std::fill(jacobian, jacobian + lower_.size() * free_.size(), 0.0);
        for (std::size_t k = 0; k < free_.size(); ++k) jacobian[k * lower_.size() + free_[k]] = 1.0;
        return true;
    }

private:
    const std::vector<double>& u_;
    std::vector<double> lower_, upper_;
    std::vector<int> free_;
    std::vector<int> at_bound_;
    std::vector<int> held_;
    ObservationEvaluator::Impl& evaluator_;
};

class Cost final : public ceres::CostFunction {
public:
    Cost(const FitRequest& request, ObservationEvaluator::Impl& evaluator)
        : request_(request), evaluator_(evaluator) {
        set_num_residuals(static_cast<int>(request.observations.size() + evaluator.measurement_.groups.size()));
        mutable_parameter_block_sizes()->push_back(static_cast<int>(evaluator.width()));
    }

    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const override {
        const std::size_t width = evaluator_.width();
        const std::vector<double> optimizer(parameters[0], parameters[0] + width);
        const bool jacobian = jacobians != nullptr && jacobians[0] != nullptr;
        ObservationEvaluation trial;
        if (!evaluator_.evaluate(optimizer, jacobian, trial)) return false;
        std::copy(trial.weighted.begin(), trial.weighted.end(), residuals);
        if (jacobian) std::copy(trial.optimizer_gradient.begin(), trial.optimizer_gradient.end(), jacobians[0]);
        return true;
    }

private:
    const FitRequest& request_;
    ObservationEvaluator::Impl& evaluator_;
};

const char* termination_name(ceres::TerminationType termination) {
    switch (termination) {
    case ceres::CONVERGENCE: return "convergence";
    case ceres::NO_CONVERGENCE: return "no_convergence";
    case ceres::FAILURE: return "failure";
    case ceres::USER_SUCCESS: return "user_success";
    case ceres::USER_FAILURE: return "user_failure";
    }
    return "failure";
}

void covariance(FitResult& result, const FitRequest& request, std::size_t n,
    const std::vector<double>& point, const std::vector<double>& lower, const std::vector<double>& upper) {
    auto& cov = result.covariance;
    const auto& set = request.coordinates;
    const auto width = static_cast<Eigen::Index>(point.size());
    cov.assumption = "unit-variance weighted residuals, independent errors and a regular local "
                     "linearization at an interior point; multiply by variance_factor for an "
                     "unknown common variance; no claim under bounds, rank deficiency or model discrepancy";
    if (request.uncertainty) cov.assumption = std::string("known absolute ")
        + (request.uncertainty->covariance || request.uncertainty->correlation ? "joint Gaussian covariance" : "independent Gaussian sigmas")
        + "; local Gaussian linearization, with reconciled inputs eliminated; no common variance estimate";
    for (std::size_t i = 0; i < n; ++i)
        if (request.weights.empty() || request.weights[i] > 0.0) cov.objective_rows.push_back(i);
    if (result.reconciliation) for (std::size_t c = 0; c < result.reconciliation->corrections.size(); ++c)
        cov.objective_rows.push_back(n + c);
    if (cov.objective_rows.empty()) { cov.status = "rank_deficient"; return; }
    Eigen::MatrixXd jacobian(static_cast<Eigen::Index>(cov.objective_rows.size()), width);
    for (Eigen::Index i = 0; i < jacobian.rows(); ++i)
        for (Eigen::Index c = 0; c < width; ++c)
            jacobian(i, c) = (result.reconciliation ? result.reconciliation->jacobian : result.optimizer_jacobian)
                [cov.objective_rows[static_cast<std::size_t>(i)] * static_cast<std::size_t>(width) + static_cast<std::size_t>(c)];
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(jacobian, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const auto& singular = svd.singularValues();
    cov.singular_values.assign(singular.data(), singular.data() + singular.size());
    if (request.singular_vectors) {
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> right = svd.matrixV(), left = svd.matrixU();
        cov.right_vectors.assign(right.data(), right.data() + right.size());
        cov.left_vectors.assign(left.data(), left.data() + left.size());
    }
    const double floor = policy::FIT_RANK_RECIPROCAL_CONDITION.value * (singular.size() > 0 ? singular(0) : 0.0);
    cov.rank = static_cast<std::size_t>(std::count_if(cov.singular_values.begin(), cov.singular_values.end(),
        [&](double value) { return value > floor && value > 0.0; }));
    result.rank_records.push_back({"weighted optimizer observation Jacobian", point.size(), cov.rank, floor});
    core::ConditioningObservation determination{"final_weighted_jacobian", "", "2", "final_SVD",
        cov.rank == point.size() ? "" : "rank_deficient", "endpoint"};
    determination.available = cov.rank == point.size();
    determination.condition = determination.available ? singular(0) / singular(singular.size() - 1) : INFINITY;
    core::observe_conditioning(result.rank_records.back(), std::move(determination), point);
    if (cov.rank < point.size()) core::retain_first(result.first_numerical_failure, [&] {
        core::Convergence record; record.residual = NAN; record.stop = core::Stop::singular;
        record.cause = "parameter Jacobian rank deficient; normal inverse unavailable";
        record.rank = std::make_shared<core::NumericalRank>(result.rank_records.back()); return record;
    });
    if (auto scope = core::Scope::detailed("fit_rank")) {
        scope->probe("singular_values", cov.singular_values);
        for (const auto row : cov.objective_rows) scope->probe("objective_rows", static_cast<double>(row));
        for (std::size_t c = 0; c < point.size(); ++c) scope->probe("optimizer_columns", static_cast<double>(c));
    }
    // A bound is active in optimizer space, within the fit's own active distance of the values it bounds.
    for (std::size_t c = 0; c < point.size(); ++c)
        if (std::min(point[c] - lower[c], upper[c] - point[c]) < policy::FIT_ACTIVE_BOUND_DISTANCE.value)
            cov.active_bounds.push_back(c);
    // The variance factor counts only rows that enter the objective (w > 0).
    if (!request.uncertainty && cov.objective_rows.size() > cov.rank) cov.variance_factor = 2.0 * result.final_cost / static_cast<double>(cov.objective_rows.size() - cov.rank);
    // An unfinished iterate keeps its rank, singular values and bound diagnostics but is not a
    // parameter-covariance estimate.
    if (result.status != "converged") { cov.status = "not_converged"; cov.variance_factor.reset(); return; }
    if (cov.rank < point.size()) { cov.status = "rank_deficient"; return; }
    if (!cov.active_bounds.empty()) { cov.status = "active_bound"; return; }
    const Eigen::MatrixXd complete = svd.matrixV() * singular.cwiseAbs2().cwiseInverse().asDiagonal()
        * svd.matrixV().transpose();
    const Eigen::MatrixXd optimizer = complete.topLeftCorner(set.coordinates.size(), set.coordinates.size());
    const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> map(
        result.chart.map.data(), static_cast<Eigen::Index>(result.physical.size()), static_cast<Eigen::Index>(set.coordinates.size()));
    const Eigen::MatrixXd physical = map * optimizer * map.transpose();
    cov.optimizer.assign(optimizer.data(), optimizer.data() + optimizer.size());
    cov.physical.assign(physical.data(), physical.data() + physical.size());
    cov.status = finite(cov.optimizer) && finite(cov.physical) ? "available" : "unavailable";
    if (cov.status == "unavailable") { cov.optimizer.clear(); cov.physical.clear(); }
}

void validate_checks(const FitRequest& request) {
    for (const auto& group : request.checks) {
        require(!group.rows.empty() && (group.maximum || group.mean || group.report_only), "check group rules");
        for (const auto threshold : {group.maximum, group.mean})
            require(!threshold || (std::isfinite(*threshold) && *threshold >= 0.0), "check threshold");
        require(!group.baseline || request.baseline.has_value(), "missing baseline parameters");
        require(group.metric == CheckMetric::AbsoluteDifference || group.metric == CheckMetric::AbsoluteLog10Ratio, "check metric");
        std::optional<bool> energy;
        for (const auto row : group.rows) {
            require(row < request.observations.size(), "check row index");
            const auto kind = request.observations[row].kind;
            require(electrolyte_observation(kind), "check requires an electrolyte observation");
            const bool current = kind == ObservationKind::SolvationGibbsEnergy || kind == ObservationKind::TransferGibbsEnergy;
            require(group.metric != CheckMetric::AbsoluteDifference || !energy || *energy == current, "mixed check units");
            energy = current;
        }
    }
}

}  // namespace

std::vector<CheckResult> evaluate_checks(const FitRequest& request, const Parameters& stored, FitResult* final,
    const std::vector<Observation>* reconciled) {
    validate_checks(request);
    std::vector<bool> selected(request.observations.size(), false);
    const auto refreshed = [&](std::size_t i) {
        return final && (final->statuses[i] == ObservationStatus::NotEvaluated || (!request.weights.empty() && request.weights[i] == 0.0));
    };
    for (const auto& group : request.checks) for (const auto i : group.rows) selected[i] = true;
    for (std::size_t i = 0; i < selected.size(); ++i) if (refreshed(i)) selected[i] = true;
    const auto predictions = [&](const HelmholtzModel& owner, const std::vector<Observation>& observations) {
        std::vector<Prediction> result(request.observations.size());
        ReferenceCache references;
        for (std::size_t i = 0; i < selected.size(); ++i) if (selected[i]) {
            validate(observations[i], owner);
            result[i] = predict(owner, observations[i], {}, {}, &references);
        }
        return result;
    };
    const auto candidate = predictions(HelmholtzModel(stored, {}, request.thermochemistry, request.controls.solver_options), request.observations);
    const auto reconciled_values = reconciled ? predictions(HelmholtzModel(stored, {}, request.thermochemistry, request.controls.solver_options), *reconciled) : std::vector<Prediction>{};
    // Checks use stored parameters; objective rows retain their accepted warm-start evaluation.
    for (std::size_t i = 0; i < candidate.size(); ++i) if (refreshed(i)) {
        const auto& prediction = reconciled ? reconciled_values[i] : candidate[i];
        final->predictions[i] = prediction.value.value_or(NAN);
        final->residuals[i] = prediction.residual.value_or(NAN);
        final->statuses[i] = prediction.status;
    }
    const auto baseline = request.baseline ? predictions(HelmholtzModel(*request.baseline, {}, request.thermochemistry, request.controls.solver_options), request.observations) : std::vector<Prediction>{};
    std::vector<CheckResult> reports;
    for (const auto& group : request.checks) {
        CheckResult out;
        out.name = group.name; out.rows = group.rows; out.report_only = group.report_only;
        out.status = "pass"; out.argmax = group.rows.front(); out.maximum_threshold = group.maximum; out.mean_threshold = group.mean;
        double baseline_sum = 0.0, baseline_max = 0.0;
        for (const auto i : group.rows) {
            const auto& row = request.observations[i];
            const auto unit = group.metric == CheckMetric::AbsoluteLog10Ratio ? "decade" : candidate[i].unit;
            out.unit = unit;
            const auto [e, slope] = check_error(group.metric, row, candidate[i].value.value_or(NAN));
            out.errors.push_back(e);
            out.slopes.push_back(slope);
            if (!candidate[i].residual || !std::isfinite(e)) { out.status = "unavailable"; out.failed_rows.push_back(i); }
            if (std::abs(e) > out.maximum) { out.maximum = std::abs(e); out.argmax = i; }
            out.mean += std::abs(e);
            if (group.baseline) {
                const double b = check_error(group.metric, row, baseline[i].value.value_or(NAN)).first;
                if (!baseline[i].residual || !std::isfinite(b)) { out.status = "unavailable"; out.failed_rows.push_back(i); }
                baseline_max = std::max(baseline_max, std::abs(b)); baseline_sum += std::abs(b);
            }
        }
        out.mean /= group.rows.size();
        if (group.baseline) {
            out.baseline_maximum = baseline_max; out.baseline_mean = baseline_sum / group.rows.size();
            if (group.maximum) out.maximum_threshold = std::max(*group.maximum, baseline_max);
            if (group.mean) out.mean_threshold = std::max(*group.mean, *out.baseline_mean);
        }
        if (out.status != "unavailable") {
            for (std::size_t r = 0; r < out.rows.size(); ++r)
                if (out.maximum_threshold && std::abs(out.errors[r]) > *out.maximum_threshold) out.failed_rows.push_back(out.rows[r]);
            if (!out.failed_rows.empty() || (out.mean_threshold && out.mean > *out.mean_threshold)) out.status = "fail";
        }
        reports.push_back(std::move(out));
    }
    return reports;
}

std::vector<Observation> observations_at(const FitRequest& request, const std::vector<double>& optimizer,
    const std::vector<double>& errors) {
    require(errors.size() == request.observations.size(), "manufactured errors width");
    std::vector<std::size_t> rows(errors.size());
    std::iota(rows.begin(), rows.end(), 0);
    FitResult failures;
    ObservationEvaluator evaluator(request, failures, rows);
    ObservationEvaluation values;
    require(evaluator.evaluate(optimizer, false, values), "manufactured prediction unavailable");
    auto observations = request.observations;
    request.coordinates.apply(observations, optimizer);
    for (std::size_t i = 0; i < errors.size(); ++i)
        observations[i].observed = observations[i].form == ResidualForm::LogRatio
            ? values.predictions[i] * std::pow(10.0, -errors[i]) : values.predictions[i] - errors[i];
    return observations;
}

std::pair<double, double> check_error(CheckMetric metric, const Observation& row, double prediction) {
    if (!std::isfinite(prediction)) return {NAN, NAN};
    if (metric == CheckMetric::AbsoluteDifference) return {prediction - row.observed, 1.0};
    if (!(prediction > 0.0 && row.observed > 0.0)) return {NAN, NAN};
    return {std::log10(prediction / row.observed), 1.0 / (std::log(10.0) * prediction)};
}
void publish_candidate(FitResult& result, ObservationEvaluation final, const std::vector<double>& point) {
    result.optimizer = point;
    result.reconciliation = std::move(final.reconciliation);
    if (result.reconciliation) {
        const auto parameters = point.size() - result.reconciliation->corrections.size(), n = final.predictions.size();
        result.optimizer.resize(parameters);
        std::vector<double> columns(n * parameters);
        for (std::size_t i = 0; i < n; ++i) std::copy_n(final.optimizer_gradient.begin() + i * point.size(), parameters, columns.begin() + i * parameters);
        final.optimizer_gradient = std::move(columns); final.weighted.resize(n);
    }
    result.chart = std::move(final.chart);
    result.physical = result.chart.physical;
    result.parameters = result.chart.parameters;
    result.saturation_predictions = std::move(final.saturation_predictions);
    result.predictions = std::move(final.predictions);
    result.residuals = std::move(final.residuals);
    result.statuses = std::move(final.statuses);
    result.physical_jacobian = std::move(final.gradient);
    result.weighted_residuals = std::move(final.weighted);
    result.optimizer_jacobian = std::move(final.optimizer_gradient);
    result.final_cost = final.cost;
}
void finish_checks(FitResult& result, const FitRequest& request, const std::vector<Observation>* reconciled) {
    auto check_request = request;
    request.coordinates.apply(check_request.observations, result.optimizer);
    result.checks = evaluate_checks(check_request, result.parameters, &result, reconciled);
    result.acceptance = "pass";
    for (const auto& check : result.checks)
        if (!check.report_only && check.status != "pass" && result.acceptance != "unavailable") result.acceptance = check.status;
}

void validate_fit_request(const FitRequest& request, bool optimizing) {
    const auto& set = request.coordinates;
    const auto n = request.observations.size();
    validate_checks(request);
    validate(set, request.parameters, request.observations, request.thermochemistry);
    const HelmholtzModel model(request.parameters, set.active_parameters(), request.thermochemistry, request.controls.solver_options);
    require(n > 0, "needs observations");
    core::Scope validation("fit_validation");
    for (const auto& observation : request.observations) {
        if (interrupted()) throw FitInterrupted{};
        validate(observation, model);
    }
    require(request.weights.empty() || request.weights.size() == n, "weights width");
    for (std::size_t i = 0; i < request.weights.size(); ++i) {
        require(std::isfinite(request.weights[i]) && request.weights[i] >= 0.0, "weight");
    }
    const auto& controls = request.controls;
    (void)controls.solver_options.resolved({});
    if (optimizing) {
        require(controls.maximum_iterations > 0 && std::isfinite(controls.maximum_elapsed_time_seconds)
            && controls.maximum_elapsed_time_seconds > 0.0,
            "controls: maximum_iterations and maximum_elapsed_time_seconds must be set and positive");
    }
}

FitResult fit(const FitRequest& request) {
    if (request.uncertainty && request.uncertainty->groups.empty() && request.uncertainty->output_sigma.empty()
        && !request.uncertainty->covariance && !request.uncertainty->correlation) {
        auto canonical = request; canonical.uncertainty.reset(); return fit(canonical);
    }
    if (request.uncertainty && request.minimax) { FitResult result; result.solver_options = request.controls.solver_options.resolved({}); result.message = "fit uncertainty with minimax is unsupported"; return result; }
    if (request.minimax) return fit_minimax(request);
    FitResult result; result.solver_options = request.controls.solver_options.resolved({});
    const auto& set = request.coordinates;
    const std::size_t n = request.observations.size();
    std::vector<double> sqrt_weights(n, 1.0);
    const auto groups = state_groups(request.observations);
    std::unique_ptr<ObservationEvaluator::Impl> owned;
    try {
        validate_fit_request(request);
        for (std::size_t i = 0; i < request.weights.size(); ++i) sqrt_weights[i] = std::sqrt(request.weights[i]);
        owned = std::make_unique<ObservationEvaluator::Impl>(request, result, groups, sqrt_weights);
    } catch (const std::invalid_argument& error) {
        core::retain_failure(result.first_numerical_failure, error, "fit setup");
        result.message = error.what();
        return result;
    }
    auto& evaluator = *owned;
    const std::size_t width = evaluator.width();
    evaluator.start_workers(request.controls.workers);
    std::vector<double> optimizer = set.start(request.parameters, request.observations, request.thermochemistry);
    optimizer.resize(width, 0.0);
    ObservationEvaluation initial;
    if (!evaluator.evaluate(optimizer, true, initial)) {
        result.status = "initial_evaluation_failed";
        result.message = result.failures.back().message;
        result.statuses = std::move(initial.statuses);
        return result;
    }
    result.initial_cost = initial.cost;

    Cost cost(request, evaluator);
    ceres::Problem::Options problem_options;
    problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    const auto [lower, upper] = evaluator.bounds();
    ceres::Solver::Options options;
    options.minimizer_type = ceres::TRUST_REGION;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.linear_solver_type = ceres::DENSE_QR;
    options.dense_linear_algebra_library_type = ceres::EIGEN;
    options.num_threads = 1;
    options.logging_type = ceres::SILENT;
    options.minimizer_progress_to_stdout = false;
    const auto& controls = request.controls;
    options.function_tolerance = controls.solver_options.value(policy::function);
    options.gradient_tolerance = controls.solver_options.value(policy::gradient);
    options.parameter_tolerance = controls.solver_options.value(policy::parameter);
    options.update_state_every_iteration = true;
    // Bounds as an active set (a projected step alone stalls against a bound it wants to cross): a
    // coordinate on a bound (within FIT_ACTIVE_BOUND_DISTANCE) whose objective gradient points out of the
    // box is held where it is for a round. A round ends on Ceres convergence, accepted new bound contact,
    // or an accepted coordinate joining the set; the set is rebuilt from the cached Jacobian, g = J^T r
    // through the same evaluation, checks and failure records as a Ceres trial. Joins are added; a held
    // coordinate whose gradient turns inward is released only after a Ceres convergence report or when none
    // would stay free (Nocedal & Wright 2006, Alg. 16.3): its sign is a multiplier only where free gradients vanish.
    std::vector<double> gradient;
    ceres::Solver::Summary summary;
    std::vector<int> held;
    int remaining = controls.maximum_iterations;
    double seconds = controls.maximum_elapsed_time_seconds;
    bool settled = false, all_held = false, exhausted = false;
    int round = 0;
    for (bool first = true;; first = false) {
        if (!evaluator.gradient_at(optimizer, gradient)) {
            summary.termination_type = ceres::FAILURE;
            summary.message = result.failures.back().message;
            break;
        }
        const auto next = select_held(optimizer, gradient, lower, upper);
        // Unchanged membership retains Ceres' stopping criteria; it does not certify the free gradient.
        all_held = next.size() == width;
        settled = (!first && summary.termination_type == ceres::CONVERGENCE && next == held) || all_held;
        exhausted = !settled && (remaining <= 0 || seconds <= 0.0);
        if (settled || exhausted) break;
        std::vector<int> kept;
        std::set_union(held.begin(), held.end(), next.begin(), next.end(), std::back_inserter(kept));
        held = summary.termination_type == ceres::CONVERGENCE || kept.size() == width ? next : kept;
        core::Scope round_scope("fit_round");
        round_scope.probe("round", round++);
        round_scope.observe([&] {
            std::vector<double> mask(width, 0.0);
            for (const int c : held) mask[c] = 1.0;
            round_scope.probe("held", mask);
        });
        ceres::Problem problem(problem_options);
        problem.AddResidualBlock(&cost, nullptr, optimizer.data());
        auto* box = new BoxChart(optimizer, lower, upper, held, evaluator);  // owned by the problem
        problem.SetManifold(optimizer.data(), box);
        options.callbacks = {box};
        options.max_num_iterations = remaining;
        options.max_solver_time_in_seconds = seconds;
        ceres::Solve(options, &problem, &summary);
        round_scope.probe("termination", static_cast<double>(summary.termination_type));
        for (std::size_t k = first ? 0 : 1; k < summary.iterations.size(); ++k) {
            result.iteration_seconds.push_back(summary.iterations[k].iteration_time_in_seconds);
            result.iteration_costs.push_back(summary.iterations[k].cost);
        }
        remaining -= static_cast<int>(summary.iterations.size()) - 1;
        seconds -= summary.total_time_in_seconds;
        if (summary.termination_type != ceres::CONVERGENCE && summary.termination_type != ceres::USER_SUCCESS) break;
    }
    result.termination = settled ? "convergence" : exhausted ? "no_convergence" : termination_name(summary.termination_type);
    result.iterations = result.iteration_costs.empty() ? 0 : result.iteration_costs.size() - 1;  // entry 0 is the start
    result.message = all_held ? "every coordinate is held at an active bound"
        : exhausted ? (remaining <= 0 ? "iteration limit reached" : "time limit reached") : summary.message;
    result.status = settled ? "converged" : summary.termination_type == ceres::FAILURE
        || summary.termination_type == ceres::USER_FAILURE ? "failed" : "limit_reached";

    result.optimizer = optimizer;
    ObservationEvaluation final;
    if (!set.within_bounds(std::vector<double>(optimizer.begin(), optimizer.begin() + set.coordinates.size()))
        || !evaluator.evaluate(optimizer, true, final)) {
        result.status = "failed";
        result.message = "the final point could not be evaluated: " + result.message;
        result.statuses = std::move(final.statuses);
        return result;
    }
    publish_candidate(result, std::move(final), optimizer);
    result.usable = result.status != "failed" && finite(result.physical);
    if (result.usable) {
        covariance(result, request, n, optimizer, lower, upper);
        if (request.uncertainty) {
            auto reconciled = request.observations;
            request.coordinates.apply(reconciled, result.optimizer);
            evaluator.measurement_.apply(request, optimizer, reconciled);
            finish_checks(result, request, &reconciled);
        } else finish_checks(result, request);
        // A limit stop leaves the round-start gradient stale; this is a cache hit on the final evaluation.
        if (!evaluator.gradient_at(optimizer, gradient)) result.covariance.projected_gradient_norm = NAN;
        else for (std::size_t c = 0; c < width; ++c)
            result.covariance.projected_gradient_norm = std::max(result.covariance.projected_gradient_norm,
                std::abs(optimizer[c] - std::clamp(optimizer[c] - gradient[c], lower[c], upper[c])));
    }
    return result;
}

std::vector<FitResult> fit_sequence(FitRequest request, const std::vector<std::vector<double>>& starts,
    int iterations, double seconds) {
    for (const auto& coordinate : request.coordinates.coordinates)
        require(!coordinate.reaction || coordinate.reaction->coefficient != ReactionCoefficient::EnthalpyShift,
            "reaction enthalpy shift is not supported in fit_sequence (starts)");
    require(iterations > 0 && std::isfinite(seconds) && seconds > 0.0, "shared fit budgets");
    const auto begin = std::chrono::steady_clock::now();
    std::vector<FitResult> results;
    const auto controls = request.controls;
    for (const auto& start : starts) {
        const double remaining = seconds - std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        if (iterations <= 0 || remaining <= 0.0) break;
        auto trial = request;
        trial.controls.maximum_iterations = std::min(controls.maximum_iterations, iterations);
        trial.controls.maximum_elapsed_time_seconds = std::min(controls.maximum_elapsed_time_seconds, remaining);
        try {
            trial.parameters = evaluate_coordinates(request.parameters, request.coordinates, start, request.holds, request.thermochemistry, request.controls.solver_options).parameters;
            request.coordinates.apply(trial.observations, start);
            results.push_back(fit(trial));
        } catch (const std::invalid_argument& error) {
            FitResult failed; failed.solver_options = request.controls.solver_options.resolved({}); failed.optimizer = start; failed.message = error.what(); results.push_back(std::move(failed));
        }
        iterations -= static_cast<int>(results.back().iterations);
    }
    return results;
}

}  // namespace epcsaft::regression
