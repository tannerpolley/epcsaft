#include "regression/fit.hpp"

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"
#include "equilibrium/stoichiometry.hpp"

#include <set>
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

// Instantiated once in equilibrium/singular_value_decompositions.cpp.
extern template class Eigen::JacobiSVD<Eigen::MatrixXd>;

namespace epcsaft::regression {
namespace {

const char* compute_cpu_directory = std::getenv("EPCSAFT_FIT_CPU_PROFILE");

void require(bool condition, const char* what) {
    if (!condition) throw std::invalid_argument(std::string("fit ") + what);
}

// The same point-local D maps raw observation columns and covariance.
void chain(const FitRequest& request, const std::vector<double>& map, double slope, const double* gradient, double* out) {
    const auto columns = request.coordinates.coordinates.size();
    if (columns == 0) return;
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

// The committed record is `optimizer`, `evaluation` and `jacobian`, written only by publish(). The rest is working
// state, which a same-point evaluation may update before its verdict.
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
// the reported index (the observation count for a refused trial). Values and the trial's first
// numerical diagnostic travel with the row evaluation.
struct Failure final {
    std::size_t row = 0, observation = 0;
    ObservationStatus status = ObservationStatus::Nonfinite;
    bool jacobian = false;
    std::string message;
    std::vector<ResidualEnergyAction> parameter_actions{};
    std::optional<core::Convergence> convergence{};
};

// One evaluated trial's disposition, decided once and only by ObservationEvaluator::Impl::evaluate. Evaluated: every
// process computed its rows and none failed. Rejected: `failure` is the first failure and no process publishes.
// Committed: entered only through publish(), in each process. A request the committed record already answers is
// reused before any of this. Workers receive the verdict as one byte: commit or not.
enum class Disposition { Evaluated, Rejected, Committed };
struct Verdict final {
    Disposition disposition = Disposition::Evaluated;
    Failure failure{};
};

FitHost fit_host;

bool interrupted() { return fit_host.interrupted && fit_host.interrupted(); }
#include "fit_worker.hpp"

}  // namespace

StabilityTransition stability_transition(std::size_t row, const std::vector<double>& previous_point, const Prediction& previous,
    const std::vector<double>& point, const Prediction& current) {
    StabilityTransition out; out.observation = row; out.previous_optimizer = previous_point; out.optimizer = point;
    const auto classify = [](const Prediction& p) { return !p.stability || !p.stability->available ? "unavailable"
        : p.stability->report.no_descent_in_trial_set ? "no_witness" : "resolved_descent"; };
    out.previous_classification = classify(previous); out.classification = classify(current);
    const auto minimum = [](const Prediction& p) -> std::optional<equilibrium::StabilityTrial> {
        if (!p.stability || !p.stability->minimizing_trial) return {}; return p.stability->report.trials.at(*p.stability->minimizing_trial); };
    out.previous_trial = minimum(previous); out.trial = minimum(current);
    const auto seeds = [](const equilibrium::StabilityMeasure& m) {
        std::set<std::string> out;
        for (auto i : m.minima) out.insert(m.report.trials.at(i).seed);
        return out;
    };
    const auto add = [&](const std::string& kind) { if (!out.kind.empty()) out.kind += "; "; out.kind += kind; };
    if (out.previous_classification != out.classification) add("classification");
    if (previous.stability && current.stability) {
        const auto& a = *previous.stability; const auto& b = *current.stability;
        if (a.report.completed && b.report.completed && a.report.no_descent_in_trial_set != b.report.no_descent_in_trial_set) add("resolution");
        if (a.available && b.available) {
            if (out.previous_trial && out.trial) {
                const auto& x = *out.previous_trial; const auto& y = *out.trial;
                if (x.phase.support != y.phase.support || x.phase.kind != y.phase.kind
                    || x.quantities.at("density_side") != y.quantities.at("density_side") || seeds(a) != seeds(b)) add("minimizer");
            }
            const auto& x = a.report.present_phases; const auto& y = b.report.present_phases;
            bool parent_changed = x.size() != y.size();
            for (std::size_t i = 0; !parent_changed && i < x.size(); ++i)
                parent_changed = x[i].at("density_side") != y[i].at("density_side");
            if (previous.compiled && current.compiled) {
                const auto& p = previous.compiled->problem.phases; const auto& q = current.compiled->problem.phases;
                parent_changed = parent_changed || p.size() != q.size();
                for (std::size_t i = 0; !parent_changed && i < p.size(); ++i)
                    parent_changed = p[i].support != q[i].support || p[i].kind != q[i].kind;
            }
            if (parent_changed) add("parent branch");
        }
        if (!b.event.empty()) add("current point: " + b.event);
    }
    return out;
}

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
    Impl(FitRequest request, FitResult& result,
              const std::vector<double>& sqrt_weights, std::optional<MeasurementChart> measurement = {})
        : needs_compiled_admission_(!measurement), request_([&] {
            if (!measurement) {
                const HelmholtzModel model(request.parameters, request.coordinates.active_parameters(), request.thermochemistry, request.controls.solver_options);
                for (const auto& observation : request.observations) {
                    (void)validate_declaration(observation, model);
                }
            }
            request.controls.solver_options = request.controls.solver_options.resolved({});
            request.controls.solver_options.nested_solve = true;
            for (auto& o : request.observations) {
                o.problem.solver_options = o.problem.solver_options.resolved(request.controls.solver_options);
                o.final_problem.solver_options = o.final_problem.solver_options.resolved(request.controls.solver_options);
            } return request; }()), result_(result), measurement_(measurement ? std::move(*measurement) : MeasurementChart(request_)),
          groups_(state_groups(request_.observations, request_.uncertainty ? &measurement_.membership : nullptr)), sqrt_weights_(sqrt_weights),
          starts_(request_.observations.size()), owner_(request_.observations.size(), 0) {
        if (request.uncertainty) for (std::size_t i = 0; i < sqrt_weights_.size(); ++i)
            if (sqrt_weights_[i] > 0.0 && request.observations[i].kind != ObservationKind::Stability) sqrt_weights_[i] = request.observations[i].scale / request.uncertainty->output_sigma[i];
    }
    std::size_t width() const { return request_.coordinates.coordinates.size() + measurement_.width(); }
    std::pair<std::vector<double>, std::vector<double>> bounds() const {
        auto box = request_.coordinates.optimizer_bounds();
        measurement_.bounds(box);
        return box;
    }

    // Fixed partition: worker k evaluates every count-th state group in row order. Electrolyte rows share
    // one reference cache and stay here, as does every row under a trace capture or in a native host.
    void start_workers(int requested) {
        interruptible_ = true;
        const std::size_t count = core::current_trace() || fit_host.worker.empty() ? 1
            : requested > 0 ? static_cast<std::size_t>(requested) : processors();
        std::size_t rank = 0;
        for (std::size_t i = 0; i < groups_.size(); ++i) {
            owner_[i] = groups_[i] != i ? owner_[groups_[i]] : electrolyte_observation(request_.observations[i].kind) ? 0 : rank++ % count;
            const auto& current = request_.observations[i];
            if (current.kind != ObservationKind::Stability || !current.pressure_reference) continue;
            auto reference = pressure_reference_row(current);
            for (std::size_t j = 0; j < i; ++j) if (request_.observations[j].kind == ObservationKind::Stability && request_.observations[j].pressure_reference) {
                auto prior = pressure_reference_row(request_.observations[j]);
                if (same_state(reference, prior)) { owner_[i] = owner_[j]; --rank; break; }
            }
        }
        for (std::size_t k = 1; k < std::min(count, rank); ++k) workers_.start();
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
        Trial fresh;
        Trial& trial = same_point ? *cached_ : fresh;
        Verdict verdict = assemble(optimizer, jacobian, out, trial, same_point, scope);
        for (const Channel fd : workers_.sockets) put(fd, verdict.disposition == Disposition::Evaluated);
        if (verdict.disposition == Disposition::Rejected) {
            if (!out.saturation_predictions.empty()) result_.saturation_predictions = out.saturation_predictions;
            if (!out.stability_predictions.empty()) result_.stability_predictions = out.stability_predictions;
            fail(optimizer, std::move(verdict.failure), out);
        } else verdict.disposition = publish(optimizer, jacobian, out, trial, fresh, same_point);
        return verdict.disposition == Disposition::Committed;
    }

    // Every process's rows of one trial merged in row order, then the weighted cost and optimizer columns: Evaluated,
    // or Rejected with the first failure. It sends no verdict and publishes nothing.
    Verdict assemble(const std::vector<double>& optimizer, bool jacobian, ObservationEvaluation& out, Trial& trial,
                     bool same_point, core::Scope& scope) {
        const std::size_t n = request_.observations.size(), width = request_.coordinates.tie().size() + request_.holds.size();
        auto failure = compute(optimizer, jacobian, out, trial, same_point, scope);
        for (std::size_t k = 1; k <= workers_.sockets.size(); ++k) {
            std::optional<Failure> theirs;
            ObservationEvaluation part;
            get(workers_.sockets[k - 1], theirs, part.predictions, part.residuals, part.statuses, part.gradient, part.input_gradient,
                part.state_densities, part.reference_values, part.saturation_predictions, part.stability_predictions, part.first_numerical_failure, part.first_numerical_failure_observation);
            // A trial refused before its rows sends none; any other reply carries every row.
            const bool rows = part.statuses.size() == n && part.predictions.size() == n && part.residuals.size() == n
                && part.gradient.size() == (jacobian ? n * width : 0)
                && part.input_gradient.size() == (jacobian ? n * measurement_.width() : 0);
            if (!rows && !(theirs && part.statuses.empty())) throw std::runtime_error("fit worker reply size");
            for (std::size_t i = 0; rows && out.statuses.size() == n && i < n; ++i) if (owner_[i] == k) {
                out.predictions[i] = part.predictions[i];
                out.residuals[i] = part.residuals[i];
                out.statuses[i] = part.statuses[i];
                if (!part.stability_predictions.empty()) { if (out.stability_predictions.empty()) out.stability_predictions.resize(n);
                    out.stability_predictions[i] = std::move(part.stability_predictions[i]); }
                if (!part.saturation_predictions.empty()) out.saturation_predictions[i] = std::move(part.saturation_predictions[i]);
                if (request_.uncertainty) { out.state_densities[i] = part.state_densities.at(i); out.reference_values[i] = part.reference_values.at(i); }
                if (jacobian) std::copy_n(part.gradient.begin() + i * width, width, out.gradient.begin() + i * width);
                if (jacobian && !measurement_.groups.empty()) std::copy_n(part.input_gradient.begin() + i * measurement_.width(),
                    measurement_.width(), out.input_gradient.begin() + i * measurement_.width());
            }
            if (theirs && (!failure || theirs->row < failure->row)) failure = std::move(theirs);
            if (part.first_numerical_failure && (!out.first_numerical_failure
                || part.first_numerical_failure_observation < out.first_numerical_failure_observation)) {
                out.first_numerical_failure = part.first_numerical_failure;
                out.first_numerical_failure_observation = part.first_numerical_failure_observation;
            }
        }
        // The serial sweep stops at its first failure: later rows stay unevaluated.
        for (std::size_t i = failure ? failure->row + 1 : n; i < out.statuses.size(); ++i) {
            out.predictions[i] = out.residuals[i] = 0.0;
            out.statuses[i] = ObservationStatus::NotEvaluated;
            if (!out.saturation_predictions.empty()) out.saturation_predictions[i].reset();
            if (!out.stability_predictions.empty()) out.stability_predictions[i].reset();
            if (jacobian) std::fill_n(out.gradient.begin() + i * width, width, 0.0);
        }
        if (!same_point && !out.stability_predictions.empty()) {
            if (stability_history_.empty()) stability_history_.resize(n);
            for (std::size_t i = 0; i < out.stability_predictions.size(); ++i) if (out.stability_predictions[i]) {
                auto& previous = stability_history_[i];
                if (!previous.first.empty()) {
                    auto event = stability_transition(i, previous.first, previous.second, optimizer, *out.stability_predictions[i]);
                    if (!event.kind.empty()) result_.stability_transitions.push_back(std::move(event));
                }
                previous = {optimizer, transition_record(*out.stability_predictions[i])};
            }
        }
        if (failure) return {Disposition::Rejected, std::move(*failure)};
        out.weighted.resize(n);
        for (std::size_t i = 0; i < n; ++i) out.weighted[i] = sqrt_weights_[i] * out.residuals[i];
        // One half of the row-ordered squared sum; reject its first nonfinite operand or sum.
        double sum = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double weighted = out.weighted[i];
            sum += weighted * weighted;
            if (!std::isfinite(weighted) || !std::isfinite(sum))
                return {Disposition::Rejected, {i, i, ObservationStatus::Nonfinite, jacobian,
                    std::isfinite(weighted) ? "accumulated weighted cost not finite" : "weighted residual not finite"}};
        }
        out.cost = 0.5 * sum;
        for (std::size_t i = 0; i < n; ++i) {
            const double cost = 0.5 * out.weighted[i] * out.weighted[i];
            (request_.observations[i].kind == ObservationKind::Stability ? out.stability_cost : out.data_cost) += cost;
        }
        const std::size_t parameters = request_.coordinates.coordinates.size(), latent = measurement_.width(), columns = this->width();
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
                return {Disposition::Rejected, {i, i, ObservationStatus::Nonfinite, true, "optimizer Jacobian not finite"}};
        }
        if (request_.uncertainty) if (const char* cause = measurement_.reconcile(out, optimizer, parameters, jacobian))
            return {Disposition::Rejected, {n, n, ObservationStatus::Nonfinite, jacobian, cause}};
        return {};
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
                out.state_densities, out.reference_values, out.saturation_predictions, out.stability_predictions, out.first_numerical_failure, out.first_numerical_failure_observation);
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
        return all_finite(out) || fail(point, {trial.weighted.size(), trial.weighted.size(), ObservationStatus::Nonfinite, true, "weighted objective gradient not finite"}, trial);
    }

    bool fail(const std::vector<double>& optimizer, Failure failure, ObservationEvaluation& out) {
        if (out.first_numerical_failure_observation > failure.row) out.first_numerical_failure.reset();
        auto first = out.first_numerical_failure;
        std::optional<double> central;
        const auto i = failure.observation;
        if (i < out.statuses.size() && out.statuses[i] != ObservationStatus::NotEvaluated && std::isfinite(out.predictions[i]))
            central = out.predictions[i];
        core::retain_failure(first, failure.message, "fit evaluation");
        core::retain_first(result_.first_numerical_failure, first);
        result_.failures.push_back({trials_, failure.observation, optimizer, failure.status, std::move(failure.message), failure.jacobian,
            std::move(failure.parameter_actions), std::move(first), central, std::move(failure.convergence)});
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

    std::vector<std::pair<std::vector<double>, Prediction>> stability_history_;
    // A separate copy of what stability_transition reads; the shared report is never edited.
    static Prediction transition_record(const Prediction& p) {
        Prediction out;
        if (p.compiled) { out.compiled.emplace(); out.compiled->problem.phases = p.compiled->problem.phases; }
        if (!p.stability) return out;
        const auto& s = *p.stability;
        auto& m = *(out.stability = std::make_shared<equilibrium::StabilityMeasure>());
        m.available = s.available; m.event = s.event; m.minimizing_trial = s.minimizing_trial; m.minima = s.minima;
        m.report.completed = s.report.completed; m.report.no_descent_in_trial_set = s.report.no_descent_in_trial_set;
        m.report.trials.resize(s.report.trials.size());
        for (const auto i : s.minima) m.report.trials.at(i).seed = s.report.trials.at(i).seed;
        if (s.minimizing_trial) m.report.trials.at(*s.minimizing_trial) = s.report.trials.at(*s.minimizing_trial);
        for (const auto& phase : s.report.present_phases) {
            auto& kept = m.report.present_phases.emplace_back();
            if (const auto side = phase.find("density_side"); side != phase.end()) kept.insert(*side);
        }
        return out;
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
            if (std::any_of(trial.observations.begin(), trial.observations.end(),
                [](const auto& o) { return o.kind == ObservationKind::Stability; })) out.stability_predictions.resize(n);
            out.residuals.assign(n, 0.0);
            out.statuses.assign(n, ObservationStatus::NotEvaluated);
            if (request_.uncertainty) { out.state_densities.resize(n); out.reference_values.resize(n); }
            if (jacobian) out.gradient.assign(n * width, 0.0);
            if (jacobian) out.input_gradient.assign(n * measurement_.width(), 0.0);
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
                    const bool unchanged = observation.kind == ObservationKind::Stability || !request_.uncertainty || (pressure_closed(observation.kind)
                        ? std::tie(observation.temperature, observation.pressure, observation.composition, observation.branch)
                            == std::tie(prior.temperature, prior.pressure, prior.composition, prior.branch)
                        : same_state(observation, prior));
                    std::vector<const Observation*> admission;
                    if (needs_compiled_admission_) for (std::size_t row = 0; row < n; ++row)
                        if (groups_[row] == g) admission.push_back(&trial.observations[row]);
                    prediction = predict(model, observation, starts_[g], std::move(prediction.compiled), &trial.references,
                        seed && unchanged ? &seed->centers[g] : nullptr, seed && unchanged ? seed->model.get() : nullptr,
                        seed && unchanged && trial.chart.physical == seed->chart.physical, deadline_, admission);
                    if (prediction.status != ObservationStatus::Available && starts_[g] && observation.kind != ObservationKind::SaturationCappedAdditionEnthalpy && observation.kind != ObservationKind::Stability) {
                        auto retry = predict(model, observation, std::nullopt, std::move(prediction.compiled), &trial.references,
                            nullptr, nullptr, false, deadline_, admission);
                        core::retain_first(retry.first_numerical_failure, prediction.first_numerical_failure, true);
                        prediction = std::move(retry);
                    }
                    // Keep the action operands of plausible trials through a same-point Jacobian request.
                    if (prediction.stability && (interruptible_ || process_)) {
                        auto measure = std::make_shared<equilibrium::StabilityMeasure>(*prediction.stability);
                        for (std::size_t t = 0; t < measure->report.trials.size(); ++t) if (t != measure->minimizing_trial) {
                            equilibrium::StabilityTrial kept; kept.seed = std::move(measure->report.trials[t].seed);
                            if (measure->available && measure->coincident_supports
                                && std::binary_search(measure->minima.begin(), measure->minima.end(), t)) {
                                auto& trial = measure->report.trials[t];
                                kept.composition = std::move(trial.composition);
                                for (const auto* name : {"volume_m3_per_mol", "numerical_plane_value_error_over_rt"})
                                    kept.quantities.insert(trial.quantities.extract(name));
                            }
                            measure->report.trials[t] = std::move(kept);
                        }
                        prediction.stability = std::move(measure);
                    }
                } else if (!pressure_closed(observation.kind)) project_prediction(model, observation, prediction);
                retain(prediction.first_numerical_failure);
                if (prediction.saturation) out.saturation_predictions[i] = prediction;
                if (observation.kind == ObservationKind::Stability) out.stability_predictions[i] = prediction;
                out.statuses[i] = prediction.status;
                out.predictions[i] = prediction.value.value_or(NAN);
                if (prediction.status != ObservationStatus::Available)
                    return Failure{i, i, prediction.status, jacobian, prediction.message, {}, prediction.refusal_convergence};
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
                    std::vector<ObservationInputDirection> directions;
                    for (std::size_t c = 0; c < measurement_.width(); ++c) if (measurement_.membership[i] == measurement_.groups[c]) {
                        const auto kind = request_.uncertainty->groups[measurement_.groups[c]].inputs[measurement_.inputs[c]].kind;
                        ObservationInputDirection d;
                        d.state.feed_amounts_mol = measurement_.composition[c];
                        d.state.active_parameters.assign(request_.coordinates.active_parameters().size() + request_.holds.size(), 0.0);
                        d.state.temperature_k = kind == MeasuredInputKind::Temperature ? 1.0 : 0.0;
                        d.state.pressure_pa = kind == MeasuredInputKind::Pressure ? 1.0 : 0.0;
                        d.molality = kind == MeasuredInputKind::Molality ? 1.0 : 0.0;
                        directions.push_back(std::move(d));
                    }
                    gradients[g] = pressure_closed(observation.kind) && directions.empty()
                        ? std::vector<PredictionGradient>{differentiate(model, observation, prediction)}
                        : differentiate_group(model, rows, prediction, request_.coordinates, directions,
                            [&](const double* gradient, double* mapped) { chain(request_, trial.chart.map, 1, gradient, mapped); });
                }
                const auto& gradient = gradients[g][position[g]++];
                if (prediction.stability && !gradient.raw_minimum.empty()) {
                    auto& mapped = out.stability_predictions[i]->stability_raw_gradient; mapped.resize(this->width());
                    chain(request_, trial.chart.map, 1, gradient.raw_minimum.data(), mapped.data());
                }
                retain(gradient.first_numerical_failure);
                if (gradient.status != ObservationStatus::Available)
                    return Failure{i, i, gradient.status, true, gradient.message, gradient.parameter_actions, gradient.refusal_convergence};
                std::copy_n(gradient.value.begin(), std::min(width, gradient.value.size()),
                          out.gradient.begin() + static_cast<std::ptrdiff_t>(i * width));
                std::size_t local = width;
                for (std::size_t c = 0; c < measurement_.width(); ++c) if (measurement_.membership[i] == measurement_.groups[c])
                    out.input_gradient[i * measurement_.width() + c] = gradient.value.at(local++);
            }
            return std::nullopt;
        } catch (const std::bad_alloc&) {
            throw;
        } catch (const std::exception& error) {
            if (!out.first_numerical_failure) out.first_numerical_failure_observation = i;
            retain_density_failure(out.first_numerical_failure, error, "fit trial");
            const auto* terminal = dynamic_cast<const core::ConvergenceFailure*>(&error);
            const auto row = i < out.statuses.size() && out.statuses[i] == ObservationStatus::Available ? i : n;
            return Failure{i, row, ObservationStatus::InadmissibleTrial, jacobian, error.what(), {},
                terminal ? std::optional(terminal->record) : std::nullopt};
        }
    }

    // Neither a failed row nor an overflowing weighted/coordinate chain publishes starts or a cache. The only commit:
    // an Evaluated trial becomes Committed here, in the fit process on evaluate's verdict and in a worker on its byte.
    Disposition publish(const std::vector<double>& optimizer, bool jacobian, ObservationEvaluation& out, Trial& trial,
                        Trial& fresh, bool same_point) {
        for (std::size_t g = 0; g < groups_.size(); ++g)
            if (trial.centers[g].solved) starts_[g] = trial.centers[g].solved;
        // After a derivative commit the cache answers every request at this point: no row here is recomputed.
        if (jacobian) {
            for (std::size_t g = 0; g < groups_.size(); ++g)
                if (request_.observations[g].kind == ObservationKind::Stability) trial.centers[g] = Prediction{};
            trial.references.pressures = {};
        }
        out.chart = trial.chart;
        trial.optimizer = optimizer;
        // A worker has no reuse path; a fit's own process keeps only its stability rows' transition fields.
        if (!process_) {
            core::retain_first(result_.first_numerical_failure, out.first_numerical_failure);
            std::vector<std::optional<Prediction>> rows;
            if (interruptible_) rows.swap(out.stability_predictions);
            trial.evaluation = out;
            for (const auto& row : rows)
                trial.evaluation.stability_predictions.push_back(row ? std::optional(transition_record(*row)) : std::nullopt);
            if (interruptible_) out.stability_predictions = std::move(rows);
        }
        trial.jacobian = jacobian;
        if (!same_point) {
            cached_.emplace(std::move(fresh));
        }
        return Disposition::Committed;
    }

    const bool needs_compiled_admission_;
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
    bool interruptible_ = false;  // a least-squares fit's own process: polls cancellation, compacts its committed rows
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
    FitResult result(request);
    ObservationEvaluator::Impl evaluator(request, result, sqrt_weights, std::move(measurement));
    get(fd, evaluator.owner_, evaluator.process_, evaluator.deadline_);
    evaluator.serve(fd);
}

ObservationEvaluator::ObservationEvaluator(const FitRequest& request, FitResult& failures, const std::vector<std::size_t>& rows)
    : ObservationEvaluator(request, failures, rows, false) {}
ObservationEvaluator::ObservationEvaluator(const FitRequest& request, FitResult& failures,
    const std::vector<std::size_t>& rows, bool admitted) {
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
    impl_ = std::make_shared<Impl>(canonical, failures, selection,
        admitted ? std::optional<MeasurementChart>(MeasurementChart(canonical)) : std::nullopt);
}
bool ObservationEvaluator::evaluate(const std::vector<double>& optimizer, bool columns, ObservationEvaluation& out) {
    out.residual_optimizer_gradient.clear();
    if (!impl_->evaluate(optimizer, columns, out)) return false;
    if (columns) out.residual_optimizer_gradient = out.optimizer_gradient;
    const auto& measurement = impl_->measurement_;
    const std::size_t parameters = impl_->request_.coordinates.coordinates.size(), latent = measurement.width();
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
        : evaluator_(evaluator) {
        set_num_residuals(static_cast<int>(request.observations.size() + evaluator.measurement_.width()));
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
    ObservationEvaluator::Impl& evaluator_;
};

FitCovarianceStatus covariance(FitResult& result, const FitRequest& request, std::size_t n, const std::vector<double>& point,
    const std::vector<double>& lower, const std::vector<double>& upper, FitTermination termination) {
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
    if (cov.objective_rows.empty()) return FitCovarianceStatus::RankDeficient;
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
    if (!request.uncertainty && cov.objective_rows.size() > cov.rank) cov.variance_factor = 2.0 * result.data_cost / static_cast<double>(cov.objective_rows.size() - cov.rank);
    // An unfinished iterate keeps its rank, singular values and bound diagnostics but is not a
    // parameter-covariance estimate.
    if (termination != FitTermination::Settled && termination != FitTermination::AllHeld) {
        cov.variance_factor.reset();
        return FitCovarianceStatus::NotConverged;
    }
    if (cov.rank < point.size()) return FitCovarianceStatus::RankDeficient;
    if (std::any_of(request.observations.begin(), request.observations.end(), [](const auto& o) { return o.kind == ObservationKind::Stability; })) {
        cov.variance_factor.reset();
        cov.assumption = "numerical rank includes requirement rows; no statistical covariance claim";
        return FitCovarianceStatus::StabilityPenalty;
    }
    if (!cov.active_bounds.empty()) return FitCovarianceStatus::ActiveBound;
    const Eigen::MatrixXd complete = svd.matrixV() * singular.cwiseAbs2().cwiseInverse().asDiagonal()
        * svd.matrixV().transpose();
    const Eigen::MatrixXd optimizer = complete.topLeftCorner(set.coordinates.size(), set.coordinates.size());
    const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> map(
        result.chart.map.data(), static_cast<Eigen::Index>(result.physical.size()), static_cast<Eigen::Index>(set.coordinates.size()));
    const Eigen::MatrixXd physical = map * optimizer * map.transpose();
    cov.optimizer.assign(optimizer.data(), optimizer.data() + optimizer.size());
    cov.physical.assign(physical.data(), physical.data() + physical.size());
    if (all_finite(cov.optimizer) && all_finite(cov.physical)) return FitCovarianceStatus::Available;
    cov.optimizer.clear(); cov.physical.clear();
    return FitCovarianceStatus::Unavailable;
}

void validate_checks(const FitRequest& request) {
    for (const auto& group : request.checks) {
        require(!group.rows.empty() && (group.maximum || group.mean || group.report_only), "check group rules");
        for (const auto threshold : {group.maximum, group.mean})
            require(!threshold || (std::isfinite(*threshold) && *threshold >= 0.0), "check threshold");
        require(!group.baseline || request.baseline.has_value(), "missing baseline parameters");
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
    const HelmholtzModel model(stored, {}, request.thermochemistry, request.controls.solver_options);
    const auto candidate = predictions(model, request.observations);
    const auto reconciled_values = reconciled ? predictions(model, *reconciled) : std::vector<Prediction>{};
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
    result.stability_predictions = std::move(final.stability_predictions);
    result.data_cost = final.data_cost; result.stability_cost = final.stability_cost;
    result.predictions = std::move(final.predictions);
    result.residuals = std::move(final.residuals);
    result.statuses = std::move(final.statuses);
    result.physical_jacobian = std::move(final.gradient);
    result.weighted_residuals = std::move(final.weighted);
    result.optimizer_jacobian = std::move(final.optimizer_gradient);
    result.final_cost = final.cost;
}
FitAcceptance finish_checks(FitResult& result, const FitRequest& request, const std::vector<Observation>* reconciled) {
    auto check_request = request;
    request.coordinates.apply(check_request.observations, result.optimizer);
    result.checks = evaluate_checks(check_request, result.parameters, &result, reconciled);
    auto acceptance = FitAcceptance::Pass;
    for (const auto& check : result.checks)
        if (!check.report_only && check.status != "pass" && acceptance != FitAcceptance::Unavailable)
            acceptance = check.status == "fail" ? FitAcceptance::Fail : FitAcceptance::Unavailable;
    return acceptance;
}

// The one place a fit's outcome becomes its labels. `message` keeps its owner's free text (a setup or Ipopt
// refusal, the Ceres summary, a trial failure; empty for PointOnly) unless the termination has a fixed text. A fit
// that ran Ceres and whose returned point has no committed evaluation is "failed", whatever its termination.
void label(FitResult& result, const FitOutcome& outcome) {
    using T = FitTermination;
    const T t = outcome.termination;
    const bool ipopt = t == T::IpoptSuccess || t == T::IpoptStop, converged = t == T::Settled || t == T::AllHeld;
    const bool ceres = !ipopt && t != T::NotRun && t != T::StartRejected && t != T::PointOnly;
    const bool absent = outcome.candidate == FitCandidate::Absent;
    result.status = t == T::NotRun ? "invalid_setup" : t == T::StartRejected ? "initial_evaluation_failed"
        : t == T::PointOnly ? "point_evaluated" : t == T::IpoptStop ? "limit_or_failure"
        : absent || t == T::Failure || t == T::CallbackFailure ? "failed"
        : ipopt || converged ? "converged" : "limit_reached";
    result.termination = ipopt ? result.minimax->solver_status : !ceres ? "not_run" : converged ? "convergence"
        : t == T::Failure ? "failure" : t == T::CallbackFailure ? "user_failure" : "no_convergence";
    if (t == T::AllHeld) result.message = "every coordinate is held at an active bound";
    if (t == T::IterationLimit) result.message = "iteration limit reached";
    if (t == T::TimeLimit) result.message = "time limit reached";
    if (ipopt) result.message = "local constrained estimate; global feasibility and optimality unresolved";
    if (ceres && absent) result.message = "the final point could not be evaluated: " + result.message;
    result.usable = outcome.candidate == FitCandidate::Usable;
    const FitAcceptance a = outcome.acceptance;
    result.acceptance = a == FitAcceptance::Pass ? "pass" : a == FitAcceptance::Fail ? "fail"
        : a == FitAcceptance::Unavailable ? "unavailable" : "not_evaluated";
}

void validate_fit_request(const FitRequest& request, bool optimizing) {
    const auto& set = request.coordinates;
    const auto n = request.observations.size();
    validate_checks(request);
    validate(set, request.parameters, request.observations, request.thermochemistry);
    const HelmholtzModel model(request.parameters, set.active_parameters(), request.thermochemistry, request.controls.solver_options);
    require(n > 0, "needs observations");
    core::Scope validation("fit_validation");
    for (const auto& observation : request.observations) try {
        if (interrupted()) throw FitInterrupted{};
        validate(observation, model);
        if (observation.kind == ObservationKind::Stability) require(!request.minimax, "stability penalty rows require least squares");
    } catch (const std::domain_error& error) {
        throw std::invalid_argument("observation " + std::to_string(&observation - request.observations.data())
            + ": central_solve_unavailable: " + error.what());
    } catch (const std::invalid_argument& error) {
        if (observation.kind != ObservationKind::AssociationSiteFraction && !observation.site) throw;
        throw std::invalid_argument("observation " + std::to_string(&observation - request.observations.data()) + ": " + error.what());
    }
    require(request.weights.empty() || request.weights.size() == n, "weights width");
    for (std::size_t i = 0; i < request.weights.size(); ++i) {
        require(std::isfinite(request.weights[i]) && request.weights[i] >= 0.0, "weight");
        if (request.observations[i].kind == ObservationKind::Stability) require(request.weights[i] > 0.0, "positive stability weight");
    }
    const auto& controls = request.controls;
    (void)controls.solver_options.resolved({});
    if (optimizing) {
        require(controls.maximum_iterations > 0 && std::isfinite(controls.maximum_elapsed_time_seconds)
            && controls.maximum_elapsed_time_seconds > 0.0,
            "controls: maximum_iterations and maximum_elapsed_time_seconds must be set and positive");
    }
}

std::vector<double> fit_start(const FitRequest& request) {
    auto point = request.coordinates.start(request.parameters, request.observations, request.thermochemistry);
    if (!request.start) return point;
    require(request.start->size() == point.size(), "fit start width");
    for (std::size_t c = 0; c < point.size(); ++c) {
        const auto& reaction = request.coordinates.coordinates[c].reaction;
        if (reaction && reaction->coefficient == ReactionCoefficient::EnthalpyShift) point[c] = (*request.start)[c];
    }
    require(request.coordinates.within_bounds(point), "start outside bounds");
    return point;
}

static FitResult fit_parameters(const FitRequest& request) {
    if (request.uncertainty && request.uncertainty->groups.empty() && request.uncertainty->output_sigma.empty()
        && !request.uncertainty->covariance && !request.uncertainty->correlation) {
        auto canonical = request; canonical.uncertainty.reset(); return fit_parameters(canonical);
    }
    if (request.uncertainty && request.minimax) { FitResult result(request); result.message = "fit uncertainty with minimax is unsupported"; return result; }
    if (request.minimax) return fit_minimax(request);
    FitResult result(request);
    const auto& set = request.coordinates;
    const std::size_t n = request.observations.size();
    std::vector<double> sqrt_weights(n, 1.0);
    std::unique_ptr<ObservationEvaluator::Impl> owned;
    try {
        validate_fit_request(request);
        for (std::size_t i = 0; i < request.weights.size(); ++i) sqrt_weights[i] = std::sqrt(request.weights[i]);
        owned = std::make_unique<ObservationEvaluator::Impl>(request, result, sqrt_weights, MeasurementChart(request));
    } catch (const std::invalid_argument& error) {
        core::retain_failure(result.first_numerical_failure, error, "fit setup");
        result.message = error.what();
        return result;
    }
    auto& evaluator = *owned;
    const std::size_t width = evaluator.width();
    evaluator.start_workers(request.controls.workers);
    std::vector<double> optimizer = fit_start(request);
    optimizer.resize(width, 0.0);
    FitOutcome outcome;
    if (ObservationEvaluation initial; evaluator.evaluate(optimizer, true, initial)) result.initial_cost = initial.cost;
    else {
        outcome.termination = FitTermination::StartRejected;
        result.message = result.failures.back().message;
        result.statuses = std::move(initial.statuses);
        label(result, outcome);
        return result;
    }
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
    int round = 0;
    using T = FitTermination;
    for (bool first = true;; first = false) {
        if (!evaluator.gradient_at(optimizer, gradient)) {
            outcome.termination = T::Failure;
            summary.message = result.failures.back().message;
            break;
        }
        const auto next = select_held(optimizer, gradient, lower, upper);
        // Unchanged membership retains Ceres' stopping criteria; it does not certify the free gradient. Until a
        // round stops the fit, the termination stays NotRun.
        outcome.termination = next.size() == width ? T::AllHeld
            : !first && summary.termination_type == ceres::CONVERGENCE && next == held ? T::Settled
            : remaining <= 0 ? T::IterationLimit : seconds <= 0.0 ? T::TimeLimit : T::NotRun;
        if (outcome.termination != T::NotRun) break;
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
        if (summary.termination_type == ceres::CONVERGENCE || summary.termination_type == ceres::USER_SUCCESS) continue;
        outcome.termination = summary.termination_type == ceres::NO_CONVERGENCE ? T::CeresLimit
            : summary.termination_type == ceres::USER_FAILURE ? T::CallbackFailure : T::Failure;
        break;
    }
    result.iterations = result.iteration_costs.empty() ? 0 : result.iteration_costs.size() - 1;  // entry 0 is the start
    result.message = summary.message;  // label() substitutes the fixed texts

    result.optimizer = optimizer;
    ObservationEvaluation final;
    if (!set.within_bounds(std::vector<double>(optimizer.begin(), optimizer.begin() + set.coordinates.size()))
        || !evaluator.evaluate(optimizer, true, final)) {
        result.statuses = std::move(final.statuses);
        label(result, outcome);
        return result;
    }
    publish_candidate(result, std::move(final), optimizer);
    outcome.candidate = outcome.termination != T::Failure && outcome.termination != T::CallbackFailure
        && all_finite(result.physical) ? FitCandidate::Usable : FitCandidate::Values;
    if (outcome.candidate == FitCandidate::Usable) {
        result.covariance.outcome = covariance(result, request, n, optimizer, lower, upper, outcome.termination);
        if (request.uncertainty) {
            auto reconciled = request.observations;
            request.coordinates.apply(reconciled, result.optimizer);
            evaluator.measurement_.apply(request, optimizer, reconciled);
            outcome.acceptance = finish_checks(result, request, &reconciled);
        } else outcome.acceptance = finish_checks(result, request);
        // A limit stop leaves the round-start gradient stale; this is a cache hit on the final evaluation.
        if (!evaluator.gradient_at(optimizer, gradient)) result.covariance.projected_gradient_norm = NAN;
        else for (std::size_t c = 0; c < width; ++c)
            result.covariance.projected_gradient_norm = std::max(result.covariance.projected_gradient_norm,
                std::abs(optimizer[c] - std::clamp(optimizer[c] - gradient[c], lower[c], upper[c])));
    }
    owned.reset();  // the replay reads only the result: release the evaluator and its workers first
    if (!result.stability_predictions.empty()) {
        auto observations = request.observations; request.coordinates.apply(observations, result.optimizer);
        const HelmholtzModel model(result.parameters, {}, request.thermochemistry, request.controls.solver_options);
        ReferenceCache references;
        for (std::size_t i = 0; i < observations.size(); ++i) if (observations[i].kind == ObservationKind::Stability) {
            auto fresh = predict(model, observations[i], {}, {}, &references);
            auto event = stability_transition(i, optimizer, *result.stability_predictions[i], optimizer, fresh);
            if (!event.kind.empty()) result.stability_transitions.push_back(std::move(event));
            result.stability_predictions[i] = fresh;
            if (!fresh.stability || !fresh.stability->available) {
                outcome.candidate = FitCandidate::Values;
                outcome.acceptance = FitAcceptance::Unavailable;
            } else if (!fresh.stability->report.no_descent_in_trial_set) outcome.acceptance = FitAcceptance::Fail;
        }
    }
    label(result, outcome);
    return result;
}

static void predict_fit(const FitRequest& request, const std::vector<Observation>& rows, FitResult& result) {
    const auto& covariance = result.covariance;
    // A published reconciliation is present for an absolute measurement model, including zero latent inputs.
    const auto scale = result.reconciliation ? std::optional{1.0} : covariance.variance_factor;
    using C = FitCovarianceStatus;
    const C outcome = covariance.outcome != C::Available ? covariance.outcome
        : covariance.optimizer.empty() ? C::Unavailable : !scale ? C::MissingVarianceFactor : C::Available;
    const C refused = outcome == C::Available ? C::Unavailable : outcome;
    for (auto row : rows) {
        auto& out = result.prediction_uncertainty.emplace_back(PredictionUncertainty{outcome});
        if (result.predictions.empty()) continue;
        row.observed = 0.0; row.scale = 1.0; row.form = ResidualForm::Difference;
        FitRequest point{request.parameters, request.thermochemistry, request.coordinates, {std::move(row)}, {}, request.controls, request.holds};
        FitResult failures; ObservationEvaluation evaluation;
        try {
            auto declarations = request.observations;
            declarations.push_back(point.observations.front());
            for (const auto& coordinate : point.coordinates.coordinates)
                if (coordinate.reaction) (void)admit_reaction_law(*coordinate.reaction, point.parameters.component_count(),
                    declarations, point.thermochemistry);
            // The unit-weight difference above makes the optimizer residual derivative the prediction derivative.
            ObservationEvaluator::Impl evaluator(std::move(point), failures, {1.0});
            if (evaluator.evaluate(result.optimizer, outcome == C::Available, evaluation)) {
                out = {outcome, {}, evaluation.predictions.front(), {}, ObservationStatus::Available};
                if (outcome != C::Available) continue;
                const Eigen::Map<const Eigen::VectorXd> j(evaluation.optimizer_gradient.data(), result.optimizer.size());
                const Eigen::Map<const Eigen::MatrixXd> c(covariance.optimizer.data(), j.size(), j.size());
                if (const double v = *scale * j.dot(c * j); std::isfinite(v) && v >= 0.0) out.variance = v;
                else out.outcome = C::Unavailable;
            } else {
                const auto& failure = failures.failures.back();
                out = {refused, failure.message, failure.central_value, {}, failure.status};
            }
        } catch (const std::logic_error& error) {
            out = {refused, error.what(), {}, {}, ObservationStatus::InadmissibleTrial};
        }
    }
}

FitResult fit(const FitRequest& request, const std::vector<Observation>& prediction_rows) {
    auto result = fit_parameters(request);
    if (!prediction_rows.empty()) predict_fit(request, prediction_rows, result);
    return result;
}

std::vector<FitResult> fit_sequence(FitRequest request, const std::vector<std::vector<double>>& starts,
    int iterations, double seconds, const std::vector<Observation>& prediction_rows) {
    require(iterations > 0 && std::isfinite(seconds) && seconds > 0.0, "shared fit budgets");
    const auto begin = std::chrono::steady_clock::now();
    std::vector<FitResult> results;
    for (const auto& start : starts) {
        const double remaining = seconds - std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        if (iterations <= 0 || remaining <= 0.0) break;
        auto trial = request;
        trial.controls.maximum_iterations = std::min(request.controls.maximum_iterations, iterations);
        trial.controls.maximum_elapsed_time_seconds = std::min(request.controls.maximum_elapsed_time_seconds, remaining);
        try {
            trial.parameters = evaluate_coordinates(request.parameters, request.coordinates, start, request.holds, request.thermochemistry, request.controls.solver_options).parameters;
            request.coordinates.apply(trial.observations, start, false);
            trial.start = start;
            results.push_back(fit_parameters(trial));
        } catch (const std::invalid_argument& error) {
            FitResult failed(request); failed.optimizer = start; failed.message = error.what(); results.push_back(std::move(failed));
        }
        iterations -= static_cast<int>(results.back().iterations);
    }
    if (!prediction_rows.empty()) for (auto& result : results) predict_fit(request, prediction_rows, result);
    return results;
}

}  // namespace epcsaft::regression
