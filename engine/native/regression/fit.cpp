#include "regression/fit.hpp"

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include <dlfcn.h>
#include <fcntl.h>
#include <poll.h>
#include <signal.h>
#include <spawn.h>
#include <sys/prctl.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <chrono>
#include <cstring>
#include <iterator>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>

namespace epcsaft::regression {
namespace {

void require(bool condition, const char* what) {
    if (!condition) throw std::invalid_argument(std::string("fit ") + what);
}

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); });
}

// The same point-local D maps raw observation columns and covariance.
void chain(const FitRequest& request, const std::vector<double>& map, std::size_t row,
           double sqrt_weight, double prediction, const double* gradient, double* out) {
    const auto columns = request.coordinates.coordinates.size();
    const double slope = sqrt_weight * residual_slope(request.observations[row], prediction);
    std::fill(out, out + columns, 0.0);
    for (std::size_t r = 0; r < map.size() / columns; ++r)
        for (std::size_t c = 0; c < columns; ++c) out[c] += slope * gradient[r] * map[r * columns + c];
}

// Each input row names the first row with the same state; this is the sole membership map.
std::vector<std::size_t> state_groups(const std::vector<Observation>& observations) {
    std::vector<std::size_t> groups;
    // ponytail: quadratic declaration comparison; index exact state records if fit datasets grow large.
    for (std::size_t i = 0; i < observations.size(); ++i) {
        const auto prior = std::find_if(observations.begin(), observations.begin() + static_cast<std::ptrdiff_t>(i),
            [&](const auto& row) { return same_state(observations[i], row); });
        const auto j = static_cast<std::size_t>(std::distance(observations.begin(), prior));
        groups.push_back(j == i ? i : groups[j]);
    }
    return groups;
}

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
};

bool interrupted() { return fit_host.interrupted && fit_host.interrupted(); }

// This binary, opened as it loads: workers run the loaded image even after a rebuild replaces its file.
const int loaded_image = [] {
    Dl_info image{};
    return dladdr(reinterpret_cast<void*>(&interrupted), &image) ? open(image.dli_fname, O_RDONLY | O_CLOEXEC) : -1;
}();

// Every member of a request aggregate in declaration order, named by position: the structured binding stops
// compiling when a member is added or removed, so a worker's request cannot silently omit one.
#define FIELDS(T, ...) auto fields(T& item) { auto& [__VA_ARGS__] = item; return std::tie(__VA_ARGS__); }
FIELDS(FixedParameters, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z)
FIELDS(IdealPolynomial, a, b, c, d)
FIELDS(IdealCorrelation, a, b)
FIELDS(IdealInterval, a, b, c, d, e)
FIELDS(ThermochemistryRecord, a, b)
FIELDS(CoordinateTarget, a, b, c)
FIELDS(ReactionCoordinate, a, b)
FIELDS(Coordinate, a, b, c, d, e, f)
FIELDS(equilibrium::PhaseSpec, a, b, c, d, e, f, g)
FIELDS(equilibrium::FeedSpec, a, b, c, d)
FIELDS(ReactionLogPolynomial, a, b, c, d, e, f, g, h)
FIELDS(equilibrium::ReactionSpec, a, b, c, d, e, f)
FIELDS(NeutralReferenceDeclaration, a, b, c, d)
FIELDS(equilibrium::CriticalDeclaration, a, b)
FIELDS(Observation, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u)
FIELDS(HydrationHold, a, b, c, d, e, f)
FIELDS(Failure, a, b, c, d, e)
// Observation validation refuses discovery options (a fit declares its phases), so none are sent.
auto fields(equilibrium::EquilibriumProblem& item) {
    auto& [a, b, c, d, e, f, g, h, discovery, j, k, l, m, n] = item;
    return std::tie(a, b, c, d, e, f, g, h, j, k, l, m, n);
}
// A derived record cannot be bound as a whole; its size pins the scalar and vector members listed here.
auto fields(Parameters& v) {
    static_assert(sizeof(Parameters) == sizeof(FixedParameters) + 2 * sizeof(double) + 12 * sizeof(std::vector<double>));
    return std::tuple_cat(fields(static_cast<FixedParameters&>(v)), std::tie(v.ionic_region_relative_permittivity,
        v.born_c_dielectric, v.packing_diameter_angstrom, v.debye_huckel_diameter_angstrom, v.born_diameter_angstrom,
        v.ion_suppression_coefficient, v.segment_count, v.sigma_constant, v.dispersion_energy_over_k, v.k_ij,
        v.k_ij_reciprocal_temperature_slope, v.l_ij, v.association_energy_over_k, v.association_volume));
}

// Whole values in one message, written or read in memory: trivially copyable bytes; a size, then the elements of
// a vector or string; a flag or alternative index, then the value; or every member of an aggregate above.
struct Wire final {
    std::string& message;
    bool out;
    std::size_t at = 0;
    template <class... T> void operator()(T&... values) { (one(values), ...); }
    template <class T> void one(T& value) {
        if constexpr (std::is_trivially_copyable_v<T>) bytes(&value, sizeof(T));
        else std::apply(*this, fields(value));
    }
    template <class T> void one(std::optional<T>& value) {
        bool present = value.has_value();
        one(present);
        if (present) one(out ? *value : value.emplace());
    }
    template <class... T> void one(std::variant<T...>& value) {
        std::size_t index = value.index(), k = 0;
        one(index);
        ((k++ == index ? one(out ? std::get<T>(value) : value.template emplace<T>()) : void()), ...);
    }
    template <class T> void one(std::vector<T>& values) { sequence(values); }
    void one(std::string& value) { sequence(value); }
    // Every element takes at least one byte, so a received count beyond the unread bytes is refused before
    // allocation, and the byte count below cannot overflow.
    template <class C> void sequence(C& values) {
        std::size_t size = values.size();
        one(size);
        if (!out && size > message.size() - at) throw std::runtime_error("fit worker message size");
        values.resize(size);
        if constexpr (std::is_trivially_copyable_v<typename C::value_type>) bytes(values.data(), size * sizeof(*values.data()));
        else for (auto& value : values) one(value);
    }
    void bytes(void* data, std::size_t size) {
        if (out) return void(message.append(static_cast<const char*>(data), size));
        if (size > message.size() - at) throw std::runtime_error("fit worker message size");
        std::memcpy(data, message.data() + at, size);
        at += size;
    }
};

// Bytes over a worker socket; a wait polls the fit's cancellation every 0.1 s.
void transfer(int fd, char* at, std::size_t left, bool out) {
    while (left > 0) {
        pollfd ready{fd, static_cast<short>(out ? POLLOUT : POLLIN), 0};
        if (poll(&ready, 1, 100) <= 0) {
            if (interrupted()) throw FitInterrupted{};
            continue;
        }
        const auto moved = out ? send(fd, at, left, MSG_NOSIGNAL) : recv(fd, at, left, 0);
        if (moved < 0 && errno == EINTR) continue;
        if (moved <= 0) throw std::runtime_error("fit worker process unavailable");
        at += moved;
        left -= static_cast<std::size_t>(moved);
    }
}
// One message: its byte count, then the values.
template <class... T> void put(int fd, T... values) {
    std::string message(sizeof(std::size_t), '\0');
    Wire{message, true}(values...);
    const std::size_t size = message.size() - sizeof(std::size_t);
    std::memcpy(message.data(), &size, sizeof size);
    transfer(fd, message.data(), message.size(), true);
}
template <class... T> void get(int fd, T&... values) {
    std::size_t size = 0;
    transfer(fd, reinterpret_cast<char*>(&size), sizeof size, false);
    std::string message(size, '\0');
    transfer(fd, message.data(), size, false);
    Wire{message, false}(values...);
}

// Ipopt 3.11.9 drives MUMPS 5.6, whose Fortran modules and C I/O state are process-global (#380), so concurrent
// states run in worker processes. They start by exec, since fork alone is unsafe in a threaded host, and inherit
// only their own socket. The fit kills and reaps them on every exit, an exception or cancellation included.
struct Workers final {
    std::vector<int> sockets;
    std::vector<pid_t> processes;
    ~Workers() {
        for (const pid_t process : processes) kill(process, SIGKILL);
        for (const int fd : sockets) close(fd);
        for (const pid_t process : processes)
            while (waitpid(process, nullptr, 0) < 0 && errno == EINTR) {}
    }
};

}  // namespace

// One successful point owns its immutable model and compiled central states through later promotion.
// Worker processes each evaluate a fixed share of the rows and keep that share's starts and centers;
// every process publishes only when the whole evaluation succeeds.
struct ObservationEvaluator::Impl final {
    Impl(const FitRequest& request, FitResult& result, const std::vector<std::size_t>& groups,
              const std::vector<double>& sqrt_weights)
        : request_(request), result_(result), groups_(groups), sqrt_weights_(sqrt_weights),
          starts_(groups.size()), owner_(groups.size(), 0) {}

    // Fixed partition: worker k evaluates every count-th state group in row order. Electrolyte rows share
    // one reference cache and stay here, as does every row under a trace capture or in a native host.
    void start_workers(int requested) {
        interruptible_ = true;
        const std::size_t count = core::active_trace || fit_host.worker.empty() || loaded_image < 3 ? 1
            : requested > 0 ? static_cast<std::size_t>(requested) : std::max(1u, std::thread::hardware_concurrency());
        std::size_t rank = 0;
        for (std::size_t i = 0; i < groups_.size(); ++i)
            owner_[i] = groups_[i] != i ? owner_[groups_[i]]
                : electrolyte_observation(request_.observations[i].kind) ? 0 : rank++ % count;
        std::vector<char*> command;
        for (auto& word : fit_host.worker) command.push_back(word.data());
        command.push_back(nullptr);
        for (std::size_t k = 1; k < std::min(count, rank); ++k) {
            int pair[2];
            if (socketpair(AF_UNIX, SOCK_STREAM | SOCK_CLOEXEC, 0, pair) != 0) throw std::runtime_error("fit worker socket unavailable");
            workers_.sockets.push_back(pair[0]);
            // The worker keeps its socket end as fd 0, this image as fd 3 and the standard streams. Its own
            // process group keeps terminal signals away: the fit, not its workers, answers an interrupt.
            posix_spawn_file_actions_t actions;
            posix_spawn_file_actions_init(&actions);
            posix_spawn_file_actions_adddup2(&actions, pair[1], 0);
            posix_spawn_file_actions_adddup2(&actions, loaded_image, 3);
            posix_spawn_file_actions_addclosefrom_np(&actions, 4);
            posix_spawnattr_t attributes;
            posix_spawnattr_init(&attributes);
            posix_spawnattr_setflags(&attributes, POSIX_SPAWN_SETPGROUP);
            pid_t process = 0;
            const int error = posix_spawn(&process, command[0], &actions, &attributes, command.data(), environ);
            posix_spawnattr_destroy(&attributes);
            posix_spawn_file_actions_destroy(&actions);
            close(pair[1]);
            if (error != 0) throw std::runtime_error("fit worker process unavailable");
            workers_.processes.push_back(process);
        }
        for (std::size_t k = 1; k <= workers_.sockets.size(); ++k) {
            put(workers_.sockets[k - 1], request_.parameters, request_.thermochemistry, request_.coordinates.coordinates,
                request_.observations, request_.holds, sqrt_weights_);
            put(workers_.sockets[k - 1], owner_, k);
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
        ++trials_;
        (jacobian ? result_.jacobian_evaluations : result_.residual_evaluations) += 1;
        for (const int fd : workers_.sockets) put(fd, optimizer, jacobian);
        const std::size_t n = request_.observations.size(), width = request_.coordinates.tie().size() + request_.holds.size();
        Trial fresh;
        Trial& trial = same_point ? *cached_ : fresh;
        auto failure = compute(optimizer, jacobian, out, trial, same_point, scope);
        for (std::size_t k = 1; k <= workers_.sockets.size(); ++k) {
            std::optional<Failure> theirs;
            ObservationEvaluation part;
            get(workers_.sockets[k - 1], theirs, part.predictions, part.residuals, part.statuses, part.gradient);
            // A trial refused before its rows sends none; any other reply carries every row.
            const bool rows = part.statuses.size() == n && part.predictions.size() == n && part.residuals.size() == n
                && part.gradient.size() == (jacobian ? n * width : 0);
            if (!rows && !(theirs && part.statuses.empty())) throw std::runtime_error("fit worker reply size");
            for (std::size_t i = 0; rows && out.statuses.size() == n && i < n; ++i) if (owner_[i] == k) {
                out.predictions[i] = part.predictions[i];
                out.residuals[i] = part.residuals[i];
                out.statuses[i] = part.statuses[i];
                if (jacobian) std::copy_n(part.gradient.begin() + i * width, width, out.gradient.begin() + i * width);
            }
            if (theirs && (!failure || theirs->row < failure->row)) failure = std::move(theirs);
        }
        // The serial sweep stops at its first failure: later rows stay unevaluated.
        for (std::size_t i = failure ? failure->row + 1 : n; i < out.statuses.size(); ++i) {
            out.predictions[i] = out.residuals[i] = 0.0;
            out.statuses[i] = ObservationStatus::NotEvaluated;
            if (jacobian) std::fill_n(out.gradient.begin() + i * width, width, 0.0);
        }
        const auto reject = [&](std::size_t index, ObservationStatus status, std::string message, bool flagged) {
            for (const int fd : workers_.sockets) put(fd, false);
            return fail(index, optimizer, status, std::move(message), flagged);
        };
        if (failure) return reject(failure->observation, failure->status, std::move(failure->message), failure->jacobian);
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
        const std::size_t columns = request_.coordinates.coordinates.size();
        if (jacobian) out.optimizer_gradient.resize(n * columns);
        for (std::size_t i = 0; jacobian && i < n; ++i) {
            if (sqrt_weights_[i] == 0.0) continue;
            chain(request_, trial.chart.map, i, sqrt_weights_[i], out.predictions[i],
                  out.gradient.data() + i * width, out.optimizer_gradient.data() + i * columns);
            if (!Eigen::Map<const Eigen::VectorXd>(out.optimizer_gradient.data() + i * columns,
                static_cast<Eigen::Index>(columns)).allFinite())
                return reject(i, ObservationStatus::Nonfinite, "optimizer Jacobian not finite", true);
        }
        for (const int fd : workers_.sockets) put(fd, true);
        publish(optimizer, jacobian, out, trial, fresh, same_point);
        return true;
    }

    // A worker answers each evaluation with its rows and first failure, then publishes on the fit's verdict.
    [[noreturn]] void serve(int fd) {
        for (std::vector<double> optimizer;;) {
            bool jacobian = false, commit = false;
            get(fd, optimizer, jacobian);
            core::Scope scope("fit_evaluation");
            const bool same_point = cached_ && cached_->optimizer == optimizer;
            Trial fresh;
            Trial& trial = same_point ? *cached_ : fresh;
            ObservationEvaluation out;
            const auto failure = compute(optimizer, jacobian, out, trial, same_point, scope);
            put(fd, failure, out.predictions, out.residuals, out.statuses, out.gradient);
            get(fd, commit);
            if (commit) publish(optimizer, jacobian, out, trial, fresh, same_point);
        }
    }

    bool gradient_at(const std::vector<double>& point, std::vector<double>& out) {
        ObservationEvaluation trial;
        if (!evaluate(point, true, trial)) return false;
        const std::size_t width = request_.coordinates.coordinates.size();
        const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> j(
            trial.optimizer_gradient.data(), static_cast<Eigen::Index>(trial.weighted.size()), static_cast<Eigen::Index>(width));
        const Eigen::Map<const Eigen::VectorXd> r(trial.weighted.data(), static_cast<Eigen::Index>(trial.weighted.size()));
        const Eigen::VectorXd gradient = j.transpose() * r;
        out.assign(gradient.data(), gradient.data() + gradient.size());
        return finite(out) || fail(trial.weighted.size(), point, ObservationStatus::Nonfinite, "weighted objective gradient not finite", true);
    }

    bool fail(std::size_t observation, const std::vector<double>& optimizer, ObservationStatus status,
              std::string message, bool jacobian) {
        result_.failures.push_back({trials_, observation, optimizer, status, std::move(message), jacobian});
        ++result_.trial_failures;
        return false;
    }

    // This process's rows of one trial; `i`, the sweep position, orders a refused trial among row failures.
    std::optional<Failure> compute(const std::vector<double>& optimizer, bool jacobian, ObservationEvaluation& out,
                                   Trial& trial, bool same_point, core::Scope& scope) {
        const std::size_t n = request_.observations.size(), width = request_.coordinates.tie().size() + request_.holds.size();
        std::size_t i = 0;
        try {
            auto observations = request_.observations;
            request_.coordinates.apply(observations, optimizer);
            if (!same_point) {
                trial.chart = evaluate_coordinates(request_.parameters, request_.coordinates, optimizer, request_.holds, request_.thermochemistry);
                trial.model = std::make_unique<HelmholtzModel>(trial.chart.parameters, trial.chart.active, request_.thermochemistry);
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
            out.residuals.assign(n, 0.0);
            out.statuses.assign(n, ObservationStatus::NotEvaluated);
            if (jacobian) out.gradient.assign(n * width, 0.0);
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
                    // Coordinates change only reaction constants, so a row's state declaration matches the seed's row.
                    prediction = predict(model, observation, starts_[g], std::move(prediction.compiled), &trial.references,
                        seed ? &seed->centers[g] : nullptr, seed ? seed->model.get() : nullptr,
                        seed && trial.chart.physical == seed->chart.physical);
                    if (prediction.status != ObservationStatus::Available && starts_[g])
                        prediction = predict(model, observation, std::nullopt, std::move(prediction.compiled), &trial.references);
                } else if (!pressure_closed(observation.kind)) project_prediction(model, observation, prediction);
                out.statuses[i] = prediction.status;
                if (prediction.status != ObservationStatus::Available)
                    return Failure{i, i, prediction.status, jacobian, prediction.message};
                out.predictions[i] = *prediction.value;
                out.residuals[i] = *prediction.residual;
                if (!jacobian) continue;
                if (gradients[g].empty()) {
                    core::Scope scope("fit_group_differentiation");
                    scope.probe("anchor", static_cast<double>(g));
                    std::vector<const Observation*> rows;
                    for (std::size_t row = 0; row < n; ++row)
                        if (groups_[row] == g && sqrt_weights_[row] > 0.0) rows.push_back(&trial.observations[row]);
                    gradients[g] = pressure_closed(observation.kind)
                        ? std::vector<PredictionGradient>{differentiate(model, observation, prediction)}
                        : differentiate_group(model, rows, prediction, request_.coordinates);
                }
                const auto& gradient = gradients[g][position[g]++];
                if (gradient.status != ObservationStatus::Available)
                    return Failure{i, i, gradient.status, true, gradient.message};
                std::copy(gradient.value.begin(), gradient.value.end(),
                          out.gradient.begin() + static_cast<std::ptrdiff_t>(i * width));
            }
            return std::nullopt;
        } catch (const std::exception& error) {
            return Failure{i, n, ObservationStatus::InadmissibleTrial, jacobian, error.what()};
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

    const FitRequest& request_;
    FitResult& result_;
    const std::vector<std::size_t> groups_;
    const std::vector<double> sqrt_weights_;
    std::vector<std::optional<equilibrium::EquilibriumResult>> starts_;
    std::optional<Trial> cached_;
    std::size_t trials_ = 0;
    std::vector<std::size_t> owner_;  // the process that evaluates each row; 0 is the fit's own
    std::size_t process_ = 0;
    bool interruptible_ = false;  // a least-squares fit's own process polls cancellation between rows
    Workers workers_;
};

// A worker never outlives the thread that started it; it reads its fit's request once, then serves it.
// The socket names the fit's process (its creator), so a fit that died before the death signal was armed
// shows as a different parent, including a subreaper's adoption, and the worker leaves without reading.
void serve_fit(int fd) {
    ucred fit{};
    socklen_t size = sizeof fit;
    if (getsockopt(fd, SOL_SOCKET, SO_PEERCRED, &fit, &size) != 0 || prctl(PR_SET_PDEATHSIG, SIGKILL) != 0
        || getppid() != fit.pid) _exit(1);
    FitRequest request;
    std::vector<double> sqrt_weights;
    get(fd, request.parameters, request.thermochemistry, request.coordinates.coordinates, request.observations,
        request.holds, sqrt_weights);
    const auto groups = state_groups(request.observations);
    FitResult result;
    ObservationEvaluator::Impl evaluator(request, result, groups, sqrt_weights);
    get(fd, evaluator.owner_, evaluator.process_);
    evaluator.serve(fd);
}

ObservationEvaluator::ObservationEvaluator(const FitRequest& request, FitResult& failures, const std::vector<std::size_t>& rows) {
    std::vector<double> selection(request.observations.size(), 0.0);
    for (const auto row : rows) { require(row < selection.size(), "evaluator row"); selection[row] = 1.0; }
    impl_ = std::make_shared<Impl>(request, failures, state_groups(request.observations), selection);
}
bool ObservationEvaluator::evaluate(const std::vector<double>& optimizer, bool columns, ObservationEvaluation& out) {
    out.residual_optimizer_gradient.clear();
    if (!impl_->evaluate(optimizer, columns, out)) return false;
    if (columns) out.residual_optimizer_gradient = out.optimizer_gradient;
    if (columns) for (std::size_t i = 0; i < out.predictions.size(); ++i) if (out.statuses[i] == ObservationStatus::Available) {
        const double slope = residual_slope(impl_->request_.observations[i], out.predictions[i]);
        for (std::size_t c = 0; c < optimizer.size(); ++c) out.optimizer_gradient[i * optimizer.size() + c] /= slope;
    }
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
        set_num_residuals(static_cast<int>(request.observations.size()));
        mutable_parameter_block_sizes()->push_back(static_cast<int>(request.coordinates.coordinates.size()));
    }

    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const override {
        const std::size_t width = request_.coordinates.coordinates.size();
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

void covariance(FitResult& result, const FitRequest& request, std::size_t n) {
    auto& cov = result.covariance;
    const auto& set = request.coordinates;
    const auto width = static_cast<Eigen::Index>(set.coordinates.size());
    cov.assumption = "unit-variance weighted residuals, independent errors and a regular local "
                     "linearization at an interior point; multiply by variance_factor for an "
                     "unknown common variance; no claim under bounds, rank deficiency or model discrepancy";
    for (std::size_t i = 0; i < n; ++i)
        if (request.weights.empty() || request.weights[i] > 0.0) cov.objective_rows.push_back(i);
    if (cov.objective_rows.empty()) { cov.status = "rank_deficient"; return; }
    Eigen::MatrixXd jacobian(static_cast<Eigen::Index>(cov.objective_rows.size()), width);
    for (Eigen::Index i = 0; i < jacobian.rows(); ++i)
        for (Eigen::Index c = 0; c < width; ++c)
            jacobian(i, c) = result.optimizer_jacobian[cov.objective_rows[static_cast<std::size_t>(i)] * static_cast<std::size_t>(width) + static_cast<std::size_t>(c)];
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
    // A bound is active in optimizer space, within the fit's own active distance of the values it bounds.
    const auto [lower, upper] = set.optimizer_bounds();
    for (std::size_t c = 0; c < set.coordinates.size(); ++c)
        if (std::min(result.optimizer[c] - lower[c], upper[c] - result.optimizer[c]) < policy::FIT_ACTIVE_BOUND_DISTANCE.value)
            cov.active_bounds.push_back(c);
    // The variance factor counts only rows that enter the objective (w > 0).
    if (cov.objective_rows.size() > cov.rank) cov.variance_factor = 2.0 * result.final_cost / static_cast<double>(cov.objective_rows.size() - cov.rank);
    // An unfinished iterate keeps its rank, singular values and bound diagnostics but is not a
    // parameter-covariance estimate.
    if (result.status != "converged") { cov.status = "not_converged"; cov.variance_factor.reset(); return; }
    if (cov.rank < set.coordinates.size()) { cov.status = "rank_deficient"; return; }
    if (!cov.active_bounds.empty()) { cov.status = "active_bound"; return; }
    const Eigen::MatrixXd optimizer = svd.matrixV() * singular.cwiseAbs2().cwiseInverse().asDiagonal()
        * svd.matrixV().transpose();
    const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> map(
        result.chart.map.data(), static_cast<Eigen::Index>(result.physical.size()), width);
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

std::vector<CheckResult> evaluate_checks(const FitRequest& request, const Parameters& stored, FitResult* final) {
    validate_checks(request);
    std::vector<bool> selected(request.observations.size(), false);
    const auto refreshed = [&](std::size_t i) {
        return final && (final->statuses[i] == ObservationStatus::NotEvaluated || (!request.weights.empty() && request.weights[i] == 0.0));
    };
    for (const auto& group : request.checks) for (const auto i : group.rows) selected[i] = true;
    for (std::size_t i = 0; i < selected.size(); ++i) if (refreshed(i)) selected[i] = true;
    const auto predictions = [&](const HelmholtzModel& owner) {
        std::vector<Prediction> result(request.observations.size());
        ReferenceCache references;
        for (std::size_t i = 0; i < selected.size(); ++i) if (selected[i]) {
            validate(request.observations[i], owner);
            result[i] = predict(owner, request.observations[i], {}, {}, &references);
        }
        return result;
    };
    const auto candidate = predictions(HelmholtzModel(stored, {}, request.thermochemistry));
    // Checks use stored parameters; objective rows retain their accepted warm-start evaluation.
    for (std::size_t i = 0; i < candidate.size(); ++i) if (refreshed(i)) {
        final->predictions[i] = candidate[i].value.value_or(NAN);
        final->residuals[i] = candidate[i].residual.value_or(NAN);
        final->statuses[i] = candidate[i].status;
    }
    const auto baseline = request.baseline ? predictions(HelmholtzModel(*request.baseline, {}, request.thermochemistry)) : std::vector<Prediction>{};
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
    result.chart = std::move(final.chart);
    result.physical = result.chart.physical;
    result.parameters = result.chart.parameters;
    result.predictions = std::move(final.predictions);
    result.residuals = std::move(final.residuals);
    result.statuses = std::move(final.statuses);
    result.physical_jacobian = std::move(final.gradient);
    result.weighted_residuals = std::move(final.weighted);
    result.optimizer_jacobian = std::move(final.optimizer_gradient);
    result.final_cost = final.cost;
}
void finish_checks(FitResult& result, const FitRequest& request) {
    auto check_request = request;
    request.coordinates.apply(check_request.observations, result.optimizer);
    result.checks = evaluate_checks(check_request, result.parameters, &result);
    result.acceptance = "pass";
    for (const auto& check : result.checks)
        if (!check.report_only && check.status != "pass" && result.acceptance != "unavailable") result.acceptance = check.status;
}

void validate_fit_request(const FitRequest& request) {
    const auto& set = request.coordinates;
    const auto n = request.observations.size();
    validate_checks(request);
    validate(set, request.parameters, request.observations, request.thermochemistry);
    const HelmholtzModel model(request.parameters, set.active_parameters(), request.thermochemistry);
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
    require(std::isfinite(controls.function_tolerance) && controls.function_tolerance >= 0.0
        && std::isfinite(controls.gradient_tolerance) && controls.gradient_tolerance >= 0.0
        && std::isfinite(controls.parameter_tolerance) && controls.parameter_tolerance >= 0.0, "controls");
    require(controls.maximum_iterations > 0 && std::isfinite(controls.maximum_elapsed_time_seconds)
        && controls.maximum_elapsed_time_seconds > 0.0,
        "controls: maximum_iterations and maximum_elapsed_time_seconds must be set and positive");
}

FitResult fit(const FitRequest& request) {
    if (request.minimax) return fit_minimax(request);
    FitResult result;
    const auto& set = request.coordinates;
    const std::size_t n = request.observations.size();
    std::vector<double> sqrt_weights(n, 1.0);
    const auto groups = state_groups(request.observations);
    try {
        validate_fit_request(request);
        for (std::size_t i = 0; i < request.weights.size(); ++i) sqrt_weights[i] = std::sqrt(request.weights[i]);
    } catch (const std::invalid_argument& error) {
        result.message = error.what();
        return result;
    }
    const std::size_t width = set.coordinates.size();
    ObservationEvaluator::Impl evaluator(request, result, groups, sqrt_weights);
    evaluator.start_workers(request.controls.workers);
    std::vector<double> optimizer = set.start(request.parameters, request.observations, request.thermochemistry);
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
    const auto [lower, upper] = set.optimizer_bounds();
    ceres::Solver::Options options;
    options.minimizer_type = ceres::TRUST_REGION;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.linear_solver_type = ceres::DENSE_QR;
    options.dense_linear_algebra_library_type = ceres::EIGEN;
    options.num_threads = 1;
    options.logging_type = ceres::SILENT;
    options.minimizer_progress_to_stdout = false;
    const auto& controls = request.controls;
    options.function_tolerance = controls.function_tolerance > 0.0 ? controls.function_tolerance : policy::FIT_FUNCTION_TOLERANCE.value;
    options.gradient_tolerance = controls.gradient_tolerance > 0.0 ? controls.gradient_tolerance : policy::FIT_GRADIENT_TOLERANCE.value;
    options.parameter_tolerance = controls.parameter_tolerance > 0.0 ? controls.parameter_tolerance : policy::FIT_PARAMETER_TOLERANCE.value;
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
    if (!set.within_bounds(optimizer) || !evaluator.evaluate(optimizer, true, final)) {
        result.status = "failed";
        result.message = "the final point could not be evaluated: " + result.message;
        result.statuses = std::move(final.statuses);
        return result;
    }
    publish_candidate(result, std::move(final), optimizer);
    result.usable = result.status != "failed" && finite(result.physical);
    if (result.usable) {
        covariance(result, request, n);
        finish_checks(result, request);
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
            trial.parameters = evaluate_coordinates(request.parameters, request.coordinates, start, request.holds, request.thermochemistry).parameters;
            request.coordinates.apply(trial.observations, start);
            results.push_back(fit(trial));
        } catch (const std::invalid_argument& error) {
            FitResult failed; failed.optimizer = start; failed.message = error.what(); results.push_back(std::move(failed));
        }
        iterations -= static_cast<int>(results.back().iterations);
    }
    return results;
}

}  // namespace epcsaft::regression
