// Included once, by fit.cpp inside its anonymous namespace: fit-worker processes and their serialization.
#ifdef _WIN32
#include "fit_windows.hpp"
#else
using Channel = int;
int extension_image = -1;

// Sources stay above inherited fd 3/4 so the child's descriptor copies cannot overwrite one another.
int open_image(const void* address) {
    Dl_info image{};
    if (!dladdr(address, &image)) return -1;
    const int opened = open(image.dli_fname, O_RDONLY | O_CLOEXEC);
    if (opened < 0) return -1;
    const int retained = fcntl(opened, F_DUPFD_CLOEXEC, 5);
    close(opened);
    return retained;
}

// Pin the mapped inode before an incremental link can replace the on-disk file.
const int native_image = open_image(reinterpret_cast<const void*>(&interrupted));
#endif

// Every member of a request aggregate in declaration order, named by position: the structured binding stops
// compiling when a member is added or removed, so a worker's request cannot silently omit one.
#define FIELDS(T, ...) auto fields(T& item) { auto& [__VA_ARGS__] = item; return std::tie(__VA_ARGS__); }
FIELDS(ResidualEnergyAction, action, status, diagnostics)
FIELDS(PropertyActionDiagnostics, site_rank, site_factorizations, site_solves, site_pivot_ratio, association_residual,
    pressure_density_derivative, pressure_residual, stability_determinant, pressure_factorizations, pressure_solves, thermal_evaluations,
    recordings, directional_paths, forward_sweeps, branch, returned_branch, stable_root_count, first_numerical_failure, site_conditioning)
FIELDS(core::ConditioningObservation, operation, basis, norm, method, cause, point_kind, criterion, dimension, rank,
    precision_bits, attempt, stage, iteration, threshold, condition, matrix_norm, error, criterion_value, coordinates,
    row_scale, column_scale, controls, phase, component, export_displacement, available)
FIELDS(core::NumericalRank, basis, dimension, rank, threshold, observation, first_conditioning_warning)
FIELDS(FixedParameters, ideal_gas, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z, A, B, C, D, E, F, G, H, I, J)
FIELDS(IdealPolynomial, a, b, c, d)
FIELDS(IdealCorrelation, a, b)
FIELDS(IdealInterval, a, b, c, d, e)
FIELDS(ThermochemistryRecord, a, b)
FIELDS(CoordinateTarget, a, b, c)
FIELDS(ReactionCoordinate, a, b, c)
FIELDS(Coordinate, a, b, c, d, e, f)
FIELDS(equilibrium::PhaseSpec, a, b, c, d, e, f, declared_composition_guess, g)
FIELDS(equilibrium::FeedSpec, a, b, c, d)
FIELDS(ReactionLogPolynomial, a, b, c, d, e, f, g, h, i, j)
FIELDS(equilibrium::ReactionSpec, a, b, c, d, e, f)
FIELDS(NeutralReferenceDeclaration, a, b, c, d)
FIELDS(equilibrium::CriticalDeclaration, a, b)
FIELDS(equilibrium::AffineLoading, a, b, c, d, e, f)
FIELDS(Observation, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z)
FIELDS(HydrationHold, a, b, c, d, e, f)
FIELDS(Failure, a, b, c, d, e, f, first_numerical_failure)
FIELDS(MeasuredInput, a, b, c, d, e, f)
FIELDS(MeasurementGroup, a, b, c, d)
FIELDS(MeasurementUncertainty, a, b, c, d)
FIELDS(MeasurementChart, a, b, c, d, e, f, g, h, i, j)
// Observation validation refuses discovery options (a fit declares its phases) and solids, so neither is sent.
auto fields(equilibrium::EquilibriumProblem& item) {
    auto& [a, b, c, d, e, f, g, h, discovery, j, k, l, m, loading, n, o, solids] = item;
    return std::tie(a, b, c, d, e, f, g, h, j, k, l, m, loading, n, o);
}
// Saturation replies carry physical records and estimates; compiled caches and optimizer memory stay in the worker.
FIELDS(core::NumericalEstimate, a, b, c, d, e)
FIELDS(core::EquilibriumComparison, a, b, c, d, e)
FIELDS(equilibrium::EquilibriumNumericalEstimates, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w)
FIELDS(equilibrium::SaturationOperand, a, b, c)
FIELDS(equilibrium::SaturationCandidate, a, b, c)
FIELDS(equilibrium::SaturationComparison, a, b, c, d)
FIELDS(equilibrium::SolvedStateActionDiagnostics, a, b, c, d, e, f, g, h, i, j)
FIELDS(equilibrium::SaturationNode, a, b, c, d, e, f, g, h, i, j)
FIELDS(equilibrium::SaturationSearch, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s)
auto fields(equilibrium::EquilibriumResult& v) {
    return std::tie(v.declaration, v.success, v.message, v.solver_status, v.temperature, v.pressure,
        v.amounts, v.mole_fractions, v.molar_densities, v.packing_fractions, v.residuals, v.classified_residuals,
        v.coordinates, v.phase_count, v.iterations, v.evaluations, v.hessian_calls, v.convergence,
        v.first_numerical_failure, v.last_rank, v.rank_records, v.requested_tolerance_met,
        v.topology_event, v.topology_phase, v.topology_pair, v.solver_settings, v.numerical_estimates);
}
auto fields(Prediction& v) {
    return std::tie(v.solver_options, v.value, v.unit, v.status, v.message, v.residual, v.first_numerical_failure,
        v.saturation, v.saturation_loading, v.endpoint_loading, v.added_amount, v.cap_branch,
        v.saturation_loading_gradient, v.endpoint_loading_gradient);
}
// A derived record cannot be bound as a whole; its size pins the scalar and vector members listed here.
auto fields(Parameters& v) {
    static_assert(sizeof(Parameters) == sizeof(FixedParameters) + 2 * sizeof(double) + 17 * sizeof(std::vector<double>));
    return std::tuple_cat(fields(static_cast<FixedParameters&>(v)), std::tie(v.ionic_region_relative_permittivity,
        v.born_c_dielectric, v.packing_diameter_angstrom, v.debye_huckel_diameter_angstrom, v.born_diameter_angstrom,
        v.ion_suppression_coefficient, v.segment_count, v.sigma_constant, v.dispersion_energy_over_k, v.k_ij,
        v.k_ij_reciprocal_temperature_slope, v.l_ij, v.association_energy_over_k, v.association_volume,
        v.dielectric_volume_m3_per_mol, v.dielectric_temperature_k, v.temperature_law_coefficient,
        v.cation_pair_mixing, v.cation_triple_mixing));
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
    void one(core::Convergence& v) {
        (*this)(v.residual, v.floor, v.amplification, v.contraction, v.iterations, v.evaluations,
            v.stop, v.cause, v.rank, v.backward_error);
        std::string name = v.criterion ? v.criterion->name : "";
        double target = v.criterion ? v.criterion->value : 0;
        (*this)(name, target);
        if (!out) {
            v.criterion = nullptr; v.criterion_owner.reset();
            for (const auto* entry : policy::ENTRIES) if (name == entry->name) {
                auto owned = std::make_shared<policy::Tolerance>(*entry); owned->value = target;
                v.criterion = owned.get(); v.criterion_owner = std::move(owned); break;
            }
            if (!name.empty() && !v.criterion) throw std::runtime_error("fit worker criterion name");
        }
    }
    void one(policy::SolverOptions& options) { for (auto& v : options.values) one(v); }
    template <class T, std::size_t N> void one(std::array<T, N>& values) { for (auto& v : values) one(v); }
    template <class T> void one(std::shared_ptr<T>& value) {
        bool present = static_cast<bool>(value); one(present);
        if (present) {
            using Value = std::remove_const_t<T>;
            if (!out) { auto owned = std::make_shared<Value>(); one(*owned); value = std::move(owned); }
            else if constexpr (std::is_const_v<T>) { auto copy = *value; one(copy); }
            else one(*value);
        } else if (!out) value.reset();
    }
    template <class A, class B> void one(std::pair<A, B>& value) { (*this)(value.first, value.second); }
    template <class K, class V> void one(std::map<K, V>& value) {
        std::size_t size = value.size(); one(size);
        if (out) for (const auto& entry : value) { auto key = entry.first; auto item = entry.second; (*this)(key, item); }
        else { value.clear(); for (std::size_t i = 0; i < size; ++i) { K key; V item; (*this)(key, item); value.emplace(std::move(key), std::move(item)); } }
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
#ifndef _WIN32
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
#endif
// One message: its byte count, then the values.
template <class... T> void put(Channel fd, T... values) {
    std::string message(sizeof(std::size_t), '\0');
    Wire{message, true}(values...);
    const std::size_t size = message.size() - sizeof(std::size_t);
    std::memcpy(message.data(), &size, sizeof size);
    transfer(fd, message.data(), message.size(), true);
}
template <class... T> void get(Channel fd, T&... values) {
    std::size_t size = 0;
    transfer(fd, reinterpret_cast<char*>(&size), sizeof size, false);
    std::string message(size, '\0');
    transfer(fd, message.data(), size, false);
    Wire{message, false}(values...);
}

// Ipopt 3.11.9 drives MUMPS 5.6, whose Fortran modules and C I/O state are process-global (#380), so concurrent
// states run in worker processes. They start by exec, since fork alone is unsafe in a threaded host, and inherit
// only their own socket and pinned images after pre-exec descriptor closure. The fit kills and reaps them on every exit,
// an exception or cancellation included.
#ifndef _WIN32
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
#endif
