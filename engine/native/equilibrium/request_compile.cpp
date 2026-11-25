#include "equilibrium/stoichiometry.hpp"

#include "core/physical_constants.hpp"
#include "core/trace.hpp"
#include "core/linear_error.hpp"
#include "contributions/helmholtz.hpp"
#include "properties/helmholtz_tape.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <cmath>
#include <Eigen/Dense>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace epcsaft::equilibrium {

struct CompiledSystem::Declaration final {
    EquilibriumProblem request;
    std::vector<std::string> species;
    CompiledSystem compiled;
};

namespace {

// Structured bindings make an added declaration field require an identity decision.
#define INPUT_FIELDS(Type, ...) auto input_fields(const Type& item) { const auto& [__VA_ARGS__] = item; return std::tie(__VA_ARGS__); }
INPUT_FIELDS(ReactionLogPolynomial, a, b, c, d, reference_temperature, temperature_min, temperature_max, standard_state, e, f)
INPUT_FIELDS(ReactionReference, reference_pressure, source_basis)
INPUT_FIELDS(ReactionSpec, name, stoichiometry, equilibrium_constant, affinity, correlation, reference)
INPUT_FIELDS(AmountSpec, kind, fraction)
INPUT_FIELDS(SolidSpec, name, stoichiometry, equilibrium_constant, correlation, reference, molar_volume, amount)
INPUT_FIELDS(NeutralReferenceDeclaration, solvent, molality_ratios, standard_molality, species_reference)
INPUT_FIELDS(CriticalDeclaration, density_guess, mode_guess)
INPUT_FIELDS(AffineLoading, initial, direction, lower, upper, seed, pressure)
#undef INPUT_FIELDS

auto input_fields(const PhaseSpec& p) {
    return std::make_tuple(p.kind, std::cref(p.name), std::cref(p.support), p.support_declared, p.declared_composition_guess,
        p.composition_guess.size(), bool(p.packing_guess), std::cref(p.amount));
}
auto input_fields(const SpinodalDeclaration& x) { return std::tuple_cat(input_fields(static_cast<const CriticalDeclaration&>(x)),
    std::tie(x.density_nodes)); }

struct InputEquality final {
    static bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(a)) == 0; }
    static bool same(const std::string& a, const std::string& b) { return a == b; }
    static bool same(const std::shared_ptr<const policy::Tolerance>& a, const std::shared_ptr<const policy::Tolerance>& b) { return a == b; }
    template<class T> static std::enable_if_t<std::is_enum_v<T> || std::is_integral_v<T>, bool> same(T a, T b) { return a == b; }
    template<class T> static auto same(const T& a, const T& b) -> decltype(input_fields(a), bool()) {
        return same(input_fields(a), input_fields(b));
    }
    template<class T> static bool same(const std::optional<T>& a, const std::optional<T>& b) {
        return bool(a) == bool(b) && (!a || same(*a, *b));
    }
    template<class T> static auto same(const T& a, const T& b) -> decltype(a.begin(), bool()) {
        return std::equal(a.begin(), a.end(), b.begin(), b.end(), [](const auto& x, const auto& y) { return same(x, y); });
    }
    template<class... T> static bool same(const std::tuple<T...>& a, const std::tuple<T...>& b) {
        return same_tuple(a, b, std::index_sequence_for<T...>{});
    }
    template<class T, std::size_t... I> static bool same_tuple(const T& a, const T& b, std::index_sequence<I...>) {
        return (same(std::get<I>(a), std::get<I>(b)) && ...);
    }
};

bool same_declaration(const EquilibriumProblem& a, const EquilibriumProblem& b) {
    const auto fields = [](const EquilibriumProblem& x) {
        return std::make_tuple(x.phase_calculation, x.topology_declared, x.azeotrope, x.temperature.free, x.temperature.absent,
            x.pressure.free, x.pressure.absent, x.feed.kind, x.feed.values.size(), std::cref(x.feed.solvent),
            x.feed.solvent_mass_kg, x.prescribed_liquid_composition ? x.prescribed_liquid_composition->size() + 1 : 0,
            std::cref(x.reactions), std::cref(x.standard_gibbs_j_per_mol), std::cref(x.source_reference), bool(x.enthalpy_j),
            bool(x.entropy_j_per_k), bool(x.volume_m3), bool(x.internal_energy_j), std::cref(x.spinodal), std::cref(x.critical),
            x.envelope ? std::optional(std::tuple(x.envelope->kind, x.envelope->component)) : std::nullopt,
            std::cref(x.loading), x.solver_options.nested_solve, std::cref(x.solver_options.values),
            std::cref(x.solver_options.owners), std::cref(x.solids));
    };
    return InputEquality::same(a.phases, b.phases) && InputEquality::same(fields(a), fields(b));
}

bool same_state_values(const EquilibriumProblem& a, const EquilibriumProblem& b) {
    const auto values = [](const EquilibriumProblem& x) {
        return std::make_tuple(x.temperature.value, x.pressure.value, std::cref(x.feed.values), std::cref(x.prescribed_liquid_composition),
            x.enthalpy_j, x.entropy_j_per_k, x.volume_m3, x.internal_energy_j, x.envelope ? x.envelope->target : 0.0);
    };
    return InputEquality::same(values(a), values(b)) && std::equal(a.phases.begin(), a.phases.end(), b.phases.begin(),
        [](const PhaseSpec& x, const PhaseSpec& y) {
            return InputEquality::same(std::tie(x.composition_guess, x.packing_guess), std::tie(y.composition_guess, y.packing_guess));
        });
}

}  // namespace

bool same_problem_state(const EquilibriumProblem& a, const EquilibriumProblem& b) {
    return !a.discovery && !b.discovery && same_declaration(a, b) && same_state_values(a, b);
}

std::vector<double> CompiledSystem::optimizer_identity() const {
    const auto& c = *this;
    if (c.ideal) throw std::invalid_argument("ideal optimizer state unavailable");
    const auto& p = c.model->parameters();
    std::vector<double> id;
    const auto append = [&](const auto& values) { id.push_back(values.size()); for (auto v : values) id.push_back(static_cast<double>(v)); };
    append(p.charge_number); append(p.molecular_ion); append(p.molar_mass_kg_per_mol); append(p.site_component);
    append(p.site_multiplicity); append(p.association_active); append(p.sigma_term_component);
    append(p.permittivity_term_component); append(p.permittivity_term_power);
    append(p.form_identity());
    append(std::vector<std::size_t>{p.component_count(), c.primal_size, c.rank, c.charge_row_count, c.mode_count,
                                   c.free_intensive_count, c.pinned_count});
    id.push_back(c.feasibility_cell()); id.push_back(c.log_packing_unbounded()); id.push_back(c.problem.critical.has_value()); id.push_back(c.problem.spinodal.has_value());
    id.push_back(c.problem.envelope ? static_cast<int>(c.problem.envelope->kind) : -1);
    id.push_back(c.problem.envelope ? c.problem.envelope->component : 0);
    for (const auto& phase : c.problem.phases) append(phase.support);
    for (const auto& row : c.rows) { id.push_back(static_cast<int>(row.kind)); id.push_back(row.phase); id.push_back(row.index); }
    return id;
}

void admit_trial_feed(CompiledSystem& compiled) {
    const auto& model = *compiled.model;
    const auto components = model.parameters().component_count();
    const auto& z = model.parameters().charge_number;
    // A neutral feed meets the state rule on its mole fractions; the ions of the excess sign are
    // then scaled so the solved feed carries no charge (formulation.md "Constraints and specification").
    std::vector<double> x_feed = compiled.n_feed;
    const double declared_total = sum(compiled.n_feed);  // resolve_feed admitted a positive finite total
    for (double& value : x_feed) value /= declared_total;
    try {
        model.validate_electroneutrality(x_feed);
    } catch (const std::invalid_argument& error) {
        throw std::invalid_argument(std::string("feed carries net z: ") + error.what());
    }
    double cation = 0.0, anion = 0.0;
    for (std::size_t i = 0; i < components; ++i) (z[i] > 0.0 ? cation : anion) += std::abs(z[i]) * compiled.n_feed[i];
    if (cation != anion && !exact_zero_dot(std::vector<double>(z.begin(), z.end()), compiled.n_feed))
        for (std::size_t i = 0; i < components; ++i)
            if (z[i] != 0.0 && (z[i] > 0.0) == (cation > anion))
                compiled.n_feed[i] *= std::min(cation, anion) / std::max(cation, anion);
    compiled.feed_total = sum(compiled.n_feed);
}

void normalize_trial_feed(CompiledSystem& compiled) {
    compiled.feed_composition = compiled.n_feed;
    for (double& value : compiled.feed_composition) value /= compiled.feed_total;
    compiled.balance = matvec(compiled.B, compiled.n_feed);
}

void admit_inventory_bound(const InventoryBound& bound) {
    if (bound.status == "infeasible_inventory") throw std::invalid_argument("infeasible_inventory: " + bound.termination);
}

namespace {


void validate_intensive_seed(const FreeValue& value, const char* name) {
    if (!value.free && !value.is_positive()) {
        throw std::invalid_argument(std::string("invalid ") + name);
    }
    if (value.free && !value.absent && !value.is_positive()) {
        throw std::invalid_argument(std::string("invalid ") + name + " guess");
    }
}








// Critical and envelope cells: one condition, neutral, nonreactive, T and P free.
void admit_conditions(const CompiledSystem& compiled) {
    const auto& p = compiled.problem;
    // An azeotrope cell: binary VLE bubble/dew on the full support of its feed (so a positive one), one of T and P free.
    if (p.azeotrope) {
        if (p.critical || p.spinodal || p.envelope || compiled.has_caloric_condition() || p.discovery || compiled.has_ions
            || compiled.has_reaction || p.source_reference || p.temperature.free == p.pressure.free || p.phases.size() != 2
            || compiled.pinned_count != 1 || !(p.phases[0].is_zero_pinned() || p.phases[1].is_zero_pinned())
            || p.phases[0].kind == p.phases[1].kind || compiled.n_feed.size() != 2 || p.phases[0].support.size() != 2
            || p.phases[1].support != p.phases[0].support)
            throw std::invalid_argument("azeotrope cell requires a binary feed, one free and one zero-pinned phase of "
                "opposite kinds on its full support, exactly one of T and P free and no other condition");
        return;
    }
    if (!p.critical && !p.spinodal && !p.envelope) return;
    if ((int(p.critical.has_value()) + int(p.spinodal.has_value()) + int(p.envelope.has_value()) != 1) || compiled.has_caloric_condition() || p.volume_m3 || p.temperature.free == p.spinodal.has_value() || !p.pressure.free
        || p.discovery || compiled.has_ions || compiled.has_reaction || p.source_reference)
        throw std::invalid_argument("critical and envelope cells require one neutral nonreactive condition with free T and P");
    const auto& support = p.phases[0].support;
    if (compiled.mode_count || p.critical || p.spinodal) {
        const auto& c = compiled.mode_seed();
        if (p.phases.size() != 1 || !p.phases[0].is_free() || support.size() != compiled.n_feed.size()
            || !(std::isfinite(c.density_guess_mol_per_m3) && c.density_guess_mol_per_m3 > 0.0)
            || (!c.mode_guess.empty() && c.mode_guess.size() != support.size())
            || std::any_of(compiled.feed_composition.begin(), compiled.feed_composition.end(), [](double z) { return !(z > 0.0); }))
            throw std::invalid_argument("critical cell requires one free phase on a positive full support and a density seed");
        if (p.spinodal) {
            if (support.size() < 2) throw std::invalid_argument("diffusive spinodal requires a mixture");
            const auto& nodes = p.spinodal->density_nodes;
            if (nodes.size() < 2 || !std::is_sorted(nodes.begin(), nodes.end())
                || std::adjacent_find(nodes.begin(), nodes.end()) != nodes.end()
                || !std::all_of(nodes.begin(), nodes.end(), [](double r) { return std::isfinite(r) && r > 0; }))
                throw std::invalid_argument("spinodal requires strictly increasing positive density nodes");
        }
        return;
    }
    const auto& e = *p.envelope;
    if (p.phases.size() != 2 || compiled.pinned_count != 1 || !std::isfinite(e.target) || p.phases[1].support != support
        || std::none_of(p.phases.begin(), p.phases.end(), [](const PhaseSpec& phase) { return phase.is_zero_pinned(); })
        || (e.kind == EnvelopeCoordinate::LogK && std::find(support.begin(), support.end(), e.component) == support.end()))
        throw std::invalid_argument("envelope requires a free feed phase and a zero-pinned incipient phase on one support");
}



// Reaction-aware feed-relative pin admission: individual [0, 1] domain always; the exact
// sum rules only when the declared reaction space conserves total moles.
void admit_pins(CompiledSystem& compiled) {
    double pin_sum = 0.0;
    std::size_t pinned = 0;
    for (const auto& phase : compiled.problem.phases) {
        if (phase.amount.kind != AmountKind::Pinned) continue;
        const double fraction = phase.amount.fraction;
        if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0) {
            throw std::invalid_argument("invalid pinned phase fraction in phase " + phase.name);
        }
        pin_sum += fraction;
        ++pinned;
    }
    if (!compiled.conserves_total_moles || pinned == 0) return;
    if (pin_sum > 1.0) {
        throw std::invalid_argument("pinned phase fractions exceed one: S_pin="
                                    + core::optional_number(pin_sum));
    }
    if (pinned == compiled.problem.phases.size() && pin_sum != 1.0) {
        throw std::invalid_argument("all-pinned phase fractions must sum to one: S_pin="
                                    + core::optional_number(pin_sum));
    }
}

// A present solid's row uses ln K^rho_s(T, P0) from solid_log_k, so the pressure slope of g_s is V_s alone. Fixed T and
// P exclude PH/PS, critical, envelope, azeotrope and pinned cells.
void admit_solids(CompiledSystem& compiled, const std::vector<bool>& covered) {
    const auto& p = compiled.problem;
    if (p.solids.empty()) return;
    for (const auto& solid : p.solids) {
        const auto refuse = [&](const std::string& cause) { return std::invalid_argument("solid " + solid.name + ": " + cause); };
        for (std::size_t i = 0; i < std::min(solid.stoichiometry.size(), covered.size()); ++i)
            if (solid.stoichiometry[i] != 0.0 && !covered[i]) throw refuse("species " + std::to_string(i) + " has no phase support");
        if (solid.amount.kind == AmountKind::Pinned && solid.amount.fraction != 0.0) throw refuse("amount must be Free or Pinned(0)");
    }
    ReferenceJets jets;
    jets.temperature_order = 1;
    const auto laws = solid_log_k(*compiled.model, p.temperature.value, p.pressure.value, p.source_reference, p.solids, jets);
    for (std::size_t s = 0; s < p.solids.size(); ++s) {
        const auto& law = laws[s];
        CompiledSolid item{std::nullopt, matvec(compiled.B, p.solids[s].stoichiometry),
                           law.ln_k + law.scale + *law.reference.terminal_contraction, std::nullopt, false};
        // A scalar K carries no temperature law: its T actions are refused.
        if (law.slope && law.reference.terminal_contraction_temperature_derivative)
            item.log_k_slope = *law.slope + *law.reference.terminal_contraction_temperature_derivative;
        item.reference = law.reference;
        compiled.solids.push_back(std::move(item));
    }
}

}  // namespace


static void append_phase_rows(CompiledSystem& compiled) {
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        for (const std::size_t index : compiled.problem.phases[phase].support)
            compiled.rows.push_back({RowKind::Chemical, phase, index});
        compiled.rows.push_back({RowKind::Pressure, phase, 0});
    }
}

static void inspect_generic_rank(CompiledSystem& compiled) {
    const auto point = compiled.generic_point();
    const bool state_only = !compiled.problem.loading
        && std::none_of(compiled.rows.begin(), compiled.rows.end(), [](const auto& row) { return condition_row(row.kind); });
    const auto prepared = compiled.prepare_values(point, state_only, nullptr, std::nullopt, {}, false, false,
        state_only && !compiled.problem.temperature.free);
    if (!compiled.problem.temperature.free && !compiled.problem.pressure.free) compiled.fixed_standard = prepared.standard;
    const auto jacobian = compiled.jacobian_at(prepared);
    Eigen::MatrixXd matrix = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        jacobian.data(), static_cast<Eigen::Index>(compiled.rows.size()), static_cast<Eigen::Index>(point.size()));
    for (Eigen::Index column = 0; column < matrix.cols(); ++column) {
        const double norm = matrix.col(column).cwiseAbs().maxCoeff();
        if (norm > 0.0) matrix.col(column) /= norm;
    }
    compiled.last_rank = core::numerical_rank(matrix.fullPivLu(), "column_scaled");
    compiled.last_rank->observation.point_kind = "generic_point";
    compiled.jacobian_rank = compiled.last_rank->rank;
    compiled.generic_nullity = point.size() - compiled.jacobian_rank;  // the rank is at most the column count
    compiled.feasible = compiled.generic_nullity == 0;
    // The compiled record is the sole source for both the driver and callers.
}

static void reuse_layout(CompiledSystem& compiled, const CompiledSystem& prepared) {
    compiled.rows = prepared.rows;
    compiled.shape = prepared.shape;
    compiled.charge_row = prepared.charge_row;
    compiled.charge_row_count = prepared.charge_row_count;
    compiled.log_charge_rows = prepared.log_charge_rows;
    compiled.primal_size = prepared.primal_size;
}

static CompiledSystem compile_prescribed_liquid(const HelmholtzModel& model, const EquilibriumProblem& request,
                                                const CompiledSystem* prepared = nullptr) {
    if (model.parameters().ideal_gas || !request.topology_declared || request.phases.size() != 2
        || request.temperature.free || !request.pressure.free || !request.feed.values.empty()
        || !request.feed.solvent.empty() || request.discovery || !request.reactions.empty() || !request.solids.empty()
        || request.source_reference || request.enthalpy_j || request.entropy_j_per_k || request.internal_energy_j
        || request.volume_m3 || request.spinodal || request.critical || request.envelope || request.loading || request.azeotrope)
        throw std::invalid_argument("prescribed composition requires two nonreacting liquids at fixed T and free P without a feed");
    validate_intensive_seed(request.temperature, "temperature");
    validate_intensive_seed(request.pressure, "pressure");
    const auto& x = *request.prescribed_liquid_composition;
    model.validate_electroneutrality(x);
    const auto& z = model.parameters().charge_number;
    std::vector<std::size_t> support;
    for (std::size_t i = 0; i < x.size(); ++i) if (x[i] > 0) {
        support.push_back(i);
    }
    const auto sign_present = [&](bool positive) {
        return std::any_of(support.begin(), support.end(), [&](auto i) { return positive ? z[i] > 0 : z[i] < 0; });
    };
    if (!sign_present(true) || !sign_present(false))
        throw std::invalid_argument("prescribed liquid support must include ions of both signs");
    if (prepared && support != prepared->declared_supports.front()) return compile_problem(model, request);
    CompiledSystem compiled;
    compiled.model = &model;
    compiled.problem = request;
    compiled.phase_calculation = request.phase_calculation;
    compiled.problem.solver_options = request.solver_options.resolved(model.solver_options());
    compiled.declared_supports.resize(2, support);
    for (auto& phase : compiled.problem.phases) {
        if (phase.kind != PhaseKind::Liquid || !phase.is_free() || phase.packing_guess
            || ((phase.support_declared || !phase.support.empty()) && phase.support != support))
            throw std::invalid_argument("prescribed liquid requires two liquid phases with one common model-order support and no amount or packing pins");
        phase.support = support;
        phase.support_declared = true;
        phase.amount = {AmountKind::Pinned, 1.0};
    }
    auto& a = compiled.problem.phases[0];
    auto& b = compiled.problem.phases[1];
    a.composition_guess.clear();
    for (auto i : support) a.composition_guess.push_back(x[i]);
    if (b.composition_guess.size() == x.size() && support.size() != x.size()) {
        const auto full_guess = std::move(b.composition_guess);
        for (std::size_t i = 0; i < x.size(); ++i) {
            if (x[i] > 0) b.composition_guess.push_back(full_guess[i]);
            else if (full_guess[i] != 0) throw std::invalid_argument("conjugate guess has nonzero fraction outside the common support");
        }
    }
    if (b.composition_guess.size() != support.size()
        || std::any_of(b.composition_guess.begin(), b.composition_guess.end(), [](double y) { return !(y > 0); }))
        throw std::invalid_argument("prescribed liquid requires a positive conjugate composition guess on its support");
    std::vector<double> guess(x.size());
    for (std::size_t j = 0; j < support.size(); ++j) guess[support[j]] = b.composition_guess[j];
    model.validate_electroneutrality(guess);
    // The feasibility driver uses a unit molar evaluation basis, with no inventory.
    compiled.feed_total = 1.0;
    compiled.has_ions = true;
    compiled.pinned_count = 2;
    compiled.free_intensive_count = 1;
    compiled.rank = support.size();
    if (prepared) {
        compiled.B = prepared->B;
        reuse_layout(compiled, *prepared);
    } else {
        compiled.B.assign(compiled.rank, std::vector<double>(x.size()));
        append_phase_rows(compiled);
        compiled.rows.push_back({RowKind::Electro, 1, 0});
        compiled.rows.push_back({RowKind::Normalization, 1, 0});
        compiled.charge_row = {std::nullopt, 0};
        compiled.charge_row_count = 1;
        std::size_t column = 0;
        for (std::size_t phase = 0; phase < 2; ++phase) {
            auto& layout = compiled.shape.emplace_back();
            for (std::size_t i = 0; i < support.size(); ++i) {
                layout.x.push_back(column++);
                if (phase == 0) compiled.rows.push_back({RowKind::PrescribedComposition, 0, i});
            }
            layout.volume = column++;
        }
        for (std::size_t i = 0; i < support.size(); ++i) compiled.B[i][support[i]] = 1;
        compiled.primal_size = column + 1;
    }
    inspect_generic_rank(compiled);
    return compiled;
}

static CompiledSystem compile_with_basis(
    const HelmholtzModel& model, const EquilibriumProblem& request,
    const std::vector<std::vector<double>>* trial_basis, const CompiledSystem* prepared = nullptr
) {
    if (model.parameters().ideal_gas && (request.phases.size() != 1 || request.phases[0].kind != PhaseKind::Vapor
        || request.temperature.free || request.pressure.free || !request.solids.empty() || request.discovery
        || request.enthalpy_j || request.entropy_j_per_k || request.internal_energy_j || request.source_reference || request.critical || request.spinodal || request.envelope
        || request.azeotrope || request.has_reaction_reference() || request.phases[0].packing_guess))
        throw std::invalid_argument("ideal-gas equilibrium requires one vapor at fixed T/P; calorics and phase searches deferred");
    if (std::any_of(request.reactions.begin(), request.reactions.end(), [](const auto& r) { return r.uses_species_thermochemistry(); })
        && (request.temperature.free || request.pressure.free || request.enthalpy_j || request.entropy_j_per_k || request.internal_energy_j || request.volume_m3
            || !request.solids.empty() || request.source_reference || request.has_reaction_reference()))
        throw std::invalid_argument("species-derived targets require ordinary fixed T/P; calorics and source-reference transfer deferred");
    if (!request.topology_declared) {
        throw std::invalid_argument("phase count unknown; declare phases");
    }
    if (request.phases.empty()) throw std::invalid_argument("no phases declared");
    validate_intensive_seed(request.temperature, "temperature");
    validate_intensive_seed(request.pressure, "pressure");
    const std::size_t components = model.parameters().component_count();
    CompiledSystem compiled;
    compiled.model = &model;
    compiled.phase_calculation = request.phase_calculation;
    compiled.problem = request;
    compiled.problem.solver_options = request.solver_options.resolved(model.solver_options());
    compiled.n_feed = resolve_feed(model, request.feed);
    compiled.has_reaction = !request.reactions.empty();
    if (request.volume_m3 && request.pressure.absent && *request.volume_m3 > 0.0)
        compiled.problem.pressure = FreeValue::variable(sum(compiled.n_feed) * GAS_CONSTANT_J_PER_MOL_K * compiled.seed_temperature() / *request.volume_m3);
    const auto& z = model.parameters().charge_number;
    admit_trial_feed(compiled);
    for (std::size_t row = 0; row < request.reactions.size(); ++row) {
        const auto& nu = request.reactions[row].stoichiometry;
        if (nu.size() != components) throw std::invalid_argument("stoichiometry size");
        // The charge defect is judged against the stoichiometric scale ||nu||_1.
        const double total = std::inner_product(nu.begin(), nu.end(), z.begin(), 0.0);
        const double scale = std::accumulate(nu.begin(), nu.end(), 0.0, [](double sum, double v) { return sum + std::abs(v); });
        if (std::abs(total) > std::sqrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, scale)) {
            throw std::invalid_argument("reaction " + reaction_name(request.reactions[row], row)
                                        + " violates z");
        }
    }
    // Supports are admitted before any reference calculation; an omitted support keeps its
    // full declared width for supplied guesses and starts.
    std::vector<std::size_t> all(components);
    std::iota(all.begin(), all.end(), 0);
    const bool reuse_support = prepared && prepared->removed_species.empty();
    for (auto& phase : compiled.problem.phases) {
        const bool omitted = phase.support.empty() && !phase.support_declared;
        phase.support = reuse_support ? prepared->problem.phases[compiled.declared_supports.size()].support
                                     : resolve_support(phase, all, z, &model.parameters());
        compiled.declared_supports.push_back(omitted ? all : phase.support);
        if (omitted && phase.composition_guess.size() == components) restrict_composition_guess(phase, all);
    }
    compiled.B = trial_basis ? *trial_basis : reuse_support && request.reactions.empty() ? prepared->B
        : balance_basis(components, request.reactions, &compiled);
    // An azeotrope cell's B' = T B makes the last row the total amount, the only balance row it keeps.
    if (request.azeotrope && request.reactions.empty()) compiled.B.back().assign(components, 1.0);
    compiled.rank = compiled.B.size();
    compiled.conserves_total_moles = conserves_total_moles(compiled.B);
    const double temperature_guess = compiled.seed_temperature();
    const double pressure_guess = compiled.seed_pressure();
    if (request.temperature.free)
        if (const auto row = thermal_reference_missing(temperature_guess, request.reactions, request.standard_gibbs_j_per_mol, model.thermochemistry() ? &*model.thermochemistry() : nullptr))
            throw std::invalid_argument("temperature-dependent standard state required for reaction "
                                        + reaction_name(request.reactions[*row], *row));
    normalize_trial_feed(compiled);
    std::vector<bool> covered(components, false);
    const auto participant = reaction_participants(components, request.reactions);
    auto inventory = compiled;
    if (!trial_basis) {
        auto& support = inventory.problem.phases.front().support;
        for (std::size_t i = 0; i < components; ++i)
            if (participant[i] && std::find(support.begin(), support.end(), i) == support.end()) support.push_back(i);
    }
    std::vector<std::size_t> indices;
    for (const auto& phase : inventory.problem.phases)
        indices.insert(indices.end(), phase.support.begin(), phase.support.end());
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    // With B = I the rows fix y = n_feed/N_feed, so u_i/N_feed is the feed fraction without a program. A fed
    // species whose box 1e-2 u_i/N_feed would fall below 1e-8 keeps HiGHS, so no box changes (#646).
    bool feed_point = !trial_basis && request.reactions.empty() && !request.azeotrope;
    for (const std::size_t i : all)
        feed_point = feed_point && (compiled.n_feed[i] == 0.0 || (std::count(indices.begin(), indices.end(), i)
            && policy::EQUILIBRIUM_BOX_INVENTORY_FRACTION.value * compiled.feed_composition[i] >= policy::EQUILIBRIUM_COMPOSITION_BOX.value));
    // A trial with one positive balance row attains each species maximum by putting the
    // entire inventory in that species. This is the same LP as inventory_bounds with one equality.
    const bool single_balance = trial_basis && compiled.rank == 1 && compiled.balance.front() > 0.0
        && std::all_of(indices.begin(), indices.end(), [&](auto i) { return compiled.B.front()[i] > 0.0; });
    const auto bounds = feed_point || single_balance ? std::vector<InventoryBound>{} : inventory_bounds(inventory, indices, false, 0.0);
    compiled.inventory_upper_bounds = feed_point ? compiled.feed_composition : std::vector<double>(components, 1.0);
    if (single_balance) for (const auto i : indices)
        compiled.inventory_upper_bounds[i] = (compiled.balance.front() / compiled.feed_total) / compiled.B.front()[i];
    std::vector<bool> absent(components, false);
    const auto certify = [&](std::size_t i, std::vector<double> candidate) {
        candidate = reconstruct_conserved_functional(std::move(candidate), request.reactions);
        return !candidate.empty() && candidate[i] > 0.0
            && std::all_of(candidate.begin(), candidate.end(), [](double value) { return value >= 0.0; })
            && exact_zero_dot(candidate, compiled.n_feed);
    };
    for (const auto& bound : bounds) {
        admit_inventory_bound(bound);
        if (bound.status == "optimal") compiled.inventory_upper_bounds[bound.species] = *bound.feed_normalized_amount;
        if (!trial_basis && bound.status == "optimal" && *bound.maximum == 0.0 && !bound.conserved_functional.empty())
            absent[bound.species] = certify(bound.species, bound.conserved_functional);
    }
    if (!trial_basis) for (std::size_t i = 0; i < components; ++i) {
        for (const auto& row : compiled.B) absent[i] = absent[i] || certify(i, row);
        if (absent[i]) compiled.removed_species.push_back(i);
    }
    for (auto& phase : compiled.problem.phases) {
        phase.support_declared = true;
        const auto support_before = phase.support;
        phase.support.erase(std::remove_if(
            phase.support.begin(), phase.support.end(), [&](std::size_t index) {
                return absent[index] || (compiled.n_feed[index] == 0.0 && !participant[index]);
            }), phase.support.end());
        if (phase.is_zero_pinned()) {
            const auto removed = std::find_if(support_before.begin(), support_before.end(),
                [&](std::size_t index) { return absent[index]; });
            if (removed != support_before.end()) {
                std::string reason;
                if (phase.support.empty()) reason = "empty effective support";
                else if (std::any_of(phase.support.begin(), phase.support.end(),
                    [&](std::size_t index) { return z[index] != 0; }))
                    reason = "effective incipient support must contain only uncharged species";
                else if (request.phases.size() != 2 || request.phases[0].kind == request.phases[1].kind)
                    reason = "requires exactly two phases of opposite liquid/vapor kinds";
                else if (request.temperature.free || !request.temperature.is_positive() || !request.pressure.free)
                    reason = "requires fixed positive T and free P";
                else if (std::count_if(request.phases.begin(), request.phases.end(),
                    [](const PhaseSpec& item) { return item.is_free(); }) != 1)
                    reason = "requires one zero pin and one free amount";
                if (!reason.empty()) throw std::invalid_argument("zero_inventory_incipient_support: component "
                    + std::to_string(*removed) + " in phase " + phase.name + ": " + reason);
            }
        }
        restrict_composition_guess(phase, support_before);
        for (const std::size_t index : support_before) {
            if (std::find(phase.support.begin(), phase.support.end(), index) == phase.support.end()
                && std::find(compiled.removed_species.begin(), compiled.removed_species.end(), index)
                    == compiled.removed_species.end()) {
                compiled.removed_species.push_back(index);
            }
        }
        if (phase.support.empty()) throw std::invalid_argument("empty support in phase " + phase.name);
        for (const std::size_t index : phase.support) covered[index] = true;
        const bool positive = std::any_of(phase.support.begin(), phase.support.end(), [&](std::size_t i) { return z[i] > 0; });
        const bool negative = std::any_of(phase.support.begin(), phase.support.end(), [&](std::size_t i) { return z[i] < 0; });
        compiled.has_ions = compiled.has_ions || positive || negative;
        if (positive != negative) throw std::invalid_argument("support z sign in phase " + phase.name);
    }
    for (std::size_t i = 0; i < components; ++i) {
        if (!covered[i] && !absent[i] && (compiled.n_feed[i] != 0.0 || participant[i])) {
            throw std::invalid_argument("species " + std::to_string(i) + " has no support");
        }
    }
    if (!trial_basis && !compiled.removed_species.empty()) restrict_inventory_balances(compiled, covered);
    const auto same_support = [](const auto& a, const auto& b) { return a.support == b.support; };
    if (prepared && (compiled.rank != prepared->rank || !std::equal(compiled.problem.phases.begin(), compiled.problem.phases.end(), prepared->problem.phases.begin(), same_support)))
        return compile_problem(model, request);
    admit_pins(compiled);
    if (!request.solids.empty() && (request.temperature.free || request.pressure.free))
        throw std::invalid_argument("solids require fixed T and P");
    const auto& phases = compiled.problem.phases;
    const auto pinned = [](const PhaseSpec& phase) { return phase.is_pinned(); };
    compiled.pinned_count = static_cast<std::size_t>(std::count_if(phases.begin(), phases.end(), pinned));
    compiled.free_intensive_count = request.temperature.free + request.pressure.free;
    const std::size_t caloric_count = compiled.has_caloric_condition();
    if (request.volume_m3 && (!std::isfinite(*request.volume_m3) || *request.volume_m3 <= 0.0
        || request.temperature.free != request.internal_energy_j.has_value() || !request.pressure.free || compiled.pinned_count
        || !request.topology_declared || request.discovery || compiled.has_reaction || (caloric_count && !request.internal_energy_j)
        || request.source_reference || !request.solids.empty() || request.critical || request.envelope || request.azeotrope
        || std::any_of(z.begin(), z.end(), [](int v) { return v != 0; })))
        throw std::invalid_argument("TV requires positive total volume, fixed T, free P and declared neutral nonreacting free phases");
    if (request.internal_energy_j && (!request.volume_m3 || phases.size() > 2))
        throw std::invalid_argument("UV requires volume and one or two declared phases");
    admit_conditions(compiled);
    // Scalar conditions e: a caloric row, two net critical equations or an envelope row.
    const auto condition_count = compiled.has_caloric_condition() + request.volume_m3.has_value()
        + 2 * request.critical.has_value() + request.spinodal.has_value() + request.envelope.has_value();
    if (compiled.free_intensive_count != compiled.pinned_count + condition_count) {
        if (caloric_count) throw std::invalid_argument("PH/PS requires declared free amounts, free T and fixed P");
        if (compiled.free_intensive_count < compiled.pinned_count)  // so at least one phase is pinned
            throw std::invalid_argument("pinned phase " + std::find_if(phases.begin(), phases.end(), pinned)->name
                                        + " requires a free T or P");
        throw std::invalid_argument(
            "free temperature or pressure requires a pinned phase amount"
        );
    }
    if (compiled.problem.phases.size() - compiled.pinned_count > compiled.rank + caloric_count + request.volume_m3.has_value()) {
        throw std::invalid_argument(
            "phase count exceeds balance rank for phase "
            + compiled.problem.phases.back().name
            + "; pin a phase or free T or P"
        );
    }
    if (!prepared) append_phase_rows(compiled);
    // Every ion-bearing phase has a charge row except the first with a free or positive amount.
    bool rowless = false;
    for (std::size_t phase = 0; !prepared && phase < phases.size(); ++phase) {
        const auto& spec = phases[phase];
        if (std::none_of(spec.support.begin(), spec.support.end(), [&](std::size_t i) { return z[i] != 0; })) continue;
        if (!rowless && (spec.is_free() || spec.amount.fraction > 0.0)) {
            if (!trial_basis && spec.kind == PhaseKind::Vapor)
                throw std::invalid_argument("charged vapor cannot be the omitted-charge phase");
            rowless = true;
        }
        else {
            if (spec.kind == PhaseKind::Vapor) compiled.log_charge_rows.emplace_back(compiled.rows.size(), phase);
            compiled.rows.push_back({RowKind::Electro, phase, 0});
            ++compiled.charge_row_count;
        }
    }
    const auto central_reactions = central_reference_reactions(
        model, request, temperature_guess, pressure_guess, compiled.reference_diagnostics, {});
    // central_reference_reactions refuses supplied standard Gibbs with a source reference.
    compiled.standard_gibbs = standard_state_series(temperature_guess, components, central_reactions, request.standard_gibbs_j_per_mol, 0, model.thermochemistry() ? &*model.thermochemistry() : nullptr, nullptr, &compiled).front();
    admit_solids(compiled, covered);
    if (prepared) {
        reuse_layout(compiled, *prepared);
        compiled.mode_count = prepared->mode_count;
        for (std::size_t s = 0; s < compiled.solids.size(); ++s) compiled.solids[s].column = prepared->solids[s].column;
    } else {
        for (std::size_t row = request.azeotrope ? compiled.B.size() - 1 : 0; row < compiled.B.size(); ++row)
            compiled.rows.push_back({RowKind::Moiety, 0, row});
        for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase)
            compiled.rows.push_back({RowKind::Normalization, phase, 0});
        for (std::size_t s = 0; s < request.solids.size(); ++s)
            if (request.solids[s].amount.kind == AmountKind::Free) compiled.rows.push_back({RowKind::SolidAffinity, 0, s});
        if (request.volume_m3) compiled.rows.push_back({RowKind::Volume, 0, 0});
        if (caloric_count) compiled.rows.push_back({request.internal_energy_j ? RowKind::InternalEnergy : request.enthalpy_j ? RowKind::Enthalpy : RowKind::Entropy, 0, 0});
        // Critical rows M u (C), u'u - 1 and c; or one envelope row. Their positions are the last rows.
        compiled.mode_count = (request.critical || request.spinodal) ? components : 0;
        for (std::size_t i = 0; i < compiled.mode_count + 2 * request.critical.has_value() + request.spinodal.has_value() + request.envelope.has_value(); ++i) {
            const auto kind = request.envelope ? RowKind::EnvelopeCondition : i < components ? RowKind::CriticalMode
                : i == components ? RowKind::CriticalNormalization : RowKind::CriticalCubic;
            compiled.rows.push_back({request.spinodal ? (i < components ? RowKind::SpinodalMode : RowKind::SpinodalNormalization) : kind, 0, i});
        }
        for (std::size_t i = 0; request.azeotrope && i + 1 < components; ++i) compiled.rows.push_back({RowKind::AzeotropeCondition, 0, i});
        if (trial_basis && phases[0].kind == PhaseKind::Vapor && compiled.has_ions)
            for (std::size_t row = 0; row < compiled.rows.size(); ++row)
                if (compiled.rows[row].kind == RowKind::Moiety && compiled.rows[row].index == 1)
                    compiled.log_charge_rows.emplace_back(row, 0);
        std::size_t offset = 0, charge = 0;
        for (const auto& phase : compiled.problem.phases) {
            auto& item = compiled.shape.emplace_back();
            for (std::size_t i = 0; i < phase.support.size(); ++i) item.x.push_back(offset++);
            if (phase.amount.kind == AmountKind::Free) item.amount = offset++;
            item.volume = offset++;
        }
        compiled.primal_size = offset;
        compiled.charge_row.resize(compiled.problem.phases.size());
        for (const auto& row : compiled.rows)
            if (row.kind == RowKind::Electro) compiled.charge_row[row.phase] = charge++;
        for (std::size_t s = 0; s < compiled.solids.size(); ++s)  // free solid amounts follow the phase coordinates
            if (request.solids[s].amount.kind == AmountKind::Free) compiled.solids[s].column = compiled.primal_size++;
        compiled.primal_size += compiled.mode_count + compiled.free_intensive_count;
    }
    if (request.loading) {
        const auto& a = *request.loading;
        if (request.temperature.free || !request.pressure.free || !request.solids.empty() || request.critical
            || request.spinodal || request.envelope || request.azeotrope || compiled.has_caloric_condition()
            || request.feed.kind != FeedKind::Amounts || a.initial.size() != components || a.direction.size() != components)
            throw std::invalid_argument("affine saturation loading needs fixed temperature, free pressure, an amounts feed, "
                "no solid or special-point condition, and initial and direction with one entry per component");
        if (compiled.pinned_count != 1 || request.phases.size() != 2
            || std::count_if(request.phases.begin(), request.phases.end(), [](const auto& p) { return p.is_zero_pinned() && p.kind == PhaseKind::Vapor; }) != 1
            || std::count_if(request.phases.begin(), request.phases.end(), [](const auto& p) { return p.is_free() && p.kind == PhaseKind::Liquid; }) != 1)
            throw std::invalid_argument("affine saturation loading needs one zero-pinned vapor and one free liquid phase");
        if (!(std::isfinite(a.lower) && std::isfinite(a.upper) && a.lower > 0 && a.upper > a.lower
              && std::isfinite(a.pressure) && a.pressure > 0))
            throw std::invalid_argument("affine saturation loading needs finite bounds 0 < lower < upper and a positive pressure");
        if (!(std::isfinite(a.seed) && a.seed >= a.lower && a.seed <= a.upper))
            throw std::invalid_argument("affine saturation loading seed lies outside [lower, upper]");
        for (std::size_t i = 0; i < components; ++i)
            if (!std::isfinite(a.initial[i]) || a.initial[i] < 0 || !std::isfinite(a.direction[i])
                || a.direction[i] < 0 || compiled.n_feed[i] != a.initial[i] + a.seed * a.direction[i])
                throw std::invalid_argument("affine saturation feed seed mismatch");
        if (!exact_zero_dot(std::vector<double>(z.begin(), z.end()), a.direction) || sum(a.direction) <= 0)
            throw std::invalid_argument("affine saturation requires neutral positive addition");
        if (!prepared) {
            ++compiled.primal_size;
            compiled.rows.push_back({RowKind::SaturationPressure, 0, 0});
        }
    }
    inspect_generic_rank(compiled);
    return compiled;
}

void restrict_composition_guess(PhaseSpec& phase, const std::vector<std::size_t>& support) {
    const auto guess = std::move(phase.composition_guess);
    phase.composition_guess.clear();
    if (guess.size() != support.size()) return;
    for (std::size_t i = 0; i < support.size(); ++i)
        if (std::find(phase.support.begin(), phase.support.end(), support[i]) != phase.support.end())
            phase.composition_guess.push_back(guess[i]);
}

IdealMoleFraction::IdealMoleFraction(std::vector<std::string> ids, std::vector<int> charges,
                                   policy::SolverOptions options, std::optional<std::string> basis,
                                   std::optional<double> volume, std::optional<std::string> solvent)
    : component_ids(std::move(ids)), charge_number(std::move(charges)), solver_options(std::move(options)),
      activity_basis(basis.value_or("mole_fraction")), solution_volume_m3(volume), unit_activity_solvent(std::move(solvent)) {
    if (component_ids.empty() || component_ids.size() != charge_number.size())
        throw std::invalid_argument("ideal species/charge size");
    auto sorted = component_ids;
    std::sort(sorted.begin(), sorted.end());
    if (sorted.front().empty() || std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
        throw std::invalid_argument("ideal species identifiers must be nonempty and unique");
    if (activity_basis != "mole_fraction" && activity_basis != "molarity")
        throw std::invalid_argument("unknown ideal activity basis");
    if (activity_basis == "mole_fraction" && (volume || unit_activity_solvent))
        throw std::invalid_argument("mole-fraction basis refuses volume and unit-activity solvent");
    if (activity_basis == "molarity" && (!volume || !std::isfinite(*volume) || *volume <= 0.0))
        throw std::invalid_argument("molarity requires finite positive solution volume in m^3");
    if (unit_activity_solvent) {
        const auto found = std::find(component_ids.begin(), component_ids.end(), *unit_activity_solvent);
        if (found == component_ids.end() || charge_number[found - component_ids.begin()] != 0)
            throw std::invalid_argument("unit-activity solvent must be a known neutral species");
        solvent_index = found - component_ids.begin();
    }
}

CompiledSystem compile_problem(const IdealMoleFraction& model, const EquilibriumProblem& request) {
    if (request.prescribed_liquid_composition) throw std::invalid_argument("prescribed liquid requires a Helmholtz model");
    core::Scope scope("compile_problem");
    if (!request.topology_declared || request.discovery || request.phases.size() != 1
        || request.phases[0].kind != PhaseKind::Liquid || !request.phases[0].is_free()
        || request.temperature.free || request.pressure.free || request.enthalpy_j || request.entropy_j_per_k || request.internal_energy_j || request.volume_m3
        || request.critical || request.spinodal || request.envelope || request.azeotrope || !request.solids.empty()
        || request.phases[0].packing_guess)
        throw std::invalid_argument("ideal mole-fraction model requires one free liquid at fixed T/P without EOS conditions");
    validate_intensive_seed(request.temperature, "temperature");
    validate_intensive_seed(request.pressure, "pressure");
    if (request.source_reference || !request.standard_gibbs_j_per_mol.empty())
        throw std::invalid_argument("ideal mole-fraction model refuses species standard Gibbs and neutral references");
    CompiledSystem c;
    c.ideal = &model; c.problem = request;
    c.problem.solver_options = request.solver_options.resolved(model.solver_options);
    c.n_feed = resolve_feed(model, request.feed);
    c.feed_total = sum(c.n_feed);
    c.feed_composition = c.n_feed;
    for (double& value : c.feed_composition) value /= c.feed_total;
    const std::vector<double> z(model.charge_number.begin(), model.charge_number.end());
    if (!exact_zero_dot(z, c.n_feed)) throw std::invalid_argument("ideal feed carries net z");
    c.B = balance_basis(c.n_feed.size(), request.reactions, &c);
    c.rank = c.B.size(); c.balance = matvec(c.B, c.n_feed);
    c.has_reaction = !request.reactions.empty();
    c.has_ions = std::any_of(z.begin(), z.end(), [](double value) { return value != 0.0; });
    c.conserves_total_moles = conserves_total_moles(c.B);
    for (const auto& reaction : request.reactions) {
        if (reaction.reference || reaction.affinity_j_per_mol)
            throw std::invalid_argument("ideal mole-fraction model refuses reaction references and affinities");
        if (!exact_zero_dot(z, reaction.stoichiometry)) throw std::invalid_argument("reaction violates z: " + reaction.name);
    }
    const auto targets = reaction_target_series(request.temperature.value, request.reactions, 0, nullptr,
                                                model.standard_state_id());
    for (double target : targets[0]) c.ideal_log_k.push_back(-target / (GAS_CONSTANT_J_PER_MOL_K * request.temperature.value));
    std::vector<std::size_t> all(c.n_feed.size()); std::iota(all.begin(), all.end(), 0);
    auto& phase = c.problem.phases[0];
    phase.support = resolve_support(phase, all, model.charge_number);
    c.declared_supports = {phase.support};
    const auto participating = reaction_participants(c.n_feed.size(), request.reactions);
    for (std::size_t i : all)
        if ((participating[i] || c.n_feed[i] > 0.0) && !std::count(phase.support.begin(), phase.support.end(), i))
            throw std::invalid_argument("ideal species has no support: " + model.component_ids[i]);
    c.inventory_upper_bounds.assign(all.size(), 0.0);
    std::vector<bool> retained(all.size(), false);
    for (const auto& bound : inventory_bounds(c, phase.support, false, 0.0)) {
        if (bound.status != "optimal" || !bound.feed_normalized_amount)
            throw std::invalid_argument("ideal inventory " + bound.status + ": " + bound.termination);
        c.inventory_upper_bounds[bound.species] = *bound.feed_normalized_amount;
        retained[bound.species] = *bound.maximum > 0.0;
        if (!retained[bound.species] && participating[bound.species])
            throw std::invalid_argument("inventory_face: reacting species forced absent: " + model.component_ids[bound.species]);
    }
    if (!phase.composition_guess.empty() && phase.composition_guess.size() != phase.support.size())
        throw std::invalid_argument("ideal composition guess size");
    const auto declared = phase.support;
    phase.support.erase(std::remove_if(phase.support.begin(), phase.support.end(), [&](auto i) { return !retained[i]; }), phase.support.end());
    restrict_composition_guess(phase, declared);
    if (phase.support.empty()) throw std::invalid_argument("empty ideal inventory support");
    if (model.solvent_index && !std::count(phase.support.begin(), phase.support.end(), *model.solvent_index))
        throw std::invalid_argument("unit-activity solvent absent from positive liquid support");
    if (phase.support.size() != all.size()) restrict_inventory_balances(c, retained);
    for (auto i : all) if (!retained[i]) c.removed_species.push_back(i);
    c.primal_size = phase.support.size();
    for (std::size_t r = 0; r < request.reactions.size(); ++r) c.rows.push_back({model.activity_basis == "molarity" ? RowKind::ReactionMolarity : RowKind::ReactionMassAction, 0, r});
    for (std::size_t k = 0; k < c.rank; ++k) c.rows.push_back({RowKind::Moiety, 0, k});
    const auto info = inspect_sensitivity(c, c.generic_point());
    c.last_rank = core::NumericalRank{"canonical", c.primal_size, info.rank, info.threshold};
    c.jacobian_rank = info.rank; c.generic_nullity = c.primal_size - info.rank; c.feasible = info.nonsingular;
    return c;
}

CompiledSystem compile_problem(const HelmholtzModel& model, const EquilibriumProblem& request) {
    if (request.prescribed_liquid_composition) return compile_prescribed_liquid(model, request);
    core::Scope scope("compile_problem");
    return compile_with_basis(model, request, nullptr);
}

CompiledSystem compile_problem(const HelmholtzModel& model, const EquilibriumProblem& request,
                               const std::vector<std::string>& species) {
    auto compiled = compile_problem(model, request);
    if (!request.discovery) compiled.original = std::make_shared<const CompiledSystem::Declaration>(
        CompiledSystem::Declaration{request, species, compiled});
    return compiled;
}

CompiledSystem refresh_problem(const HelmholtzModel& model, const EquilibriumProblem& request,
                               const std::vector<std::string>& species, const CompiledSystem& prepared) {
    const auto& original = prepared.original;
    if (request.discovery || !original || original->compiled.model != &model || original->species != species
        || !same_declaration(original->request, request)) return compile_problem(model, request);
    if (same_state_values(original->request, request)) return original->compiled;
    core::Scope scope("refresh_problem");
    return request.prescribed_liquid_composition ? compile_prescribed_liquid(model, request, &original->compiled)
        : compile_with_basis(model, request, nullptr, &original->compiled);
}

// The reduced trial space is internal: the public no-reaction compiler retains species balances.
CompiledSystem compile_stability_trial(const HelmholtzModel& model, double T, double P,
    const PhaseSpec& phase, const std::vector<double>& w, const std::vector<double>& c,
    std::vector<std::vector<double>>& transformation) {
    if (model.parameters().ideal_gas) throw std::invalid_argument("ideal-gas stability trial compilation is deferred");
    core::Scope scope("compile_problem");
    const auto& support = phase.support;
    const auto& z = model.parameters().charge_number;
    const bool has_charge = std::any_of(support.begin(), support.end(), [&](auto i) { return z[i] != 0; });
    Eigen::MatrixXd physical(has_charge ? 2 : 1, support.size());
    physical.row(0).setOnes();
    if (has_charge) for (std::size_t i = 0; i < support.size(); ++i) physical(1, i) = z[support[i]];
    std::vector<std::vector<double>> B;
    std::vector<core::NumericalRank> ranks;
    transformation.clear();
    if (has_charge && phase.kind == PhaseKind::Vapor) {
        B.assign(2, std::vector<double>(w.size(), 0.0));
        for (std::size_t i = 0; i < support.size(); ++i) {
            B[0][support[i]] = 1.0; B[1][support[i]] = z[support[i]];
        }
        transformation = {{1.0, 0.0}, {0.0, 1.0}};
    } else {
        const auto factor = physical.fullPivLu();
        core::retain_rank(ranks, factor, "trial_normalization_charge");
        const Eigen::MatrixXd kernel = factor.kernel();
        std::vector<ReactionSpec> directions;
        if (support.size() > static_cast<std::size_t>(physical.rows())) {
            for (Eigen::Index j = 0; j < kernel.cols(); ++j) {
                ReactionSpec direction;
                direction.stoichiometry.assign(kernel.col(j).data(), kernel.col(j).data() + kernel.rows());
                // Only the basis owner's input validation consumes this zero target.
                // These directions never enter the trial's reaction declaration.
                direction.affinity_j_per_mol = 0.0;
                directions.push_back(std::move(direction));
            }
        }
        const auto reduced = balance_basis(support.size(), directions);
        B.assign(reduced.size(), std::vector<double>(w.size(), 0.0));
        for (std::size_t k = 0; k < B.size(); ++k) {
            for (std::size_t i = 0; i < support.size(); ++i) B[k][support[i]] = reduced[k][i];
            const auto transform = physical.transpose().colPivHouseholderQr();
            core::retain_rank(ranks, transform, "trial_basis_transform");
            const Eigen::VectorXd coefficients = transform.solve(
                Eigen::Map<const Eigen::VectorXd>(reduced[k].data(), reduced[k].size()));
            transformation.emplace_back(coefficients.data(), coefficients.data() + coefficients.size());
        }
    }
    EquilibriumProblem request;
    request.temperature = FreeValue::fixed(T); request.pressure = FreeValue::fixed(P);
    request.feed.values = w;
    request.phase_calculation = PhaseCalculation::Discovery;
    request.phases = {phase};
    request.phases[0].amount = AmountSpec::free();
    request.solver_options.values[policy::row] = policy::EQUILIBRIUM_ROW_RESIDUAL.value;
    for (double coefficient : c) request.standard_gibbs_j_per_mol.push_back(coefficient * GAS_CONSTANT_J_PER_MOL_K * T);
    auto compiled = compile_with_basis(model, request, &B);
    for (const auto& r : ranks) core::retain_rank(compiled.rank_records, r);
    return compiled;
}

std::vector<double> CompiledSystem::generic_point() const {
    core::Scope scope("equilibrium_generic_point");
    if (ideal) {
        std::vector<double> point;
        const auto& phase = problem.phases[0];
        const auto& guess = phase.composition_guess;
        if (!guess.empty() && (guess.size() != phase.support.size() || !all_numeric(guess)
            || *std::min_element(guess.begin(), guess.end()) <= 0.0 || std::abs(sum(guess) - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value))
            throw std::invalid_argument("invalid ideal composition guess");
        for (auto i : phase.support) point.push_back(std::log(inventory_upper_bounds[i] / phase.support.size()));
        if (!guess.empty()) for (std::size_t j = 0; j < point.size(); ++j) point[j] = std::log(guess[j]);
        return point;
    }
    std::vector<double> point(primal_size + rank + charge_row_count, 0.0);
    const auto& parameters = model->parameters();
    for (std::size_t phase = 0; phase < shape.size(); ++phase) {
        const auto& support = problem.phases[phase].support;
        double total = 0.0, solute = 0.0;
        std::vector<double> fractions(parameters.component_count(), 0.0);
        for (std::size_t i = 0; i < shape[phase].x.size(); ++i) {
            const double value = 1.0 + 0.37 * static_cast<double>(i) + 0.11 * static_cast<double>(phase)
                + (phase > 1 ? 0.05 * std::pow(static_cast<double>(i), phase) : 0.0);  // 3+ phases not coplanar
            point[shape[phase].x[i]] = fractions[support[i]] = value;
            total += value;
            if (parameters.charge_number[support[i]] != 0) solute += value;
        }
        for (auto& x : fractions) x /= total;
        const double D = ion_suppression_denominator(fractions, parameters);
        if (const double margin = policy::RANK_PROBE_PERMITTIVITY_DENOMINATOR.value; total > solute && D < margin) {
            const double scale = (1 - margin) * (total - solute) / ((1 - D) * total - (1 - margin) * solute);
            for (std::size_t i = 0; i < support.size(); ++i)
                if (parameters.charge_number[support[i]] != 0) point[shape[phase].x[i]] *= scale;
            total = total - solute + scale * solute;
        }
        for (const std::size_t coordinate : shape[phase].x)
            point[coordinate] = std::log(point[coordinate] / total);
        if (shape[phase].amount.has_value()) point[*shape[phase].amount] = 1.0;
        point[shape[phase].volume] = model->parameters().ideal_gas
            ? std::log(seed_pressure() / (GAS_CONSTANT_J_PER_MOL_K * seed_temperature()))
            : std::log(0.25 + 0.05 * static_cast<double>(phase));
    }
    for (std::size_t i = 0; i < mode_count; ++i)  // a generic mode direction, aligned with no axis
        point[mode_offset() + i] = (1.0 + 0.29 * static_cast<double>(i)) / static_cast<double>(mode_count + 1);
    if (problem.loading) point[loading_column()] = problem.loading->seed;
    std::size_t offset = primal_size - free_intensive_count;
    if (problem.temperature.free) point[offset++] = seed_temperature();
    if (problem.pressure.free) point[offset++] = seed_pressure();
    return point;
}
}  // namespace epcsaft::equilibrium
