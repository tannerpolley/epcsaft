// Numerical falsifiers of the regression layer (#46 R1, #86 density actions, #47 R2, #48 R3)
// and the native rows the public R4 transport reproduces (#49).
#include "regression/fit.hpp"

#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_support.hpp"
#include "properties_support.hpp"
#include "recording_support.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace epcsaft;
using namespace epcsaft::regression;
using equilibrium::AmountSpec;
using equilibrium::EquilibriumProblem;
using equilibrium::EquilibriumResult;
using equilibrium::FreeValue;
using equilibrium::PhaseKind;
using recording_support::field;
using recording_support::fields;
using recording_support::number;

constexpr double EPS = std::numeric_limits<double>::epsilon();
constexpr double FIRST_ORDER = 5.0e-6;  // P4's frozen first-order action criterion

// |model - reference| <= max(atol, rtol max(|model|, |reference|)); records and returns the verdict.
bool agree(const std::string& id, const char* quantity, const char* unit, double model, double reference,
           double atol, double rtol, const std::string& inputs, const char* claim, const char* group) {
    const bool verdict = std::isfinite(model) && std::isfinite(reference)
        && std::abs(model - reference) <= std::max(atol, rtol * std::max(std::abs(model), std::abs(reference)));
    recording_support::compare(id.c_str(), quantity, unit, model, reference, inputs, "Relationship", "regression",
                               claim, "symmetric_close", atol, rtol, verdict, group, id.c_str());
    return verdict;
}

bool truth(const std::string& id, const char* claim, bool ok, const std::string& inputs, const char* group) {
    recording_support::compare(id.c_str(), "verdict", "bool", ok ? 1.0 : 0.0, 1.0, inputs, "Relationship",
                               "regression", claim, "absolute", 0.0, 0.0, ok, group, id.c_str());
    return ok;
}

std::string rejection(const std::function<void()>& call) {
    try { call(); } catch (const std::invalid_argument& error) { return error.what(); }
    return "";
}

// ---- declared specifications

EquilibriumProblem pure_saturation(double temperature, double guess, bool dew = false) {
    EquilibriumProblem problem;
    problem.phases = {equilibrium_support::phase("L", PhaseKind::Liquid, dew ? AmountSpec::pinned(0.0) : AmountSpec::free()),
                      equilibrium_support::phase("V", PhaseKind::Vapor, dew ? AmountSpec::free() : AmountSpec::pinned(0.0))};
    for (auto& phase : problem.phases) phase.composition_guess = {1.0};
    problem.temperature = FreeValue::fixed(temperature);
    problem.pressure = FreeValue::variable(guess);
    problem.feed.kind = equilibrium::FeedKind::Amounts;
    problem.feed.values = {1.0};
    problem.row_tolerance = 1.0e-10;
    return problem;
}

EquilibriumProblem bubble_at(double temperature) {
    auto problem = equilibrium_support::boundary();
    problem.feed.values = {0.5, 0.5};
    problem.temperature = FreeValue::fixed(temperature);
    return problem;
}

Observation solved(ObservationKind kind, const EquilibriumProblem& problem, std::size_t phase = 0,
                   std::vector<std::size_t> species = {}, bool mass = false) {
    Observation result;
    result.kind = kind;
    result.problem = problem;
    result.phase = phase;
    result.species = std::move(species);
    result.mass_basis = mass;
    result.observed = 1.0;
    result.form = ResidualForm::LogRatio;
    return result;
}

Observation density(double temperature, double pressure, std::vector<double> composition, DensityBranch branch,
                    bool mass = false) {
    Observation result;
    result.kind = ObservationKind::Density;
    result.temperature = temperature;
    result.pressure = pressure;
    result.composition = std::move(composition);
    result.branch = branch;
    result.mass_basis = mass;
    result.observed = 1.0;
    result.form = ResidualForm::LogRatio;
    return result;
}

// The E4-B interior tie line: neutral binary at 180 K, feed at the mean of the bubble compositions,
// both amounts free at half the bubble pressure; the recorded seed only if the cold solve is refused.
struct Interior final { EquilibriumProblem problem; EquilibriumResult start; double bubble_pressure = 0.0; };
Interior interior(const HelmholtzModel& model) {
    Interior result;
    const auto bubble = equilibrium::solve_equilibrium(model, bubble_at(180.0));
    CHECK(bubble.success);
    result.bubble_pressure = bubble.pressure.front();
    result.problem = bubble_at(180.0);
    result.problem.feed.values = {0.5 * (bubble.mole_fractions[0] + bubble.mole_fractions[2]),
                                  0.5 * (bubble.mole_fractions[1] + bubble.mole_fractions[3])};
    result.problem.phases[1].amount = AmountSpec::free();
    result.problem.pressure = FreeValue::fixed(0.5 * result.bubble_pressure);
    result.start = equilibrium::solve_equilibrium(model, result.problem);
    if (!(result.start.success && result.start.requested_tolerance_met.value_or(false))) {
        const std::vector<double> seed = {-1.4805904749125276, -0.25812755981945695, 0.25979825064587636,
            -0.9468581761373848, -0.10102638484899645, -2.342461524948757, 0.7402017493541235,
            -4.6620373689869785, 6.106829374848107, 3.7251308960023377};
        result.start = equilibrium::solve_equilibrium_from_coordinates(model, result.problem, seed);
    }
    CHECK(result.start.success);
    return result;
}

// ---- R1 F1: units and basis on the retained NIST SRD 69 ethane rows (packet v1,
// source b01333e827933c0a…, parameters 7dd1899712aacfbc…, E6 probe inputs).
struct NistRow final { double T, p_source, rho_source, p_model, rho_model, rho_molar; };
const NistRow NIST[10] = {
    {100.0, 11.080787, 640.94852, 11.123831887357404, 645.7410638010867, 21474.59473897861},
    {120.0, 352.30167, 618.94997, 351.9423497001842, 622.0818292190638, 20687.78946521662},
    {140.0, 3813.5564, 596.58156, 3814.3604750250333, 598.8924184547142, 19916.60852859043},
    {160.0, 21405.224, 573.55144, 21443.595352485714, 575.8127513131142, 19149.077196977527},
    {180.0, 78638.137, 549.50874, 78838.3060997493, 552.2582967935026, 18365.75646137355},
    {200.0, 217232.94, 523.97698, 217800.6949469248, 527.4654263154658, 17541.25129083691},
    {220.0, 492046.39, 496.27145, 493283.9144780184, 500.41506993891124, 16641.671763848062},
    {240.0, 966787.79, 465.30887, 969152.4807171348, 469.5944555701788, 15616.70953010239},
    {260.0, 1711835.4, 429.07617, 1715342.640248023, 432.36903897862464, 14378.750880566167},
    {280.0, 2806735.8, 382.72712, 2807103.3101879596, 382.94892442644885, 12735.248567557328},
};

void check_f1() {
    const HelmholtzModel model(properties_support::ethane());
    const double mass = 0.03007;
    for (const auto& row : NIST) {
        const std::string point = "T" + std::to_string(static_cast<int>(row.T));
        const auto problem = pure_saturation(row.T, row.p_source);
        auto pressure = solved(ObservationKind::EquilibriumPressure, problem);
        pressure.observed = row.p_source;
        const auto p = predict(model, pressure);
        CHECK(p.status == ObservationStatus::Available && p.unit == "Pa" && p.residual.has_value());
        const auto inputs = fields({field("T_K", row.T), field("source_p_sat_Pa", row.p_source),
                                    field("source_rho_liq_kg_m3", row.rho_source), field("branch", "L free / V pinned 0")});
        CHECK(agree("r1-f1/" + point + "/p_sat", "pressure", "Pa", *p.value, row.p_model, 0.0, 1.0e-9, inputs,
                    "the O2 value reproduces the retained E6 model pressure", "r1-f1"));
        CHECK(agree("r1-f1/" + point + "/p_deviation", "relative_deviation", "1", (*p.value - row.p_source) / row.p_source,
                    (row.p_model - row.p_source) / row.p_source, 1.0e-9, 0.0, inputs,
                    "the retained NIST deviation is reproduced, observation only", "r1-f1"));
        auto rho_mass = solved(ObservationKind::PhaseDensity, problem, 0, {}, true);
        rho_mass.observed = row.rho_source;
        auto rho_molar = solved(ObservationKind::PhaseDensity, problem, 0, {}, false);
        const auto kg = predict(model, rho_mass), mol = predict(model, rho_molar);
        CHECK(kg.status == ObservationStatus::Available && kg.unit == "kg/m^3" && mol.unit == "mol/m^3");
        CHECK(agree("r1-f1/" + point + "/rho_liq", "saturated_liquid_density", "kg/m^3", *kg.value, row.rho_model,
                    0.0, 1.0e-9, inputs, "the O1b mass-basis value reproduces the retained E6 row", "r1-f1"));
        CHECK(agree("r1-f1/" + point + "/mass_chain", "saturated_liquid_density", "kg/m^3", *kg.value, *mol.value * mass,
                    0.0, 1.0e-12, inputs, "kg/m^3 equals mol/m^3 times the fixed molar mass", "r1-f1"));
    }
    const auto value = predict(model, density(200.0, 217800.6949469248, {1.0}, DensityBranch::Liquid));
    CHECK(value.status == ObservationStatus::Available);
    CHECK(agree("r1-f1/T200/o1a_density", "molar_density", "mol/m^3", *value.value, 17541.25129083691, 0.0, 1.0e-9,
                fields({field("T_K", 200.0), field("P_Pa", 217800.6949469248), field("branch", "liquid")}),
                "the pressure-closed O1a density equals the retained saturated-liquid root", "r1-f1"));
    // A neutral species always carries a positive molar mass by Properties admission; the missing-mass
    // path is reachable only for an ion without a recorded mass.
    auto ionic = properties_support::aqueous_electrolyte();
    ionic.molar_mass_kg_per_mol[1] = 0.0;
    const HelmholtzModel ions(ionic);
    const auto missing = predict(ions, density(298.15, 1.0e5, {0.96, 0.02, 0.02}, DensityBranch::Liquid, true));
    CHECK(truth("r1-f1/missing_molar_mass", "a mass-basis density with an unrecorded ionic mass is typed, not 0 kg/m^3",
                missing.status == ObservationStatus::MissingMolarMass && !missing.value.has_value(),
                fields({field("status", observation_status_name(missing.status))}), "r1-f1"));
}

// ---- R1 F2: affine coordinates, bounds and admission on all seven families.
CoordinateSet seven_families() {
    CoordinateSet set;
    set.coordinates = {
        {ParameterFamily::DispersionEnergyOverK, {{0}}, 221.0, 10.0, 150.0, 260.0},
        {ParameterFamily::SegmentCount, {{0}}, 1.9489, 0.1, 1.5, 6.0},
        {ParameterFamily::SigmaConstant, {{0}}, 3.7209, 0.1, 2.0, 4.2},
        {ParameterFamily::Kij, {{1, 2}}, 0.0, 1.0, -1.0, 1.0},
        {ParameterFamily::Lij, {{1, 2}}, 0.0, 1.0, -1.0, 1.0},
        {ParameterFamily::AssociationEnergyOverK, {{1, 4}}, 0.0, 1.0, 0.0, 5000.0},
        {ParameterFamily::AssociationVolume, {{1, 4}}, 0.0, 1.0, 0.0, 1.0},
    };
    return set;
}

const std::vector<double>& family_of(const Parameters& p, ParameterFamily family) {
    switch (family) {
    case ParameterFamily::SegmentCount: return p.segment_count;
    case ParameterFamily::SigmaConstant: return p.sigma_constant;
    case ParameterFamily::DispersionEnergyOverK: return p.dispersion_energy_over_k;
    case ParameterFamily::Kij: return p.k_ij;
    case ParameterFamily::Lij: return p.l_ij;
    case ParameterFamily::AssociationEnergyOverK: return p.association_energy_over_k;
    case ParameterFamily::AssociationVolume: return p.association_volume;
    }
    throw std::logic_error("family");
}

void check_f2() {
    const auto parameters = properties_support::associating_binary();
    const auto set = seven_families();
    validate(set, parameters);
    const auto start = set.start(parameters);
    const auto physical = set.physical(start);
    const auto active = set.active_parameters();
    CHECK(physical.size() == 7 && active.size() == 7);
    for (std::size_t r = 0; r < 7; ++r) {
        const double truth_value = family_of(parameters, active[r].family)[active[r].indices.front()];
        const double ulp = 2.0 * EPS * std::max(std::abs(truth_value), std::abs(set.coordinates[r].origin));
        CHECK(agree("r1-f2/roundtrip/" + std::to_string(r), "theta", "native", physical[r], truth_value, ulp, 0.0,
                    fields({field("origin", set.coordinates[r].origin), field("scale", set.coordinates[r].scale),
                            field("u", start[r])}), "u -> theta -> u returns the start within 2 ulp of max(|theta|, |origin|)", "r1-f2"));
    }
    const auto applied = set.apply(parameters, start);
    bool identical = true;
    for (const auto& item : active) {
        const auto& lhs = family_of(applied, item.family);
        const auto& rhs = family_of(parameters, item.family);
        identical = identical && lhs == rhs;
    }
    CHECK(truth("r1-f2/apply_identity", "applying the start reproduces every family array bit for bit", identical, "", "r1-f2"));
    // Bounds map through the same affine transform and are inclusive.
    auto lower = start;
    lower[0] = (150.0 - 221.0) / 10.0;
    CHECK(truth("r1-f2/bound_inclusive", "theta = lower is admitted", set.within_bounds(lower), field("u0", lower[0]), "r1-f2"));
    lower[0] -= 1.0e-9;
    CHECK(truth("r1-f2/bound_exclusive", "theta below lower is rejected", !set.within_bounds(lower), field("u0", lower[0]), "r1-f2"));
    struct Rejection final { const char* phrase; std::function<void(CoordinateSet&)> mutate; };
    const std::vector<Rejection> rejections = {
        {"start outside bounds", [](CoordinateSet& s) { s.coordinates[0].lower = 230.0; }},
        {"scale", [](CoordinateSet& s) { s.coordinates[1].scale = 0.0; }},
        {"bounds", [](CoordinateSet& s) { s.coordinates[2].lower = 5.0; }},
        {"slot index", [](CoordinateSet& s) { s.coordinates[1].slots = {{99}}; }},
        {"slots", [](CoordinateSet& s) { s.coordinates[1].slots.clear(); }},
        {"set is empty", [](CoordinateSet& s) { s.coordinates.clear(); }},
        {"diagonal parameter coordinate", [](CoordinateSet& s) { s.coordinates[3].slots = {{0}}; }},
        {"reciprocal parameter coordinate", [](CoordinateSet& s) { s.coordinates[3].slots = {{1}}; }},
    };
    for (const auto& item : rejections) {
        auto broken = set;
        item.mutate(broken);
        const auto message = rejection([&] { validate(broken, parameters); });
        CHECK(truth(std::string("r1-f2/reject/") + item.phrase, "rejected by name", message.find(item.phrase) != std::string::npos,
                    field("message", message), "r1-f2"));
    }
}

// ---- R1 F3 and F4: same-family sharing, reciprocal ownership, fixing and the active-list order.
void check_f3_f4() {
    auto parameters = properties_support::associating_binary();
    parameters.sigma_constant = {2.5, 2.5};
    CoordinateSet tied;
    tied.coordinates = {{ParameterFamily::SigmaConstant, {{0}, {1}}, 2.5, 0.1, 2.0, 4.0},
                        {ParameterFamily::Kij, {{1, 2}}, 0.0, 1.0, -1.0, 1.0}};
    validate(tied, parameters);
    const std::vector<std::size_t> tie = tied.tie();
    CHECK(truth("r1-f3/tie_map", "the tie map has one column per physical row: [0, 0, 1]",
                tie == std::vector<std::size_t>{0, 0, 1}, "", "r1-f3"));
    const std::vector<double> u = {0.3, -0.02};
    const auto applied = tied.apply(parameters, u);
    auto hand = parameters;
    hand.sigma_constant = {2.5 + 0.1 * 0.3, 2.5 + 0.1 * 0.3};
    hand.k_ij[1] = hand.k_ij[2] = -0.02;
    CHECK(truth("r1-f3/shared_slots", "a tied coordinate writes both slots and a pair writes both mirrored entries",
                applied.sigma_constant == hand.sigma_constant && applied.k_ij == hand.k_ij, field("u", u), "r1-f3"));
    const auto lhs = HelmholtzModel(applied).state(300.0, 5000.0, {0.4, 0.6});
    const auto rhs = HelmholtzModel(hand).state(300.0, 5000.0, {0.4, 0.6});
    CHECK(agree("r1-f3/state_identity", "residual_helmholtz", "1", lhs.residual_helmholtz, rhs.residual_helmholtz, 0.0, 0.0,
                fields({field("T", 300.0), field("rho", 5000.0)}), "the applied and hand-built models agree bit for bit", "r1-f3"));
    auto half = parameters;
    half.k_ij[1] = -0.02;
    const auto message = rejection([&] { (void)HelmholtzModel(half); });
    CHECK(truth("r1-f3/one_slot_rejected", "writing one mirrored slot only fails the model's own symmetry admission",
                message.find("k_ij symmetry") != std::string::npos, field("message", message), "r1-f3"));
    auto disagree = parameters;
    disagree.sigma_constant[1] = 2.6;
    const auto start_message = rejection([&] { (void)tied.start(disagree); });
    CHECK(truth("r1-f3/tied_start_disagrees", "tied slots with different values are rejected, never averaged",
                start_message.find("tied values disagree") != std::string::npos, field("message", start_message), "r1-f3"));
    // F4: the model's active list is the coordinate set's list in order; a fixed family projects to zero.
    const HelmholtzModel model(applied, tied.active_parameters());
    const auto& active = model.active_parameters();
    bool ordered = active.size() == 3;
    for (std::size_t r = 0; ordered && r < 3; ++r)
        ordered = active[r].family == tied.active_parameters()[r].family && active[r].indices == tied.active_parameters()[r].indices;
    CHECK(truth("r1-f4/active_order", "HelmholtzModel::active_parameters equals CoordinateSet::active in order", ordered, "", "r1-f4"));
    const auto fixed = model.pressure_actions(300.0, 1.0e5, {0.4, 0.6}, DensityBranch::Liquid, {},
        {{0.0, 0.0, {0.0, 0.0}, {{ParameterFamily::DispersionEnergyOverK, {0}, 1.0}}}});
    CHECK(fixed.density_status == PropertyActionStatus::Available && fixed.density_action.size() == 1);
    CHECK(agree("r1-f4/fixed_projects_to_zero", "density_action", "mol/m^3 per K", fixed.density_action[0], 0.0, 1.0e-12, 0.0,
                field("family", "DispersionEnergyOverK fixed"), "a well-formed coordinate absent from the active list projects to zero", "r1-f4"));
    auto outside = tied;
    outside.coordinates[0].slots = {{0}, {5}};
    const auto slot_message = rejection([&] { validate(outside, parameters); });
    CHECK(truth("r1-f4/out_of_range_rejected", "an out-of-range slot raises by name",
                slot_message.find("slot index") != std::string::npos, field("message", slot_message), "r1-f4"));
}

// ---- R1 F5 and F6: pressure/temperature specifications and the pure saturation identity.
void check_f5_f6() {
    const HelmholtzModel model(properties_support::neutral_binary());
    const auto bubble = bubble_at(180.0);
    const auto direct = equilibrium::solve_equilibrium(model, bubble);
    const auto pressure = predict(model, solved(ObservationKind::EquilibriumPressure, bubble));
    CHECK(pressure.status == ObservationStatus::Available && direct.success);
    CHECK(agree("r1-f5/bubble_pressure", "pressure", "Pa", *pressure.value, direct.pressure.front(), 0.0, 0.0,
                field("T_K", 180.0), "the O2 value is the driver's pressure bit for bit", "r1-f5"));
    // The pinned-zero free-T solve is solved cold from a guess above the answer (a fit's first trial
    // is always cold); a warm start from the O2 result would return the seed bit for bit and prove nothing.
    auto free_t = bubble;
    free_t.temperature = FreeValue::variable(190.0);
    free_t.pressure = FreeValue::fixed(*pressure.value);
    const auto temperature = predict(model, solved(ObservationKind::EquilibriumTemperature, free_t));
    CHECK(truth("r1-f5/bubble_temperature_status", "the cold pinned-zero free-T solve from above is accepted (new route)",
                temperature.status == ObservationStatus::Available, field("status", observation_status_name(temperature.status)), "r1-f5"));
    if (temperature.status == ObservationStatus::Available) {
        CHECK(agree("r1-f5/bubble_temperature", "temperature", "K", *temperature.value, 180.0, 1.0e-6, 0.0,
                    fields({field("P_Pa", *pressure.value), field("T_guess", 190.0)}),
                    "the round trip P(T) -> T(P) returns 180 K within 1e-6 K (observed 1e-11 K)", "r1-f5"));
    }
    free_t.temperature = FreeValue::variable(170.0);
    const auto below = predict(model, solved(ObservationKind::EquilibriumTemperature, free_t));
    if (below.status == ObservationStatus::Available) {
        CHECK(agree("r1-f5/bubble_temperature_from_below", "temperature", "K", *below.value, 180.0, 1e-6, 0.0,
                    field("T_guess", 170.0), "an available cold result recovers the declared bubble state", "r1-f5"));
    } else {
        CHECK(truth("r1-f5/bubble_temperature_from_below", "an unrecovered cold start retains its cause and no value",
                    (below.status == ObservationStatus::TopologyEvent
                        || below.status == ObservationStatus::CentralSolveUnavailable) && !below.value.has_value(),
                    fields({field("status", observation_status_name(below.status)), field("message", below.message),
                            field("T_guess", 170.0)}), "r1-f5"));
    }
    auto dew = bubble;
    dew.phases[0].amount = AmountSpec::pinned(0.0);
    dew.phases[1].amount = AmountSpec::free();
    const auto dew_pressure = predict(model, solved(ObservationKind::EquilibriumPressure, dew));
    CHECK(dew_pressure.status == ObservationStatus::Available);
    CHECK(truth("r1-f5/dew_distinct", "the dew specification gives a distinct pressure",
                std::abs(*dew_pressure.value / *pressure.value - 1.0) > 1.0e-3,
                fields({field("P_bubble", *pressure.value), field("P_dew", *dew_pressure.value)}), "r1-f5"));
    // F6
    const HelmholtzModel ethane(properties_support::ethane());
    const auto up = predict(ethane, solved(ObservationKind::EquilibriumPressure, pure_saturation(200.0, 2.0e5)));
    const auto down = predict(ethane, solved(ObservationKind::EquilibriumPressure, pure_saturation(200.0, 2.0e5, true)));
    CHECK(up.status == ObservationStatus::Available && down.status == ObservationStatus::Available);
    CHECK(agree("r1-f6/pure_bubble_dew", "pressure", "Pa", *up.value, *down.value, 0.0, 1.0e-8, field("T_K", 200.0),
                "a one-component bubble and dew specification give the same saturation pressure", "r1-f6"));
    for (std::size_t phase = 0; phase < 2; ++phase) {
        const auto x = predict(ethane, solved(ObservationKind::PhaseMoleFraction, pure_saturation(200.0, 2.0e5), phase, {0}));
        CHECK(agree("r1-f6/pure_fraction/" + std::to_string(phase), "mole_fraction", "1", *x.value, 1.0, 0.0, 0.0,
                    field("phase", phase), "a pure phase composition is exactly one", "r1-f6"));
    }
}

// ---- R1 F7 and F8: phase fractions on the E4-B tie line; aggregate compositions on May row 001.
void check_f7_f8() {
    const HelmholtzModel model(properties_support::neutral_binary());
    const auto tie = interior(model);
    std::vector<double> beta(2), fraction(4);
    for (std::size_t phase = 0; phase < 2; ++phase) {
        const auto b = predict(model, solved(ObservationKind::PhaseFraction, tie.problem, phase), tie.start);
        CHECK(b.status == ObservationStatus::Available);
        beta[phase] = *b.value;
        for (std::size_t i = 0; i < 2; ++i)
            fraction[phase * 2 + i] = *predict(model, solved(ObservationKind::PhaseMoleFraction, tie.problem, phase, {i}), tie.start).value;
    }
    const auto inputs = fields({field("T_K", 180.0), field("P_Pa", 0.5 * tie.bubble_pressure), field("z", tie.problem.feed.values)});
    CHECK(agree("r1-f7/fraction_sum", "sum_beta", "1", beta[0] + beta[1], 1.0, 1.0e-14, 0.0, inputs,
                "mole-basis phase fractions sum to one", "r1-f7"));
    for (std::size_t i = 0; i < 2; ++i)
        CHECK(agree("r1-f7/lever/" + std::to_string(i), "feed_fraction", "1", beta[0] * fraction[i] + beta[1] * fraction[2 + i],
                    tie.problem.feed.values[i], 1.0e-12, 0.0, inputs, "the lever rule z_i = sum_p beta_p x_pi holds", "r1-f7"));
    auto scaled = tie.problem;
    for (auto& value : scaled.feed.values) value *= 3.0;
    const auto b3 = predict(model, solved(ObservationKind::PhaseFraction, scaled, 0), tie.start);
    CHECK(b3.status == ObservationStatus::Available);
    CHECK(agree("r1-f7/fraction_homogeneous", "beta_L", "1", *b3.value, beta[0], 0.0, 1.0e-12, inputs,
                "tripling the feed leaves the phase fraction unchanged", "r1-f7"));
    CHECK(agree("r1-f7/amount_scales", "n_L", "mol", b3.solved->amounts[0], 3.0 * tie.start.amounts[0], 0.0, 1.0e-12, inputs,
                "tripling the feed triples the absolute phase amount", "r1-f7"));
    // F8: May 2015 row 001 (ref-may/may2015-ch4-c2h6-001/solve, E6 retained compositions).
    const std::vector<double> retained = {0.3025223259394309, 0.697477674060569, 0.6703563353246006, 0.3296436646753994};
    auto may = equilibrium_support::flash(2);
    may.phases = {equilibrium_support::phase("L", PhaseKind::Liquid), equilibrium_support::phase("V", PhaseKind::Vapor)};
    may.phases[0].composition_guess = {0.3099, 0.6901};
    may.phases[1].composition_guess = {0.6664, 0.3336};
    may.temperature = FreeValue::fixed(243.58);
    may.pressure = FreeValue::fixed(3.949e6);
    may.feed.values = {0.5 * (retained[0] + retained[2]), 0.5 * (retained[1] + retained[3])};
    may.row_tolerance = 1.0e-10;
    const auto may_inputs = fields({field("T_K", 243.58), field("P_Pa", 3.949e6), field("z", may.feed.values)});
    const auto both = predict(model, solved(ObservationKind::PhaseMoleFraction, may, 0, {0, 1}));
    CHECK(both.status == ObservationStatus::Available);
    CHECK(agree("r1-f8/aggregate_sum", "x_CH4+x_C2H6", "1", *both.value, 1.0, 1.0e-15, 0.0, may_inputs,
                "the full species set sums to one", "r1-f8"));
    const auto x0 = predict(model, solved(ObservationKind::PhaseMoleFraction, may, 0, {0}), both.solved);
    // The start re-solves through the driver; two accepted solves of one problem agree to rounding.
    CHECK(agree("r1-f8/single_species", "x_CH4", "1", *x0.value, both.solved->mole_fractions[0], 0.0, 1.0e-13, may_inputs,
                "a one-element set is the driver's mole fraction to solve reproducibility", "r1-f8"));
    CHECK(agree("r1-f8/retained_x", "x_CH4", "1", *x0.value, retained[0], 1.0e-9, 0.0, may_inputs,
                "the retained E6 liquid composition is reproduced", "r1-f8"));
    const auto y0 = predict(model, solved(ObservationKind::PhaseMoleFraction, may, 1, {0}), both.solved);
    CHECK(agree("r1-f8/retained_y", "y_CH4", "1", *y0.value, retained[2], 1.0e-9, 0.0, may_inputs,
                "the retained E6 vapor composition is reproduced", "r1-f8"));
    // A mass basis is defined only for the two density kinds: a two-phase flash whose mole-basis
    // phase fraction differs from the physical mass fraction is refused by name, never reinterpreted.
    auto two = bubble_at(180.0);
    two.phases[1].amount = AmountSpec::free();
    two.pressure = FreeValue::fixed(1.0e6);
    const auto liquid_fraction = predict(model, solved(ObservationKind::PhaseFraction, two, 0));
    CHECK(liquid_fraction.status == ObservationStatus::Available && liquid_fraction.solved.has_value());
    const auto& state = *liquid_fraction.solved;
    const auto& masses = model.parameters().molar_mass_kg_per_mol;
    const double mass_liquid = state.amounts[0] * (state.mole_fractions[0] * masses[0] + state.mole_fractions[1] * masses[1]);
    const double mass_vapor = state.amounts[1] * (state.mole_fractions[2] * masses[0] + state.mole_fractions[3] * masses[1]);
    const auto two_inputs = fields({field("T_K", 180.0), field("P_Pa", 1.0e6), field("amounts", state.amounts)});
    CHECK(agree("r1-f8/phase_fraction_mole_basis", "beta_L", "1", *liquid_fraction.value, 0.6738453209601837, 0.0, 1.0e-9, two_inputs,
                "the mole-basis liquid fraction of the two-phase flash", "r1-f8"));
    CHECK(agree("r1-f8/phase_fraction_mass_control", "w_L", "1", mass_liquid / (mass_liquid + mass_vapor), 0.7580747775119783, 0.0, 1.0e-9,
                two_inputs, "the physical mass fraction, which differs from the mole basis by more than 0.05", "r1-f8"));
    CHECK(std::abs(mass_liquid / (mass_liquid + mass_vapor) - *liquid_fraction.value) > 0.05);
    struct Rejection final { const char* phrase; Observation observation; };
    auto supported = may;
    supported.phases[0].support = {0};
    supported.phases[0].support_declared = true;
    const std::vector<Rejection> rejections = {
        {"species set is empty", solved(ObservationKind::PhaseMoleFraction, may, 0, {})},
        {"species repeated", solved(ObservationKind::PhaseMoleFraction, may, 0, {0, 0})},
        {"species outside the phase support", solved(ObservationKind::PhaseMoleFraction, supported, 0, {1})},
        {"phase index", solved(ObservationKind::PhaseFraction, may, 2)},
        {"mass basis is defined only for a density", solved(ObservationKind::PhaseFraction, two, 0, {}, true)},
        {"mass basis is defined only for a density", solved(ObservationKind::PhaseMoleFraction, two, 0, {0}, true)},
    };
    std::size_t rejection_index = 0;
    for (const auto& item : rejections) {
        const auto message = rejection([&] { validate(item.observation, model); });
        CHECK(truth("r1-f8/reject/" + std::to_string(rejection_index++) + "/" + item.phrase, "rejected by name",
                    message.find(item.phrase) != std::string::npos, field("message", message), "r1-f8"));
    }
}

std::string summary(const FitResult& result);

// ---- R1 F9 and F10: preserved failure causes and the residual forms.
void check_f9_f10() {
    const HelmholtzModel model(properties_support::neutral_binary());
    auto unattainable = bubble_at(180.0);
    unattainable.row_tolerance = 1.0e-18;
    const auto floor = predict(model, solved(ObservationKind::EquilibriumPressure, unattainable));
    CHECK(truth("r1-f9/unmet_tolerance", "a floor-limited rejected solve keeps its cause and returns no value",
                floor.status == ObservationStatus::CentralSolveUnavailable && !floor.value.has_value()
                    && !floor.residual.has_value() && floor.requested_tolerance_met.has_value()
                    && !*floor.requested_tolerance_met && floor.convergence.has_value()
                    && std::isfinite(floor.convergence->residual),
                fields({field("status", observation_status_name(floor.status)), field("message", floor.message),
                        field("residual", floor.convergence ? floor.convergence->residual : -1.0)}), "r1-f9"));
    auto coincident = equilibrium_support::boundary();
    coincident.feed.values = {0.5, 0.5};
    coincident.temperature.value = 280.0;
    const auto topology = predict(model, solved(ObservationKind::EquilibriumPressure, coincident));
    CHECK(truth("r1-f9/topology_event", "the s15c coincident state is typed topology_event whatever success says",
                topology.status == ObservationStatus::TopologyEvent && !topology.value.has_value()
                    && topology.message.find("phases_coincide") != std::string::npos,
                fields({field("status", observation_status_name(topology.status)), field("message", topology.message)}), "r1-f9"));
    const HelmholtzModel ethane(properties_support::ethane());
    // A declared branch selects among two stable roots; one stable root is returned as Unspecified with
    // its count, so the caller can see a branch that was not honoured. No root at all is the failure.
    const auto single = predict(ethane, density(200.0, 1.0e8, {1.0}, DensityBranch::Vapor));
    CHECK(truth("r1-f9/single_root_branch", "a declared vapor branch at a one-root state reports Unspecified and count 1",
                single.status == ObservationStatus::Available && single.branch == DensityBranch::Unspecified
                    && single.stable_root_count == std::size_t{1},
                fields({field("status", observation_status_name(single.status)), field("rho", single.value.value_or(-1.0))}), "r1-f9"));
    const auto vapor = predict(ethane, density(200.0, 1.0e10, {1.0}, DensityBranch::Liquid));
    CHECK(truth("r1-f9/density_unavailable", "a state without a stable root is typed with the thrown cause",
                vapor.status == ObservationStatus::DensityUnavailable && !vapor.value.has_value() && !vapor.message.empty(),
                fields({field("status", observation_status_name(vapor.status)), field("message", vapor.message),
                        field("exception_shape", vapor.convergence.has_value() ? "ConvergenceFailure with record" : "runtime_error without record")}),
                "r1-f9"));
    // F10
    Observation difference = density(200.0, 1.0e5, {1.0}, DensityBranch::Liquid);
    difference.form = ResidualForm::Difference;
    difference.observed = 1.0;
    difference.scale = 0.5;
    ObservationStatus status = ObservationStatus::Available;
    CHECK(agree("r1-f10/difference", "residual", "1", *residual(difference, 3.0, status), 4.0, 0.0, 0.0,
                fields({field("y", 3.0), field("y_obs", 1.0), field("scale", 0.5)}), "(3 - 1) / 0.5 is exactly 4", "r1-f10"));
    Observation ratio = difference;
    ratio.form = ResidualForm::LogRatio;
    CHECK(agree("r1-f10/log_ratio", "residual", "1", *residual(ratio, std::exp(1.0), status), 1.0, 4.0 * EPS, 0.0,
                fields({field("y", std::exp(1.0)), field("y_obs", 1.0)}), "ln(e / 1) is one within 4 ulp", "r1-f10"));
    ratio.observed = std::exp(1.0);
    CHECK(agree("r1-f10/log_ratio_swapped", "residual", "1", *residual(ratio, 1.0, status), -1.0, 4.0 * EPS, 0.0,
                fields({field("y", 1.0), field("y_obs", std::exp(1.0))}), "ln(1 / e) is minus one", "r1-f10"));
    ratio.observed = 0.0;
    const auto rejected = residual(ratio, 1.0, status);
    CHECK(truth("r1-f10/nonpositive_log", "a nonpositive log argument is typed, not NaN",
                !rejected.has_value() && status == ObservationStatus::NonpositiveLogArgument,
                field("status", observation_status_name(status)), "r1-f10"));
    difference.scale = 0.0;
    const auto message = rejection([&] { validate(difference, ethane); });
    CHECK(truth("r1-f10/scale_rejected", "a nonpositive Difference scale is rejected by name",
                message.find("scale") != std::string::npos, field("message", message), "r1-f10"));
    // The log ratio is finite across the dynamic range and keeps its accuracy near unity.
    ratio.observed = 1.0e-308;
    const double wide = *residual(ratio, 19936.48319859253, status);
    CHECK(agree("r1-f10/log_ratio_wide_range", "residual", "1", wide, std::log(19936.48319859253) - std::log(1.0e-308), 0.0, 1.0e-12,
                fields({field("y", 19936.48319859253), field("y_obs", 1.0e-308)}),
                "a ratio that overflows still gives the finite log residual (finiteness and branch selection; the formula is the reference)", "r1-f10"));
    // Off an exactly divisible base the naive log(y / y_obs) loses five digits; the reference is the
    // 60-digit decimal ln y - ln y_obs of the rounded doubles.
    ratio.observed = 17541.251290836932;
    const double near = ratio.observed * (1.0 + 1.0e-12), naive = std::log(near / ratio.observed);
    CHECK(agree("r1-f10/log_ratio_near_unity", "residual", "1", *residual(ratio, near, status), 1.0000617126413953e-12, 0.0, 1.0e-9,
                fields({field("y", near), field("y_obs", ratio.observed), field("naive_log", naive)}),
                "the residual near unity keeps its leading digits (log1p), where log(y / y_obs) is off by 3e-5", "r1-f10"));
    CHECK(std::abs(naive - 1.0000617126413953e-12) > 1.0e-6 * 1.0000617126413953e-12);
    Observation overflow = density(200.0, 1.0e5, {1.0}, DensityBranch::Liquid);
    overflow.form = ResidualForm::Difference;
    overflow.observed = -1.0e308;
    overflow.scale = 1.0e-308;
    status = ObservationStatus::Available;
    const auto overflowed = residual(overflow, 1.0, status);
    CHECK(truth("r1-f10/difference_overflow", "an overflowing Difference residual is typed nonfinite, not inf",
                !overflowed.has_value() && status == ObservationStatus::Nonfinite, field("status", observation_status_name(status)), "r1-f10"));
    const auto kept = predict(ethane, overflow);
    CHECK(truth("r1-f10/value_survives_nonfinite_residual", "the prediction keeps its finite value beside a nonfinite residual",
                kept.status == ObservationStatus::Nonfinite && kept.value.has_value() && !kept.residual.has_value() && !kept.message.empty(),
                fields({field("status", observation_status_name(kept.status)), field("value", kept.value.value_or(-1.0))}), "r1-f10"));
    FitRequest request;
    request.parameters = properties_support::ethane();
    request.coordinates.coordinates = {{ParameterFamily::DispersionEnergyOverK, {{0}}, 191.42, 10.0, 150.0, 260.0}};
    request.observations = {overflow};
    const auto refused = fit(request);
    CHECK(truth("r1-f10/fit_nonfinite_typed", "a nonfinite residual at the start is a typed initial failure with its row",
                refused.status == "initial_evaluation_failed" && refused.failures.size() == 1 && refused.failures[0].observation == 0
                    && refused.failures[0].status == ObservationStatus::Nonfinite && !refused.usable && std::isfinite(refused.initial_cost),
                summary(refused), "r1-f10"));
    // A finite residual whose weighted form overflows is recorded by the initial-evaluation recorder with its row.
    Observation steep = density(200.0, 1.0e5, {1.0}, DensityBranch::Liquid);
    steep.form = ResidualForm::Difference;
    steep.observed = *kept.value - 1.0e4;
    steep.scale = 1.0e-300;
    request.observations = {steep};
    request.weights = {1.0e40};
    const auto weighted = fit(request);
    CHECK(truth("r1-f10/weighted_overflow_typed", "an overflowing weighted residual is a typed initial failure naming its row",
                weighted.status == "initial_evaluation_failed" && weighted.failures.size() == 1 && weighted.failures[0].observation == 0
                    && weighted.failures[0].status == ObservationStatus::Nonfinite
                    && weighted.failures[0].message == "weighted residual not finite" && weighted.initial_cost == 0.0,
                summary(weighted), "r1-f10"));
    // A finite weighted residual whose optimizer Jacobian overflows is recorded by the trial and final recorders.
    steep.observed = *kept.value - 1.0;
    steep.scale = 1.0e-150;
    request.observations = {steep};
    request.weights.clear();
    request.coordinates.coordinates = {{ParameterFamily::DispersionEnergyOverK, {{0}}, 191.42, 1.0e200, 150.0, 260.0}};
    const auto steep_fit = fit(request);
    const bool trial_recorded = std::any_of(steep_fit.failures.begin(), steep_fit.failures.end(), [](const TrialFailure& item) {
        return item.status == ObservationStatus::Nonfinite && item.observation == 0 && item.message == "optimizer Jacobian not finite"; });
    const bool final_recorded = std::any_of(steep_fit.failures.begin(), steep_fit.failures.end(), [](const TrialFailure& item) {
        return item.observation == 0 && item.message.find("final optimizer Jacobian not finite") != std::string::npos; });
    CHECK(truth("r1-f10/jacobian_overflow_typed", "an overflowing optimizer Jacobian is recorded by the trial and final recorders and the fit fails typed",
                trial_recorded && final_recorded && steep_fit.status == "failed" && !steep_fit.usable, summary(steep_fit), "r1-f10"));
    // Three finite rows whose squares are representable but whose accumulated objective is not: the
    // failure is typed at the row where the running sum stops being finite, never an infinite cost.
    Observation steep_row = density(200.0, 1.0e5, {1.0}, DensityBranch::Liquid);
    steep_row.form = ResidualForm::Difference;
    steep_row.observed = *kept.value - 1.0;
    steep_row.scale = 1.0 / 1.2e154;
    request.observations = {steep_row, steep_row, steep_row};
    request.weights.clear();
    request.coordinates.coordinates = {{ParameterFamily::DispersionEnergyOverK, {{0}}, 191.42, 10.0, 150.0, 260.0}};
    const auto accumulated = fit(request);
    CHECK(truth("r1-f10/accumulated_cost_typed", "an accumulated objective that overflows is a typed initial failure naming the row where the sum stops being finite",
                accumulated.status == "initial_evaluation_failed" && !accumulated.usable && accumulated.failures.size() == 1
                    && accumulated.failures[0].observation == 1 && accumulated.failures[0].status == ObservationStatus::Nonfinite
                    && accumulated.failures[0].message == "accumulated weighted cost not finite" && accumulated.initial_cost == 0.0
                    && accumulated.statuses.size() == 3,
                fields({summary(accumulated), field("row_residual", *predict(ethane, steep_row).residual)}), "r1-f10"));
}

// ---- #86: first-order pressure-chart density actions against rebuilt-model closures.
struct Ladder final { double coarse = 0.0, fine = 0.0; };
template <class Value>
Ladder ladder(const Value& value, double step) {
    Ladder result;
    result.coarse = (value(step) - value(-step)) / (2.0 * step);
    result.fine = (value(0.5 * step) - value(-0.5 * step)) / step;
    return result;
}

bool ladder_agrees(const std::string& id, const char* quantity, double action, const Ladder& rungs, const std::string& inputs,
                   const char* claim, const char* group) {
    // The refined rung meets the criterion; the coarse rung may carry up to the fourfold
    // truncation of a second-order central difference.
    const double scale = std::max(1.0, std::abs(rungs.fine));
    const bool fine = agree(id + "/fine", quantity, "native", action, rungs.fine, FIRST_ORDER * scale, 0.0, inputs, claim, group);
    const bool coarse = agree(id + "/coarse", quantity, "native", action, rungs.coarse, 4.0 * FIRST_ORDER * scale, 0.0, inputs, claim, group);
    return fine && coarse;
}

void check_d86() {
    const std::vector<ActiveParameter> active = {{ParameterFamily::DispersionEnergyOverK, {0}},
                                                 {ParameterFamily::SegmentCount, {0}},
                                                 {ParameterFamily::SigmaConstant, {0}}};
    const double T = 200.0, P = 217800.6949469248;
    const HelmholtzModel model(properties_support::ethane(), active);
    std::vector<PressureActionDirection> directions;
    for (const auto& item : active) directions.push_back({0.0, 0.0, {0.0}, {{item.family, item.indices, 1.0}}});
    const auto batch = model.pressure_actions(T, P, {1.0}, DensityBranch::Liquid, {}, directions);
    CHECK(batch.density_status == PropertyActionStatus::Available && batch.density_action.size() == 3
          && batch.results.empty() && batch.molar_density.has_value());
    const auto central = model.state_at_pressure(T, P, {1.0}, DensityBranch::Liquid);
    CHECK(agree("d86/density_value", "molar_density", "mol/m^3", *batch.molar_density, central.molar_density, 0.0, 0.0,
                fields({field("T_K", T), field("P_Pa", P)}), "the density-only batch carries the resolved root", "d86"));
    const double steps[3] = {0.05, 5.0e-4, 5.0e-4};
    const char* names[3] = {"epsilon", "m", "sigma"};
    for (std::size_t k = 0; k < 3; ++k) {
        const auto rebuilt = [&](double shift) {
            auto p = properties_support::ethane();
            (k == 0 ? p.dispersion_energy_over_k : k == 1 ? p.segment_count : p.sigma_constant)[0] += shift;
            return p;
        };
        const auto closure = [&](double shift) {
            return HelmholtzModel(rebuilt(shift)).state_at_pressure(T, P, {1.0}, DensityBranch::Liquid).molar_density;
        };
        const auto inputs = fields({field("T_K", T), field("P_Pa", P), field("step", steps[k])});
        CHECK(ladder_agrees(std::string("d86/parameter/") + names[k], "drho_dtheta", batch.density_action[k], ladder(closure, steps[k]),
                            inputs, "the pressure-chart density action equals the rebuilt-model closure ladder", "d86"));
        // The implicit identity: drho/dtheta|P = -(dP/dtheta|rho) / (dP/drho), the explicit pressure slope from
        // rebuilt models at the resolved density; an explicit-only stand-in (zero) is off by the whole action.
        const auto explicit_pressure = [&](double shift) {
            return HelmholtzModel(rebuilt(shift)).state(T, central.molar_density, {1.0}).pressure;
        };
        const auto slope = ladder(explicit_pressure, steps[k]);
        CHECK(agree(std::string("d86/implicit/") + names[k], "drho_dtheta", "native", batch.density_action[k],
                    -slope.fine / central.pressure_density_derivative, FIRST_ORDER * std::max(1.0, std::abs(batch.density_action[k])), 0.0,
                    inputs, "the implicit lift equals minus the explicit pressure slope over dP/drho", "d86"));
        CHECK(truth(std::string("d86/omitted_term/") + names[k], "omitting the implicit term leaves an order-unity error",
                    std::abs(batch.density_action[k]) > 1.0e3 * FIRST_ORDER * std::max(1.0, std::abs(batch.density_action[k])),
                    field("action", batch.density_action[k]), "d86"));
    }
    // State directions: temperature, pressure and a zero-sum composition direction on the neutral binary.
    const auto thermal = model.pressure_actions(T, P, {1.0}, DensityBranch::Liquid, {}, {{1.0, 0.0, {0.0}, {}}, {0.0, 1.0, {0.0}, {}}});
    CHECK(thermal.density_status == PropertyActionStatus::Available && thermal.density_action.size() == 2);
    const auto over_t = [&](double shift) { return model.state_at_pressure(T + shift, P, {1.0}, DensityBranch::Liquid).molar_density; };
    CHECK(ladder_agrees("d86/temperature", "drho_dT", thermal.density_action[0], ladder(over_t, 0.05),
                        fields({field("T_K", T), field("P_Pa", P)}), "drho/dT at fixed P equals the closure ladder", "d86"));
    CHECK(agree("d86/pressure", "drho_dP", "mol/(m^3 Pa)", thermal.density_action[1], 1.0 / central.pressure_density_derivative, 0.0, 1.0e-12,
                field("dP_drho", central.pressure_density_derivative), "drho/dP is the reciprocal pressure slope", "d86"));
    const HelmholtzModel binary(properties_support::neutral_binary());
    const double Tb = 180.0, Pb = 2.0 * equilibrium::solve_equilibrium(binary, bubble_at(180.0)).pressure.front();
    const auto mixed = binary.pressure_actions(Tb, Pb, {0.5, 0.5}, DensityBranch::Liquid, {}, {{0.0, 0.0, {0.01, -0.01}, {}}});
    CHECK(mixed.density_status == PropertyActionStatus::Available);
    const auto over_x = [&](double shift) {
        return binary.state_at_pressure(Tb, Pb, {0.5 + 0.01 * shift, 0.5 - 0.01 * shift}, DensityBranch::Liquid).molar_density;
    };
    CHECK(ladder_agrees("d86/composition", "drho_dx", mixed.density_action[0], ladder(over_x, 0.1),
                        fields({field("T_K", Tb), field("P_Pa", Pb), field("direction", std::vector<double>{0.01, -0.01})}),
                        "a zero-sum composition direction at fixed P equals the closure ladder", "d86"));
    // Failure meaning and the association route.
    const auto absent = model.pressure_actions(T, 1.0e10, {1.0}, DensityBranch::Liquid, {}, directions);
    CHECK(truth("d86/density_unavailable", "a state without a stable root keeps density_unavailable and no action",
                absent.density_status == PropertyActionStatus::DensityUnavailable && absent.density_action.empty()
                    && !absent.molar_density.has_value(), "", "d86"));
    // The root search's cause travels with the batch; a record-less throw is an absent record, distinct
    // from the measured (finite, possibly zero) pressure residual of a resolved root.
    // An anchor outside the density domain throws with a record: the batch carries it, with its stop.
    const auto anchored = model.pressure_actions(T, P, {1.0}, DensityBranch::Liquid, {}, directions, 1.0e-6);
    CHECK(truth("d86/density_failure_record", "a root search that throws with a convergence record carries it (stop left_domain)",
                anchored.density_status == PropertyActionStatus::DensityUnavailable && anchored.density_convergence.has_value()
                    && anchored.density_convergence->stop == core::Stop::left_domain && !anchored.density_message.empty(),
                fields({field("message", anchored.density_message), field("anchor", 1.0e-6)}), "d86"));
    CHECK(truth("d86/density_failure_cause", "the failed root search keeps its message and its absent record",
                absent.density_message.find("no stable density root") != std::string::npos && !absent.density_convergence.has_value()
                    && batch.density_message.empty() && !batch.density_convergence.has_value()
                    && std::isfinite(batch.diagnostics.pressure_residual),
                fields({field("message", absent.density_message), field("resolved_pressure_residual", batch.diagnostics.pressure_residual)}), "d86"));
    const auto message = rejection([&] { (void)model.density_actions(T, 1000.0, {1.0}, {}, {{0.0, 0.0, {0.0}, {}}}); });
    CHECK(truth("d86/density_chart_rejected", "a density-chart request without observables is rejected by name",
                message.find("observables") != std::string::npos, field("message", message), "d86"));
    const auto support = binary.pressure_actions(Tb, Pb, {1.0, 0.0}, DensityBranch::Liquid, {}, {{0.0, 0.0, {-0.01, 0.01}, {}}});
    CHECK(truth("d86/support_change", "a direction into an absent species is typed support_change with no action",
                support.density_status == PropertyActionStatus::SupportChange && support.density_action.empty(),
                field("status", static_cast<int>(support.density_status)), "d86"));
    // Unfalsified arms of the density status: AssociationResidual and the ill-conditioned/singular site causes
    // share the caloric gates but the P4 stress fixture throws in state() before a pressure-chart request forms;
    // PressureSingular is unreachable here because the root search admits only positive pressure slopes.
    const std::vector<ActiveParameter> site = {{ParameterFamily::AssociationEnergyOverK, {1, 4}}};
    const HelmholtzModel associating(properties_support::associating_binary(), site);
    const auto pair = associating.pressure_actions(300.0, 1.0e5, {0.4, 0.6}, DensityBranch::Liquid, {}, {{0.0, 0.0, {0.0, 0.0}, {{site[0].family, site[0].indices, 1.0}}}});
    CHECK(pair.density_status == PropertyActionStatus::Available);
    const auto over_site = [&](double shift) {
        auto p = properties_support::associating_binary();
        p.association_energy_over_k[1] += shift;
        p.association_energy_over_k[4] += shift;
        return HelmholtzModel(p).state_at_pressure(300.0, 1.0e5, {0.4, 0.6}, DensityBranch::Liquid).molar_density;
    };
    CHECK(ladder_agrees("d86/association", "drho_depsAB", pair.density_action[0], ladder(over_site, 1.0),
                        fields({field("T_K", 300.0), field("P_Pa", 1.0e5)}), "the association route carries the implicit site chain", "d86"));
}

// (J^T J)_{ab} of a result's n x 2 optimizer Jacobian.
double gram_entry(const FitResult& result, std::size_t a, std::size_t b) {
    double gram = 0.0;
    const std::size_t rows = result.optimizer_jacobian.size() / 2;
    for (std::size_t i = 0; i < rows; ++i) gram += result.optimizer_jacobian[i * 2 + a] * result.optimizer_jacobian[i * 2 + b];
    return gram;
}

// ---- R2 J1: prediction gradients against rebuilt-model central differences, then the residual chain.
CoordinateSet binary_coordinates() {
    CoordinateSet set;
    set.coordinates = {{ParameterFamily::DispersionEnergyOverK, {{0}}, 150.0, 10.0, 130.0, 170.0},
                       {ParameterFamily::Kij, {{1, 2}}, 0.0, 0.1, -0.2, 0.2},
                       {ParameterFamily::SegmentCount, {{1}}, 1.6, 0.1, 1.0, 2.5}};
    return set;
}

void check_j1() {
    const auto parameters = properties_support::neutral_binary();
    const auto set = binary_coordinates();
    validate(set, parameters);
    const HelmholtzModel model(parameters, set.active_parameters());
    const auto tie = interior(model);
    const auto bubble = bubble_at(180.0);
    const auto bubble_solve = equilibrium::solve_equilibrium(model, bubble);
    const double P_bubble = bubble_solve.pressure.front();
    auto free_t = bubble;
    free_t.temperature = FreeValue::variable(180.0);
    free_t.pressure = FreeValue::fixed(P_bubble);
    struct Case final { const char* name; Observation observation; };
    std::vector<Case> cases = {
        {"bubble_pressure", solved(ObservationKind::EquilibriumPressure, bubble)},
        {"bubble_temperature", solved(ObservationKind::EquilibriumTemperature, free_t)},
        {"vapor_methane", solved(ObservationKind::PhaseMoleFraction, bubble, 1, {0})},
        {"phase_fraction", solved(ObservationKind::PhaseFraction, tie.problem, 0)},
        {"liquid_mass_density", solved(ObservationKind::PhaseDensity, bubble, 0, {}, true)},
        {"density", density(180.0, 2.0 * P_bubble, {0.5, 0.5}, DensityBranch::Liquid, true)},
    };
    const auto u0 = set.start(parameters);
    for (const auto& item : cases) {
        const std::optional<EquilibriumResult> start = item.observation.kind == ObservationKind::PhaseFraction
            ? std::optional<EquilibriumResult>(tie.start)
            : item.observation.kind == ObservationKind::EquilibriumTemperature ? std::optional<EquilibriumResult>(bubble_solve)
            : std::nullopt;
        const auto prediction = predict(model, item.observation, start);
        CHECK(prediction.status == ObservationStatus::Available);
        const auto gradient = differentiate(model, item.observation, prediction);
        CHECK(gradient.status == ObservationStatus::Available && gradient.value.size() == 3);
        for (std::size_t c = 0; c < 3; ++c) {
            const auto value = [&](double shift) {
                auto u = u0;
                u[c] += shift;
                const HelmholtzModel rebuilt(set.apply(parameters, u), set.active_parameters());
                const auto p = predict(rebuilt, item.observation, start);
                if (p.status != ObservationStatus::Available) throw std::runtime_error("rebuilt prediction failed");
                return *p.value;
            };
            const double scale = set.coordinates[c].scale;
            auto rungs = ladder(value, 0.005);
            rungs.coarse /= scale;
            rungs.fine /= scale;
            CHECK(ladder_agrees(std::string("r2-j1/") + item.name + "/" + std::to_string(c), "dy_dtheta", gradient.value[c], rungs,
                                fields({field("value", *prediction.value), field("unit", prediction.unit), field("u_step", 0.005)}),
                                "the exact prediction gradient equals the rebuilt-model central difference", "r2-j1"));
        }
    }
    // The residual chain: optimizer_jacobian = sqrt(w) (dr/dy) (dy/dtheta) S diag(scale) with a tie and w != 1.
    auto tied = parameters;
    tied.sigma_constant = {3.6, 3.6};
    CoordinateSet shared;
    shared.coordinates = {{ParameterFamily::SigmaConstant, {{0}, {1}}, 3.6, 0.1, 3.0, 4.0},
                          {ParameterFamily::DispersionEnergyOverK, {{0}}, 150.0, 10.0, 130.0, 170.0}};
    FitRequest request;
    request.parameters = tied;
    request.coordinates = shared;
    request.observations = {solved(ObservationKind::EquilibriumPressure, bubble),
                            solved(ObservationKind::PhaseMoleFraction, bubble, 1, {0})};
    request.observations[0].observed = 1.0e6;
    request.observations[1].form = ResidualForm::Difference;
    request.observations[1].observed = 0.9;
    request.observations[1].scale = 0.01;
    request.weights = {2.0, 0.5};
    request.controls.maximum_iterations = 1;
    const auto result = fit(request);
    CHECK(truth("r2-j1/chain_setup", "the one-iteration tied fit is usable", result.usable,
                fields({summary(result), field("message", result.message)}), "r2-j1"));
    CHECK(result.optimizer_jacobian.size() == 4 && result.physical_jacobian.size() == 6);
    for (std::size_t i = 0; i < 2; ++i) {
        const double slope = std::sqrt(request.weights[i]) * residual_slope(request.observations[i], result.predictions[i]);
        const double* j = result.physical_jacobian.data() + i * 3;
        const double hand[2] = {slope * (j[0] + j[1]) * 0.1, slope * j[2] * 10.0};
        // A tied sum may cancel (row 1: y_CH4 against a shared sigma), so the floor is the sum of the terms.
        const double magnitude[2] = {std::abs(slope) * (std::abs(j[0]) + std::abs(j[1])) * 0.1, std::abs(slope * j[2]) * 10.0};
        for (std::size_t c = 0; c < 2; ++c)
            CHECK(agree("r2-j1/chain/" + std::to_string(i) + "/" + std::to_string(c), "dr_du", "1", result.optimizer_jacobian[i * 2 + c],
                        hand[c], 1.0e-13 * magnitude[c], 1.0e-13, fields({field("weight", request.weights[i]), field("form", i == 0 ? "LogRatio" : "Difference")}),
                        "the optimizer Jacobian is the weighted residual slope times the tied, scaled physical gradient", "r2-j1"));
    }
}

// ---- R2 J2 and R3: manufactured recovery, bounds, failures, weights, rank and covariance.
struct Truth final { Parameters parameters; std::vector<Observation> observations; };
Truth manufactured() {
    Truth result;
    result.parameters = properties_support::neutral_binary();
    result.parameters.dispersion_energy_over_k[0] = 152.0;
    result.parameters.k_ij[1] = result.parameters.k_ij[2] = -0.02;
    const HelmholtzModel model(result.parameters);
    for (const double T : {170.0, 180.0, 190.0}) result.observations.push_back(solved(ObservationKind::EquilibriumPressure, bubble_at(T)));
    auto y = solved(ObservationKind::PhaseMoleFraction, bubble_at(180.0), 1, {0});
    y.form = ResidualForm::Difference;
    y.scale = 0.01;
    result.observations.push_back(y);
    const double P_bubble = equilibrium::solve_equilibrium(model, bubble_at(180.0)).pressure.front();
    result.observations.push_back(density(180.0, 2.0 * P_bubble, {0.5, 0.5}, DensityBranch::Liquid, true));
    for (auto& observation : result.observations) observation.observed = *predict(model, observation).value;
    return result;
}

FitRequest recovery_request(const Truth& truth_case, double epsilon, double kij, double upper = 170.0) {
    FitRequest request;
    request.parameters = properties_support::neutral_binary();
    request.parameters.dispersion_energy_over_k[0] = epsilon;
    request.parameters.k_ij[1] = request.parameters.k_ij[2] = kij;
    request.coordinates.coordinates = {{ParameterFamily::DispersionEnergyOverK, {{0}}, 150.0, 10.0, 140.0, upper},
                                       {ParameterFamily::Kij, {{1, 2}}, 0.0, 0.1, -0.2, 0.2}};
    request.observations = truth_case.observations;
    return request;
}

std::string summary(const FitResult& result) {
    return fields({field("status", result.status), field("termination", result.termination), field("iterations", result.iterations),
                   field("residual_evaluations", result.residual_evaluations), field("jacobian_evaluations", result.jacobian_evaluations),
                   field("trial_failures", result.trial_failures), field("final_cost", result.final_cost),
                   field("covariance", result.covariance.status)});
}

void check_j2() {
    const auto truth_case = manufactured();
    const auto result = fit(recovery_request(truth_case, 148.0, 0.0));
    CHECK(truth("r2-j2/converged", "the noise-free recovery converges with no trial failure",
                result.status == "converged" && result.usable && result.trial_failures == 0, summary(result), "r2-j2"));
    CHECK(agree("r2-j2/epsilon", "dispersion_energy_over_k", "K", result.physical[0], 152.0, 0.0, 1.0e-6, summary(result),
                "the manufactured epsilon is recovered", "r2-j2"));
    CHECK(agree("r2-j2/kij", "k_ij", "1", result.physical[1], -0.02, 1.0e-6, 0.0, summary(result), "the manufactured k_ij is recovered", "r2-j2"));
    CHECK(agree("r2-j2/final_cost", "cost", "1", result.final_cost, 0.0, 1.0e-14, 0.0, summary(result), "the final cost vanishes", "r2-j2"));
    CHECK(truth("r2-j2/identifiable", "the interior full-rank fit carries a covariance",
                result.covariance.status == "available" && result.covariance.rank == 2, summary(result), "r2-j2"));
    auto unfinished = recovery_request(truth_case, 148.0, 0.0);
    unfinished.controls.maximum_iterations = 1;
    const auto iterate = fit(unfinished);
    CHECK(truth("r2-j2/not_converged", "an iteration-limited fit keeps its candidate and rank but no covariance",
                iterate.status == "limit_reached" && iterate.usable && iterate.covariance.status == "not_converged"
                    && iterate.covariance.rank == 2 && iterate.covariance.optimizer.empty() && iterate.covariance.physical.empty()
                    && !iterate.covariance.variance_factor.has_value() && iterate.final_cost > 1.0e-6, summary(iterate), "r2-j2"));
    const auto other = fit(recovery_request(truth_case, 156.0, 0.05));
    CHECK(agree("r2-j2/second_start/epsilon", "dispersion_energy_over_k", "K", other.physical[0], 152.0, 0.0, 1.0e-6, summary(other),
                "a second start recovers the same parameters", "r2-j2"));
    CHECK(agree("r2-j2/second_start/kij", "k_ij", "1", other.physical[1], -0.02, 1.0e-6, 0.0, summary(other),
                "a second start recovers the same parameters", "r2-j2"));
    const auto bounded = fit(recovery_request(truth_case, 148.0, 0.0, 151.0));
    CHECK(truth("r2-j2/active_bound", "a bound below the truth ends at the bound with covariance withheld",
                bounded.usable && bounded.physical[0] == 151.0 && bounded.covariance.status == "active_bound"
                    && bounded.covariance.active_bounds == std::vector<std::size_t>{0}, summary(bounded), "r2-j2"));
    auto failing = recovery_request(truth_case, 148.0, 0.0);
    failing.observations.push_back(density(200.0, 1.0e10, {0.5, 0.5}, DensityBranch::Liquid));
    const auto failed = fit(failing);
    const std::vector<ObservationStatus> expected_statuses = {ObservationStatus::Available, ObservationStatus::Available,
        ObservationStatus::Available, ObservationStatus::Available, ObservationStatus::Available, ObservationStatus::DensityUnavailable};
    CHECK(truth("r2-j2/initial_failure", "a failing observation at the start is retained with its cause and index",
                failed.status == "initial_evaluation_failed" && !failed.usable && failed.failures.size() == 1
                    && failed.failures[0].observation == 5 && failed.failures[0].status == ObservationStatus::DensityUnavailable
                    && failed.predictions.empty() && failed.statuses == expected_statuses, summary(failed), "r2-j2"));
    // R4 rows: the public transport reproduces these native values for the same explicit inputs.
    const auto row = [&](const char* id, const char* unit, double value) {
        recording_support::observe(id, "native", unit, value, std::nullopt, summary(result), "Relationship", "regression",
                                   "native fit value for the public R4 comparison", "r4");
    };
    row("r4/recovery/epsilon", "K", result.physical[0]);
    row("r4/recovery/kij", "1", result.physical[1]);
    row("r4/recovery/final_cost", "1", result.final_cost);
    row("r4/recovery/initial_cost", "1", result.initial_cost);
    for (std::size_t i = 0; i < result.predictions.size(); ++i)
        row((std::string("r4/recovery/prediction/") + std::to_string(i)).c_str(), "native", result.predictions[i]);
    for (std::size_t k = 0; k < 4; ++k) {
        row((std::string("r4/recovery/covariance/") + std::to_string(k)).c_str(), "1", result.covariance.optimizer[k]);
        row((std::string("r4/recovery/physical_covariance/") + std::to_string(k)).c_str(), "1", result.covariance.physical[k]);
    }
    row("r4/recovery/variance_factor", "1", result.covariance.variance_factor.value_or(-1.0));
    for (std::size_t k = 0; k < result.optimizer_jacobian.size(); ++k)
        row((std::string("r4/recovery/optimizer_jacobian/") + std::to_string(k)).c_str(), "1", result.optimizer_jacobian[k]);
}

void check_r3() {
    const auto truth_case = manufactured();
    auto request = recovery_request(truth_case, 148.0, 0.0);
    const auto base = fit(request);
    request.weights.assign(request.observations.size(), 2.0);
    const auto doubled = fit(request);
    CHECK(base.usable && doubled.usable && base.covariance.status == "available" && doubled.covariance.status == "available");
    for (std::size_t c = 0; c < 2; ++c)
        CHECK(agree("r3/weight_scale/parameter/" + std::to_string(c), "theta", "native", doubled.physical[c], base.physical[c], 0.0, 1.0e-9,
                    summary(doubled), "doubling every weight leaves the fitted parameters unchanged", "r3"));
    for (std::size_t k = 0; k < 4; ++k)
        CHECK(agree("r3/weight_scale/covariance/" + std::to_string(k), "covariance", "1", doubled.covariance.optimizer[k],
                    0.5 * base.covariance.optimizer[k], 0.0, 1.0e-6, summary(doubled), "doubling every weight halves the covariance", "r3"));
    CHECK(agree("r3/weight_scale/cost", "cost", "1", doubled.initial_cost, 2.0 * base.initial_cost, 0.0, 1.0e-12, summary(doubled),
                "doubling every weight doubles the cost", "r3"));
    // A zero weight removes a row from the objective but not from the report.
    auto zero = recovery_request(truth_case, 148.0, 0.0);
    auto inconsistent = solved(ObservationKind::EquilibriumPressure, bubble_at(175.0));
    inconsistent.observed = 1.5 * *predict(HelmholtzModel(truth_case.parameters), inconsistent).value;
    zero.observations.push_back(inconsistent);
    zero.weights.assign(zero.observations.size(), 1.0);
    zero.weights.back() = 0.0;
    const auto ignored = fit(zero);
    CHECK(ignored.usable && ignored.residuals.size() == 6 && ignored.statuses.size() == 6);
    CHECK(agree("r3/zero_weight/epsilon", "dispersion_energy_over_k", "K", ignored.physical[0], base.physical[0], 0.0, 1.0e-9, summary(ignored),
                "a zero-weighted inconsistent row does not move the fit", "r3"));
    CHECK(agree("r3/zero_weight/reported_residual", "residual", "1", ignored.residuals[5], std::log(1.0 / 1.5), 0.0, 1.0e-6, summary(ignored),
                "the zero-weighted row keeps its residual in the report", "r3"));
    CHECK(agree("r3/zero_weight/optimizer_residual", "weighted_residual", "1", ignored.weighted_residuals[5], 0.0, 0.0, 0.0, summary(ignored),
                "the zero-weighted row contributes nothing to the objective", "r3"));
    double ignored_sum = 0.0;
    for (const double r : ignored.weighted_residuals) ignored_sum += r * r;
    CHECK(agree("r3/zero_weight/variance_factor", "variance_factor", "1", *ignored.covariance.variance_factor, ignored_sum / (5.0 - 2.0), 0.0, 1.0e-12,
                summary(ignored), "the variance factor divides by the five weighted rows minus the rank, not six", "r3"));
    // Rank deficiency: one observation cannot identify two coordinates.
    auto single = recovery_request(truth_case, 148.0, 0.0);
    single.observations = {truth_case.observations[0]};
    const auto deficient = fit(single);
    CHECK(truth("r3/rank_deficient", "one row and two coordinates withhold covariance as rank deficient",
                deficient.usable && deficient.covariance.status == "rank_deficient" && deficient.covariance.rank == 1
                    && deficient.covariance.optimizer.empty() && !deficient.covariance.variance_factor.has_value(), summary(deficient), "r3"));
    // Covariance transformation through a tie and scales, and the variance factor.
    auto shared = recovery_request(truth_case, 148.0, 0.0);
    shared.parameters.sigma_constant = {3.6, 3.6};
    shared.coordinates.coordinates = {{ParameterFamily::SigmaConstant, {{0}, {1}}, 3.6, 0.1, 3.0, 4.0},
                                      {ParameterFamily::DispersionEnergyOverK, {{0}}, 150.0, 10.0, 130.0, 170.0}};
    const auto tied = fit(shared);  // converged: an unfinished iterate carries no covariance (r2-j2/not_converged)
    CHECK(tied.usable && tied.covariance.status == "available" && tied.covariance.physical.size() == 9);
    const auto& cu = tied.covariance.optimizer;
    const double hand[9] = {0.01 * cu[0], 0.01 * cu[0], 1.0 * cu[1], 0.01 * cu[0], 0.01 * cu[0], 1.0 * cu[1],
                            1.0 * cu[2], 1.0 * cu[2], 100.0 * cu[3]};
    for (std::size_t k = 0; k < 9; ++k)
        CHECK(agree("r3/physical_covariance/" + std::to_string(k), "covariance", "native", tied.covariance.physical[k], hand[k], 0.0, 1.0e-12,
                    summary(tied), "Cov_theta = D Cov_u D^T with D = S diag(scale)", "r3"));
    // Cov_u is the inverse of J^T J of the independently qualified optimizer Jacobian.
    double worst = 0.0;
    for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t b = 0; b < 2; ++b)
            worst = std::max(worst, std::abs(cu[a * 2] * gram_entry(tied, 0, b) + cu[a * 2 + 1] * gram_entry(tied, 1, b) - (a == b ? 1.0 : 0.0)));
    CHECK(agree("r3/covariance_inverse", "max_abs(Cov_u J^T J - I)", "1", worst, 0.0, 1.0e-10, 0.0, summary(tied),
                "the optimizer covariance inverts the Gram matrix of the qualified Jacobian", "r3"));
    double weighted = 0.0;
    for (const double r : tied.weighted_residuals) weighted += r * r;
    CHECK(agree("r3/variance_factor", "variance_factor", "1", *tied.covariance.variance_factor, weighted / (5.0 - 2.0), 0.0, 1.0e-12,
                summary(tied), "the variance factor is sum(w r^2) / (n - rank)", "r3"));
    // A declared bound the model itself refuses (m = 0) is a typed trial failure, not an exception: an
    // unattainable target (pressures fifty times the truth) drives the optimizer onto that bound.
    FitRequest inadmissible;
    inadmissible.parameters = properties_support::neutral_binary();
    inadmissible.parameters.dispersion_energy_over_k[0] = 148.0;
    inadmissible.coordinates.coordinates = {{ParameterFamily::DispersionEnergyOverK, {{0}}, 150.0, 10.0, 140.0, 170.0},
                                            {ParameterFamily::SegmentCount, {{0}}, 1.0, 2.0, 0.0, 5.0}};
    for (std::size_t i = 0; i < 3; ++i) {
        inadmissible.observations.push_back(truth_case.observations[i]);
        inadmissible.observations.back().observed *= 50.0;
    }
    const auto refused = fit(inadmissible);
    const bool typed = std::any_of(refused.failures.begin(), refused.failures.end(), [&](const TrialFailure& item) {
        return item.status == ObservationStatus::InadmissibleTrial && item.observation == 3
            && item.message.find("segment_count") != std::string::npos; });
    // The optimizer may back off the refused bound and still finish or give up; either way the refusal is
    // retained by name and the result returns typed.
    CHECK(truth("r3/inadmissible_trial", "a trial the model refuses by name is retained as inadmissible_trial and the fit converges at the bound",
                typed && refused.trial_failures > 0 && refused.status == "converged" && refused.usable
                    && refused.covariance.status == "active_bound", summary(refused), "r3"));
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.regression.r1-f1", "Regression_R1_F1", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f1/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F1", "r1-f1/*:42"}, check_f1, "R1"},
        {{"native.regression.r1-f2", "Regression_R1_F2", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f2/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F2", "r1-f2/*:18"}, check_f2, "R1"},
        {{"native.regression.r1-f3", "Regression_R1_F3_F4", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f3/*,r1-f4/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F3", "r1-f3/*:5,r1-f4/*:3"}, check_f3_f4, "R1"},
        {{"native.regression.r1-f5", "Regression_R1_F5_F6", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f5/*,r1-f6/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F5", "r1-f5/*:5,r1-f6/*:3"}, check_f5_f6, "R1"},
        {{"native.regression.r1-f7", "Regression_R1_F7_F8", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f7/*,r1-f8/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F7", "r1-f7/*:5,r1-f8/*:12"}, check_f7_f8, "R1"},
        {{"native.regression.r1-f9", "Regression_R1_F9_F10", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f9/*,r1-f10/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F9", "r1-f9/*:4,r1-f10/*:13"}, check_f9_f10, "R1"},
        {{"native.regression.d86", "Regression_D86_density_actions", "Numerical", "Properties", "Explicit", "Numerical", nullptr, "D86", "d86/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "D86", "d86/*:25"}, check_d86, "R1"},
        {{"native.regression.r2-j1", "Regression_R2_J1", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "r2-j1/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R2-J1", "r2-j1/*:40"}, check_j1, "R2"},
        {{"native.regression.r2-j2", "Regression_R2_J2", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "r2-j2/*,r4/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R2-J2", "r2-j2/*:9"}, check_j2, "R2"},
        {{"native.regression.r3", "Regression_R3", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R3", "r3/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R3", "r3/*:22"}, check_r3, "R2"},
    });
}
