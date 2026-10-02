// Numerical falsifiers of the regression layer (#46 R1, #86 density actions, #47 R2, #48 R3)
// and the native rows the public R4 transport reproduces (#49).
#include "regression/fit.hpp"

#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_support.hpp"
#include "properties_support.hpp"
#include "recording_support.hpp"
#include "surface/parameter_builder.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
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

bool reference_scale_agrees(const std::string& id, const char* quantity, const char* unit, double model,
                            double reference, double atol, double rtol, const std::string& inputs,
                            const char* claim, const char* group) {
    const bool verdict = std::isfinite(model) && std::isfinite(reference)
        && std::abs(model - reference) <= atol + rtol * std::abs(reference);
    recording_support::compare(id.c_str(), quantity, unit, model, reference, inputs, "Identity", "direct_reference",
        claim, "reference_scale", atol, rtol, verdict, group, id.c_str());
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

void check_f2() {
    const auto parameters = properties_support::associating_binary();
    const auto set = seven_families();
    validate(set, parameters);
    const auto start = set.start(parameters);
    const auto physical = set.physical(start);
    const auto active = set.active_parameters();
    CHECK(physical.size() == 7 && active.size() == 7);
    for (std::size_t r = 0; r < 7; ++r) {
        const double truth_value = detail::family_values(parameters, active[r].family).first[active[r].indices.front()];
        const double ulp = 2.0 * EPS * std::max(std::abs(truth_value), std::abs(set.coordinates[r].origin));
        CHECK(agree("r1-f2/roundtrip/" + std::to_string(r), "theta", "native", physical[r], truth_value, ulp, 0.0,
                    fields({field("origin", set.coordinates[r].origin), field("scale", set.coordinates[r].scale),
                            field("u", start[r])}), "u -> theta -> u returns the start within 2 ulp of max(|theta|, |origin|)", "r1-f2"));
    }
    const auto applied = set.apply(parameters, start);
    bool identical = true;
    for (const auto& item : active) {
        const auto& lhs = detail::family_values(applied, item.family);
        const auto& rhs = detail::family_values(parameters, item.family);
        identical = identical && lhs.second == rhs.second && std::equal(lhs.first, lhs.first + lhs.second, rhs.first);
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
    request.observations = {overflow, overflow, overflow};
    request.weights = {0.0, 1.0, 1.0};
    const auto refused = fit(request);
    CHECK(truth("r1-f10/fit_nonfinite_typed", "a nonfinite residual at the start is a typed initial failure with its row",
                refused.status == "initial_evaluation_failed" && refused.failures.size() == 1 && refused.failures[0].observation == 1
                    && refused.failures[0].status == ObservationStatus::Nonfinite && !refused.usable && std::isfinite(refused.initial_cost)
                    && refused.statuses == std::vector<ObservationStatus>{ObservationStatus::NotEvaluated, ObservationStatus::Nonfinite, ObservationStatus::NotEvaluated},
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
    // The coordinate chain can overflow with finite physical gradients and a finite objective;
    // refusal precedes publication of the initial successful-point cache and warm starts.
    steep.observed = *kept.value - 1.0;
    steep.scale = 1.0e-150;
    request.observations = {steep};
    request.weights.clear();
    request.coordinates.coordinates = {{ParameterFamily::DispersionEnergyOverK, {{0}}, 191.42, 1.0e200, 150.0, 260.0}};
    const auto steep_fit = fit(request);
    CHECK(truth("r1-f10/jacobian_overflow_typed", "an overflowing coordinate chain is refused at its row before publishing the initial trial",
                steep_fit.status == "initial_evaluation_failed" && !steep_fit.usable && steep_fit.failures.size() == 1
                    && steep_fit.failures[0].status == ObservationStatus::Nonfinite && steep_fit.failures[0].observation == 0
                    && steep_fit.failures[0].message == "optimizer Jacobian not finite" && steep_fit.failures[0].jacobian,
                summary(steep_fit), "r1-f10"));
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
                       {ParameterFamily::SegmentCount, {{1}}, 1.6, 0.1, 1.0, 2.5},
                       {ParameterFamily::KijReciprocalTemperatureSlope, {{1, 2}}, 0.0, 10.0, -100.0, 100.0}};
    return set;
}

void check_j1() {
    // k_ij(T) = k_ij + b (1/T - 1/T_ref) with b = 0 at T_ref = 200 K: the model is unchanged, dk/db is not.
    auto parameters = properties_support::neutral_binary();
    parameters.k_ij_reference_temperature_k[1] = parameters.k_ij_reference_temperature_k[2] = 200.0;
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
        {"methane_partial_pressure", solved(ObservationKind::PartialPressure, bubble, 1, {0})},
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
        CHECK(gradient.status == ObservationStatus::Available && gradient.value.size() == 4);
        for (std::size_t c = 0; c < 4; ++c) {
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
    // Interleaved physical states, pooled species, fraction and both density bases keep their own rows.
    auto declared = tie.problem;
    for (std::size_t q = 0; q < declared.phases.size(); ++q) {
        declared.phases[q].composition_guess.assign(tie.start.mole_fractions.begin() + q * 2,
                                                    tie.start.mole_fractions.begin() + (q + 1) * 2);
        declared.phases[q].packing_guess = tie.start.packing_fractions[q];
    }
    FitRequest grouped;
    grouped.parameters = parameters;
    grouped.coordinates.coordinates = {set.coordinates.front()};
    grouped.observations = {solved(ObservationKind::PhaseFraction, declared, 0),
        solved(ObservationKind::EquilibriumPressure, bubble_at(190.0)),
        solved(ObservationKind::PhaseMoleFraction, declared, 0, {0, 1}),
        solved(ObservationKind::PartialPressure, declared, 1, {0}),
        solved(ObservationKind::PhaseDensity, declared, 0),
        solved(ObservationKind::PhaseDensity, declared, 0, {}, true),
        solved(ObservationKind::EquilibriumPressure, declared),
        solved(ObservationKind::PhaseFraction, declared, 1)};
    for (std::size_t i = 0; i < grouped.observations.size(); ++i) {
        auto& row = grouped.observations[i];
        row.observed = *predict(model, row).value;
        if (i % 2 == 0) { row.form = ResidualForm::Difference; row.scale = row.observed; }
        grouped.weights.push_back(1.0 + i);
    }
    grouped.controls.maximum_iterations = 1;
    grouped.controls.maximum_elapsed_time_seconds = 1.0e-6;
    const auto batch = fit(grouped);
    CHECK(batch.usable && batch.predictions.size() == grouped.observations.size());
    const auto point = grouped.coordinates.start(parameters);
    for (std::size_t i = 0; i < grouped.observations.size(); ++i) {
        const auto& row = grouped.observations[i];
        CHECK(agree("r2-j1/group_value/" + std::to_string(i), "prediction", "native", batch.predictions[i], row.observed,
                    1.0e-10 * row.observed, 1.0e-10, summary(batch),
                    "interleaved declarations and independent projections retain the scalar-owner central value", "r2-j1"));
        const auto value = [&](double shift) {
            auto u = point;
            u[0] += shift;
            const HelmholtzModel rebuilt(grouped.coordinates.apply(parameters, u));
            return *predict(rebuilt, row).value;
        };
        auto rungs = ladder(value, 0.005);
        rungs.coarse /= grouped.coordinates.coordinates.front().scale;
        rungs.fine /= grouped.coordinates.coordinates.front().scale;
        CHECK(ladder_agrees("r2-j1/group_action/" + std::to_string(i), "dy_dtheta", batch.physical_jacobian[i], rungs,
                            fields({field("value", batch.predictions[i]), field("u_step", 0.005)}),
                            "shared-state actions agree with independently rebuilt-model central differences", "r2-j1"));
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

// Held, Cameretti and Sadowski (2008), Table 1 / Eq. (7), packet v3 neutral water subset.
Parameters water_parameters() {
    auto p = properties_support::parameters(1);
    p.segment_count = {1.204659};
    p.sigma_constant = {2.7927};
    p.dispersion_energy_over_k = {353.9449};
    p.molar_mass_kg_per_mol = {0.01801528};
    p.sigma_term_component = {0, 0};
    p.sigma_term_amplitude = {10.11, -1.417};
    p.sigma_term_exponent = {-0.01775, -0.01146};
    p.site_component = {0, 0};
    p.site_multiplicity = {1, 1};
    p.association_active = {0, 1, 1, 0};
    p.association_energy_over_k = {0.0, 2425.6714, 2425.6714, 0.0};
    p.association_volume = {0.0, 0.0450989, 0.0450989, 0.0};
    return p;
}

ThermochemistryRecord constant_cp(std::size_t components) {
    ThermochemistryRecord record;
    record.standard_pressure = 1.0e5;
    record.component_intervals.assign(components, {{150.0, 500.0, true, true,
        {IdealPolynomial{{4.0 * GAS_CONSTANT_J_PER_MOL_K}, 300.0, 0.0, 0.0}, 1.0e5}}});
    return record;
}

CoordinateSet water_coordinates() {
    CoordinateSet set;
    set.coordinates = {
        {ParameterFamily::SegmentCount, {{0}}, 1.204659, 0.1, 0.8, 2.0},
        {ParameterFamily::SigmaConstant, {{0}}, 2.7927, 0.05, 2.4, 3.2},
        {ParameterFamily::DispersionEnergyOverK, {{0}}, 353.9449, 10.0, 200.0, 500.0},
        {ParameterFamily::AssociationEnergyOverK, {{1, 2}}, 2425.6714, 100.0, 1500.0, 3500.0},
        {ParameterFamily::AssociationVolume, {{1, 2}}, 0.0450989, 0.005, 0.005, 0.2}};
    return set;
}

void check_cp_derivatives() {
    struct Case { const char* name; Parameters parameters; double T, P; std::vector<double> x;
                  std::vector<ActiveParameter> active; std::vector<double> steps; };
    const auto water_active = water_coordinates().active_parameters();
    const std::vector<Case> cases = {
        {"a1", properties_support::ethane(), 200.0, 217800.6949469248, {1.0},
         {{ParameterFamily::SegmentCount, {0}}, {ParameterFamily::SigmaConstant, {0}},
          {ParameterFamily::DispersionEnergyOverK, {0}}}, {5e-4, 5e-4, 0.05}},
        {"a2-298", water_parameters(), 298.15, 1e5, {1.0}, water_active, {5e-4, 5e-4, 0.05, 1.0, 1e-5}},
        {"a2-400", water_parameters(), 400.0, 1e6, {1.0}, water_active, {5e-4, 5e-4, 0.05, 1.0, 1e-5}},
        {"a3", properties_support::associating_binary(), 313.15, 1e5, {0.25, 0.75},
         {{ParameterFamily::Kij, {1, 2}}, {ParameterFamily::SegmentCount, {0}},
          {ParameterFamily::AssociationEnergyOverK, {3, 12}},
          {ParameterFamily::AssociationVolume, {11, 14}}}, {1e-4, 5e-4, 1.0, 1e-5}}};
    for (const auto& c : cases) {
        const auto thermal = constant_cp(c.x.size());
        const HelmholtzModel model(c.parameters, c.active, thermal);
        auto item = density(c.T, c.P, c.x, DensityBranch::Liquid);
        item.kind = ObservationKind::IsobaricHeatCapacity;
        const auto value = predict(model, item);
        const auto gradient = differentiate(model, item, value);
        CHECK(value.status == ObservationStatus::Available && gradient.status == ObservationStatus::Available);
        const auto state = model.state_at_pressure(c.T, c.P, c.x, item.branch);
        const std::string id = std::string("cp-a/") + c.name;
        const auto inputs = fields({field("T_K", c.T), field("P_Pa", c.P), field("x", c.x)});
        CHECK(agree(id + "/ideal-value", "ideal_cp", "J/(mol K)",
            *value.value - *state.residual_isobaric_heat_capacity, 4.0 * GAS_CONSTANT_J_PER_MOL_K,
            1e-12 * std::max(1.0, *value.value), 0.0, inputs,
            "the total prediction contains the known constant ideal contribution", "cp-a"));
        double lift_ratio = 0.0;
        for (std::size_t k = 0; k < c.active.size(); ++k) {
            const auto& active = c.active[k];
            const auto closure = [&](double shift) {
                auto p = c.parameters;
                for (const auto index : active.indices) detail::family_values(p, active.family).first[index] += shift;
                return *predict(HelmholtzModel(p, {}, thermal), item).value;
            };
            const auto rungs = ladder(closure, c.steps[k]);
            const double richardson = (4.0 * rungs.fine - rungs.coarse) / 3.0;
            const auto column_inputs = fields({inputs, field("step", c.steps[k]), field("fd_signal", std::abs(rungs.fine)),
                                               field("richardson", richardson), field("fine_minus_coarse", rungs.fine - rungs.coarse),
                                               field("richardson_correction", richardson - rungs.fine),
                                               field("richardson_error", std::abs(gradient.value[k] - richardson))});
            const auto column = id + "/" + std::to_string(k);
            // Owner amendment: these two physical sensitivities lie below the native-unit signal floor.
            const bool weak = k == 2 && (std::string(c.name) == "a2-298" || std::string(c.name) == "a3");
            if (weak) {
                CHECK(agree(column + "/fine", "dcp_dtheta", "native", gradient.value[k], rungs.fine,
                            FIRST_ORDER * std::abs(rungs.fine), 0.0, column_inputs,
                            "the weak physical action agrees within the stricter relative criterion", "cp-a"));
                CHECK(agree(column + "/coarse", "dcp_dtheta", "native", gradient.value[k], rungs.coarse,
                            4.0 * FIRST_ORDER * std::abs(rungs.fine), 0.0, column_inputs,
                            "the coarse relative rung carries at most four times the criterion", "cp-a"));
                CHECK(truth(column + "/resolved", "finite-difference refinement resolves the relative criterion",
                    std::isfinite(rungs.fine) && rungs.fine != 0.0
                    && std::abs(rungs.fine - rungs.coarse) < FIRST_ORDER * std::abs(rungs.fine)
                    && std::abs(richardson - rungs.fine) < FIRST_ORDER * std::abs(rungs.fine), column_inputs, "cp-a"));
            } else {
                CHECK(ladder_agrees(column, "dcp_dtheta", gradient.value[k], rungs, column_inputs,
                                   "the fixed-pressure action agrees with rebuilt-model central differences", "cp-a"));
            }
            CHECK(agree(column + "/richardson", "dcp_dtheta", "native", gradient.value[k], richardson,
                        FIRST_ORDER * (weak ? std::abs(richardson) : std::max(1.0, std::abs(richardson))), 0.0, column_inputs,
                        "the action agrees with Richardson extrapolation", "cp-a"));
            CHECK(truth(column + "/signal", "the difference carries a nonzero parameter signal",
                        weak ? std::isfinite(rungs.fine) && rungs.fine != 0.0 : std::abs(rungs.fine) >= 1e-2,
                        column_inputs, "cp-a"));
            const ParameterActionDirection parameter{active.family, active.indices, 1.0};
            const auto pressure = model.pressure_actions(c.T, c.P, c.x, item.branch,
                {PropertyObservable::TotalIsobaricHeatCapacity, PropertyObservable::ResidualIsobaricHeatCapacity},
                {{0.0, 0.0, std::vector<double>(c.x.size(), 0.0), {parameter}}});
            const auto fixed = model.density_actions(c.T, state.molar_density, c.x,
                {PropertyObservable::TotalIsobaricHeatCapacity},
                {{0.0, 0.0, std::vector<double>(c.x.size(), 0.0), {parameter}}});
            CHECK(pressure.results[0].action && pressure.results[1].action && fixed.results[0].action);
            CHECK(agree(column + "/ideal-action", "dcp_dtheta", "native", *pressure.results[0].action,
                        *pressure.results[1].action, 0.0, 1e-12, column_inputs,
                        "fixed ideal thermochemistry carries no parameter action", "cp-a"));
            lift_ratio = std::max(lift_ratio, std::abs(gradient.value[k] - *fixed.results[0].action)
                / (FIRST_ORDER * std::max(1.0, std::abs(rungs.fine))));
        }
        CHECK(truth(id + "/density-lift", "at least one omitted density lift exceeds 1000 criteria",
                    lift_ratio > 1e3, fields({inputs, field("lift_criterion_ratio", lift_ratio)}), "cp-a"));
    }
    auto item = density(200.0, 217800.6949469248, {1.0}, DensityBranch::Liquid);
    item.kind = ObservationKind::IsobaricHeatCapacity;
    const auto missing = predict(HelmholtzModel(properties_support::ethane()), item);
    CHECK(truth("cp-a/missing", "missing thermochemistry is typed and never a value",
        missing.status == ObservationStatus::ThermochemistryUnavailable && !missing.value && missing.message == "missing_record",
        field("message", missing.message), "cp-a"));
    item.pressure = 1e10;
    const auto absent = predict(HelmholtzModel(properties_support::ethane(), {}, constant_cp(1)), item);
    CHECK(truth("cp-a/root", "an absent density root retains its typed cause",
        absent.status == ObservationStatus::DensityUnavailable && !absent.value && !absent.message.empty(),
        field("message", absent.message), "cp-a"));
}

void check_cp_recovery() {
    FitRequest request;
    request.parameters = water_parameters();
    request.coordinates = water_coordinates();
    request.thermochemistry = constant_cp(1);
    const HelmholtzModel known(request.parameters, {}, request.thermochemistry);
    for (double T : {300.0, 340.0, 380.0, 420.0}) {
        auto rho = density(T, 1e6, {1.0}, DensityBranch::Liquid);
        auto cp = rho;
        cp.kind = ObservationKind::IsobaricHeatCapacity;
        auto sat = solved(ObservationKind::EquilibriumPressure, pure_saturation(T, 1e5));
        for (auto item : {rho, cp, sat}) {
            const auto prediction = predict(known, item);
            CHECK(prediction.status == ObservationStatus::Available);
            item.observed = *prediction.value;
            request.observations.push_back(item);
        }
    }
    for (const auto& active : request.coordinates.active_parameters())
        for (const auto index : active.indices)
            detail::family_values(request.parameters, active.family).first[index] *=
                active.family == ParameterFamily::AssociationVolume ? 1.05 : 1.02;
    const auto result = fit(request);
    const auto inputs = fields({summary(result), field("singular_values", result.covariance.singular_values),
                               field("rank", result.covariance.rank)});
    CHECK(truth("cp-b/converged", "five-coordinate recovery converges with no trial failures",
        result.status == "converged" && result.usable && result.trial_failures == 0, inputs, "cp-b"));
    for (std::size_t k = 0; k < 5; ++k)
        CHECK(agree("cp-b/value/" + std::to_string(k), "parameter", "native", result.physical[k],
            request.coordinates.coordinates[k].origin, 0.0, 1e-6, inputs,
            "the known physical coordinate is recovered", "cp-b"));
    CHECK(agree("cp-b/cost", "cost", "1", result.final_cost, 0.0, 1e-14, 0.0, inputs,
                "the manufactured residual cost vanishes", "cp-b"));
    CHECK(truth("cp-b/rank", "the conditional covariance is available with rank five",
                result.covariance.status == "available" && result.covariance.rank == 5, inputs, "cp-b"));
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
    // The same manufactured targets, starting on the upper epsilon bound with descent initially
    // inward. Coupling to k_ij turns that gradient outward without a new bound contact (#193).
    auto bound_start = recovery_request(truth_case, 151.0, -0.1, 151.0);
    bound_start.controls.maximum_iterations = 12;  // shared work limit; no iteration-count pass criterion
    const auto initial_epsilon_gradient = [](const FitRequest& start) {
        const HelmholtzModel model(start.parameters, start.coordinates.active_parameters());
        double gradient = 0.0;
        for (std::size_t row = 0; row < start.observations.size(); ++row) {
            const auto& observation = start.observations[row];
            const auto prediction = predict(model, observation);
            const auto action = differentiate(model, observation, prediction);
            if (!(prediction.value && prediction.residual && action.status == ObservationStatus::Available))
                throw std::runtime_error("initial epsilon gradient row " + std::to_string(row) + ": prediction "
                    + observation_status_name(prediction.status) + " " + prediction.message + "; action "
                    + observation_status_name(action.status) + " " + action.message);
            gradient += *prediction.residual * residual_slope(observation, *prediction.value)
                        * action.value[0] * start.coordinates.coordinates[0].scale;
        }
        return gradient;
    };
    CHECK(initial_epsilon_gradient(bound_start) > 0.0);  // descent from the upper bound is initially inward
    // Removing epsilon from the fit gives the independent scalar minimization at epsilon = 151 K.
    auto reduced = bound_start;
    reduced.coordinates.coordinates.erase(reduced.coordinates.coordinates.begin());
    const auto reference = fit(reduced);
    const auto refreshed = fit(bound_start);
    CHECK(reference.usable && refreshed.usable && refreshed.optimizer.size() == 2);
    CHECK(agree("r2-j2/bound_start/cost", "cost", "1", refreshed.final_cost, reference.final_cost,
                1.0e-12, 1.0e-8, summary(refreshed),
                "the coupled fit reaches the independently reduced bounded minimum within the work limit", "r2-j2"));
    CHECK(agree("r2-j2/bound_start/kij", "k_ij", "1", refreshed.physical[1], reference.physical[0],
                1.0e-7, 0.0, summary(refreshed), "the free coordinate agrees with the scalar minimum", "r2-j2"));
    double projected = 0.0;
    const auto [lower, upper] = bound_start.coordinates.optimizer_bounds();
    for (std::size_t c = 0; c < 2; ++c) {
        double gradient = 0.0;
        for (std::size_t i = 0; i < refreshed.weighted_residuals.size(); ++i)
            gradient += refreshed.optimizer_jacobian[i * 2 + c] * refreshed.weighted_residuals[i];
        CHECK(std::isfinite(gradient));
        projected = std::max(projected, std::abs(refreshed.optimizer[c]
            - std::clamp(refreshed.optimizer[c] - gradient, lower[c], upper[c])));
    }
    CHECK(truth("r2-j2/bound_start/feasible", "the accepted optimum is feasible with epsilon at its upper bound",
                bound_start.coordinates.within_bounds(refreshed.optimizer) && refreshed.physical[0] == 151.0,
                summary(refreshed), "r2-j2"));
    CHECK(agree("r2-j2/bound_start/projected_gradient", "projected_gradient_max", "1", projected, 0.0,
                1.0e-6, 0.0, summary(refreshed), "the full weighted gradient satisfies bounded stationarity", "r2-j2"));
    // Initially held epsilon must be released as k_ij moves and the gradient turns inward.
    auto interior_truth = truth_case;
    interior_truth.parameters.dispersion_energy_over_k[0] = 150.0;
    const HelmholtzModel interior_model(interior_truth.parameters);
    for (auto& observation : interior_truth.observations) {
        const auto prediction = predict(interior_model, observation);
        if (!prediction.value)
            throw std::runtime_error("interior truth row " + std::to_string(&observation - interior_truth.observations.data())
                + ": " + observation_status_name(prediction.status) + " " + prediction.message);
        observation.observed = *prediction.value;
        // The default cold seed coalesces at the release start (+0.1); these measured two-phase
        // guesses admit the same T/feed/bounds/parameter start, without a topology-event bypass.
        if (prediction.solved) {
            const auto& state = *prediction.solved;
            for (std::size_t q = 0; q < observation.problem.phases.size(); ++q) {
                observation.problem.phases[q].composition_guess.assign(state.mole_fractions.begin() + q * 2,
                                                                        state.mole_fractions.begin() + (q + 1) * 2);
                observation.problem.phases[q].packing_guess = state.packing_fractions[q];
            }
            observation.problem.pressure = FreeValue::variable(state.pressure.front());
        }
    }
    auto inward = recovery_request(interior_truth, 151.0, 0.1, 151.0);
    CHECK(initial_epsilon_gradient(inward) < 0.0);  // outward descent at the upper bound: initially held
    const auto released = fit(inward);
    CHECK(released.usable);
    CHECK(agree("r2-j2/bound_start/inward", "dispersion_energy_over_k", "K", released.physical[0], 150.0,
                0.0, 1.0e-6, summary(released), "inward descent releases the held coordinate and recovers the interior truth", "r2-j2"));
    CHECK(agree("r2-j2/bound_start/inward_cost", "cost", "1", released.final_cost, 0.0,
                1.0e-14, 0.0, summary(released), "the released interior fit recovers the manufactured targets", "r2-j2"));
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

// #187: source-basis reaction constants through the existing solved-state first actions.
Coordinate reaction_coordinate(const char* name, ReactionCoefficient coefficient = ReactionCoefficient::A,
                               double scale = 1.0) {
    Coordinate result;
    result.reaction = ReactionCoordinate{name, coefficient};
    result.scale = scale;
    return result;
}

Observation amount_ratio(const EquilibriumProblem& problem, std::size_t numerator, std::size_t denominator,
                         std::vector<std::size_t> species, std::vector<std::size_t> denominator_species) {
    auto result = solved(ObservationKind::PhaseAmountRatio, problem, numerator, species);
    result.denominator_phase = denominator;
    result.denominator_species = std::move(denominator_species);
    return result;
}

// A state is used only after the same owner's full-rank, residual and topology gates admit a nonzero action.
void regular_reaction_state(const HelmholtzModel& model, const Observation& observation, const std::string& id) {
    const auto prediction = predict(model, observation);
    CHECK(prediction.status == ObservationStatus::Available);
    equilibrium::SolvedStateActionDirection direction;
    direction.feed_amounts_mol.assign(model.parameters().component_count(), 0.0);
    direction.active_parameters.assign(model.active_parameters().size(), 0.0);
    direction.reaction_log_constants.assign(observation.problem.reactions.size(), 0.0);
    direction.reaction_log_constants[0] = 1.0;
    equilibrium::SolvedStateObservable observable;
    observable.kind = equilibrium::SolvedStateObservableKind::PhaseMoleFraction;
    observable.component = observation.problem.phases[0].support_declared ? observation.problem.phases[0].support.front() : 0;
    const auto batch = equilibrium::solved_state_first_actions(*prediction.compiled, *prediction.solved, {observable}, {direction}).front();
    const auto& result = batch.results.front();
    CHECK(truth(id + "/regular", "the central state admits a finite nonzero first reaction action",
        result.status == equilibrium::SolvedStateActionStatus::Available && result.action && std::abs(*result.action) > 0.0
        && prediction.solved->topology_event == equilibrium::TopologyEvent::None,
        fields({field("central_residual", batch.diagnostics.central_residual.value_or(std::numeric_limits<double>::quiet_NaN())),
                field("dx_dlnK", result.action.value_or(std::numeric_limits<double>::quiet_NaN())),
                field("condition", batch.diagnostics.condition_estimate.value_or(std::numeric_limits<double>::quiet_NaN()))}), "reaction-fit"));
}

void reaction_data(FitRequest& request) {
    const HelmholtzModel model(request.parameters, {}, request.thermochemistry);
    for (auto& row : request.observations) {
        const auto prediction = predict(model, row);
        CHECK(prediction.status == ObservationStatus::Available);
        row.observed = *prediction.value;
    }
}

FitResult reaction_columns(const std::string& id, FitRequest request) {
    reaction_data(request);  // zero residuals keep the returned Jacobian at this exact centre
    const auto result = fit(request);
    CHECK(result.usable && result.trial_failures == 0 && result.status == "converged");
    const auto u0 = request.coordinates.start(request.parameters, request.observations);
    const auto columns = u0.size();
    for (std::size_t row = 0; row < request.observations.size(); ++row) {
        for (std::size_t c = 0; c < columns; ++c) {
            const auto value = [&](double shift) {
                auto u = u0;
                u[c] += shift;
                auto observations = request.observations;
                request.coordinates.apply(observations, u);
                const HelmholtzModel rebuilt(request.coordinates.apply(request.parameters, u));
                const auto prediction = predict(rebuilt, observations[row]);
                CHECK(prediction.status == ObservationStatus::Available);
                return std::log(*prediction.value);
            };
            const auto rungs = ladder(value, 0.005);
            const auto key = id + "/" + std::to_string(row) + "/" + std::to_string(c);
            const auto inputs = fields({field("fd_signal", std::abs(rungs.fine)), field("u_step", 0.005),
                                      field("prediction", result.predictions[row])});
            CHECK(truth(key + "/signal", "the dimensionless finite-difference column has signal at least 1e-2",
                        std::abs(rungs.fine) >= 1.0e-2, inputs, "reaction-fit"));
            CHECK(ladder_agrees(key, "dlny_du", result.optimizer_jacobian[row * columns + c], rungs, inputs,
                               "the exact log-residual column agrees with independently rebuilt central differences", "reaction-fit"));
        }
    }
    return result;
}

EquilibriumProblem reacting_bubble(const HelmholtzModel& model) {
    auto problem = bubble_at(180.0);
    const auto center = equilibrium::solve_equilibrium(model, problem);
    CHECK(center.success && center.topology_event == equilibrium::TopologyEvent::None);
    const auto compiled = equilibrium::compile_problem(model, problem);
    equilibrium::SolvedStateObservable first, second;
    first.kind = second.kind = equilibrium::SolvedStateObservableKind::PhaseLogActivity;
    second.component = 1;
    const auto batch = equilibrium::solved_state_actions(compiled, center,
        {{first, second}, {{0.0, 0.0, {0.0, 0.0}, {}}}});
    CHECK(batch.results[0].value && batch.results[1].value);
    const double L = *batch.results[1].value - *batch.results[0].value;
    problem.reactions = {{"iso", {-1.0, 1.0}, std::nullopt, std::nullopt,
        ReactionLogPolynomial{L - 300.0 / 180.0 - 0.5 * std::log(180.0 / 200.0) - 0.001 * 180.0,
                            300.0, 0.5, 0.001, 200.0, 150.0, 220.0, equilibrium::EOS_STANDARD_STATE_ID}}};
    const auto reacting = equilibrium::solve_equilibrium(model, problem);
    CHECK(reacting.success && reacting.topology_event == equilibrium::TopologyEvent::None);
    CHECK(agree("reaction-a3/pinned_pressure", "P", "Pa", reacting.pressure.front(), center.pressure.front(), 0.0, 1.0e-10,
                "T=180 K", "the reacting manufactured bubble retains the original state", "reaction-fit"));
    for (std::size_t i = 0; i < center.mole_fractions.size(); ++i)
        CHECK(agree("reaction-a3/pinned_x/" + std::to_string(i), "x", "1", reacting.mole_fractions[i], center.mole_fractions[i],
                    1.0e-10, 0.0, "T=180 K", "the reacting manufactured bubble retains the original compositions", "reaction-fit"));
    return problem;
}

EquilibriumProblem referenced_formula(std::optional<double> source_pressure) {
    auto problem = equilibrium_support::flash(3);
    problem.temperature = FreeValue::fixed(298.15);
    problem.feed.values = {1.0, 1.0e-6, 1.0e-6};
    problem.source_reference = NeutralReferenceDeclaration{0, {0.0, 1.0, 1.0}, 1.0};
    problem.reactions = {{"H2O=Na++Cl-", {-1.0, 1.0, 1.0}, std::nullopt, std::nullopt,
        ReactionLogPolynomial{std::log(1.0e-8), 0.0, 0.0, 0.0, 298.15, 250.0, 350.0, "analytic-source"},
        ReactionReference{source_pressure, ReferenceSourceBasis::RawMoleFractionInfiniteDilution}}};
    return problem;
}

struct AnalyticMinimaxProblem final { FitRequest request; std::vector<CheckResult> frozen; };

HydrationHold analytic_minimax_hold(double upper) {
    HydrationHold hold;
    hold.reference.kind = ObservationKind::SolvationGibbsEnergy;
    hold.reference.temperature = 298.15; hold.reference.pressure = 1.0e5;
    hold.reference.composition = {1.0}; hold.reference.ion_indices = {0, 1};
    hold.reference.stoichiometry = {1, 1}; hold.reference.branch = DensityBranch::Liquid;
    hold.reference.species = {0}; hold.reference.form = ResidualForm::Difference;
    hold.lower = 0.0; hold.upper = upper;
    hold.energy_tolerance = 1.0e-6; hold.diameter_tolerance = 1.0e-10; hold.maximum_iterations = 100;
    return hold;
}

AnalyticMinimaxProblem analytic_minimax_problem(bool b1, double q, double sigma) {
    AnalyticMinimaxProblem problem;
    auto& request = problem.request;
    request.coordinates.coordinates.push_back({ParameterFamily::Kij, {{1, 2}}, 0.0, 1.0, 0.0, b1 ? 2.0 : 3.0});
    if (b1) request.coordinates.coordinates.push_back({ParameterFamily::SegmentCount, {{0}}, 0.0, 1.0, 0.0, 2.0});
    request.observations.resize(5);
    for (auto& row : request.observations) {
        row.kind = ObservationKind::MeanIonicActivityCoefficient;
        row.form = ResidualForm::LogRatio; row.observed = 1.0; row.scale = 1.0;
    }
    request.weights.assign(5, 1.0);
    request.holds = {analytic_minimax_hold(b1 ? 4.0 : 6.0)};
    request.controls.maximum_iterations = 100; request.controls.maximum_elapsed_time_seconds = 5.0;
    request.minimax = MinimaxRequest{};
    request.minimax->targets = b1 ? std::vector<std::size_t>{0, 1} : std::vector<std::size_t>{0, 1, 2, 3, 4};
    request.minimax->group_scales = b1 ? std::vector<double>{sigma, sigma} : std::vector<double>{};
    request.minimax->q = q; request.minimax->tolerance = 1.0e-10;
    request.minimax->activity_distance = 1.0e-6; request.minimax->margin = 0.1;
    if (b1) {
        CheckGroup maximum; maximum.name = "maximum"; maximum.rows = {2};
        maximum.metric = CheckMetric::AbsoluteLog10Ratio; maximum.maximum = 0.0075;
        request.checks.push_back(maximum);
        CheckResult frozen_maximum; frozen_maximum.name = maximum.name; frozen_maximum.unit = "decade";
        frozen_maximum.status = "pass"; frozen_maximum.rows = maximum.rows;
        frozen_maximum.maximum_threshold = maximum.maximum; problem.frozen.push_back(frozen_maximum);
        CheckGroup mean; mean.name = "mean"; mean.rows = {3, 4};
        mean.metric = CheckMetric::AbsoluteLog10Ratio; mean.mean = 0.0075;
        request.checks.push_back(mean);
        CheckResult frozen_mean; frozen_mean.name = mean.name; frozen_mean.unit = "decade";
        frozen_mean.status = "pass"; frozen_mean.rows = mean.rows;
        frozen_mean.mean_threshold = mean.mean; problem.frozen.push_back(frozen_mean);
    }
    return problem;
}

ObservationFunction analytic_minimax_evaluator(const FitRequest& request) {
    return [&request](const std::vector<double>& u, bool, ObservationEvaluation& out) {
        constexpr double a = 0.01;
        const double ln10 = std::log(10.0);
        const bool b1 = u.size() == 2;
        const std::size_t physical = b1 ? 3 : 2, n = request.observations.size();
        const double x = u[0], z = b1 ? u[1] : 0.0, d = 2.0 * x;
        out = ObservationEvaluation{};
        out.predictions.resize(n); out.residuals.resize(n); out.weighted.resize(n);
        out.statuses.assign(n, ObservationStatus::Available);
        out.gradient.assign(n * physical, 0.0); out.optimizer_gradient.assign(n * u.size(), 0.0);
        out.chart.parameters = request.parameters; out.chart.active = request.coordinates.active_parameters();
        out.chart.active.push_back({ParameterFamily::BornDiameter, {0}});
        out.chart.physical = b1 ? std::vector<double>{x, z, d} : std::vector<double>{x, d};
        out.chart.map = b1 ? std::vector<double>{1.0, 0.0, 0.0, 1.0, 2.0, 0.0}
                           : std::vector<double>{1.0, 2.0};
        HoldResult hold; hold.derivative = 1.0;
        hold.parameter_derivative = b1 ? std::vector<double>{-2.0, 0.0, 1.0} : std::vector<double>{-2.0, 1.0};
        hold.cause = "converged"; hold.iterations = 1; out.chart.holds = {hold};
        for (std::size_t i = 0; i < n; ++i) {
            double e = 0.0; std::vector<double> de(physical, 0.0);
            if (b1 && i == 0) { e = a * (x + z); de = {a, a, 0.0}; }
            else if (b1 && i == 1) { e = a * (3.0 - d / 2.0 - z); de = {0.0, -a, -a / 2.0}; }
            else if (b1 && i == 2) { e = a * x; de = {a, 0.0, 0.0}; }
            else if (b1 && i == 3) { e = a * z; de = {0.0, a, 0.0}; }
            else if (b1) { e = a * (z + 0.5); de = {0.0, a, 0.0}; }
            else if (i < 4) { e = a * x; de = {a, 0.0}; }
            else { e = a * (d / 2.0 - 3.0); de = {0.0, a / 2.0}; }
            const double y = std::pow(10.0, e), weight = request.weights.empty() ? 1.0 : request.weights[i];
            out.predictions[i] = y; out.residuals[i] = std::log(y / request.observations[i].observed);
            out.weighted[i] = std::sqrt(weight) * out.residuals[i]; out.cost += 0.5 * out.weighted[i] * out.weighted[i];
            for (std::size_t p = 0; p < physical; ++p) {
                out.gradient[i * physical + p] = ln10 * y * de[p];
                for (std::size_t c = 0; c < u.size(); ++c)
                    out.optimizer_gradient[i * u.size() + c] += out.gradient[i * physical + p] * out.chart.map[p * u.size() + c];
            }
        }
        return true;
    };
}

void check_minimax_kkt(const std::string& id, const MinimaxResult& m, const std::string& inputs) {
    CHECK(agree(id + "/primal", "scaled_primal_norm", "1", m.primal, 0.0, 1.0e-8, 0.0,
        inputs, "scaled primal feasibility is within tolerance", "r2-minimax"));
    CHECK(agree(id + "/stationarity", "scaled_stationarity_norm", "1", m.stationarity, 0.0, 1.0e-8, 0.0,
        inputs, "scaled stationarity is within tolerance", "r2-minimax"));
    CHECK(agree(id + "/dual_sign", "scaled_dual_sign_norm", "1", m.dual_sign, 0.0, 1.0e-8, 0.0,
        inputs, "scaled dual signs are feasible", "r2-minimax"));
    CHECK(agree(id + "/complementarity", "scaled_complementarity_norm", "1", m.complementarity, 0.0,
        1.0e-8, 0.0, inputs, "scaled complementarity is within tolerance", "r2-minimax"));
}

void check_b1_minimax(const std::string& id, const FitRequest& request, const FitResult& result,
                      const std::vector<double>& start) {
    const auto inputs = fields({field("q", request.minimax->q), field("sigma", request.minimax->group_scales[0]),
        field("start_x", start[0]), field("start_z", start[1]), field("a_per_decade", 0.01)});
    CHECK(result.minimax.has_value());
    const auto& m = *result.minimax;
    CHECK(truth(id + "/available", "the analytic constrained candidate and KKT diagnostics are available",
        result.usable && m.availability == "available" && m.feasible && m.stationary, inputs, "r2-minimax"));
    const std::vector<double> coordinates = {0.75, 0.5, 1.5}, errors = {0.0125, 0.0175, 0.0075, 0.005, 0.010};
    for (std::size_t i = 0; i < coordinates.size(); ++i)
        CHECK(agree(id + "/coordinate/" + std::to_string(i), "physical_coordinate", "1", result.chart.physical[i],
            coordinates[i], 1.0e-7, 0.0, inputs, "the unique guarded minimax point is recovered", "r2-minimax"));
    CHECK(agree(id + "/H", "H", "1", result.chart.physical[2] - 2.0 * result.chart.physical[0], 0.0,
        1.0e-12, 0.0, inputs, "the held equality is eliminated", "r2-minimax"));
    for (std::size_t i = 0; i < errors.size(); ++i)
        CHECK(agree(id + "/error/" + std::to_string(i), "log10_prediction_ratio", "decade",
            check_error(CheckMetric::AbsoluteLog10Ratio, request.observations[i], result.predictions[i]).first,
            errors[i], 1.0e-9, 0.0, inputs, "the analytic observation error is recovered", "r2-minimax"));
    CHECK(agree(id + "/t", "minimax_t", "decade", m.t, 0.0175, 1.0e-9, 0.0,
        inputs, "the active target sets the epigraph", "r2-minimax"));
    CHECK(agree(id + "/target_maximum", "target_maximum", "decade", m.target_maximum, 0.0175, 1.0e-9, 0.0,
        inputs, "the reported maximum matches the active target", "r2-minimax"));
    CHECK(agree(id + "/gap", "epigraph_gap", "decade", m.gap, 0.0, 1.0e-9, 0.0,
        inputs, "the epigraph is tight", "r2-minimax"));
    CHECK(m.rows.size() == m.physical_values.size() && m.rows.size() == m.physical_multipliers.size());
    double mu_e2 = 0.0, mu_max = 0.0, mu_e4 = 0.0, mu_e5 = 0.0;
    for (std::size_t i = 0; i < m.rows.size(); ++i) {
        const auto& r = m.rows[i]; double value = 0.0, multiplier = 0.0;
        if (r.kind == "target") {
            value = r.sign * errors[r.row] - m.t;
            if (r.row == 1 && r.sign == 1.0) { multiplier = 1.0; mu_e2 = m.physical_multipliers[i]; }
        } else if (r.kind == "maximum") {
            value = r.sign * errors[r.row] - 0.0075;
            if (r.group == 0 && r.sign == 1.0) { multiplier = 1.0; mu_max = m.physical_multipliers[i]; }
        } else if (r.kind == "link") {
            const double slack = errors[r.row]; value = r.sign * errors[r.row] - slack;
            if (r.sign == 1.0) {
                multiplier = 0.5;
                if (r.row == 3) mu_e4 = m.physical_multipliers[i]; else mu_e5 = m.physical_multipliers[i];
                CHECK(r.auxiliary.has_value());
                CHECK(agree(id + "/slack/" + std::to_string(r.row), "mean_slack", "decade",
                    m.point.at(*r.auxiliary), slack, 1.0e-9, 0.0, inputs,
                    "each mean slack equals its absolute row error", "r2-minimax"));
            }
        } else { value = 0.0; multiplier = 1.0; }
        CHECK(agree(id + "/row/" + std::to_string(i), "physical_constraint", "decade", m.physical_values[i],
            value, 1.0e-9, 0.0, inputs, "the physical target and guard rows match the analytic values", "r2-minimax"));
        CHECK(agree(id + "/multiplier/" + std::to_string(i), "physical_row_multiplier", "1", m.physical_multipliers[i],
            multiplier, 1.0e-6, 0.0, inputs, "physical row multipliers match the analytic active set", "r2-minimax"));
    }
    for (const double z : m.physical_z_lower)
        CHECK(agree(id + "/lower_multiplier", "physical_bound_multiplier", "1", z, 0.0, 1.0e-6, 0.0,
            inputs, "the interior point has no lower-bound multiplier", "r2-minimax"));
    for (const double z : m.physical_z_upper)
        CHECK(agree(id + "/upper_multiplier", "physical_bound_multiplier", "1", z, 0.0, 1.0e-6, 0.0,
            inputs, "the interior point has no upper-bound multiplier", "r2-minimax"));
    CHECK(agree(id + "/nu", "equality_hold_multiplier", "decade_per_H", m.hold_multipliers.at(0), 0.005,
        1.0e-6, 0.0, inputs, "the eliminated equality reconstructs nu=a/2", "r2-minimax"));
    CHECK(m.full_stationarity.size() == 3);
    for (std::size_t i = 0; i < m.full_stationarity.size(); ++i)
        CHECK(agree(id + "/full_stationarity/" + std::to_string(i), "full_space_stationarity", "decade_per_coordinate",
            m.full_stationarity[i], 0.0, 1.0e-8, 0.0, inputs, "full-space stationarity vanishes", "r2-minimax"));
    CHECK(agree(id + "/stationarity_x", "stationarity_x", "decade_per_x", 0.01 * mu_max - 2.0 * m.hold_multipliers[0],
        0.0, 1.0e-8, 0.0, inputs, "a*mu_max-2*nu=0", "r2-minimax"));
    CHECK(agree(id + "/stationarity_z", "stationarity_z", "decade_per_z", -0.01 * mu_e2 + 0.01 * (mu_e4 + mu_e5),
        0.0, 1.0e-8, 0.0, inputs, "-a*mu_e2+a*(mu_e4+mu_e5)=0", "r2-minimax"));
    CHECK(agree(id + "/stationarity_d", "stationarity_d", "decade_per_d", -0.01 * mu_e2 / 2.0 + m.hold_multipliers[0],
        0.0, 1.0e-8, 0.0, inputs, "-a*mu_e2/2+nu=0", "r2-minimax"));
    CHECK(result.chart.map == std::vector<double>({1.0, 0.0, 0.0, 1.0, 2.0, 0.0}));
    CHECK(result.chart.holds[0].parameter_derivative == std::vector<double>({-2.0, 0.0, 1.0}));
    CHECK(agree(id + "/H_d", "H_d", "1", result.chart.holds[0].derivative, 1.0, 0.0, 0.0,
        inputs, "the analytic hold has H_d=1", "r2-minimax"));
    CHECK(agree(id + "/hold_residual", "H", "1", result.chart.holds[0].energy_error, 0.0, 1.0e-12, 0.0,
        inputs, "the analytic equality has zero residual", "r2-minimax"));
    check_minimax_kkt(id, m, inputs);
}

void check_d7_minimax() {
    for (const auto& start : std::vector<std::vector<double>>{{0.25, 0.25}, {1.0, 1.0}}) {
        const auto problem = analytic_minimax_problem(true, 0.02, 0.02);
        const auto evaluate = analytic_minimax_evaluator(problem.request);
        const auto result = solve_minimax(problem.request, start, problem.frozen, evaluate);
        check_b1_minimax("r2-d7-b1/start-" + number(start[0]), problem.request, result, start);
    }
    {
        const std::vector<double> start = {0.25, 0.25};
        const auto problem = analytic_minimax_problem(true, 0.01, 0.04);
        const auto evaluate = analytic_minimax_evaluator(problem.request);
        const auto result = solve_minimax(problem.request, start, problem.frozen, evaluate);
        check_b1_minimax("r2-d7-b1/rescaled", problem.request, result, start);
    }
    const auto problem = analytic_minimax_problem(false, 0.02, 0.02);
    const auto evaluate = analytic_minimax_evaluator(problem.request);
    ObservationEvaluation ls;
    const std::vector<double> ls_point = {3.0 / 5.0};
    CHECK(evaluate(ls_point, true, ls));
    double ls_max = 0.0, ls_gradient = 0.0;
    for (std::size_t i = 0; i < problem.request.observations.size(); ++i) {
        const auto [error, slope] = check_error(CheckMetric::AbsoluteLog10Ratio, problem.request.observations[i], ls.predictions[i]);
        ls_max = std::max(ls_max, std::abs(error));
        ls_gradient += error * slope * ls.optimizer_gradient[i];
    }
    const auto ls_inputs = fields({field("a_per_decade", 0.01), field("equal_LogRatio_weights", 1.0)});
    CHECK(agree("r2-d7-b2/least_squares_x", "x", "1", ls.chart.physical[0], 0.6, 1.0e-12, 0.0,
        ls_inputs, "four ax rows and one a(x-3) row have their analytic least-squares minimum", "r2-minimax"));
    CHECK(agree("r2-d7-b2/least_squares_d", "d", "1", ls.chart.physical[1], 1.2, 1.0e-12, 0.0,
        ls_inputs, "the least-squares point obeys d=2x", "r2-minimax"));
    CHECK(agree("r2-d7-b2/least_squares_gradient", "least_squares_gradient", "decade_squared_per_x", ls_gradient, 0.0,
        1.0e-14, 0.0, ls_inputs, "4x+(x-3)=0 is the equal-weight least-squares condition", "r2-minimax"));
    CHECK(agree("r2-d7-b2/least_squares_maximum", "maximum_target_error", "decade", ls_max, 0.024, 1.0e-12, 0.0,
        ls_inputs, "the least-squares maximum is 0.024 decade", "r2-minimax"));
    CHECK(truth("r2-d7-b2/least_squares_guard", "the least-squares point misses the 0.02 decade maximum guard",
        ls_max > 0.02, fields({ls_inputs, field("guard_decade", 0.02)}), "r2-minimax"));
    for (const double x : {0.5, 2.5}) {
        const std::vector<double> start = {x};
        const auto result = solve_minimax(problem.request, start, problem.frozen, evaluate);
        const std::string id = "r2-d7-b2/start-" + number(x);
        const auto inputs = fields({field("q", problem.request.minimax->q), field("start_x", x), field("a_per_decade", 0.01)});
        CHECK(result.minimax.has_value());
        const auto& m = *result.minimax;
        CHECK(truth(id + "/available", "the analytic minimax candidate and KKT diagnostics are available",
            result.usable && m.availability == "available" && m.feasible && m.stationary, inputs, "r2-minimax"));
        CHECK(agree(id + "/x", "x", "1", result.chart.physical[0], 1.5, 1.0e-7, 0.0,
            inputs, "the opposing maximum errors balance at x=1.5", "r2-minimax"));
        CHECK(agree(id + "/d", "d", "1", result.chart.physical[1], 3.0, 1.0e-7, 0.0,
            inputs, "the eliminated equality gives d=2x", "r2-minimax"));
        CHECK(agree(id + "/H", "H", "1", result.chart.physical[1] - 2.0 * result.chart.physical[0], 0.0,
            1.0e-12, 0.0, inputs, "the held equality is eliminated", "r2-minimax"));
        CHECK(agree(id + "/t", "minimax_t", "decade", m.t, 0.015, 1.0e-9, 0.0,
            inputs, "the epigraph is the maximum target error", "r2-minimax"));
        CHECK(agree(id + "/target_maximum", "target_maximum", "decade", m.target_maximum, 0.015, 1.0e-9, 0.0,
            inputs, "the reported maximum matches the analytic target", "r2-minimax"));
        CHECK(agree(id + "/gap", "epigraph_gap", "decade", m.gap, 0.0, 1.0e-9, 0.0,
            inputs, "the epigraph is tight", "r2-minimax"));
        CHECK(m.rows.size() == m.physical_values.size());
        for (std::size_t i = 0; i < m.rows.size(); ++i) {
            const auto& row = m.rows[i];
            const double e = check_error(row.metric, problem.request.observations[row.row], result.predictions[row.row]).first;
            CHECK(agree(id + "/row/" + std::to_string(i), "physical_constraint", "decade", m.physical_values[i],
                row.sign * e - m.t, 1.0e-9, 0.0, inputs,
                "the physical max rows retain their analytic errors", "r2-minimax"));
        }
        check_minimax_kkt(id, m, inputs);
    }
}

void check_reaction_derivatives() {
    FitRequest a1;
    a1.parameters = equilibrium_support::synthetic_reactive_ternary();
    auto liquids = equilibrium_support::synthetic_reactive_problem();
    liquids.standard_gibbs_j_per_mol.clear();
    a1.observations = {solved(ObservationKind::PhaseMoleFraction, liquids, 0, {2}),
                       solved(ObservationKind::PhaseFraction, liquids, 0), amount_ratio(liquids, 0, 1, {2}, {2})};
    a1.coordinates.coordinates = {reaction_coordinate("A+B=C")};
    regular_reaction_state(HelmholtzModel(a1.parameters), a1.observations[0], "reaction-a1");
    reaction_columns("reaction-a1/scalar", a1);
    a1.coordinates.coordinates = {{ParameterFamily::Kij, {{2, 6}}, 0.045, 0.05}, reaction_coordinate("A+B=C"),
                                 {ParameterFamily::SigmaConstant, {{1}}, 3.0, 0.5}};
    const auto interleaved = reaction_columns("reaction-a1/joint", a1);
    const auto& j = interleaved.physical_jacobian;
    for (std::size_t row = 0; row < a1.observations.size(); ++row) {
        const double slope = 1.0 / interleaved.predictions[row];
        for (std::size_t c = 0; c < 3; ++c) {
            const std::size_t physical = c == 1 ? 2 : c == 2 ? 1 : 0;
            CHECK(agree("reaction-a1/order/" + std::to_string(row) + "/" + std::to_string(c), "dlny_du", "1",
                interleaved.optimizer_jacobian[row * 3 + c], slope * j[row * 3 + physical] * a1.coordinates.coordinates[c].scale, 0.0, 1.0e-12, "",
                "optimizer declaration order and EOS-first physical-row order carry the same derivatives", "reaction-fit"));
        }
    }
    for (const auto source_pressure : {std::optional<double>{}, std::optional<double>{2.0e5}}) {
        FitRequest a2;
        a2.parameters = properties_support::aqueous_electrolyte();
        a2.observations = {solved(ObservationKind::PhaseMoleFraction, referenced_formula(source_pressure), 0, {1})};
        a2.coordinates.coordinates = {reaction_coordinate("H2O=Na++Cl-"),
                                     reaction_coordinate("H2O=Na++Cl-", ReactionCoefficient::B, 298.15)};
        const auto direct = predict(HelmholtzModel(a2.parameters), a2.observations.front());
        regular_reaction_state(HelmholtzModel(a2.parameters), a2.observations.front(), source_pressure ? "reaction-a2/fixed" : "reaction-a2/system");
        const auto result = reaction_columns(source_pressure ? "reaction-a2/fixed" : "reaction-a2/system", a2);
        CHECK(agree(source_pressure ? "reaction-a2/fixed_basis" : "reaction-a2/system_basis", "x_Na", "1",
                    result.predictions[0], *direct.value, 0.0, 1.0e-12, "",
                    "rewriting the native coefficient retains the declared source reference", "reaction-fit"));
    }
    FitRequest a3;
    a3.parameters = properties_support::neutral_binary();
    const auto bubble = reacting_bubble(HelmholtzModel(a3.parameters));
    a3.observations = {solved(ObservationKind::EquilibriumPressure, bubble), solved(ObservationKind::PartialPressure, bubble, 1, {0})};
    a3.coordinates.coordinates = {reaction_coordinate("iso"), reaction_coordinate("iso", ReactionCoefficient::B, 180.0),
        reaction_coordinate("iso", ReactionCoefficient::C), reaction_coordinate("iso", ReactionCoefficient::D, 1.0 / 180.0)};
    regular_reaction_state(HelmholtzModel(a3.parameters), a3.observations.front(), "reaction-a3");
    const auto result = reaction_columns("reaction-a3/free_pressure", a3);
    const double basis[] = {1.0, 1.0 / 180.0, std::log(180.0 / 200.0), 180.0};
    for (std::size_t row = 0; row < 2; ++row)
        for (std::size_t c = 1; c < 4; ++c)
            CHECK(agree("reaction-a3/basis/" + std::to_string(row) + "/" + std::to_string(c), "dy_dtheta", "native",
                        result.physical_jacobian[row * 4 + c], result.physical_jacobian[row * 4] * basis[c], 0.0, 1.0e-12,
                        "T=180 K; Tref=200 K", "physical columns obey the source ln K basis vector", "reaction-fit"));
    auto free_t = bubble;
    free_t.pressure = FreeValue::fixed(result.predictions.front());
    free_t.temperature = FreeValue::variable(180.0);
    const HelmholtzModel model(a3.parameters);
    const auto j1 = equilibrium::solve_equilibrium(model, bubble);
    for (std::size_t q = 0; q < free_t.phases.size(); ++q) {
        free_t.phases[q].composition_guess = {j1.mole_fractions[q * 2], j1.mole_fractions[q * 2 + 1]};
        free_t.phases[q].packing_guess = j1.packing_fractions[q];
    }
    const auto recovered = equilibrium::solve_equilibrium(model, free_t);
    CHECK(recovered.success && recovered.topology_event == equilibrium::TopologyEvent::None);
    CHECK(agree("reaction-a4/pinned_temperature", "T", "K", recovered.temperature.front(), 180.0, 0.0, 1.0e-10,
                "declared J1 compositions and packing fractions", "the free-temperature solve returns the specified J1 state", "reaction-fit"));
    CHECK(agree("reaction-a4/pinned_pressure", "P", "Pa", recovered.pressure.front(), j1.pressure.front(), 0.0, 1.0e-10,
                "", "the free-temperature solve retains the J1 pressure", "reaction-fit"));
    for (std::size_t i = 0; i < j1.mole_fractions.size(); ++i)
        CHECK(agree("reaction-a4/pinned_x/" + std::to_string(i), "x", "1", recovered.mole_fractions[i], j1.mole_fractions[i],
                    1.0e-10, 0.0, "", "the free-temperature solve retains the J1 compositions", "reaction-fit"));
    for (std::size_t q = 0; q < j1.packing_fractions.size(); ++q)
        CHECK(agree("reaction-a4/pinned_packing/" + std::to_string(q), "eta", "1", recovered.packing_fractions[q], j1.packing_fractions[q],
                    0.0, 1.0e-10, "", "the free-temperature solve retains the J1 packing fractions", "reaction-fit"));
    a3.observations = {solved(ObservationKind::EquilibriumTemperature, free_t)};
    a3.coordinates.coordinates = {reaction_coordinate("iso"), reaction_coordinate("iso", ReactionCoefficient::D, 1.0 / 180.0)};
    regular_reaction_state(HelmholtzModel(a3.parameters), a3.observations.front(), "reaction-a4");
    reaction_columns("reaction-a4/free_temperature", a3);
}

void recover_reactions(const std::string& id, const FitRequest& request, const std::vector<double>& expected) {
    const auto result = fit(request);
    CHECK(truth(id + "/converged", "the manufactured fit converges without failed trials",
                result.status == "converged" && result.usable && result.trial_failures == 0, summary(result), "reaction-fit"));
    for (std::size_t row = 0; row < expected.size(); ++row) {
        const bool b = request.coordinates.coordinates[request.coordinates.tie()[row]].reaction
            && request.coordinates.coordinates[request.coordinates.tie()[row]].reaction->coefficient == ReactionCoefficient::B;
        CHECK(agree(id + "/value/" + std::to_string(row), "theta", b ? "K" : "1", result.physical[row], expected[row],
                    b ? 0.0 : 1.0e-6, b ? 1.0e-6 : 0.0, summary(result), "the known physical coefficient is recovered", "reaction-fit"));
    }
    CHECK(agree(id + "/cost", "cost", "1", result.final_cost, 0.0, 1.0e-14, 0.0, summary(result),
                "the noise-free final objective vanishes", "reaction-fit"));
    CHECK(truth(id + "/covariance", "the interior fit has full rank and an available conditional covariance",
                result.covariance.status == "available" && result.covariance.rank == expected.size(), summary(result), "reaction-fit"));
}

void check_reaction_recovery() {
    FitRequest b1;
    b1.parameters = equilibrium_support::synthetic_reactive_ternary();
    auto problem = equilibrium_support::synthetic_reactive_problem();
    problem.standard_gibbs_j_per_mol.clear();
    problem.reactions[0].equilibrium_constant = std::exp(0.1);
    b1.parameters.k_ij[2] = b1.parameters.k_ij[6] = 0.05;
    b1.observations = {solved(ObservationKind::PhaseMoleFraction, problem, 0, {2}),
        solved(ObservationKind::PhaseMoleFraction, problem, 1, {0}), solved(ObservationKind::PhaseFraction, problem, 0),
        amount_ratio(problem, 0, 1, {2}, {2})};
    b1.coordinates.coordinates = {{ParameterFamily::Kij, {{2, 6}}, 0.045, 0.1}, reaction_coordinate("A+B=C")};
    reaction_data(b1);
    b1.parameters.k_ij[2] = b1.parameters.k_ij[6] = 0.045;
    for (auto& row : b1.observations) row.problem.reactions[0].equilibrium_constant = 1.0;
    recover_reactions("reaction-b1", b1, {0.05, 0.1});
    FitRequest b2;
    b2.parameters = properties_support::neutral_binary();
    const auto bubble = reacting_bubble(HelmholtzModel(b2.parameters));
    const double a_truth = bubble.reactions[0].correlation->a;
    for (const double T : {170.0, 180.0, 190.0}) {
        auto cell = bubble;
        cell.temperature = FreeValue::fixed(T);
        b2.observations.push_back(solved(ObservationKind::EquilibriumPressure, cell));
        b2.observations.push_back(solved(ObservationKind::PhaseMoleFraction, cell, 1, {0}));
    }
    b2.coordinates.coordinates = {reaction_coordinate("iso"), reaction_coordinate("iso", ReactionCoefficient::B, 180.0)};
    reaction_data(b2);
    for (auto& row : b2.observations) { row.problem.reactions[0].correlation->a += 0.1; row.problem.reactions[0].correlation->b -= 20.0; }
    recover_reactions("reaction-b2", b2, {a_truth, 300.0});
    FitRequest b3;
    b3.parameters = properties_support::neutral_ternary();
    auto vapor = equilibrium_support::lattice_problem(equilibrium_support::PhaseSet::OneVapor,
        equilibrium_support::AmountPattern::AllFree, equilibrium_support::IntensivePattern::BothFixed,
        equilibrium_support::Chemistry::TwoReactions);
    b3.coordinates.coordinates = {reaction_coordinate("reaction_a"), reaction_coordinate("reaction_b")};
    for (const auto& feed : {std::vector<double>{0.6, 0.3, 0.1}, std::vector<double>{0.2, 0.3, 0.5}}) {
        vapor.feed.values = feed;
        if (!b3.observations.empty()) std::reverse(vapor.reactions.begin(), vapor.reactions.end());
        for (auto& r : vapor.reactions) r.equilibrium_constant = 1.0;
        regular_reaction_state(HelmholtzModel(b3.parameters), solved(ObservationKind::PhaseMoleFraction, vapor, 0, {0}), "reaction-b3/feed_" + std::to_string(b3.observations.size() / 2));
        for (auto& r : vapor.reactions) r.equilibrium_constant = std::exp(r.name == "reaction_a" ? 0.2 : -0.1);
        b3.observations.push_back(solved(ObservationKind::PhaseMoleFraction, vapor, 0, {0}));
        b3.observations.push_back(solved(ObservationKind::PhaseMoleFraction, vapor, 0, {2}));
    }
    reaction_data(b3);
    for (auto& row : b3.observations) for (auto& r : row.problem.reactions) r.equilibrium_constant = 1.0;
    recover_reactions("reaction-b3/name_tying", b3, {0.2, -0.1});
    // A large affine offset loses the supplied 0.1 ln K during resolution. Validation's
    // compiled problem is therefore stale even though its EOS active-value key is empty.
    FitRequest preparation;
    preparation.parameters = equilibrium_support::synthetic_reactive_ternary();
    auto source = equilibrium_support::synthetic_reactive_problem();
    source.standard_gibbs_j_per_mol.clear();
    preparation.observations = {solved(ObservationKind::PhaseMoleFraction, source, 0, {2})};
    preparation.coordinates.coordinates = {reaction_coordinate("A+B=C")};
    preparation.coordinates.coordinates[0].origin = std::ldexp(1.0, 54);
    reaction_data(preparation);  // data at the independently solved resolved value ln K = 0
    preparation.observations[0].problem.reactions[0].equilibrium_constant = std::exp(0.1);
    const auto u = preparation.coordinates.start(preparation.parameters, preparation.observations);
    const double resolved = preparation.coordinates.physical(u).front();
    CHECK(resolved == 0.0);
    const auto supplied = predict(HelmholtzModel(preparation.parameters), preparation.observations[0]);
    CHECK(supplied.status == ObservationStatus::Available);
    const double signal = std::abs(std::log(*supplied.value / preparation.observations[0].observed));
    CHECK(signal >= 1.0e-2);
    const auto rebuilt = fit(preparation);
    const auto inputs = fields({field("supplied_lnK", 0.1), field("resolved_lnK", resolved), field("log_prediction_change", signal)});
    CHECK(agree("reaction-preparation/prediction", "x_C", "1", rebuilt.predictions.front(), preparation.observations[0].observed,
                0.0, 1.0e-12, inputs, "a changed resolved reaction constant invalidates validation preparation", "reaction-fit"));
    CHECK(agree("reaction-preparation/initial_cost", "cost", "1", rebuilt.initial_cost, 0.0, 1.0e-14, 0.0, inputs,
                "initial prediction uses the resolved constant rather than the supplied compiled problem", "reaction-fit"));
}

void check_reaction_rejections() {
    FitRequest request;
    request.parameters = properties_support::neutral_ternary();
    auto problem = equilibrium_support::lattice_problem(equilibrium_support::PhaseSet::OneVapor,
        equilibrium_support::AmountPattern::AllFree, equilibrium_support::IntensivePattern::BothFixed,
        equilibrium_support::Chemistry::TwoReactions);
    request.coordinates.coordinates = {reaction_coordinate("reaction_a")};
    request.observations = {solved(ObservationKind::PhaseMoleFraction, problem, 0, {0})};
    const auto reject = [&](FitRequest candidate, const std::string& cause) {
        const auto result = fit(candidate);
        CHECK(truth("reaction-rejection/" + cause, "an overdetermined or fixed reaction input is refused by name",
            result.status == "invalid_setup" && result.message.find(cause) != std::string::npos
            && result.message.find("reaction_a") != std::string::npos,
            fields({field("message", result.message)}), "reaction-fit"));
    };
    auto affinity = request;
    affinity.observations[0].problem.reactions[0].equilibrium_constant.reset();
    affinity.observations[0].problem.reactions[0].affinity_j_per_mol = 0.0;
    reject(affinity, "affinity target");
    auto gibbs = request;
    gibbs.observations[0].problem.standard_gibbs_j_per_mol = {0.0, 0.0, 0.0};
    reject(gibbs, "supplied standard Gibbs");
    auto dependent = request;
    dependent.observations[0].problem.reactions[1].stoichiometry = problem.reactions[0].stoichiometry;
    reject(dependent, "dependent reaction");
    auto fixed = request;
    ThermochemistryRecord record;
    record.standard_pressure = 1.0e5;
    record.component_intervals.assign(3, {{200.0, 400.0, true, true, {IdealPolynomial{{30.0}, 300.0, 0.0, 0.0}, 1.0e5}}});
    fixed.thermochemistry = record;
    reject(fixed, "K fixed by thermochemistry");
    auto free_t = request;
    free_t.observations[0].problem.reactions[0].reference = ReactionReference{};
    free_t.observations[0].problem.temperature = FreeValue::variable(300.0);
    reject(free_t, "referenced reaction at free temperature");
    // Two named source constants with different declarations are physically different coordinates.
    auto conflicting = request;
    conflicting.observations.push_back(request.observations[0]);
    conflicting.observations[1].problem.reactions[0].equilibrium_constant = std::exp(0.01);
    reject(conflicting, "tied declarations or starts disagree");
    auto scalar_b = request;
    scalar_b.coordinates.coordinates[0].reaction->coefficient = ReactionCoefficient::B;
    reject(scalar_b, "scalar K admits only A");
    auto reference = request;
    reference.parameters = properties_support::aqueous_electrolyte();
    auto source = referenced_formula(std::nullopt);
    source.reactions[0].name = "reaction_a";
    source.reactions[0].reference->source_basis = ReferenceSourceBasis::CommonMolalityInfiniteDilution;
    reference.observations = {solved(ObservationKind::PhaseMoleFraction, source, 0, {1}),
                              solved(ObservationKind::PhaseMoleFraction, source, 0, {1})};
    reference.observations[1].problem.source_reference->standard_molality_mol_per_kg = 2.0;
    reject(reference, "tied source bases disagree");
    reference.observations[1].problem.source_reference = source.source_reference;
    reference.observations[1].problem.reactions[0].correlation->reference_temperature = 300.0;
    reject(reference, "tied declarations or starts disagree");
}

// #211: the source constant is measured directly; no EOS reference transfer enters it.
Observation log_constant(const EquilibriumProblem& problem, const char* name) {
    Observation row;
    row.kind = ObservationKind::LogEquilibriumConstant;
    row.problem = problem;
    row.reaction = name;
    return row;
}

void check_log_constant() {
    const auto parameters = properties_support::neutral_ternary();
    auto problem = equilibrium_support::lattice_problem(equilibrium_support::PhaseSet::OneVapor,
        equilibrium_support::AmountPattern::AllFree, equilibrium_support::IntensivePattern::BothFixed,
        equilibrium_support::Chemistry::TwoReactions);
    problem.reactions[0].equilibrium_constant.reset();
    problem.reactions[0].correlation = ReactionLogPolynomial{0.4, 80.0, 0.5, 0.001, 300.0, 150.0, 400.0,
                                                           equilibrium::EOS_STANDARD_STATE_ID};
    CoordinateSet coordinates{{reaction_coordinate("reaction_a", ReactionCoefficient::B, 300.0),
        {ParameterFamily::SegmentCount, {{0}}, 0.0, 1.0}, reaction_coordinate("reaction_a"),
        reaction_coordinate("reaction_b"), reaction_coordinate("reaction_a", ReactionCoefficient::C),
        reaction_coordinate("reaction_a", ReactionCoefficient::D, 1.0 / 300.0)}};
    const HelmholtzModel model(parameters, coordinates.active_parameters());
    for (const double T : {200.0, 250.0, 300.0, 350.0}) {
        problem.temperature = FreeValue::fixed(T);
        std::reverse(problem.reactions.begin(), problem.reactions.end());
        auto row = log_constant(problem, "reaction_a");
        validate(row, model);
        const auto prediction = predict(model, row);
        const double value = 0.4 + 80.0 / T + 0.5 * std::log(T / 300.0) + 0.001 * T;
        const auto id = "ln-k/" + std::to_string(T);
        CHECK(prediction.status == ObservationStatus::Available);
        CHECK(agree(id + "/value", "ln_K_src", "1", *prediction.value, value, 0.0, 0.0,
                    fields({field("T_K", T), field("Tref_K", 300.0)}), "declared source polynomial", "ln-k"));
        const auto gradient = differentiate_group(model, {&row}, prediction, coordinates).front();
        CHECK(gradient.status == ObservationStatus::Available);
        // Physical rows: EOS, b, a, other reaction a, c, d; optimizer order is interleaved.
        const std::vector<double> expected{0.0, 1.0 / T, 1.0, 0.0, std::log(T / 300.0), T};
        for (std::size_t k = 0; k < expected.size(); ++k)
            CHECK(agree(id + "/column/" + std::to_string(k), "dln_K_src_dtheta", "native",
                gradient.value[k], expected[k], 0.0, 0.0, fields({field("T_K", T), field("Tref_K", 300.0)}),
                "exact coefficient basis; EOS and other reaction columns vanish", "ln-k"));
    }
    // A referenced source value remains the source value even though its density-basis target differs.
    auto reference = log_constant(referenced_formula(std::nullopt), "H2O=Na++Cl-");
    const HelmholtzModel aqueous(properties_support::aqueous_electrolyte());
    validate(reference, aqueous);
    CHECK(agree("ln-k/source", "ln_K_src", "1", *predict(aqueous, reference).value, std::log(1.0e-8),
                0.0, 0.0, "T=298.15 K; raw mole-fraction source", "no EOS reference transfer in direct data", "ln-k"));
    auto scalar = log_constant(problem, "reaction_b");
    CHECK(agree("ln-k/scalar", "ln_K_src", "1", *predict(model, scalar).value,
                std::log(*problem.reactions[1].equilibrium_constant), 0.0, 0.0,
                "T=350 K", "scalar K is observed as ln K", "ln-k"));

    FitRequest request;
    request.parameters = properties_support::neutral_binary();
    const auto bubble = reacting_bubble(HelmholtzModel(request.parameters));
    const auto law = *bubble.reactions[0].correlation;
    request.coordinates.coordinates = {reaction_coordinate("iso"), reaction_coordinate("iso", ReactionCoefficient::B, 180.0)};
    for (const double T : {170.0, 180.0, 190.0}) {
        auto cell = bubble;
        cell.temperature = FreeValue::fixed(T);
        auto row = log_constant(cell, "iso");
        row.observed = law.a + law.b / T + law.c * std::log(T / law.reference_temperature) + law.d * T;
        request.observations.push_back(row);
    }
    auto equilibrium = solved(ObservationKind::EquilibriumPressure, bubble);
    equilibrium.observed = *predict(HelmholtzModel(request.parameters), equilibrium).value;
    request.observations.push_back(equilibrium);
    for (auto& row : request.observations) { row.problem.reactions[0].correlation->a += 0.1; row.problem.reactions[0].correlation->b -= 20.0; }
    recover_reactions("ln-k/recovery", request, {law.a, law.b});

    // Observation admission must reject fixed targets even when only an EOS coordinate is fitted.
    FitRequest rejected;
    rejected.parameters = parameters;
    rejected.coordinates.coordinates = {{ParameterFamily::SegmentCount, {{0}}, 0.0, 1.0}};
    rejected.observations = {log_constant(problem, "reaction_a")};
    const auto reject = [&](FitRequest candidate, const std::string& cause) {
        const auto result = fit(candidate);
        CHECK(truth("ln-k/rejection/" + cause, "the direct source target is refused by reaction name",
            result.status == "invalid_setup" && result.message.find(cause) != std::string::npos
            && result.message.find("reaction_a") != std::string::npos,
            fields({field("message", result.message)}), "ln-k"));
    };
    // At the last reversal reaction_a is first.
    auto affinity = rejected;
    affinity.observations[0].problem.reactions[0].correlation.reset();
    affinity.observations[0].problem.reactions[0].affinity_j_per_mol = 0.0;
    reject(affinity, "affinity target");
    auto dependent = rejected;
    dependent.observations[0].problem.reactions[1].stoichiometry = problem.reactions[0].stoichiometry;
    reject(dependent, "dependent reaction");
    auto gibbs = rejected;
    gibbs.observations[0].problem.standard_gibbs_j_per_mol = {0.0, 0.0, 0.0};
    reject(gibbs, "supplied standard Gibbs");
    auto fixed = rejected;
    ThermochemistryRecord record;
    record.standard_pressure = 1.0e5;
    record.component_intervals.assign(3, {{200.0, 400.0, true, true, {IdealPolynomial{{30.0}, 300.0, 0.0, 0.0}, 1.0e5}}});
    fixed.thermochemistry = record;
    reject(fixed, "K fixed by thermochemistry");
}

// #210: finite enthalpy increments on a referenced, reaction-completed ethane pair.
Parameters addition_parameters() {
    auto result = properties_support::parameters(2);
    const auto ethane = properties_support::ethane();
    for (std::size_t i = 0; i < 2; ++i) {
        result.segment_count[i] = ethane.segment_count[0];
        result.sigma_constant[i] = ethane.sigma_constant[0];
        result.dispersion_energy_over_k[i] = ethane.dispersion_energy_over_k[0];
        result.molar_mass_kg_per_mol[i] = ethane.molar_mass_kg_per_mol[0];
    }
    return result;
}

ThermochemistryRecord addition_record() {
    ThermochemistryRecord result;
    result.standard_pressure = 1.0e5;
    result.component_intervals = {
        {{200.0, 400.0, true, true, {IdealPolynomial{{52.0, 0.12}, 300.0, 0.0, 0.0}, 1.0e5}}},
        {}};
    return result;
}

ReactionLogPolynomial addition_law(double reference_temperature = 300.0) {
    const double a = std::log(0.3) + 3.0 * std::log(reference_temperature / 300.0);
    return {a, 200.0, 3.0, -0.004, reference_temperature, 250.0, 350.0, "analytic-source"};
}

EquilibriumProblem addition_problem(double temperature, double pressure, std::vector<double> feed,
                                    std::optional<double> reference_pressure = {}, bool free_pressure = false,
                                    std::optional<double> vapor_fraction = {}) {
    EquilibriumProblem result;
    result.phases = {equilibrium_support::phase("L", PhaseKind::Liquid)};
    if (vapor_fraction) {
        auto vapor = equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(*vapor_fraction));
        result.phases.push_back(std::move(vapor));
    } else if (free_pressure) {
        result.phases.push_back(equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.0)));
    }
    result.temperature = FreeValue::fixed(temperature);
    result.pressure = free_pressure ? FreeValue::variable(pressure) : FreeValue::fixed(pressure);
    result.feed.kind = equilibrium::FeedKind::Amounts;
    result.feed.values = std::move(feed);
    result.source_reference = NeutralReferenceDeclaration{0, {0.0, 1.0}, 1.0};
    result.reactions = {{"2A=B", {-2.0, 1.0}, std::nullopt, std::nullopt, addition_law(),
        ReactionReference{reference_pressure, ReferenceSourceBasis::RawMoleFractionInfiniteDilution}}};
    result.row_tolerance = 1.0e-10;
    return result;
}

Observation addition_observation(const EquilibriumProblem& initial, const EquilibriumProblem& final,
                                 double amount, double feed_temperature = 260.0,
                                 double feed_pressure = 1.0e5, std::vector<double> feed = {1.0, 0.0}) {
    Observation result;
    result.kind = ObservationKind::AdditionEnthalpy;
    result.problem = initial;
    result.final_problem = final;
    result.feed_amount = amount;
    result.temperature = feed_temperature;
    result.pressure = feed_pressure;
    result.composition = std::move(feed);
    result.branch = DensityBranch::Vapor;
    result.form = ResidualForm::Difference;
    result.scale = 1.0;
    return result;
}

CoordinateSet addition_coordinates() {
    CoordinateSet result;
    result.coordinates = {
        {ParameterFamily::Kij, {{1, 2}}, 0.0, 1.0, -0.2, 0.2},
        {ParameterFamily::DispersionEnergyOverK, {{0}, {1}}, 191.42, 1.0, 150.0, 250.0},
        reaction_coordinate("2A=B"),
        reaction_coordinate("2A=B", ReactionCoefficient::B),
        reaction_coordinate("2A=B", ReactionCoefficient::C),
        reaction_coordinate("2A=B", ReactionCoefficient::D)};
    return result;
}

struct AdditionEndpoint final {
    double enthalpy = 0.0;
    std::vector<double> phase_amount, phase_enthalpy;
};

AdditionEndpoint endpoint_enthalpy(const HelmholtzModel& model, const EquilibriumProblem& problem,
                                   const equilibrium::CompiledSystem& compiled,
                                   const EquilibriumResult& solved_state) {
    AdditionEndpoint result;
    result.phase_amount.assign(problem.phases.size(), 0.0);
    result.phase_enthalpy.assign(problem.phases.size(), 0.0);
    std::vector<equilibrium::SolvedStateObservable> observables;
    std::vector<std::size_t> phases;
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        if (problem.phases[phase].is_pinned() && problem.phases[phase].amount.fraction == 0.0) continue;
        equilibrium::SolvedStateObservable amount, enthalpy;
        amount.kind = equilibrium::SolvedStateObservableKind::PhaseAmount;
        amount.phase = phase;
        enthalpy.kind = equilibrium::SolvedStateObservableKind::PhaseProperty;
        enthalpy.phase = phase;
        enthalpy.property = PropertyObservable::TotalEnthalpy;
        observables.push_back(amount);
        observables.push_back(enthalpy);
        phases.push_back(phase);
    }
    equilibrium::SolvedStateActionDirection zero;
    zero.feed_amounts_mol.assign(model.parameters().component_count(), 0.0);
    zero.active_parameters.assign(model.active_parameters().size(), 0.0);
    zero.reaction_log_constants.assign(problem.reactions.size(), 0.0);
    const auto batch = equilibrium::solved_state_actions(compiled, solved_state, {observables, {zero}});
    CHECK(batch.results.size() == 2 * phases.size());
    for (std::size_t i = 0; i < phases.size(); ++i) {
        const auto& amount = batch.results[2 * i];
        const auto& enthalpy = batch.results[2 * i + 1];
        CHECK(amount.status == equilibrium::SolvedStateActionStatus::Available && amount.value);
        CHECK(enthalpy.status == equilibrium::SolvedStateActionStatus::Available && enthalpy.value);
        result.phase_amount[phases[i]] = *amount.value;
        result.phase_enthalpy[phases[i]] = *enthalpy.value;
        result.enthalpy += *amount.value * *enthalpy.value;
    }
    return result;
}

double addition_oracle_phase_enthalpy(const HelmholtzModel& pure, const EquilibriumResult& solved_state,
                                     const ReactionLogPolynomial& law, std::size_t phase) {
    const double T = solved_state.temperature.front();
    const double P_ref = solved_state.pressure.front();
    const double dt = T - 300.0;
    const double h_a = 52.0 * dt + 0.06 * dt * dt;
    const double ell = -law.b / (T * T) + law.c / T + law.d;
    const double h_b = 2.0 * h_a + GAS_CONSTANT_J_PER_MOL_K * T * T * ell
        + pure.state_at_pressure(T, P_ref, {1.0}, DensityBranch::Liquid).residual_enthalpy;
    const double x_a = solved_state.mole_fractions[phase * 2];
    const double x_b = solved_state.mole_fractions[phase * 2 + 1];
    const double residual = pure.state(T, solved_state.molar_densities[phase], {1.0}).residual_enthalpy;
    return x_a * h_a + x_b * h_b + residual;
}

double addition_oracle_endpoint(const HelmholtzModel& pure, const EquilibriumResult& solved_state,
                                const ReactionLogPolynomial& law) {
    double result = 0.0;
    for (std::size_t phase = 0; phase < solved_state.phase_count; ++phase)
        result += solved_state.amounts[phase] * addition_oracle_phase_enthalpy(pure, solved_state, law, phase);
    return result;
}

double total_b_amount(const EquilibriumResult& solved_state) {
    double result = 0.0;
    for (std::size_t phase = 0; phase < solved_state.phase_count; ++phase)
        result += solved_state.amounts[phase] * solved_state.mole_fractions[phase * 2 + 1];
    return result;
}

void addition_timing(const std::string& id, double seconds, std::size_t c_ttheta_e3_count,
                     const std::string& inputs) {
    recording_support::observe((id + "/jacobian_seconds").c_str(), "differentiate_group_wall_time", "s", seconds,
        std::nullopt, inputs, "Observation", "steady_clock", "wall time for one native observation Jacobian",
        "addition-enthalpy");
    recording_support::observe((id + "/new_c_ttheta_e3_count").c_str(), "C_Ttheta_E3_contractions", "count",
        static_cast<double>(c_ttheta_e3_count), std::nullopt, inputs, "Count", "reference_path_polarization",
        "new C_Ttheta contractions only; excludes pressure and other E3 work", "addition-enthalpy");
}

struct AdditionVariant final {
    std::string name;
    Observation observation;
    Prediction prediction;
    std::vector<double> gradient;
};

AdditionVariant addition_variant(const std::string& name, double amount, std::optional<double> reference_pressure,
                                 bool free_pressure = false, std::optional<double> vapor_fraction = {},
                                 bool unequal_temperature = false) {
    const double pressure = free_pressure ? 1.5e6 : 3.0e6;  // original Case N regularity input
    auto initial = addition_problem(260.0, pressure, {2.0, 0.2}, reference_pressure, free_pressure, vapor_fraction);
    auto final = addition_problem(unequal_temperature ? 265.0 : 260.0, pressure, {2.0 + amount, 0.2},
                                  reference_pressure, free_pressure, vapor_fraction);
    if (unequal_temperature) {
        initial.reactions[0].correlation = addition_law(260.0);
        final.reactions[0].correlation = addition_law(260.0);
    }
    auto row = addition_observation(initial, final, amount);
    return {name, std::move(row), {}, {}};
}

double addition_coordinate_gradient(const std::vector<double>& gradient, std::size_t column) {
    return column == 1 ? gradient[1] + gradient[2] : gradient[column < 2 ? column : column + 1];
}

void check_addition_enthalpy() {
    const auto parameters = addition_parameters();
    const auto record = addition_record();
    const auto coordinates = addition_coordinates();
    const HelmholtzModel model(parameters, coordinates.active_parameters(), record);
    const double amount = 0.2;
    std::vector<AdditionVariant> variants;
    variants.push_back(addition_variant("fixed_system", amount, {}));
    variants.push_back(addition_variant("fixed_source", amount, 5.0e6));
    variants.push_back(addition_variant("free_pressure", amount, {}, true));
    variants.push_back(addition_variant("two_phase", amount, {}, true, 0.3));
    variants.push_back(addition_variant("unequal_temperature", amount, {}, false, {}, true));

    const double steps[] = {1.0e-4, 0.05, 5.0e-3, 1.0, 5.0e-3, 2.0e-5};
    const char* names[] = {"kij", "epsilon_tied", "reaction_a", "reaction_b", "reaction_c", "reaction_d"};
    for (auto& variant : variants) {
        auto center = predict(model, variant.observation);
        CHECK(center.status == ObservationStatus::Available && center.value && center.solved
            && center.final_solved && center.compiled && center.final_compiled);
        CHECK(center.solved->success && center.final_solved->success
            && center.solved->topology_event == equilibrium::TopologyEvent::None
            && center.final_solved->topology_event == equilibrium::TopologyEvent::None
            && center.solved->requested_tolerance_met.value_or(false)
            && center.final_solved->requested_tolerance_met.value_or(false)
            && center.solved->jacobian_rank == center.solved->unknown_count
            && center.final_solved->jacobian_rank == center.final_solved->unknown_count);
        variant.prediction = center;
        const auto start = coordinates.start(parameters, {variant.observation}, record);
        const auto t0 = std::chrono::steady_clock::now();
        const auto exact = differentiate_group(model, {&variant.observation}, center, coordinates).front();
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        CHECK(exact.status == ObservationStatus::Available && exact.value.size() == 7);
        variant.gradient = exact.value;
        const std::size_t eos_parameters = model.active_parameters().size();
        addition_timing("addition-f3/" + variant.name, seconds,
            4 * eos_parameters * 3 * variant.observation.problem.reactions.size() * 2,
            fields({field("variant", variant.name), field("n_add_mol", amount), field("active_parameters", eos_parameters)}));

        for (std::size_t column = 0; column < 6; ++column) {
            const auto value = [&](double shift) {
                auto u = start;
                u[column] += shift;
                auto changed_observations = std::vector<Observation>{variant.observation};
                coordinates.apply(changed_observations, u);
                const HelmholtzModel rebuilt(coordinates.apply(parameters, u), coordinates.active_parameters(), record);
                const auto rebuilt_prediction = predict(rebuilt, changed_observations.front());
                if (rebuilt_prediction.status != ObservationStatus::Available || !rebuilt_prediction.value)
                    throw std::runtime_error("addition-f3 " + variant.name + " " + names[column] + ": "
                        + observation_status_name(rebuilt_prediction.status) + " " + rebuilt_prediction.message);
                return *rebuilt_prediction.value;
            };
            auto rungs = ladder(value, steps[column]);
            const double action = addition_coordinate_gradient(exact.value, column);
            const auto inputs = fields({field("variant", variant.name), field("T_initial_K", 260.0),
                field("T_final_K", variant.observation.final_problem.temperature.value), field("step", steps[column]),
                field("gradient", exact.value)});
            const std::string id = "addition-f3/" + variant.name + "/" + names[column];
            CHECK(ladder_agrees(id, "d_delta_h_add_dtheta", action, rungs, inputs,
                "the exact finite-addition gradient agrees with rebuilt-model prediction re-solves", "addition-f3"));
            if (column == 1 || column == 3) {
                const double fine_tolerance = FIRST_ORDER * std::max(1.0, std::abs(rungs.fine));
                double omitted = 0.0;
                if (column == 1) {
                    const HelmholtzModel pure(properties_support::ethane(),
                        {{ParameterFamily::DispersionEnergyOverK, {0}}});
                    const auto pure_action = [&](double T, double P) {
                        const PressureActionDirection direction{0.0, 0.0, {0.0},
                            {{ParameterFamily::DispersionEnergyOverK, {0}, 1.0}}};
                        const auto batch = pure.pressure_actions(T, P, {1.0}, DensityBranch::Liquid,
                            {PropertyObservable::ResidualEnthalpy}, {direction});
                        CHECK(batch.density_status == PropertyActionStatus::Available && batch.results.size() == 1
                            && batch.results[0].status == PropertyActionStatus::Available && batch.results[0].action);
                        return *batch.results[0].action;
                    };
                    const double P1 = variant.observation.problem.reactions[0].reference->reference_pressure_pa
                        .value_or(center.solved->pressure.front());
                    const double P2 = variant.observation.final_problem.reactions[0].reference->reference_pressure_pa
                        .value_or(center.final_solved->pressure.front());
                    omitted = (total_b_amount(*center.final_solved) * pure_action(center.final_solved->temperature.front(), P2)
                        - total_b_amount(*center.solved) * pure_action(center.solved->temperature.front(), P1)) / amount;
                } else {
                    omitted = -GAS_CONSTANT_J_PER_MOL_K
                        * (total_b_amount(*center.final_solved) - total_b_amount(*center.solved)) / amount;
                }
                const std::string direct_name = column == 1 ? "C_Ttheta_borne_term" : "reaction_b_endpoint_term";
                CHECK(truth(id + "/omission_discrimination", "the independently calculated omitted endpoint term exceeds 1e3 times the fine criterion",
                    std::isfinite(omitted) && std::abs(omitted) > 1.0e3 * fine_tolerance,
                    fields({field("direct_term", omitted), field("fine_tolerance", fine_tolerance), field("term", direct_name)}),
                    "addition-f3"));
            }
        }
        if (variant.observation.problem.temperature.value == variant.observation.final_problem.temperature.value) {
            const double temperature = variant.observation.problem.temperature.value;
            const double tref = variant.observation.problem.reactions[0].correlation->reference_temperature;
            const double G = exact.value[3];
            // g_b = G/T - D/T^2, so recover the source-slope coefficient D from a and b.
            const double D = temperature * G - temperature * temperature * exact.value[4];
            CHECK(agree("addition-f3/" + variant.name + "/basis_c", "d_delta_h_add_dc", "J/mol",
                exact.value[5], G * std::log(temperature / tref) + D / temperature, 0.0, 1.0e-10,
                fields({field("T_K", temperature), field("Tref_K", tref), field("G", G), field("D", D)}),
                "the c column follows the source reaction basis identity before coordinate scaling", "addition-f3"));
            CHECK(agree("addition-f3/" + variant.name + "/basis_d", "d_delta_h_add_dd", "J/mol",
                exact.value[6], G * temperature + D, 0.0, 1.0e-10,
                fields({field("T_K", temperature), field("G", G), field("D", D)}),
                "the d column follows the source reaction basis identity before coordinate scaling", "addition-f3"));
        } else {
            const double initial_b = total_b_amount(*center.solved);
            const double initial_c_term = -GAS_CONSTANT_J_PER_MOL_K * 260.0 * initial_b / amount;
            CHECK(truth("addition-f3/unequal_temperature/initial_c_term", "c has a nonzero initial-endpoint caloric term at Tref with no state motion",
                std::isfinite(initial_c_term) && initial_c_term != 0.0,
                fields({field("phi_c", 0.0), field("phi_c_prime_per_K", 1.0 / 260.0),
                    field("initial_caloric_term_J_per_mol", initial_c_term), field("initial_n_B_mol", initial_b)}),
                "addition-f3"));
        }
    }

    auto& reference = variants.front();
    const auto first_h = endpoint_enthalpy(model, reference.observation.problem,
        *reference.prediction.compiled, *reference.prediction.solved);
    const auto second_h = endpoint_enthalpy(model, reference.observation.final_problem,
        *reference.prediction.final_compiled, *reference.prediction.final_solved);
    auto shifted_record = record;
    std::get<IdealPolynomial>(shifted_record.component_intervals[0][0].correlation.form).reference_enthalpy = 1.0e4;
    const HelmholtzModel shifted_model(parameters, coordinates.active_parameters(), shifted_record);
    const auto shifted_prediction = predict(shifted_model, reference.observation);
    CHECK(shifted_prediction.status == ObservationStatus::Available && shifted_prediction.value);
    CHECK(agree("addition-f4/value", "delta_h_add", "J/mol", *shifted_prediction.value, *reference.prediction.value,
        0.0, 1.0e-9, fields({field("anchor_h_J_per_mol", 1.0e4)}),
        "a common record enthalpy shift cancels from the finite addition", "addition-f4"));
    const auto base_gradient = differentiate_group(model, {&reference.observation}, reference.prediction, coordinates).front();
    const auto shifted_gradient = differentiate_group(shifted_model, {&reference.observation}, shifted_prediction, coordinates).front();
    CHECK(base_gradient.status == ObservationStatus::Available && shifted_gradient.status == ObservationStatus::Available
        && base_gradient.value.size() == shifted_gradient.value.size());
    for (std::size_t column = 0; column < coordinates.coordinates.size(); ++column)
        CHECK(agree("addition-f4/gradient/" + std::to_string(column), "d_delta_h_add_dtheta", "native",
            addition_coordinate_gradient(shifted_gradient.value, column),
            addition_coordinate_gradient(base_gradient.value, column), 0.0, 1.0e-9,
            fields({field("anchor_h_J_per_mol", 1.0e4), field("column", column)}),
            "a common record enthalpy shift preserves every gradient column", "addition-f4"));
    const double condition = (std::abs(first_h.enthalpy) + std::abs(second_h.enthalpy))
        / (amount * std::abs(*reference.prediction.value));
    recording_support::observe("addition-f4/conditioning", "enthalpy_cancellation_condition", "1", condition,
        std::nullopt, fields({field("H_initial_J", first_h.enthalpy), field("H_final_J", second_h.enthalpy),
            field("n_add_mol", amount), field("delta_h_add_J_per_mol", *reference.prediction.value)}),
        "Observation", "endpoint_enthalpy_ratio", "the reference cancellation conditioning ratio", "addition-f4");
    CHECK(truth("addition-f4/conditioning_limit", "the endpoint enthalpy cancellation ratio is at most 1e4",
        std::isfinite(condition) && condition <= 1.0e4, field("condition", condition), "addition-f4"));

    for (const double increment : {0.2, 2.0e-3}) {
        auto initial = addition_problem(260.0, 3.0e6, {2.0, 0.2});
        auto final = addition_problem(260.0, 3.0e6, {2.0 + increment, 0.2});
        auto row = addition_observation(initial, final, increment);
        const auto prediction = predict(model, row);
        CHECK(prediction.status == ObservationStatus::Available && prediction.value && prediction.solved
            && prediction.final_solved && prediction.compiled && prediction.final_compiled);
        equilibrium::SolvedStateObservable activity;
        activity.kind = equilibrium::SolvedStateObservableKind::PhaseLogActivity;
        activity.phase = 0;
        activity.component = 0;
        equilibrium::SolvedStateActionDirection thermal;
        thermal.temperature_k = 1.0;
        thermal.feed_amounts_mol.assign(2, 0.0);
        thermal.active_parameters.assign(model.active_parameters().size(), 0.0);
        thermal.reaction_log_constants.assign(1, 0.0);
        const auto action = equilibrium::solved_state_first_actions(*prediction.compiled, *prediction.solved,
            {activity}, {thermal}).front();
        CHECK(action.results.size() == 1 && action.results[0].status == equilibrium::SolvedStateActionStatus::Available
            && action.results[0].action);
        const HelmholtzModel pure(properties_support::ethane());
        const double h_feed = pure.state_at_pressure(260.0, 1.0e5, {1.0}, DensityBranch::Vapor).residual_enthalpy;
        const double q0 = -GAS_CONSTANT_J_PER_MOL_K * 260.0 * 260.0 * *action.results[0].action
            - GAS_CONSTANT_J_PER_MOL_K * 260.0 - h_feed;
        const auto h_initial = endpoint_enthalpy(model, row.problem, *prediction.compiled, *prediction.solved);
        const auto h_final = endpoint_enthalpy(model, row.final_problem, *prediction.final_compiled, *prediction.final_solved);
        const double roundoff = 64.0 * EPS * (std::abs(h_initial.enthalpy) + std::abs(h_final.enthalpy)) / increment;
        CHECK(reference_scale_agrees("addition-f5/" + number(increment), "finite_addition_limit", "J/mol", *prediction.value, q0,
            roundoff, 1.0e-8, fields({field("n_add_mol", increment), field("D_n", *prediction.value),
                field("q0", q0), field("H_initial_J", h_initial.enthalpy), field("H_final_J", h_final.enthalpy),
                field("dln_activity_dT_per_K", *action.results[0].action), field("h_feed_residual_J_per_mol", h_feed)}),
            "each finite quotient equals the Gibbs-Helmholtz limit for the linear case", "addition-f5"));
    }

    FitRequest recovery;
    auto truth_parameters = parameters;
    truth_parameters.dispersion_energy_over_k = {191.42, 191.42};
    const HelmholtzModel truth_model(truth_parameters, {}, record);
    recovery.parameters = parameters;
    recovery.parameters.dispersion_energy_over_k = {195.2484, 195.2484};
    recovery.thermochemistry = record;
    recovery.coordinates.coordinates = {
        {ParameterFamily::DispersionEnergyOverK, {{0}, {1}}, 191.42, 10.0, 150.0, 250.0},
        reaction_coordinate("2A=B", ReactionCoefficient::B, 260.0)};
    recovery.coordinates.coordinates[1].origin = 200.0;
    recovery.coordinates.coordinates[1].lower = -2000.0;
    recovery.coordinates.coordinates[1].upper = 2000.0;
    for (const double T : {255.0, 260.0, 265.0, 270.0}) {
        auto initial = addition_problem(T, 3.0e6, {2.0, 0.2});
        auto final = addition_problem(T, 3.0e6, {2.2, 0.2});
        auto row = addition_observation(initial, final, 0.2, T);
        row.scale = 100.0;
        const auto truth_prediction = predict(truth_model, row);
        CHECK(truth_prediction.status == ObservationStatus::Available && truth_prediction.value);
        row.observed = *truth_prediction.value;
        row.problem.reactions[0].correlation->b = 180.0;
        row.final_problem.reactions[0].correlation->b = 180.0;
        recovery.observations.push_back(std::move(row));
    }
    const HelmholtzModel start_model(recovery.parameters, recovery.coordinates.active_parameters(), record);
    for (std::size_t row_index = 0; row_index < recovery.observations.size(); ++row_index) {
        const auto& row = recovery.observations[row_index];
        const auto prediction = predict(start_model, row);
        CHECK(prediction.status == ObservationStatus::Available);
        const auto begin = std::chrono::steady_clock::now();
        const auto gradient = differentiate_group(start_model, {&row}, prediction, recovery.coordinates).front();
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        CHECK(gradient.status == ObservationStatus::Available && gradient.value.size() == 3);
        addition_timing("addition-f6/row_" + std::to_string(row_index), seconds,
            4 * start_model.active_parameters().size() * 3 * row.problem.reactions.size() * 2,
            fields({field("T_K", row.temperature), field("active_parameters", start_model.active_parameters().size()),
                field("jacobian_row", row_index)}));
    }
    const auto recovered = fit(recovery);
    CHECK(truth("addition-f6/converged", "the noise-free known-parameter fit converges without rejected trials",
        recovered.status == "converged" && recovered.usable && recovered.trial_failures == 0,
        fields({summary(recovered), field("trial_failures", recovered.trial_failures)}), "addition-f6"));
    for (std::size_t row = 0; row < recovery.observations.size(); ++row)
        recording_support::observe(("addition-f6/data/" + std::to_string(row)).c_str(), "delta_h_add", "J/mol",
            recovered.predictions[row], recovery.observations[row].observed,
            fields({field("T_K", recovery.observations[row].temperature), field("P_Pa", 3.0e6),
                field("scale_J_per_mol", recovery.observations[row].scale)}),
            "Observation", "known_parameter_recovery", "synthetic datum and fitted prediction", "addition-f6");
    CHECK(recovered.physical.size() == 3);
    CHECK(agree("addition-f6/epsilon_s0", "epsilon_over_k", "K", recovered.physical[0], 191.42, 0.0, 1.0e-6,
        summary(recovered), "the tied epsilon coordinate recovers the known value on s0", "addition-f6"));
    CHECK(agree("addition-f6/epsilon_s1", "epsilon_over_k", "K", recovered.physical[1], 191.42, 0.0, 1.0e-6,
        summary(recovered), "the tied epsilon coordinate recovers the known value on s1", "addition-f6"));
    CHECK(agree("addition-f6/reaction_b", "reaction_b", "K", recovered.physical[2], 200.0, 0.0, 1.0e-6,
        summary(recovered), "the reaction b coordinate recovers the known value", "addition-f6"));
    CHECK(agree("addition-f6/cost", "cost", "1", recovered.final_cost, 0.0, 1.0e-14, 0.0,
        summary(recovered), "the noise-free residual cost is at most 1e-14", "addition-f6"));
    CHECK(truth("addition-f6/covariance", "the fit has an available full-rank covariance",
        recovered.covariance.status == "available" && recovered.covariance.rank == 2,
        fields({field("covariance", recovered.covariance.status), field("rank", recovered.covariance.rank),
            field("singular_values", recovered.covariance.singular_values)}), "addition-f6"));

    auto two_phase = std::find_if(variants.begin(), variants.end(), [](const AdditionVariant& item) {
        return item.name == "two_phase";
    });
    CHECK(two_phase != variants.end());
    const auto pure = HelmholtzModel(properties_support::ethane());
    const double oracle_initial = addition_oracle_endpoint(pure, *two_phase->prediction.solved,
        *two_phase->observation.problem.reactions[0].correlation);
    const double oracle_final = addition_oracle_endpoint(pure, *two_phase->prediction.final_solved,
        *two_phase->observation.final_problem.reactions[0].correlation);
    const auto public_initial = endpoint_enthalpy(model, two_phase->observation.problem,
        *two_phase->prediction.compiled, *two_phase->prediction.solved);
    const auto public_final = endpoint_enthalpy(model, two_phase->observation.final_problem,
        *two_phase->prediction.final_compiled, *two_phase->prediction.final_solved);
    for (const auto& [label, oracle, public_value] : {
            std::tuple<const char*, double, double>{"initial", oracle_initial, public_initial.enthalpy},
            {"final", oracle_final, public_final.enthalpy}}) {
        CHECK(reference_scale_agrees(std::string("addition-f8/") + label + "/endpoint_enthalpy", "H", "J", public_value, oracle,
            1.0e-6, 1.0e-9, fields({field("oracle_H_J", oracle), field("public_H_J", public_value)}),
            "each endpoint total enthalpy agrees with the independent pure-fluid oracle", "addition-f8"));
    }
    const auto& solved_initial = *two_phase->prediction.solved;
    const auto& solved_final = *two_phase->prediction.final_solved;
    const double hfeed = pure.state_at_pressure(260.0, 1.0e5, {1.0}, DensityBranch::Vapor).residual_enthalpy
        + 52.0 * (260.0 - 300.0) + 0.06 * (260.0 - 300.0) * (260.0 - 300.0);
    const double oracle_delta = (oracle_final - oracle_initial) / amount - hfeed;
    const double value_tolerance = 1.0e-9 * (std::abs(oracle_initial) + std::abs(oracle_final)) / amount + 1.0e-6;
    CHECK(agree("addition-f8/value", "delta_h_add", "J/mol", *two_phase->prediction.value, oracle_delta,
        value_tolerance, 0.0, fields({field("H_initial_J", oracle_initial), field("H_final_J", oracle_final),
            field("h_feed_J_per_mol", hfeed), field("delta_h_add_J_per_mol", *two_phase->prediction.value)}),
        "the total-system finite addition agrees with the independent endpoint oracle", "addition-f8"));
    const double h_initial_vapor = addition_oracle_phase_enthalpy(pure, solved_initial,
        *two_phase->observation.problem.reactions[0].correlation, 1);
    const double h_final_vapor = addition_oracle_phase_enthalpy(pure, solved_final,
        *two_phase->observation.final_problem.reactions[0].correlation, 1);
    const double vapor_contribution = std::abs(solved_final.amounts[1] * h_final_vapor
        - solved_initial.amounts[1] * h_initial_vapor) / amount;
    CHECK(truth("addition-f8/vapor_discrimination", "the changed vapor contribution exceeds 1e3 times the value tolerance",
        std::isfinite(vapor_contribution) && vapor_contribution > 1.0e3 * value_tolerance,
        fields({field("vapor_contribution_J_per_mol", vapor_contribution), field("delta_h_tolerance_J_per_mol", value_tolerance),
            field("nV_initial_mol", solved_initial.amounts[1]), field("nV_final_mol", solved_final.amounts[1]),
            field("hV_initial_J_per_mol", h_initial_vapor), field("hV_final_J_per_mol", h_final_vapor)}), "addition-f8"));

    auto failure = variants.front().observation;
    failure.temperature = 200.0;
    failure.pressure = 1.0e10;
    failure.branch = DensityBranch::Liquid;
    auto failed_feed = predict(model, failure);
    CHECK(truth("addition-f7/feed_density", "the declared feed root failure is typed and names feed",
        failed_feed.status == ObservationStatus::DensityUnavailable && !failed_feed.value
            && failed_feed.message.find("feed:") != std::string::npos,
        fields({field("status", observation_status_name(failed_feed.status)), field("message", failed_feed.message)}),
        "addition-f7"));
    auto missing_feed = variants.front().observation;
    missing_feed.composition = {0.0, 1.0};
    missing_feed.final_problem.feed.values = {2.0, 0.4};
    const auto missing = predict(model, missing_feed);
    CHECK(truth("addition-f7/feed_record", "the balanced feed species without an ideal record is typed and names feed",
        missing.status == ObservationStatus::ThermochemistryUnavailable && !missing.value
            && missing.message.find("feed:") != std::string::npos,
        fields({field("status", observation_status_name(missing.status)), field("message", missing.message),
            field("feed_B_mol", 1.0), field("n_add_mol", amount)}), "addition-f7"));
    ThermochemistryRecord no_records;
    no_records.standard_pressure = 1.0e5;
    no_records.component_intervals = {{}, {}};
    const HelmholtzModel undetermined_model(parameters, coordinates.active_parameters(), no_records);
    const auto undetermined = predict(undetermined_model, variants.front().observation);
    CHECK(truth("addition-f7/initial_record", "undetermined endpoint calorics are typed and name initial",
        undetermined.status == ObservationStatus::ThermochemistryUnavailable && !undetermined.value
            && undetermined.message.find("initial:") != std::string::npos,
        fields({field("status", observation_status_name(undetermined.status)), field("message", undetermined.message)}),
        "addition-f7"));
    auto mismatch = variants.front().observation;
    mismatch.final_problem.feed.values[0] += 0.01;
    const auto mismatch_message = rejection([&] { validate(mismatch, model); });
    CHECK(truth("addition-f7/material_balance", "a mismatched component balance is rejected by component name",
        mismatch_message.find("component") != std::string::npos,
        fields({field("message", mismatch_message), field("delta_A_mol", 0.01)}), "addition-f7"));
}

}  // namespace

// Native rejected holds: a cached point followed by no root, and a Born-free singular relation.
void check_hydration_failures() {
    auto parameters = properties_support::aqueous_electrolyte();
    Observation energy;
    energy.kind = ObservationKind::SolvationGibbsEnergy; energy.temperature = 298.15; energy.pressure = 1e5;
    energy.composition = {1, 0, 0}; energy.ion_indices = {1, 2}; energy.stoichiometry = {1, 1};
    energy.species = {1}; energy.branch = DensityBranch::Liquid; energy.observed = -381930;
    FitRequest request; request.parameters = parameters;
    request.coordinates.coordinates = {{ParameterFamily::Kij, {{1, 3}}, 0, .2, -1, 1}};
    request.observations = {energy}; request.holds = {{energy, 1, 8, 1e-6, 1e-10, 50}};
    const auto start = request.coordinates.start(parameters);
    FitResult failures;
    ObservationEvaluator good(request, failures, {0});
    ObservationEvaluation values;
    CHECK(good.evaluate(start, true, values));
    const auto saved = values.chart.physical;
    for (const bool singular : {false, true}) {
        auto bad = request;
        bad.holds[0].lower = saved.back() - 1e-6;
        bad.holds[0].upper = saved.back() + 1e-6;
        if (singular) {
            bad.parameters.electrolyte_form = ElectrolyteForm::DebyeHuckel; // manufactured H independent of d_Born
            HelmholtzModel constant(bad.parameters);
            bad.holds[0].reference.observed = *predict(constant, energy).value;
        }
        ObservationEvaluator evaluator(bad, failures, {0});
        auto rejected = start;
        if (!singular) { CHECK(evaluator.evaluate(start, true, values)); rejected[0] += 1.0; }
        const auto admitted = values.chart.physical;
        const bool accepted = evaluator.evaluate(rejected, true, values);
        const auto& cause = failures.failures.back().message;
        CHECK(truth(singular ? "hydration-failure/singular" : "hydration-failure/no-root",
            "A failed native held-coordinate trial publishes no diameter or Jacobian", !accepted && values.chart.physical.empty()
            && values.chart.map.empty() && values.gradient.empty() && cause.find(singular ? "singular" : "outside bracket") != std::string::npos,
            fields({field("target_J_per_mol", number(bad.holds[0].reference.observed)), field("cause", cause), field("optimizer", rejected)}), "issue-219"));
        CHECK((singular ? good : evaluator).evaluate(start, true, values));
        CHECK(truth(singular ? "hydration-failure/reuse-singular" : "hydration-failure/reuse-no-root",
            "A rejected trial preserves the admitted chart", values.chart.physical == admitted,
            "{}", "issue-219"));
        CHECK(good.evaluate(start, true, values));
    }
}

// B3 selects the actual native electrolyte evaluator and complete nonunit tied chart.
void check_minimax_consistency() {
    auto parameters = properties_support::aqueous_electrolyte();
    Observation energy;
    energy.kind = ObservationKind::SolvationGibbsEnergy; energy.temperature = 298.15; energy.pressure = 1e5;
    energy.composition = {1, 0, 0}; energy.ion_indices = {1, 2}; energy.stoichiometry = {1, 1};
    energy.species = {1}; energy.branch = DensityBranch::Liquid; energy.observed = -381930;
    FitRequest request; request.parameters = parameters;
    request.coordinates.coordinates = {{ParameterFamily::Kij, {{1, 3}}, 0, .2, -1, 1},
        {ParameterFamily::IonSuppressionCoefficient, {{1}, {2}}, 0, .5, 0, 20}};
    request.holds = {{energy, 1, 8, 1e-6, 1e-10, 50}};
    auto gamma = energy; gamma.kind = ObservationKind::MeanIonicActivityCoefficient;
    gamma.species.clear(); gamma.molality = .1; gamma.form = ResidualForm::LogRatio; gamma.observed = 1;
    auto phi = gamma; phi.kind = ObservationKind::OsmoticCoefficient;
    auto epsilon = gamma; epsilon.kind = ObservationKind::BulkRelativePermittivity;
    epsilon.form = ResidualForm::Difference;
    epsilon.composition = formula_composition({1, 0, 0}, parameters, {1, 2}, {1, 1}, .1);
    request.observations = {gamma, phi, epsilon, energy};
    const auto u = request.coordinates.start(parameters);
    request.observations = observations_at(request, u, {0, 0, .1, 0});
    request.weights = {1, 1, 0, 0};
    request.checks = {{"decades", {0, 1}, CheckMetric::AbsoluteLog10Ratio, .02, {}, false, false},
        {"permittivity", {2}, CheckMetric::AbsoluteDifference, .2, .2, false, false},
        {"energy", {3}, CheckMetric::AbsoluteDifference, 10, {}, false, false}};
    request.minimax = MinimaxRequest{{0, 1}, {.02, 1, 1000}, .02, 1e-10, 1e-6, .1};
    FitResult failures; ObservationEvaluator evaluator(request, failures, {0, 1, 2, 3}); ObservationEvaluation a;
    CHECK(evaluator.evaluate(u, true, a));
    const auto b = fit_minimax(request, u);
    CHECK(b.minimax && b.minimax->availability == "values_only");
    std::size_t quantity = 0;
    for (const auto& pair : {std::pair{b.predictions, a.predictions},
            std::pair{b.physical_jacobian, a.gradient}, std::pair{b.optimizer_jacobian, a.optimizer_gradient}}) {
        CHECK(pair.first.size() == pair.second.size());
        const auto columns = quantity == 0 ? 1 : quantity == 1 ? a.chart.physical.size() : u.size();
        for (std::size_t i = 0; i < pair.first.size(); ++i) {
            const bool energy_row = request.observations[i / columns].kind == ObservationKind::SolvationGibbsEnergy;
            const bool diameter_column = quantity == 1 && i % columns == columns - 1;
            const char* unit = energy_row ? (diameter_column ? "J/(mol Angstrom)" : "J/mol")
                : (diameter_column ? "1/Angstrom" : "1");
            const double tolerance = pair.second[i] == 0.0 ? 1e-12 : 1e-12 * std::abs(pair.second[i]);
            CHECK(agree("minimax-consistency/" + std::to_string(quantity) + "/" + std::to_string(i),
                quantity == 0 ? "prediction" : quantity == 1 ? "physical_column" : "optimizer_column",
                unit, pair.first[i], pair.second[i], tolerance, 0.0,
                "T=298.15 K,P=100000 Pa,m=.1 mol/kg", "B preserves A at the same held point", "r2-minimax"));
        }
        ++quantity;
    }
    const auto checks = evaluate_checks(request, b.parameters);
    for (std::size_t g = 0; g < checks.size(); ++g) {
        CHECK(agree("minimax-consistency/check/" + std::to_string(g), "maximum", checks[g].unit.c_str(),
            b.checks[g].maximum, checks[g].maximum, checks[g].maximum == 0 ? 1e-12 : 1e-12 * std::abs(checks[g].maximum), 0.0, "same stored diameters",
            "B final checks use the A check owner", "r2-minimax"));
        CHECK(agree("minimax-consistency/threshold/" + std::to_string(g), "maximum_threshold", checks[g].unit.c_str(),
            *b.minimax->frozen_checks[g].maximum_threshold, *checks[g].maximum_threshold,
            *checks[g].maximum_threshold == 0 ? 1e-12 : 1e-12 * std::abs(*checks[g].maximum_threshold), 0.0,
            "same explicit thresholds", "B freezes declared unweighted checks", "r2-minimax"));
    }
}

// Independent bisection uses the Properties value owner, with no chart or cached Newton root.
int electrolyte_reference_file(const char* filename) {
    std::ifstream in(filename);
    epcsaft::ParameterDocument doc;
    std::size_t n, sites, entries;
    in >> n >> sites >> entries;
    doc.component_ids.resize(n); doc.site_component.resize(sites); doc.site_multiplicity.resize(sites);
    for (auto& id : doc.component_ids) in >> std::quoted(id);
    for (std::size_t i = 0; i < sites; ++i) in >> doc.site_component[i] >> doc.site_multiplicity[i];
    in >> std::quoted(doc.base_choice) >> std::quoted(doc.association_choice) >> std::quoted(doc.electrolyte_choice)
       >> std::quoted(doc.ionic_dispersion_choice) >> std::quoted(doc.permittivity_choice) >> std::quoted(doc.polar_choice);
    for (std::size_t i = 0; i < entries; ++i) {
        epcsaft::ParameterEntry entry;
        std::size_t width;
        in >> std::quoted(entry.identity) >> std::quoted(entry.field) >> width;
        entry.index.resize(width);
        for (auto& index : entry.index) in >> index;
        in >> entry.magnitude >> std::quoted(entry.unit);
        doc.entries.push_back(entry);
    }
    auto parameters = epcsaft::build_parameters(doc);
    std::size_t count;
    in >> count;
    std::cout << "ion,stored_diameter,independent_diameter,energy_error,bracket_width\n" << std::setprecision(17);
    for (std::size_t h = 0; h < count; ++h) {
        std::size_t ion, counter, solvent;
        double t, p, target;
        in >> ion >> counter >> solvent >> t >> p >> target;
        if (!in) throw std::runtime_error("incomplete reference inputs");
        const double stored = parameters.born_diameter_angstrom.at(ion);
        std::vector<double> x(n, 0.0); x.at(solvent) = 1.0;
        const auto energy = [&](double d) {
            parameters.born_diameter_angstrom[ion] = d;
            epcsaft::HelmholtzModel model(parameters);
            return epcsaft::GAS_CONSTANT_J_PER_MOL_K * t
                * epcsaft::infinite_dilution_reference(model, t, p, x, {ion, counter}, {1, 1}).log_fugacity[0] - target;
        };
        double lower = 1.0, upper = 8.0, low = energy(lower), mid = 0.0, error = 0.0;
        if (std::signbit(low) == std::signbit(energy(upper))) throw std::runtime_error("independent root outside bracket");
        for (int iteration = 0; iteration < 100; ++iteration) {
            mid = 0.5 * (lower + upper); error = energy(mid);
            if (std::abs(error) <= 1e-7 && upper - lower <= 1e-12) break;
            if (std::signbit(error) == std::signbit(low)) { lower = mid; low = error; } else upper = mid;
        }
        std::cout << ion << ',' << stored << ',' << mid << ',' << error << ',' << upper-lower << '\n';
        parameters.born_diameter_angstrom[ion] = stored;
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "--electrolyte-reference") return electrolyte_reference_file(argv[2]);
    return recording_support::run(argc, argv, {
        {{"native.regression.hydration-failures", "Dependent_hydration_failure", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "issue-219", "hydration-failure/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "hydration-failures"}, check_hydration_failures, "R1"},
        {{"native.regression.addition-enthalpy", "Finite_addition_enthalpy_and_exact_actions", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "issue-210", "addition-f3/*,addition-f4/*,addition-f5/*,addition-f6/*,addition-f7/*,addition-f8/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "addition-enthalpy"}, check_addition_enthalpy, "addition-enthalpy"},
        {{"native.regression.ln-k", "Source_log_equilibrium_constant", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "ln-k/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "ln-k"}, check_log_constant, "R2"},
        {{"native.regression.r1-f1", "Regression_R1_F1", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f1/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F1", "r1-f1/*:42"}, check_f1, "R1"},
        {{"native.regression.r1-f2", "Regression_R1_F2", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f2/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F2", "r1-f2/*:18"}, check_f2, "R1"},
        {{"native.regression.r1-f3", "Regression_R1_F3_F4", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f3/*,r1-f4/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F3", "r1-f3/*:5,r1-f4/*:3"}, check_f3_f4, "R1"},
        {{"native.regression.r1-f5", "Regression_R1_F5_F6", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f5/*,r1-f6/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F5", "r1-f5/*:5,r1-f6/*:3"}, check_f5_f6, "R1"},
        {{"native.regression.r1-f7", "Regression_R1_F7_F8", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f7/*,r1-f8/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F7", "r1-f7/*:5,r1-f8/*:12"}, check_f7_f8, "R1"},
        {{"native.regression.r1-f9", "Regression_R1_F9_F10", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R1", "r1-f9/*,r1-f10/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R1-F9", "r1-f9/*:4,r1-f10/*:13"}, check_f9_f10, "R1"},
        {{"native.regression.d86", "Regression_D86_density_actions", "Numerical", "Properties", "Explicit", "Numerical", nullptr, "D86", "d86/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "D86", "d86/*:25"}, check_d86, "R1"},
        {{"native.regression.r2-j1", "Regression_R2_J1", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "r2-j1/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R2-J1", "r2-j1/*:84"}, check_j1, "R2"},
        {{"native.regression.r2-j2", "Regression_R2_J2", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "r2-j2/*,r4/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R2-J2", "r2-j2/*:15"}, check_j2, "R2"},
        {{"native.regression.r2-minimax-consistency", "Minimax_same_point_columns", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "minimax-consistency/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "minimax-consistency"}, check_minimax_consistency, "R2"},
        {{"native.regression.r2-d7-minimax-analytic", "Regression_R2_D7_Analytic_Minimax", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "r2-d7-b1/*,r2-d7-b2/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R2-D7", "r2-d7/*"}, check_d7_minimax, "R2"},
        {{"native.regression.reaction-derivatives", "Reaction_constant_derivatives", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "reaction-a*", "engine/tests/native/regression_checks.cpp", "regression_checks", "reaction-derivatives"}, check_reaction_derivatives, "R2"},
        {{"native.regression.reaction-recovery", "Reaction_constant_recovery", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "reaction-b*", "engine/tests/native/regression_checks.cpp", "regression_checks", "reaction-recovery"}, check_reaction_recovery, "R2"},
        {{"native.regression.reaction-rejections", "Reaction_constant_admission", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "reaction-rejection/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "reaction-rejections"}, check_reaction_rejections, "R2"},
        {{"native.regression.cp-derivatives", "Isobaric_heat_capacity_derivatives", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "cp-a/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "cp-derivatives"}, check_cp_derivatives, "R2"},
        {{"native.regression.cp-recovery", "Isobaric_heat_capacity_recovery", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R2", "cp-b/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "cp-recovery"}, check_cp_recovery, "R2"},
        {{"native.regression.r3", "Regression_R3", "Numerical", "Regression", "Explicit", "Numerical", nullptr, "R3", "r3/*", "engine/tests/native/regression_checks.cpp", "regression_checks", "R3", "r3/*:22"}, check_r3, "R2"},
    });
}
