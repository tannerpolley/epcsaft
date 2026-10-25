// A Python-free caller of the native thermodynamic interface. It constructs physical inputs,
// evaluates full and selected properties, solves declared equilibria, requests exact actions and
// inspects partial and failed results, using only the two headers below. The associating binary
// repeats the constants of the retained native fixture properties_support::associating_binary()
// (an amine-water-like pair) and the ideal-gas record is manufactured: they exercise the calling
// surface, units and identities, with no source or validity claim. SI molar units throughout.
// Build, run and debug: engine/README.md, "C++ usage and debugging".
#include "equilibrium/problem.hpp"
#include "surface/parameter_builder.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

using namespace epcsaft;
using namespace epcsaft::equilibrium;

namespace {

int failures = 0;

void expect(bool verdict, const char* claim) {
    std::printf("  %s  %s\n", verdict ? "ok  " : "FAIL", claim);
    failures += verdict ? 0 : 1;
}

// The accepted symmetric criterion for values of one owner reached by two routes.
bool same(double a, double b, double tolerance = 1.0e-12) {
    return std::abs(a - b) <= std::max(tolerance, tolerance * std::max(std::abs(a), std::abs(b)));
}

// Full-precision machine-readable rows: the C++ reference of the installed-wheel Python checks.
void row(const std::string& id, const char* unit, double value) { std::printf("EPCSAFT_ROW\t%s\t%s\t%.17g\n", id.c_str(), unit, value); }
void rows(const std::string& id, const char* unit, const std::vector<double>& values) {
    for (std::size_t i = 0; i < values.size(); ++i) row(id + "/" + std::to_string(i), unit, values[i]);
}

// Species index = position in component_ids, in every vector below. Every entry states its unit.
ParameterDocument associating_binary() {
    ParameterDocument doc;
    doc.component_ids = {"associating-a", "associating-b"};
    doc.site_component = {0, 0, 1, 1};  // two sites on each species
    doc.site_multiplicity = {2, 1, 2, 2};
    doc.base_choice = "pc-saft";
    doc.association_choice = "general-site";
    doc.electrolyte_choice = "none";
    doc.permittivity_choice = "none";
    const double segment[2] = {4.5354, 2.1945}, sigma[2] = {2.6019, 2.229}, epsilon[2] = {204.0438, 141.66};
    const double molar_mass[2] = {61.08, 18.01528};
    for (std::size_t i = 0; i < 2; ++i) {
        doc.entries.push_back({"molar mass", "molar_mass", {i}, molar_mass[i], "gram / mole"});
        doc.entries.push_back({"charge", "charge_number", {i}, 0.0, "elementary-charge"});
        doc.entries.push_back({"segments", "segment_count", {i}, segment[i], "dimensionless"});
        doc.entries.push_back({"diameter", "sigma_constant", {i}, sigma[i], "angstrom"});
        doc.entries.push_back({"dispersion", "dispersion_energy_over_k", {i}, epsilon[i], "kelvin"});
    }
    doc.entries.push_back({"k_ij", "k_ij", {0, 1}, -0.052, "dimensionless"});
    const std::size_t first[4] = {0, 0, 1, 2}, second[4] = {1, 3, 2, 3};  // site pairs
    const double energy[4] = {2383.4744, 2093.8222, 2093.8222, 1804.17};
    const double volume[4] = {0.118488, 0.1540467874409226, 0.1540467874409226, 0.2039};
    for (std::size_t k = 0; k < 4; ++k) {
        doc.entries.push_back({"association energy", "association_energy_over_k", {first[k], second[k]}, energy[k], "kelvin"});
        doc.entries.push_back({"association volume", "association_volume", {first[k], second[k]}, volume[k], "dimensionless"});
    }
    return doc;
}

// Constant ideal-gas cp with an enthalpy and entropy anchor at 300 K and 1 bar, on [250, 450] K.
ThermochemistryRecord ideal_record() {
    ThermochemistryRecord record;
    record.standard_pressure = 1.0e5;
    const double cp[2] = {80.0, 35.0}, h[2] = {1000.0, 2000.0}, s[2] = {50.0, 70.0};
    for (std::size_t i = 0; i < 2; ++i) {
        record.component_intervals.push_back({IdealInterval{250.0, 450.0, true, true,
            IdealCorrelation{IdealPolynomial{{cp[i], 0.0}, 300.0, h[i], s[i]}, 1.0e5}}});
    }
    return record;
}

SolvedStateObservable observable(SolvedStateObservableKind kind, std::size_t phase, std::size_t component = 0,
                                 PropertyObservable property = PropertyObservable::TotalEnthalpy) {
    SolvedStateObservable result;
    result.kind = kind; result.phase = phase; result.component = component; result.property = property;
    return result;
}

}  // namespace

int main() {
    const double T = 353.15, R = GAS_CONSTANT_J_PER_MOL_K;
    const std::vector<double> feed = {0.25, 0.75};  // mol, absolute

    std::printf("1. physical inputs -> one reusable model\n");
    const HelmholtzModel model(build_parameters(associating_binary()), {}, ideal_record());
    const HelmholtzModel bare(build_parameters(associating_binary()));  // no thermochemistry

    std::printf("2. full state in the T-P chart (liquid branch) and in the T-rho chart\n");
    const State at_pressure = model.state_at_pressure(T, 1.0e5, feed, DensityBranch::Liquid);
    const State at_density = model.state(T, at_pressure.molar_density, feed);
    std::printf("  rho = %.6f mol/m^3, P = %.3f Pa, h_res = %.4f J/mol, cp_total = %.4f J/(mol K)\n", at_pressure.molar_density,
                at_pressure.pressure, at_pressure.residual_enthalpy, at_pressure.total_isobaric_heat_capacity.value_or(std::nan("")));
    const core::Convergence& root = *at_pressure.density_convergence;  // the density owner's own verdict
    std::printf("  density root: %s, residual %.1e against %s = %.1e; P(state) - P = %.1e Pa\n", core::stop_name(root.stop), root.residual,
                root.criterion->name, root.criterion->value, at_pressure.pressure - 1.0e5);
    rows("state/calorics", "J/mol or J/(mol K)", {*at_pressure.ideal_enthalpy, *at_pressure.ideal_internal_energy, *at_pressure.ideal_entropy,
        *at_pressure.ideal_isobaric_heat_capacity, *at_pressure.ideal_isochoric_heat_capacity, *at_pressure.total_enthalpy, *at_pressure.total_internal_energy,
        *at_pressure.total_entropy, *at_pressure.total_isobaric_heat_capacity, *at_pressure.total_isochoric_heat_capacity});  // ideal, then total: h, u, s, cp, cv
    rows("state/root", "mixed", {at_pressure.molar_density, at_pressure.pressure, at_pressure.pressure_density_derivative, root.residual,
                                 static_cast<double>(at_pressure.stable_root_count.value_or(0))});
    expect(at_pressure.density_branch == DensityBranch::Liquid && root.residual <= root.criterion->value, "the liquid density root meets its owner's criterion");
    expect(same(at_density.pressure, at_pressure.pressure) && same(at_density.residual_enthalpy, at_pressure.residual_enthalpy), "both charts give the same state");

    std::printf("3. selected observables with exact actions, first and mixed second order\n");
    const std::vector<PropertyObservable> selected = {PropertyObservable::ResidualEnthalpy, PropertyObservable::TotalIsobaricHeatCapacity};
    const PressureActionDirection heat{1.0, 0.0, {0.0, 0.0}, {}}, squeeze{0.0, 1.0e5, {0.0, 0.0}, {}};
    const auto first = model.pressure_actions(T, 1.0e5, feed, DensityBranch::Liquid, selected, {heat});
    const auto mixed = model.pressure_actions(T, 1.0e5, feed, DensityBranch::Liquid, selected, {heat, squeeze});
    const DensityActionDirection heat_at_density{1.0, 0.0, {0.0, 0.0}, {}}, compress{0.0, 1.0, {0.0, 0.0}, {}};
    const auto fixed_density = model.density_actions(T, at_pressure.molar_density, feed, selected, {heat_at_density});
    const auto mixed_density = model.density_actions(T, at_pressure.molar_density, feed, selected, {heat_at_density, compress});
    std::printf("  d h_res/dT at P = %.6f J/(mol K); d2 h_res/(dT dP) * 1 bar = %.6e J/(mol K)\n", first.results[0].action.value_or(std::nan("")),
                mixed.results[0].action.value_or(std::nan("")));
    std::printf("  d h_res/dT at rho = %.6f J/(mol K); d2 h_res/(dT d rho) * 1 mol/m^3 = %.6e J/(mol K)\n", fixed_density.results[0].action.value_or(std::nan("")),
                mixed_density.results[0].action.value_or(std::nan("")));
    expect(first.results[0].status == PropertyActionStatus::Available && same(*first.results[0].value, at_pressure.residual_enthalpy)
           && same(*first.results[1].value, *at_pressure.total_isobaric_heat_capacity)
           && same(*mixed_density.results[0].value, at_pressure.residual_enthalpy), "selected central values equal the full state in both charts");
    // The T-direction Taylor chain cancels O(1e4) J/(mol K) terms to ~40; its rounding differs from the direct Hessian cp by up to ~2e-12.
    expect(same(*first.results[0].action, *at_pressure.residual_isobaric_heat_capacity, 1.0e-11), "d h_res/dT at constant P is the residual cp");

    std::printf("4. declared equilibria: bubble point (pinned incipient vapor, free P), then a two-phase flash\n");
    EquilibriumProblem bubble;
    PhaseSpec liquid, vapor;
    liquid.name = "L"; liquid.kind = PhaseKind::Liquid;
    vapor.name = "V"; vapor.kind = PhaseKind::Vapor; vapor.amount = AmountSpec::pinned(0.0);  // phase fraction of the feed
    bubble.phases = {liquid, vapor};
    bubble.temperature = FreeValue::fixed(T);
    bubble.pressure = FreeValue::variable(3.0e4);
    bubble.feed.kind = FeedKind::Amounts;
    bubble.feed.values = feed;
    const EquilibriumResult boiling = solve_equilibrium(model, bubble);
    EquilibriumProblem flash = bubble;
    flash.phases[1].amount = AmountSpec::free();
    flash.pressure = FreeValue::fixed(1.8e4);
    const EquilibriumResult split = solve_equilibrium(model, flash);
    expect(boiling.success && split.success && *split.requested_tolerance_met, "both solves meet the caller row tolerance");
    std::printf("  bubble P = %.2f Pa; flash: N_L = %.6f mol, N_V = %.6f mol, x_L = (%.6f, %.6f), y_V = (%.6f, %.6f), |r| = %.2e\n", *boiling.pressure,
                split.amounts[0], split.amounts[1], split.mole_fractions[0], split.mole_fractions[1], split.mole_fractions[2], split.mole_fractions[3], split.residual_norm());
    std::printf("  rho_L = %.4f mol/m^3, rho_V = %.6f mol/m^3\n", split.molar_densities[0], split.molar_densities[1]);
    row("bubble/pressure", "Pa", *boiling.pressure);
    rows("flash/amounts", "mol", split.amounts);
    rows("flash/mole_fractions", "1", split.mole_fractions);
    rows("flash/multipliers", "1", split.multipliers);  // reduced-energy balance coordinates
    rows("flash/constraint_multipliers", "1", split.constraint_multipliers);  // extensive-objective duals at returned state
    expect(split.molar_densities[0] > split.molar_densities[1], "the declared liquid is the dense phase and the declared vapor the dilute one");
    double balance = std::abs(split.amounts[0] + split.amounts[1] - 1.0);
    for (std::size_t i = 0; i < 2; ++i) balance = std::max(balance, std::abs(split.amounts[0] * split.mole_fractions[i] + split.amounts[1] * split.mole_fractions[2 + i] - feed[i]));
    std::printf("  largest balance defect = %.1e mol\n", balance);
    expect(balance <= 1.0e-12, "absolute phase amounts close the species balances");

    std::printf("5. exact solved-state actions and the caller-formed fugacity\n");
    const CompiledSystem compiled = compile_problem(model, flash);  // the model must outlive it
    SolvedStateActionRequest request;
    using Kind = SolvedStateObservableKind;
    request.observables = {observable(Kind::PhaseComponentAmount, 0, 1), observable(Kind::PhaseLogActivity, 0, 1),
                           observable(Kind::PhaseMoleFraction, 0, 1), observable(Kind::PhaseProperty, 0), observable(Kind::PhaseMolarDensity, 0)};
    SolvedStateActionDirection warm, scale;
    warm.temperature_k = 1.0; warm.feed_amounts_mol = {0.0, 0.0};
    scale.feed_amounts_mol = feed;
    request.directions = {warm};
    const auto thermal = solved_state_actions(compiled, split, request);
    request.directions = {scale};
    const auto radial = solved_state_actions(compiled, split, request);
    std::printf("  along the two-phase state at fixed P and feed: d x_L[b]/dT = %.6e 1/K, d h_L/dT = %.4f J/(mol K) (composition shifts),\n"
                "  d rho_L/dT = %.4f mol/(m^3 K), condition = %.2e\n", thermal.results[2].action.value_or(std::nan("")), thermal.results[3].action.value_or(std::nan("")),
                thermal.results[4].action.value_or(std::nan("")), thermal.diagnostics.condition_estimate.value_or(std::nan("")));
    for (std::size_t i = 0; i < request.observables.size(); ++i) {
        row("solved/thermal/" + std::to_string(i), thermal.results[i].unit.c_str(), *thermal.results[i].action);
        row("solved/radial/" + std::to_string(i), radial.results[i].unit.c_str(), *radial.results[i].action);
    }
    row("solved/condition", "1", *thermal.diagnostics.condition_estimate);
    expect(same(*thermal.results[4].value, split.molar_densities[0]), "the solved-state density value is the result's liquid density");
    expect(same(*radial.results[0].action, *radial.results[0].value, 5.0e-6), "Euler: D n_i[n_feed] = n_i");
    expect(std::abs(*radial.results[1].action) <= 1.0e-12 * std::max(1.0, std::abs(*radial.results[1].value)), "Euler: D ln a_i[n_feed] = 0");
    const std::vector<double> x_liquid = {split.mole_fractions[0], split.mole_fractions[1]};
    const State solved = model.state(*split.temperature, split.molar_densities[0], x_liquid);
    const double fugacity = std::exp(*radial.results[1].value) * R * *split.temperature * 1.0;  // a_i R T (1 mol/m^3)
    const double from_state = x_liquid[1] * std::exp((*solved.log_fugacity_coefficient)[1].value()) * solved.pressure;
    std::printf("  f_b = %.6f Pa (from ln a), %.6f Pa (x phi P of the state), sum x - 1 = %.1e\n", fugacity, from_state, x_liquid[0] + x_liquid[1] - 1.0);
    expect(same(fugacity, from_state), "f_i = a_i R T rho0 equals x_i phi_i P");

    std::printf("6. a partial result: no thermochemistry\n");
    const EquilibriumResult bare_split = solve_equilibrium(bare, flash);
    const auto partial = solved_state_actions(compile_problem(bare, flash), bare_split, request);
    std::printf("  total enthalpy: %s; mole fraction: %s\n", solved_state_status_name(partial.results[3].status), solved_state_status_name(partial.results[2].status));
    expect(partial.results[3].status == SolvedStateActionStatus::MissingThermochemistry && !partial.results[3].value.has_value()
           && partial.results[2].status == SolvedStateActionStatus::Available && partial.results[2].action.has_value(), "the caloric row is typed missing and its siblings stay available");

    std::printf("7. failures a caller can inspect\n");
    std::string message;
    EquilibriumProblem unnormalized = flash;
    unnormalized.feed.kind = FeedKind::MoleFractions;
    unnormalized.feed.values = {0.25, 0.74};
    try { (void)solve_equilibrium(model, unnormalized); } catch (const std::invalid_argument& error) { message = error.what(); }
    std::printf("  admission: %s\n", message.c_str());
    expect(message.find("mole fractions must sum to one") != std::string::npos, "a fraction feed that does not sum to one is rejected, not renormalized");
    EquilibriumProblem unattainable = bubble;
    unattainable.solver_options.values[epcsaft::policy::row] = 1.0e-18;
    const EquilibriumResult floor = solve_equilibrium(model, unattainable);
    std::printf("  unmet tolerance: success = %d, stop = %s, achieved |r| = %.2e\n  %s\n", floor.success, core::stop_name(floor.convergence.stop), floor.convergence.residual, floor.message.c_str());
    row("floor/residual", "1", floor.convergence.residual);
    expect(!floor.success && floor.requested_tolerance_met.has_value() && !*floor.requested_tolerance_met && std::isfinite(floor.convergence.residual),
           "an unattainable tolerance is reported as unmet with the achieved residual");

    std::printf("%s\n", failures == 0 ? "surface example passed" : "surface example FAILED");
    return failures == 0 ? 0 : 1;
}
