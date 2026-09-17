#include "equilibrium/problem.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "core/physical_constants.hpp"
#include "equilibrium_support.hpp"
#include "recording_support.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
namespace {

const char* phase_set_name(equilibrium_support::PhaseSet value) {
    switch (value) {
    case equilibrium_support::PhaseSet::OneLiquid: return "one-liquid";
    case equilibrium_support::PhaseSet::OneVapor: return "one-vapor";
    case equilibrium_support::PhaseSet::Two: return "two";
    case equilibrium_support::PhaseSet::TwoVaporLiquid: return "two-vapor-liquid";
    case equilibrium_support::PhaseSet::Three: return "three";
    case equilibrium_support::PhaseSet::ThreeVaporLiquid: return "three-vapor-liquid";
    }
    return "unknown";
}

const char* amount_pattern_name(equilibrium_support::AmountPattern value) {
    switch (value) {
    case equilibrium_support::AmountPattern::AllFree: return "all-free";
    case equilibrium_support::AmountPattern::FirstPinned: return "first-pinned";
    case equilibrium_support::AmountPattern::LastPinned: return "last-pinned";
    case equilibrium_support::AmountPattern::SplitPinned: return "split-pinned";
    case equilibrium_support::AmountPattern::AllPinned: return "all-pinned";
    }
    return "unknown";
}

const char* intensive_pattern_name(equilibrium_support::IntensivePattern value) {
    switch (value) {
    case equilibrium_support::IntensivePattern::BothFixed: return "both-fixed";
    case equilibrium_support::IntensivePattern::TemperatureFree: return "temperature-free";
    case equilibrium_support::IntensivePattern::PressureFree: return "pressure-free";
    case equilibrium_support::IntensivePattern::BothFree: return "both-free";
    }
    return "unknown";
}

const char* chemistry_name(equilibrium_support::Chemistry value) {
    switch (value) {
    case equilibrium_support::Chemistry::Neutral: return "neutral";
    case equilibrium_support::Chemistry::Isomerisation: return "isomerisation";
    case equilibrium_support::Chemistry::Electrolyte: return "electrolyte";
    case equilibrium_support::Chemistry::TwoReactions: return "two-reactions";
    }
    return "unknown";
}

void check_c1() {
    constexpr std::size_t expected_cases = 6 * 5 * 4 * 4;
    recording_support::expect_group(
        recording_support::active_test_id(), "c1-structural", "structural/",
        expected_cases,
        recording_support::fields({
            recording_support::field("phase_sets", equilibrium_support::phase_sets().size()),
            recording_support::field("amount_patterns", equilibrium_support::amount_patterns().size()),
            recording_support::field("intensive_patterns", equilibrium_support::intensive_patterns().size()),
            recording_support::field("chemistries", equilibrium_support::chemistry_values().size()),
        }), "phase_set,amount,intensive,chemistry");
    for (const auto phase_set : equilibrium_support::phase_sets()) {
        for (const auto amount_pattern : equilibrium_support::amount_patterns()) {
            for (const auto intensive_pattern : equilibrium_support::intensive_patterns()) {
                for (const auto chemistry : equilibrium_support::chemistry_values()) {
                    auto problem = equilibrium_support::lattice_problem(
                        phase_set, amount_pattern, intensive_pattern, chemistry);
                    const auto parameters = chemistry == equilibrium_support::Chemistry::Neutral ? properties_support::neutral_binary() : chemistry == equilibrium_support::Chemistry::Electrolyte ? properties_support::aqueous_electrolyte() : equilibrium_support::synthetic_reactive_ternary();
                    epcsaft::HelmholtzModel model(parameters);
                    const int k = equilibrium_support::pinned_phases(phase_set, amount_pattern);
                    const int f = equilibrium_support::free_intensives(intensive_pattern);
                    const int reaction_count = static_cast<int>(epcsaft::equilibrium::reaction_rank(problem.feed.values.size(), problem.reactions));
                    const int rank = static_cast<int>(problem.feed.values.size()) - reaction_count;
                    bool accepted = false;
                    std::size_t jacobian_rank = 0, generic_nullity = 0, row_count = 0, unknown_count = 0;
                    std::size_t support_count = 0, charge_rows = 0;
                    std::string solver_status = "not-compiled";
                    try {
                        const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
                        accepted = compiled.rows_minus_unknowns == 0;
                        jacobian_rank = compiled.jacobian_rank;
                        generic_nullity = compiled.generic_nullity;
                        row_count = compiled.rows.size();
                        unknown_count = compiled.unknown_count();
                        charge_rows = compiled.charge_row_count;
                        for (const auto& item : compiled.problem.phases) support_count += item.support.size();
                        // Structural degeneracy is reported before any nonlinear work; every
                        // lattice all-pinned cell is rank deficient and returns the exact NOT_RUN record.
                        if (!compiled.feasible) {
                            const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
                            solver_status = result.solver_status;
                            CHECK(!result.success && result.solver_status == "NOT_RUN"
                                  && result.residuals.empty() && !result.has_state()
                                  && result.convergence.stop == epcsaft::core::Stop::infeasible
                                  && result.generic_nullity == generic_nullity);
                        } else {
                            solver_status = "full-rank";
                        }
                    } catch (const std::invalid_argument& error) {
                        accepted = false;
                        solver_status = error.what();
                    }
                    const int phases = equilibrium_support::pinned_phases(phase_set, equilibrium_support::AmountPattern::AllFree) == 0 ? (phase_set == equilibrium_support::PhaseSet::OneLiquid || phase_set == equilibrium_support::PhaseSet::OneVapor ? 1 : phase_set == equilibrium_support::PhaseSet::Two || phase_set == equilibrium_support::PhaseSet::TwoVaporLiquid ? 2 : 3) : 1;
                    const bool structurally_square = f == k && phases - k <= rank;
                    const auto case_id = std::string(phase_set_name(phase_set)) + "/"
                        + amount_pattern_name(amount_pattern) + "/"
                        + intensive_pattern_name(intensive_pattern) + "/"
                        + chemistry_name(chemistry);
                    recording_support::compare(
                        (std::string("structural/") + case_id).c_str(),
                        "structurally_square", "bool", accepted ? 1.0 : 0.0,
                        structurally_square ? 1.0 : 0.0,
                        recording_support::fields({
                            recording_support::field("phase_set", phase_set_name(phase_set)),
                            recording_support::field("amount_pattern", amount_pattern_name(amount_pattern)),
                            recording_support::field("intensive_pattern", intensive_pattern_name(intensive_pattern)),
                            recording_support::field("chemistry", chemistry_name(chemistry)),
                            recording_support::field("free_intensives", f),
                            recording_support::field("pinned_phases", k),
                            recording_support::field("phases", phases),
                            recording_support::field("rank", rank),
                            recording_support::field("accepted", accepted ? 1U : 0U),
                            recording_support::field("support_count", support_count),
                            recording_support::field("charge_rows", charge_rows),
                            recording_support::field("rows", row_count),
                            recording_support::field("unknowns", unknown_count),
                            recording_support::field("jacobian_rank", jacobian_rank),
                            recording_support::field("generic_nullity", generic_nullity),
                            recording_support::field("solve_outcome", solver_status),
                        }), "Relationship", "rank-count",
                        "accepted iff structural rows are square", "absolute", 0.0, 0.0,
                        accepted == structurally_square, "c1-structural", case_id.c_str(),
                        "accepted", "phase_set,amount,intensive,chemistry", "accepted=structural",
                        "case");
                    CHECK((f - k) == 0 || !accepted);
                    CHECK(accepted == structurally_square);
                    CHECK(!(accepted && k == phases) || generic_nullity > 0);
                }
            }
        }
    }
}
void check_c4() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto amounts = equilibrium_support::flash(2);
    amounts.feed.values = {0.5, 0.5};
    auto fractions = amounts;
    fractions.feed.kind = epcsaft::equilibrium::FeedKind::MoleFractions;
    auto mass = amounts;
    mass.feed.kind = epcsaft::equilibrium::FeedKind::MassFractions;
    mass.feed.values = {0.5 * model.parameters().molar_mass_kg_per_mol[0],
                        0.5 * model.parameters().molar_mass_kg_per_mol[1]};
    auto molality = amounts;
    molality.feed.kind = epcsaft::equilibrium::FeedKind::Molalities;
    molality.feed.values = {0.0, 1.0 / model.parameters().molar_mass_kg_per_mol[0]};
    molality.feed.solvent = {0};
    const auto reference = epcsaft::equilibrium::compile_problem(model, amounts);
    recording_support::expect_group(
        recording_support::active_test_id(), "c4-feed", "feed-composition/", 6,
        recording_support::fields({
            recording_support::field("reference_feed", std::vector<double>{0.5, 0.5}),
            recording_support::field("feed_forms", 3U),
        }), "basis,index");
    std::size_t basis_index = 0;
    for (const auto& item : {fractions, mass, molality}) {
        const auto compiled = epcsaft::equilibrium::compile_problem(model, item);
        CHECK(compiled.feed_composition.size() == reference.feed_composition.size());
        for (std::size_t i = 0; i < reference.feed_composition.size(); ++i) {
            const double difference = compiled.feed_composition[i]
                - reference.feed_composition[i];
            const auto point = std::to_string(basis_index) + "/" + std::to_string(i);
            recording_support::compare(
                (std::string("feed-composition/") + point).c_str(),
                "feed_composition", "1", compiled.feed_composition[i],
                reference.feed_composition[i],
                recording_support::fields({
                    recording_support::field("T", item.temperature.value),
                    recording_support::field("P", item.pressure.value),
                    recording_support::field("feed", item.feed.values),
                    recording_support::field("feed_kind", epcsaft::equilibrium::feed_kind_name(item.feed.kind)),
                    recording_support::field("reference_feed", reference.feed_composition),
                    recording_support::field("feed_basis", basis_index),
                    recording_support::field("component", i),
                }), "Relationship", "feed-normalization",
                "compiled feed composition matches amount basis", "strict_absolute",
                1e-12, 0.0, std::abs(difference) < 1e-12,
                "c4-feed", point.c_str(), "x_i", "T,P,feed,basis,index",
                "x_i^{compiled}=x_i^{reference}", "basis,index");
            CHECK(std::abs(compiled.feed_composition[i] - reference.feed_composition[i]) < 1e-12);
        }
        CHECK(compiled.rows.size() == reference.rows.size());
        for (std::size_t row = 0; row < reference.rows.size(); ++row) {
            CHECK(compiled.rows[row].kind == reference.rows[row].kind);
            CHECK(compiled.rows[row].phase == reference.rows[row].phase);
            CHECK(compiled.rows[row].index == reference.rows[row].index);
            CHECK(std::string(compiled.rows[row].eqid) == std::string(reference.rows[row].eqid));
        }
        ++basis_index;
    }
}
void check_c5() {
    auto parameters = properties_support::parameters(3);
    parameters.segment_count = {1.0, 1.0, 1.0};
    parameters.sigma_constant = {3.0, 3.1, 3.2};
    parameters.dispersion_energy_over_k = {100.0, 110.0, 120.0};
    parameters.molar_mass_kg_per_mol = {0.01, 0.02, 0.03};
    epcsaft::HelmholtzModel model(parameters);
    auto problem = equilibrium_support::flash(3);
    problem.feed.values = {1.0, 1.0, 1.0};
    problem.reactions = {{"exchange", {1.0, -1.0, 0.0}, 1.0, std::nullopt}};
    const double temperature = 300.0;
    const std::size_t components = problem.feed.values.size();
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    recording_support::compare(
        "rank", "reaction_rank", "count", static_cast<double>(compiled.rank), 2.0,
        recording_support::fields({recording_support::field("T", temperature),
                                   recording_support::field("components", components),
                                   recording_support::field("reactions", problem.reactions.size()),
                                   recording_support::field("feed", problem.feed.values)}),
        "Relationship", "stoichiometric-rank", "independent reaction basis has rank two",
        "reference_or_one", 0.0, 0.0, compiled.rank == 2,
        "c5-rank", "basis", "rank", "T,components,reactions,feed", "rank(B)=2", "rank");
    const std::vector<epcsaft::equilibrium::ReactionSpec> basis = {
        {"r1", {1.0, -1.0, 0.0}, 1.0, std::nullopt},
        {"r2", {0.0, 1.0, -1.0}, 1.0, std::nullopt},
    };
    const std::vector<epcsaft::equilibrium::ReactionSpec> mixed = {
        {"sum", {1.0, 0.0, -1.0}, 1.0, std::nullopt},
        {"difference", {1.0, -2.0, 1.0}, 1.0, std::nullopt},
    };
    const auto basis_coordinates = epcsaft::equilibrium::balance_basis(components, basis);
    const auto mixed_coordinates = epcsaft::equilibrium::balance_basis(components, mixed);
    const auto standard = epcsaft::equilibrium::standard_state(temperature, components, basis, {});
    const auto basis_affinities = epcsaft::equilibrium::reaction_affinities(basis, standard);
    const auto mixed_affinities = epcsaft::equilibrium::reaction_affinities(mixed, standard);
    recording_support::expect_group(
        recording_support::active_test_id(), "c5-basis", "basis-entry/", 3,
        recording_support::fields({
            recording_support::field("T", temperature),
            recording_support::field("components", components),
            recording_support::field("reactions", basis.size()),
            recording_support::field("basis", "independent"),
        }), "row,column");
    recording_support::expect_group(
        recording_support::active_test_id(), "c5-mixed-basis", "mixed-basis-entry/", 3,
        recording_support::fields({
            recording_support::field("T", temperature),
            recording_support::field("components", components),
            recording_support::field("reactions", mixed.size()),
            recording_support::field("basis", "mixed"),
        }), "row,column");
    recording_support::expect_group(
        recording_support::active_test_id(), "c5-affinity", "affinity/", 2,
        recording_support::fields({
            recording_support::field("T", temperature),
            recording_support::field("components", components),
            recording_support::field("reactions", basis.size()),
        }), "reaction");
    recording_support::expect_group(
        recording_support::active_test_id(), "c5-basis-invariance", "basis-subspace", 1,
        recording_support::fields({
            recording_support::field("T", temperature),
            recording_support::field("components", components),
            recording_support::field("reactions", basis.size()),
        }), "basis");
    for (std::size_t row = 0; row < basis_coordinates.size(); ++row) {
        for (std::size_t column = 0; column < basis_coordinates[row].size(); ++column) {
            const auto point = std::to_string(row) + "/" + std::to_string(column);
            recording_support::observe(
                (std::string("basis-entry/") + point).c_str(), "balance_basis", "1",
                basis_coordinates[row][column], std::nullopt,
                recording_support::fields({
                    recording_support::field("T", temperature),
                    recording_support::field("basis_row", basis_coordinates[row]),
                    recording_support::field("row", row),
                    recording_support::field("column", column),
                }), "Relationship", "balance-basis", "basis entry is finite",
                "c5-basis", point.c_str(), "B_{ij}", "T,row,column", "B_{ij}", "row,column");
            recording_support::observe(
                (std::string("mixed-basis-entry/") + point).c_str(), "balance_basis", "1",
                mixed_coordinates[row][column], std::nullopt,
                recording_support::fields({
                    recording_support::field("T", temperature),
                    recording_support::field("basis_row", mixed_coordinates[row]),
                    recording_support::field("row", row),
                    recording_support::field("column", column),
                }), "Relationship", "balance-basis", "mixed basis entry is finite",
                "c5-mixed-basis", point.c_str(), "B_{ij}", "T,row,column", "B_{ij}", "row,column");
        }
    }
    for (std::size_t index = 0; index < basis_affinities.size(); ++index) {
        const auto point = std::to_string(index);
        recording_support::compare(
            (std::string("affinity/") + point).c_str(), "reaction_affinity", "J/mol",
            basis_affinities[index], mixed_affinities[index],
            recording_support::fields({
                recording_support::field("T", temperature),
                recording_support::field("basis_affinities", basis_affinities),
                recording_support::field("mixed_affinities", mixed_affinities),
                recording_support::field("reaction", index),
            }), "Relationship", "reaction-basis", "reaction affinity is invariant under basis change",
            "absolute", 0.0, 0.0, basis_affinities[index] == mixed_affinities[index],
            "c5-affinity", point.c_str(), "A", "T,reaction", "A_{basis}=A_{mixed}", "reaction");
    }
    const bool equal_subspace = epcsaft::equilibrium::same_subspace(
        basis_coordinates, mixed_coordinates, 1.0e-12);
    const bool equal_affinities = basis_affinities == mixed_affinities;
    recording_support::compare(
        "basis-subspace", "same_subspace", "bool", equal_subspace ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("T", temperature),
            recording_support::field("basis", "independent"),
            recording_support::field("mixed_basis", "mixed"),
            recording_support::field("tolerance", 1.0e-12),
        }), "Relationship", "basis-invariance", "reaction bases span the same balance subspace",
        "absolute", 0.0, 0.0, equal_subspace,
        "c5-basis", "subspace", "same_subspace", "T,basis,mixed_basis", "B_{basis}=B_{mixed}", "basis");
    CHECK(compiled.rank == 2);
    CHECK(compiled.generic_nullity == 0);
    CHECK(equal_subspace);
    CHECK(equal_affinities);
    auto dependent = problem;
    dependent.reactions.push_back({"dependent", {2.0, -2.0, 0.0}, 1.0, std::nullopt});
    bool rejected = false;
    try { (void)epcsaft::equilibrium::compile_problem(model, dependent); }
    catch (const std::invalid_argument& error) { rejected = std::string(error.what()).find("dependent") != std::string::npos; }
    CHECK(rejected);
    auto three_phase = equilibrium_support::flash(2);
    three_phase.phases = {equilibrium_support::phase("L1", epcsaft::equilibrium::PhaseKind::Liquid), equilibrium_support::phase("L2", epcsaft::equilibrium::PhaseKind::Liquid), equilibrium_support::phase("L3", epcsaft::equilibrium::PhaseKind::Liquid)};
    rejected = false;
    epcsaft::HelmholtzModel binary(properties_support::neutral_binary());
    try { (void)epcsaft::equilibrium::compile_problem(binary, three_phase); }
    catch (const std::invalid_argument& error) { rejected = std::string(error.what()).find("phase count") != std::string::npos; }
    CHECK(rejected);
}
void check_c8() {
    epcsaft::HelmholtzModel model(properties_support::aqueous_electrolyte());
    auto problem = equilibrium_support::flash(3);
    problem.feed.values = {0.98, 0.01, 0.01};
    problem.phases.front().support = {0, 1, 2};
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    auto second = problem;
    second.phases.push_back(equilibrium_support::phase(
        "V", epcsaft::equilibrium::PhaseKind::Vapor,
        epcsaft::equilibrium::AmountSpec::pinned(0.0)));
    second.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    const auto boundary = epcsaft::equilibrium::compile_problem(model, second);
    recording_support::compare(
        "charge-rows", "charge_row_count", "count", static_cast<double>(compiled.charge_row_count), 0.0,
        recording_support::fields({recording_support::field("T", problem.temperature.value),
                                   recording_support::field("P", problem.pressure.value),
                                   recording_support::field("x", problem.feed.values),
                                   recording_support::field("boundary_charge_rows", boundary.charge_row_count)}),
        "Relationship", "charge-balance", "single electrolyte phase has no charge row",
        "reference_or_one", 0.0, 0.0, compiled.charge_row_count == 0,
        "c8-charge", "single", "N_q", "T,P,\\mathbf{x},topology", "N_q=0", "topology");
    recording_support::compare(
        "boundary-charge-row", "charge_row_count", "count", static_cast<double>(boundary.charge_row_count), 1.0,
        recording_support::fields({recording_support::field("T", second.temperature.value),
                                   recording_support::field("P", second.pressure.value),
                                   recording_support::field("x", second.feed.values),
                                   recording_support::field("phase_count", second.phases.size()),
                                   recording_support::field("vapor_amount", second.phases[1].amount.fraction)}),
        "Relationship", "charge-balance", "boundary electrolyte topology adds one charge row",
        "reference_or_one", 0.0, 0.0, boundary.charge_row_count == 1,
        "c8-boundary-charge", "boundary", "N_q", "T,P,\\mathbf{x},topology", "N_q=1", "topology");
    CHECK(compiled.charge_row_count == 0);
    CHECK(boundary.charge_row_count == 1);
}
void check_c2() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto under = equilibrium_support::boundary();
    under.temperature = epcsaft::equilibrium::FreeValue::fixed(300.0);
    under.pressure = epcsaft::equilibrium::FreeValue::fixed(1.0e5);
    under.phases[1].name = "vapour-offender";
    bool named = false;
    try { (void)epcsaft::equilibrium::compile_problem(model, under); }
    catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        named = message.find("vapour-offender") != std::string::npos
            && message.find("pinned") != std::string::npos;
    }
    recording_support::compare(
        "named-rejection", "rejection", "bool", named ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", under.temperature.value),
                                   recording_support::field("P", under.pressure.value),
                                   recording_support::field("phase", under.phases[1].name),
                                   recording_support::field("amount", under.phases[1].amount.fraction)}),
        "Relationship", "admission", "pinned phase error names the phase",
        "reference_or_one", 0.0, 0.0, named, "c2-named-rejection", "phase",
        "accepted", "T,P,phase,amount", "error names phase and pinned amount", "phase");
    CHECK(named);
    auto omitted = equilibrium_support::flash(2);
    omitted.topology_declared = false;
    bool phase_count = false;
    try { (void)epcsaft::equilibrium::compile_problem(model, omitted); }
    catch (const std::invalid_argument& error) {
        phase_count = std::string(error.what()).find("phase count") != std::string::npos;
    }
    recording_support::compare(
        "topology-rejection", "rejection", "bool", phase_count ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", omitted.temperature.value),
                                   recording_support::field("P", omitted.pressure.value),
                                   recording_support::field("topology_declared", omitted.topology_declared ? 1U : 0U),
                                   recording_support::field("phase_count", omitted.phases.size())}),
        "Relationship", "admission", "undeclared topology is rejected",
        "reference_or_one", 0.0, 0.0, phase_count, "c2-topology-rejection", "topology",
        "accepted", "T,P,topology,phase_count", "undeclared topology raises phase-count error", "topology");
    CHECK(phase_count);
}
void check_c3() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    const auto problem = equilibrium_support::flash(2);
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    const auto point = compiled.generic_point();
    const auto values = compiled.evaluate(point).values;
    const auto jacobian = compiled.jacobian(point);
    recording_support::compare(
        "jacobian-shape", "jacobian_entries", "count", static_cast<double>(jacobian.size()),
        static_cast<double>(values.size() * point.size()),
        recording_support::fields({recording_support::field("T", problem.temperature.value),
                                   recording_support::field("P", problem.pressure.value),
                                   recording_support::field("x", problem.feed.values),
                                   recording_support::field("rows", values.size()),
                                   recording_support::field("variables", point.size())}),
        "Relationship", "jacobian-shape", "Jacobian has one entry per row and coordinate",
        "reference_or_one", 0.0, 0.0,
        jacobian.size() == values.size() * point.size(), "c3-jacobian", "shape",
        "J", "T,P,\\mathbf{x},shape", "|J|=N_{rows}N_{vars}", "shape");
    CHECK(values.size() == compiled.rows.size() && jacobian.size() == values.size() * point.size());
    CHECK(compiled.generic_nullity == 0);
    auto singular = equilibrium_support::flash(2);
    singular.phases.push_back(equilibrium_support::phase("third", epcsaft::equilibrium::PhaseKind::Liquid));
    singular.phases.push_back(equilibrium_support::phase("fourth", epcsaft::equilibrium::PhaseKind::Liquid));
    bool named = false;
    try { (void)epcsaft::equilibrium::compile_problem(model, singular); }
    catch (const std::invalid_argument& error) {
        named = std::string(error.what()).find("phase count") != std::string::npos;
    }
    CHECK(named);
}
void check_c6() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto flash = equilibrium_support::flash(2);
    auto boundary = equilibrium_support::boundary();
    flash.phases = boundary.phases;
    flash.phases[1].amount = epcsaft::equilibrium::AmountSpec::free();
    flash.pressure = epcsaft::equilibrium::FreeValue::fixed(1.0e5);
    const auto a = epcsaft::equilibrium::compile_problem(model, flash);
    const auto b = epcsaft::equilibrium::compile_problem(model, boundary);
    auto flash_point = a.generic_point();
    flash_point[6] = 0.0;
    auto boundary_point = b.generic_point();
    for (std::size_t i = 0; i < 6; ++i) boundary_point[i] = flash_point[i];
    boundary_point[6] = flash_point[7];
    boundary_point[7] = 1.0e5;
    for (std::size_t i = 8; i < boundary_point.size(); ++i) boundary_point[i] = flash_point[i];
    const auto flash_values = a.evaluate(flash_point).values;
    const auto boundary_values = b.evaluate(boundary_point).values;
    const auto comparison_inputs = recording_support::fields({
        recording_support::field("T", flash.temperature.value),
        recording_support::field("P", flash.pressure.value),
        recording_support::field("x", flash.feed.values),
        recording_support::field("phase_count", flash.phases.size()),
        recording_support::field("boundary_phase_count", boundary.phases.size()),
        recording_support::field("flash_rows", flash_values.size()),
        recording_support::field("boundary_rows", boundary_values.size()),
    });
    recording_support::expect_group(
        recording_support::active_test_id(), "c6-boundary-limit", "flash-boundary/",
        flash_values.size(), comparison_inputs, "row");
    for (std::size_t row = 0; row < flash_values.size(); ++row) {
        const bool has_reference = row < boundary_values.size();
        const double reference = has_reference
            ? boundary_values[row] : std::numeric_limits<double>::quiet_NaN();
        const bool equal = has_reference
            && std::abs(flash_values[row] - reference) < 1.0e-12;
        const auto point = std::to_string(row);
        recording_support::compare(
            (std::string("flash-boundary/") + point).c_str(), "residual", "1",
            flash_values[row], reference,
            recording_support::fields({comparison_inputs,
                                       recording_support::field("row", row)}),
            "Relationship", "topology-limit", "flash and boundary residuals agree at zero vapor amount",
            "strict_absolute", 1.0e-12, 0.0, equal,
            "c6-boundary-limit", point.c_str(), "r", "T,P,\\mathbf{x},topology,row",
            "\\Delta r=0", "row");
        CHECK(equal);
    }
    CHECK(flash_values.size() == boundary_values.size());
}
void check_c7() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto first = equilibrium_support::boundary();
    auto second = first;
    std::swap(second.phases[0], second.phases[1]);
    const auto left = epcsaft::equilibrium::compile_problem(model, first);
    const auto right = epcsaft::equilibrium::compile_problem(model, second);
    auto point = left.generic_point();
    auto permuted = right.generic_point();
    for (std::size_t i = 0; i < 3; ++i) permuted[i] = point[4 + i];
    for (std::size_t i = 0; i < 4; ++i) permuted[3 + i] = point[i];
    permuted[7] = point[7];
    auto a = left.evaluate(point).values;
    auto b = right.evaluate(permuted).values;
    const auto phase_inputs = recording_support::fields({
        recording_support::field("T", first.temperature.value),
        recording_support::field("P", first.pressure.value),
        recording_support::field("x", first.feed.values),
        recording_support::field("left_phase_count", first.phases.size()),
        recording_support::field("right_phase_count", second.phases.size()),
    });
    const double phase_residual_difference = std::abs(
        left.evaluate(point).residual_norm() - right.evaluate(permuted).residual_norm());
    recording_support::expect_group(
        recording_support::active_test_id(), "c7-phase-residual", "phase-residual", 1,
        phase_inputs, "phase_set");
    recording_support::expect_group(
        recording_support::active_test_id(), "c7-phase-coordinate", "phase-coordinate/",
        a.size(), phase_inputs, "row");
    recording_support::compare(
        "phase-residual", "residual_norm_difference", "1", phase_residual_difference, 0.0,
        phase_inputs, "Relationship", "permutation", "phase permutation preserves residual norm",
        "absolute", 1.0e-14, 0.0, phase_residual_difference <= 1.0e-14,
        "c7-phase-residual", "phase_set", "\\Delta ||r||", "T,P,\\mathbf{x},phase_set",
        "\\Delta||r||=0", "phase_set");
    CHECK(left.generic_nullity == right.generic_nullity && a.size() == b.size());
    std::sort(a.begin(), a.end()); std::sort(b.begin(), b.end());
    for (std::size_t i = 0; i < a.size(); ++i) {
        const bool equal = i < b.size() && std::abs(std::abs(a[i]) - std::abs(b[i])) <= 1.0e-14;
        const double reference = i < b.size()
            ? std::abs(b[i]) : std::numeric_limits<double>::quiet_NaN();
        const auto row = std::to_string(i);
        recording_support::compare(
            (std::string("phase-coordinate/") + row).c_str(), "residual_coordinate", "1",
            std::abs(a[i]), reference,
            recording_support::fields({phase_inputs, recording_support::field("row", i)}),
            "Relationship", "permutation", "phase permutation preserves residual coordinates",
            "absolute", 1.0e-14, 0.0, equal,
            "c7-phase-coordinate", row.c_str(), "r_i", "T,P,\\mathbf{x},phase_set,row",
            "|r_i|_{left}=|r_i|_{right}", "row");
        CHECK(equal);
    }
    epcsaft::HelmholtzModel relabelled(equilibrium_support::relabel_neutral_binary());
    auto relabelled_problem = first;
    relabelled_problem.feed.values = {first.feed.values[1], first.feed.values[0]};
    const auto relabelled_system = epcsaft::equilibrium::compile_problem(relabelled, relabelled_problem);
    auto relabelled_point = point;
    std::swap(relabelled_point[0], relabelled_point[1]); std::swap(relabelled_point[4], relabelled_point[5]);
    const auto relabelled_values = relabelled_system.evaluate(relabelled_point).values;
    const double reference_norm = left.evaluate(point).residual_norm();
    const double relabelled_error = std::abs(
        relabelled_system.evaluate(relabelled_point).residual_norm()
        - reference_norm);
    const bool relabelled_ok = std::isfinite(relabelled_error)
        && std::isfinite(reference_norm)
        && relabelled_error <= 1.0e-14 * reference_norm;
    recording_support::compare_custom(
        "phase-permutation", "residual_norm_difference", "1", relabelled_error, reference_norm,
        recording_support::fields({recording_support::field("T", first.temperature.value),
                                   recording_support::field("P", first.pressure.value),
                                   recording_support::field("x", first.feed.values),
                                   recording_support::field("permutation", "L,V->V,L"),
                                   recording_support::field("reference_norm", reference_norm)}),
        "Relationship", "permutation", "phase permutation preserves residual norm",
        "isfinite(model) && isfinite(reference) && model <= reference * 1e-14", relabelled_ok,
        "c7-phase-permutation", "phase", "\\Delta ||r||", "T,P,\\mathbf{x},phase",
        "\\Delta||r||\\leq10^{-14}||r||", "phase");
    CHECK(relabelled_values.size() == a.size());
    CHECK(relabelled_ok);
}
void check_c9() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    std::size_t case_index = 0;
    for (const auto phase_set : {equilibrium_support::PhaseSet::OneLiquid,
                                 equilibrium_support::PhaseSet::OneVapor}) {
        auto problem = equilibrium_support::lattice_problem(
            phase_set, equilibrium_support::AmountPattern::AllFree,
            equilibrium_support::IntensivePattern::BothFixed,
            equilibrium_support::Chemistry::Neutral
        );
        const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
        recording_support::compare(
            (std::string("square/") + std::to_string(case_index)).c_str(),
            "rows_minus_unknowns", "count",
            static_cast<double>(compiled.rows_minus_unknowns), 0.0,
            recording_support::fields({recording_support::field("phase_set", phase_set_name(phase_set)),
                                       recording_support::field("T", problem.temperature.value),
                                       recording_support::field("P", problem.pressure.value),
                                       recording_support::field("x", problem.feed.values),
                                       recording_support::field("phase_count", problem.phases.size())}),
            "Relationship", "rank-count", "single-phase specification is square",
            "reference_or_one", 0.0, 0.0, compiled.rows_minus_unknowns == 0,
            "c9-square", std::to_string(case_index).c_str(), "N_{rows}-N_{unknowns}", "phase,T,P,\\mathbf{x}",
            "N_{rows}=N_{unknowns}", "phase");
        ++case_index;
        CHECK(compiled.rows_minus_unknowns == 0);
        CHECK(compiled.generic_nullity == 0);
    }
}

// ---------------------------------------------------------------- E2 admission

namespace e2 {

using epcsaft::equilibrium::AmountSpec;
using epcsaft::equilibrium::FreeValue;
using epcsaft::equilibrium::PhaseKind;

std::string rejection(const epcsaft::HelmholtzModel& model,
                      const epcsaft::equilibrium::EquilibriumProblem& problem) {
    try { (void)epcsaft::equilibrium::compile_problem(model, problem); }
    catch (const std::invalid_argument& error) { return error.what(); }
    return "";
}

bool record_bool(const char* id, const char* group, const char* claim, bool verdict,
                 const std::string& inputs) {
    recording_support::compare(id, "verdict", "bool", verdict ? 1.0 : 0.0, 1.0, inputs,
        "Relationship", "admission", claim, "absolute", 0.0, 0.0, verdict, group, id,
        "verdict", "inputs", claim, "case");
    return verdict;
}

bool record_rejection(const char* id, const char* group, const epcsaft::HelmholtzModel& model,
                      const epcsaft::equilibrium::EquilibriumProblem& problem,
                      const char* phrase, const std::string& inputs) {
    const std::string message = rejection(model, problem);
    const bool verdict = message.find(phrase) != std::string::npos;
    return record_bool(id, group, phrase,
                       verdict, recording_support::fields({inputs, recording_support::field("message", message)}));
}

std::string pin_inputs(const epcsaft::equilibrium::EquilibriumProblem& problem) {
    std::vector<double> pins;
    for (const auto& phase : problem.phases) pins.push_back(phase.is_pinned() ? phase.amount.fraction : -1.0);
    double total = 0.0, sum = 0.0;
    for (const double value : problem.feed.values) total += value;
    for (const double value : pins) if (value >= 0.0) sum += value;
    return recording_support::fields({
        recording_support::field("feed", problem.feed.values), recording_support::field("N_feed", total),
        recording_support::field("pins", pins), recording_support::field("S_pin", sum),
        recording_support::field("A_pin", sum * total),
        recording_support::field("T_free", problem.temperature.free ? 1U : 0U),
        recording_support::field("P_free", problem.pressure.free ? 1U : 0U)});
}

epcsaft::equilibrium::EquilibriumProblem two_phase(std::vector<double> feed, double first, double second) {
    auto problem = equilibrium_support::flash(feed.size());
    problem.feed.values = std::move(feed);
    problem.phases = {equilibrium_support::phase("L", PhaseKind::Liquid, AmountSpec::pinned(first)),
                      equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(second))};
    problem.temperature = FreeValue::variable(300.0);
    problem.pressure = FreeValue::variable(1.0e5);
    return problem;
}

}  // namespace e2

void check_e2_pins() {
    using namespace e2;
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto over = two_phase({1.0, 1.0}, 0.6, 0.5);
    {
        epcsaft::core::Capture capture;
        CHECK(record_rejection("e2-0/over-sum", "e2-0-pins", model, over,
            "pinned phase fractions exceed one", pin_inputs(over)));
        bool generic = false;
        for (const auto& span : capture.trace().spans)
            generic = generic || std::string(span.name) == "equilibrium_generic_point"
                || std::string(span.name) == "boundary_start";
        CHECK(record_bool("e2-0/over-sum-before-rank", "e2-0-pins",
            "rejection precedes generic-point and start work", !generic, pin_inputs(over)));
    }
    auto under = two_phase({1.0, 1.0}, 0.4, 0.5);
    CHECK(record_rejection("e2-0/all-pinned-under-sum", "e2-0-pins", model, under,
        "all-pinned phase fractions must sum to one", pin_inputs(under)));
    auto bad = two_phase({1.0, 1.0}, 1.2, 0.0);
    bad.phases[1].amount = AmountSpec::free();
    bad.pressure = FreeValue::fixed(1.0e5);
    CHECK(record_rejection("e2-0/above-one", "e2-0-pins", model, bad,
        "invalid pinned phase fraction in phase L", pin_inputs(bad)));
    // Sum-consistent but rank-deficient: pin admission passes and the pre-solve record is returned.
    auto exact = two_phase({1.0, 1.0}, 0.5, 0.5);
    const auto compiled = epcsaft::equilibrium::compile_problem(model, exact);
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, exact);
    CHECK(record_bool("e2-0/exact-sum-not-run", "e2-0-pins",
        "exact all-pinned sum passes admission and returns the pre-solve infeasible record",
        compiled.is_square() && !compiled.feasible && result.solver_status == "NOT_RUN"
            && result.convergence.stop == epcsaft::core::Stop::infeasible
            && result.message.find("invalid pinned") == std::string::npos,
        recording_support::fields({pin_inputs(exact),
            recording_support::field("generic_nullity", compiled.generic_nullity)})));
    // Reactive counterexamples: A <-> B conserves moles (rules apply, both orientations).
    for (const double sign : {1.0, -1.0}) {
        auto reactive_over = two_phase({1.0, 0.0}, 0.6, 0.5);
        reactive_over.reactions = {{"A=B", {-sign, sign}, sign > 0.0 ? 2.0 : 0.5, std::nullopt}};
        reactive_over.phases[0].support = {0, 1};
        reactive_over.phases[1].support = {0, 1};
        const auto id = std::string("e2-0/isomerisation-over-sum/") + (sign > 0.0 ? "forward" : "reversed");
        CHECK(record_rejection(id.c_str(), "e2-0-pins", model, reactive_over,
            "pinned phase fractions exceed one", pin_inputs(reactive_over)));
        auto reactive_under = two_phase({1.0, 0.0}, 0.4, 0.5);
        reactive_under.reactions = reactive_over.reactions;
        reactive_under.phases[0].support = {0, 1};
        reactive_under.phases[1].support = {0, 1};
        const auto under_id = std::string("e2-0/isomerisation-under-sum/") + (sign > 0.0 ? "forward" : "reversed");
        CHECK(record_rejection(under_id.c_str(), "e2-0-pins", model, reactive_under,
            "all-pinned phase fractions must sum to one", pin_inputs(reactive_under)));
    }
    // A <-> 2B: amounts (0.75, 0.5) satisfy B n = B feed and sum to 1.25; no sum rejection.
    auto doubling = two_phase({1.0, 0.0}, 0.625, 0.625);
    doubling.reactions = {{"A=2B", {-1.0, 2.0}, 2.0, std::nullopt}};
    for (auto& phase : doubling.phases) { phase.support = {0, 1}; phase.composition_guess = {0.6, 0.4}; }
    const auto doubling_message = rejection(model, doubling);
    const auto doubling_basis = epcsaft::equilibrium::balance_basis(2, doubling.reactions);
    const auto balance_gap = [](const std::vector<std::vector<double>>& basis, const std::vector<double>& amounts,
                                const std::vector<double>& feed) {
        const auto left = epcsaft::equilibrium::matvec(basis, amounts);
        const auto right = epcsaft::equilibrium::matvec(basis, feed);
        double result = 0.0;
        for (std::size_t i = 0; i < left.size(); ++i) result = std::max(result, std::abs(left[i] - right[i]));
        return result;
    };
    const double doubling_balance = balance_gap(doubling_basis, {0.75, 0.5}, doubling.feed.values);
    CHECK(record_bool("e2-0/doubling-no-sum-rejection", "e2-0-pins",
        "nonconserving A=2B pins are not rejected by a total-mole sum rule",
        doubling_message.find("pinned phase fractions") == std::string::npos
            && doubling_message.find("all-pinned") == std::string::npos && doubling_balance <= 1.0e-12,
        recording_support::fields({pin_inputs(doubling), recording_support::field("balance_residual", doubling_balance),
            recording_support::field("message", doubling_message)})));
    if (doubling_message.empty()) {
        const auto doubling_compiled = epcsaft::equilibrium::compile_problem(model, doubling);
        recording_support::observe("e2-0/doubling-generic-nullity", "generic_nullity", "count",
            static_cast<double>(doubling_compiled.generic_nullity), std::nullopt, pin_inputs(doubling),
            "Relationship", "rank-count", "retained rank admission of the A=2B all-pinned declaration",
            "e2-0-pins", "doubling");
    }
    // 2A <-> B: one phase, pin 0.75, one free intensive passes the sum stage.
    auto halving = equilibrium_support::flash(2);
    halving.feed.values = {1.0, 0.0};
    halving.phases = {equilibrium_support::phase("L", PhaseKind::Liquid, AmountSpec::pinned(0.75))};
    halving.phases[0].support = {0, 1};
    halving.phases[0].composition_guess = {2.0 / 3.0, 1.0 / 3.0};
    halving.temperature = FreeValue::variable(300.0);
    halving.reactions = {{"2A=B", {-2.0, 1.0}, 2.0, std::nullopt}};
    const auto halving_message = rejection(model, halving);
    const auto halving_basis = epcsaft::equilibrium::balance_basis(2, halving.reactions);
    const double halving_balance = balance_gap(halving_basis, {0.5, 0.25}, halving.feed.values);
    CHECK(record_bool("e2-0/halving-no-sum-rejection", "e2-0-pins",
        "nonconserving 2A=B pin passes the sum stage",
        halving_message.find("pinned phase fractions") == std::string::npos
            && halving_message.find("all-pinned") == std::string::npos && halving_balance <= 1.0e-12,
        recording_support::fields({pin_inputs(halving), recording_support::field("balance_residual", halving_balance),
            recording_support::field("message", halving_message)})));
}

void check_e2_support() {
    using namespace e2;
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto base = equilibrium_support::flash(2);
    base.feed.values = {0.5, 0.5};
    const auto inputs = recording_support::fields({recording_support::field("feed", base.feed.values)});
    const auto omitted = epcsaft::equilibrium::compile_problem(model, base);
    auto full = base;
    full.phases[0].support = {0, 1};
    const auto explicit_full = epcsaft::equilibrium::compile_problem(model, full);
    CHECK(record_bool("e2-3/omitted-and-full", "e2-3-support", "omitted and explicit full support agree",
        omitted.problem.phases[0].support == std::vector<std::size_t>{0, 1}
            && explicit_full.problem.phases[0].support == std::vector<std::size_t>{0, 1}
            && omitted.removed_species.empty() && explicit_full.removed_species.empty(), inputs));
    auto empty = base;
    empty.phases[0].support_declared = true;
    CHECK(record_rejection("e2-3/explicit-empty", "e2-3-support", model, empty, "empty support in phase L", inputs));
    auto partial = base;
    partial.phases[0].support = {0};
    CHECK(record_rejection("e2-3/missing-species", "e2-3-support", model, partial, "species 1 has no support", inputs));
    auto duplicate = base;
    duplicate.phases[0].support = {0, 0};
    CHECK(record_rejection("e2-3/duplicate", "e2-3-support", model, duplicate, "duplicate support index", inputs));
    auto range = base;
    range.phases[0].support = {0, 5};
    CHECK(record_rejection("e2-3/out-of-range", "e2-3-support", model, range, "support index out of range", inputs));
    auto pruned = full;
    pruned.feed.values = {0.5, 0.0};
    const auto pruned_compiled = epcsaft::equilibrium::compile_problem(model, pruned);
    auto pruned_omitted = base;
    pruned_omitted.feed.values = {0.5, 0.0};
    const auto pruned_omitted_compiled = epcsaft::equilibrium::compile_problem(model, pruned_omitted);
    CHECK(record_bool("e2-3/zero-feed-pruning", "e2-3-support",
        "zero-feed nonparticipant is removed exactly once from explicit and omitted support",
        pruned_compiled.problem.phases[0].support == std::vector<std::size_t>{0}
            && pruned_compiled.removed_species == std::vector<std::size_t>{1}
            && pruned_omitted_compiled.problem.phases[0].support == std::vector<std::size_t>{0}
            && pruned_omitted_compiled.removed_species == std::vector<std::size_t>{1},
        recording_support::fields({recording_support::field("feed", pruned.feed.values)})));
}

void check_e2_reaction() {
    using namespace e2;
    epcsaft::HelmholtzModel model(properties_support::neutral_ternary());
    auto forward = equilibrium_support::flash(3);
    forward.feed.values = {1.0, 0.0, 0.0};
    forward.reactions = {{"A=B", {-1.0, 1.0, 0.0}, 2.0, std::nullopt}};
    auto reversed = forward;
    reversed.reactions = {{"B=A", {1.0, -1.0, 0.0}, 0.5, std::nullopt}};
    const auto left = epcsaft::equilibrium::compile_problem(model, forward);
    const auto right = epcsaft::equilibrium::compile_problem(model, reversed);
    bool same_rows = left.rows.size() == right.rows.size();
    for (std::size_t row = 0; same_rows && row < left.rows.size(); ++row)
        same_rows = left.rows[row].kind == right.rows[row].kind && left.rows[row].phase == right.rows[row].phase
            && left.rows[row].index == right.rows[row].index;
    const auto inputs = recording_support::fields({recording_support::field("feed", forward.feed.values),
        recording_support::field("nu", forward.reactions[0].stoichiometry), recording_support::field("K", 2.0)});
    CHECK(record_bool("e2-2/orientation", "e2-2-reaction",
        "support, rank, balance subspace and rows are independent of reaction orientation",
        left.problem.phases[0].support == std::vector<std::size_t>{0, 1}
            && right.problem.phases[0].support == std::vector<std::size_t>{0, 1}
            && left.removed_species == std::vector<std::size_t>{2} && right.removed_species == std::vector<std::size_t>{2}
            && left.rank == 2 && right.rank == 2 && epcsaft::equilibrium::same_subspace(left.B, right.B, 1.0e-12)
            && same_rows, inputs));
    auto dependent = forward;
    dependent.reactions.push_back({"twice", {-2.0, 2.0, 0.0}, 4.0, std::nullopt});
    CHECK(record_rejection("e2-2/dependent", "e2-2-reaction", model, dependent, "dependent reaction twice", inputs));
    // Standard-state falsifier at T = 300 K, nu = (-1, +1, 0), K = 2.
    const double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const double g_b = -R * 300.0 * std::log(2.0);
    auto supplied = equilibrium_support::flash(3);
    supplied.reactions = forward.reactions;
    supplied.standard_gibbs_j_per_mol = {0.0, g_b, 0.0};
    auto gauge = supplied;
    gauge.standard_gibbs_j_per_mol = {100.0, g_b + 100.0, 0.0};
    const auto supplied_result = epcsaft::equilibrium::solve_equilibrium(model, supplied);
    const auto gauge_result = epcsaft::equilibrium::solve_equilibrium(model, gauge);
    const auto supplied_compiled = epcsaft::equilibrium::compile_problem(model, supplied);
    const auto gauge_compiled = epcsaft::equilibrium::compile_problem(model, gauge);
    double gauge_difference = 0.0;
    if (supplied_result.success && gauge_result.success) {
        for (std::size_t i = 0; i < supplied_result.mole_fractions.size(); ++i)
            gauge_difference = std::max(gauge_difference,
                std::abs(supplied_result.mole_fractions[i] - gauge_result.mole_fractions[i]));
        gauge_difference = std::max(gauge_difference,
            std::abs(supplied_result.molar_densities[0] / gauge_result.molar_densities[0] - 1.0));
    }
    const double gauge_limit = 1.0e-9;
    recording_support::compare("e2-2/gauge-shift", "state_difference", "1", gauge_difference, 0.0,
        recording_support::fields({recording_support::field("T", 300.0), recording_support::field("g_supplied", supplied.standard_gibbs_j_per_mol),
            recording_support::field("g_gauge", gauge.standard_gibbs_j_per_mol),
            recording_support::field("x_supplied", supplied_result.mole_fractions),
            recording_support::field("x_gauge", gauge_result.mole_fractions)}),
        "Relationship", "gauge-invariance", "a conserved constant Gibbs shift leaves the solved state unchanged",
        "strict_absolute", gauge_limit, 0.0, supplied_result.success && gauge_result.success && gauge_difference <= gauge_limit,
        "e2-2-reaction", "gauge", "\\Delta x", "T,g", "\\Delta x=0", "case");
    CHECK(supplied_result.success && gauge_result.success && gauge_difference <= gauge_limit);
    CHECK(supplied_compiled.rank == gauge_compiled.rank && supplied_compiled.rows.size() == gauge_compiled.rows.size());
    auto zero = supplied;
    zero.standard_gibbs_j_per_mol = {0.0, 0.0, 0.0};
    CHECK(record_rejection("e2-2/zero-vector-with-K", "e2-2-reaction", model, zero, "inconsistent standard state", inputs));
    auto both = supplied;
    both.reactions[0].affinity_j_per_mol = 12.0;
    CHECK(record_rejection("e2-2/two-targets", "e2-2-reaction", model, both, "exactly one target", inputs));
    auto free_t = supplied;
    free_t.temperature = FreeValue::variable(300.0);
    free_t.pressure = FreeValue::fixed(1.0e5);
    free_t.phases[0].amount = AmountSpec::pinned(1.0);
    CHECK(record_rejection("e2-2/free-T-point-consistent", "e2-2-reaction", model, free_t,
        "temperature-dependent standard state required", inputs));
    auto affinity_free_t = free_t;
    affinity_free_t.reactions = {{"A=B", {-1.0, 1.0, 0.0}, std::nullopt, g_b}};
    CHECK(record_bool("e2-2/free-T-affinity-admitted", "e2-2-reaction",
        "a constant affinity target admits a free temperature with a consistent supplied vector",
        rejection(model, affinity_free_t).empty(), inputs));
    // Charged compiler/permutation algebra: (water, Na+, Cl-), feed (0,1,1), nu = (2,-1,-1), K = 4.
    epcsaft::HelmholtzModel salt(properties_support::aqueous_electrolyte());
    auto charged = equilibrium_support::flash(3);
    charged.feed.values = {0.0, 1.0, 1.0};
    charged.temperature.value = 298.15;
    charged.phases[0].support = {0, 1, 2};
    charged.reactions = {{"Na+Cl=2W", {2.0, -1.0, -1.0}, 4.0, std::nullopt}};
    const auto charged_compiled = epcsaft::equilibrium::compile_problem(salt, charged);
    const std::vector<std::size_t> order = {2, 0, 1};
    epcsaft::HelmholtzModel permuted_salt(equilibrium_support::permute_parameters(properties_support::aqueous_electrolyte(), order));
    auto permuted = charged;
    permuted.feed.values = equilibrium_support::permute_vector(charged.feed.values, order);
    permuted.reactions[0].stoichiometry = equilibrium_support::permute_vector(charged.reactions[0].stoichiometry, order);
    const auto permuted_compiled = epcsaft::equilibrium::compile_problem(permuted_salt, permuted);
    std::vector<std::vector<double>> mapped_basis;
    for (const auto& row : permuted_compiled.B) {
        std::vector<double> original(3, 0.0);
        for (std::size_t i = 0; i < 3; ++i) original[order[i]] = row[i];
        mapped_basis.push_back(original);
    }
    bool same_topology = charged_compiled.rows.size() == permuted_compiled.rows.size();
    for (std::size_t row = 0; same_topology && row < charged_compiled.rows.size(); ++row)
        same_topology = charged_compiled.rows[row].kind == permuted_compiled.rows[row].kind
            && charged_compiled.rows[row].phase == permuted_compiled.rows[row].phase;
    auto charge_feed = charged;
    charge_feed.feed.values = {0.98, 0.02, 0.0};
    auto charge_reaction = charged;
    charge_reaction.reactions = {{"bad", {0.0, 1.0, 0.0}, 1.0, std::nullopt}};
    CHECK(record_bool("e2-2/charged-permutation", "e2-2-reaction",
        "charged reaction rank, balance subspace, participant support and rows follow species identity",
        charged_compiled.rank == 2 && permuted_compiled.rank == 2
            && charged_compiled.problem.phases[0].support == std::vector<std::size_t>{0, 1, 2}
            && permuted_compiled.problem.phases[0].support == std::vector<std::size_t>{0, 1, 2}
            && same_topology
            && charged_compiled.charge_row_count == permuted_compiled.charge_row_count
            && epcsaft::equilibrium::same_subspace(charged_compiled.B, mapped_basis, 1.0e-12)
            && rejection(salt, charge_feed).find("feed carries net z") != std::string::npos
            && rejection(salt, charge_reaction).find("violates z") != std::string::npos,
        recording_support::fields({recording_support::field("feed", charged.feed.values),
            recording_support::field("nu", charged.reactions[0].stoichiometry), recording_support::field("K", 4.0),
            recording_support::field("permutation", std::vector<double>{2, 0, 1})})));
}

void check_e2_charge() {
    using namespace e2;
    epcsaft::HelmholtzModel model(properties_support::aqueous_electrolyte());
    auto base = equilibrium_support::flash(3);
    base.feed.values = {0.98, 0.01, 0.01};
    base.temperature.value = 298.15;
    base.phases[0].support = {0, 1, 2};
    const auto inputs = recording_support::fields({recording_support::field("feed", base.feed.values),
        recording_support::field("T", 298.15), recording_support::field("P", 1.0e5)});
    const auto record_rows = [&](const char* id, const epcsaft::equilibrium::CompiledSystem& compiled,
                                 std::size_t expected_count, const std::vector<std::size_t>& expected_phases) {
        std::vector<double> phases, values;
        const auto generic = compiled.evaluate(compiled.generic_point()).values;
        bool finite = true;
        for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
            if (compiled.rows[row].kind != epcsaft::equilibrium::RowKind::Electro) continue;
            phases.push_back(static_cast<double>(compiled.rows[row].phase));
            values.push_back(generic[row]);
            finite = finite && std::isfinite(generic[row]);
        }
        std::vector<double> expected(expected_phases.begin(), expected_phases.end());
        const bool verdict = compiled.charge_row_count == expected_count && phases == expected && finite;
        recording_support::compare(id, "charge_row_count", "count", static_cast<double>(compiled.charge_row_count),
            static_cast<double>(expected_count),
            recording_support::fields({inputs, recording_support::field("charge_row_phases", phases),
                recording_support::field("generic_row_values", values)}),
            "Relationship", "charge-balance", "retained charge rows follow the first-eligible carrier rule",
            "absolute", 0.0, 0.0, verdict, "e2-4-charge", id, "N_q", "topology", "N_q=expected", "case");
        CHECK(verdict);
    };
    record_rows("e2-4/one-free", epcsaft::equilibrium::compile_problem(model, base), 0, {});
    auto free_plus_zero = base;
    free_plus_zero.phases.push_back(equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.0)));
    free_plus_zero.phases[1].support = {0, 1, 2};
    free_plus_zero.pressure = FreeValue::variable(1.0e5);
    record_rows("e2-4/free-plus-incipient", epcsaft::equilibrium::compile_problem(model, free_plus_zero), 1, {1});
    auto all_pinned = free_plus_zero;
    all_pinned.phases[0].amount = AmountSpec::pinned(1.0);
    all_pinned.temperature = FreeValue::variable(298.15);
    record_rows("e2-4/all-pinned-one-zero", epcsaft::equilibrium::compile_problem(model, all_pinned), 1, {1});
    auto two_zero = base;
    two_zero.phases[0].amount = AmountSpec::pinned(0.0);
    two_zero.phases.push_back(equilibrium_support::phase("L2", PhaseKind::Liquid, AmountSpec::pinned(0.0)));
    two_zero.phases[1].support = {0, 1, 2};
    two_zero.phases.push_back(equilibrium_support::phase("V", PhaseKind::Vapor));
    two_zero.phases[2].support = {0};
    two_zero.temperature = FreeValue::variable(298.15);
    two_zero.pressure = FreeValue::variable(1.0e5);
    const auto two_zero_compiled = epcsaft::equilibrium::compile_problem(model, two_zero);
    record_rows("e2-4/two-incipient-ionic", two_zero_compiled, 2, {0, 1});
    const auto two_zero_result = epcsaft::equilibrium::solve_equilibrium(model, two_zero);
    recording_support::observe("e2-4/two-incipient-ionic-solve", "success", "bool",
        two_zero_result.success ? 1.0 : 0.0, std::nullopt,
        recording_support::fields({inputs, recording_support::field("solver_status", two_zero_result.solver_status),
            recording_support::field("stop", epcsaft::core::stop_name(two_zero_result.convergence.stop))}),
        "Relationship", "charge-balance", "no free phase can carry the ionic feed; inspectable rows, expected infeasible",
        "e2-4-charge", "two-incipient-ionic");
    CHECK(!two_zero_result.success);
    auto one_sign = base;
    one_sign.phases[0].support = {0, 1};
    CHECK(record_rejection("e2-4/one-sign", "e2-4-charge", model, one_sign, "support z sign in phase L", inputs));
    const auto solved = epcsaft::equilibrium::solve_equilibrium(model, base);
    double charge = 0.0;
    if (solved.success) for (std::size_t i = 0; i < 3; ++i) charge += model.parameters().charge_number[i] * solved.mole_fractions[i];
    recording_support::compare("e2-4/solved-neutrality", "phase_charge", "1", charge, 0.0, inputs,
        "Relationship", "charge-balance", "solved rowless ionic phase is electroneutral", "strict_absolute",
        epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value, 0.0,
        solved.success && std::abs(charge) <= epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value,
        "e2-4-charge", "solved", "q", "T,P,feed", "q=0", "case");
    CHECK(solved.success && std::abs(charge) <= epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value);
}

}  // namespace
int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.specification.c1", "Specification_C1", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "structural/*", "engine/tests/native/specification_checks.cpp", "specification_checks", "C1", "structural/*:480"}, check_c1},
        {{"native.specification.c2", "Specification_C2", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "named-rejection,topology-rejection", "engine/tests/native/specification_checks.cpp", "specification_checks", "C2"}, check_c2},
        {{"native.specification.c3", "Specification_C3", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "jacobian-shape", "engine/tests/native/specification_checks.cpp", "specification_checks", "C3"}, check_c3},
        {{"native.specification.c4", "Specification_C4", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "feed-composition/*", "engine/tests/native/specification_checks.cpp", "specification_checks", "C4", "feed-composition/*:6"}, check_c4},
        {{"native.specification.c5", "Specification_C5", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "rank,basis-entry/*,mixed-basis-entry/*,affinity/*,basis-subspace", "engine/tests/native/specification_checks.cpp", "specification_checks", "C5", "basis-entry/*:3,mixed-basis-entry/*:3,affinity/*:2,basis-subspace:1"}, check_c5},
        {{"native.specification.c6", "Specification_C6", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "flash-boundary/*", "engine/tests/native/specification_checks.cpp", "specification_checks", "C6"}, check_c6},
        {{"native.specification.c7", "Specification_C7", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "phase-residual,phase-coordinate/*,phase-permutation", "engine/tests/native/specification_checks.cpp", "specification_checks", "C7"}, check_c7},
        {{"native.specification.c8", "Specification_C8", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "charge-rows,boundary-charge-row", "engine/tests/native/specification_checks.cpp", "specification_checks", "C8"}, check_c8},
        {{"native.specification.c9", "Specification_C9", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "square/*", "engine/tests/native/specification_checks.cpp", "specification_checks", "C9"}, check_c9},
        {{"native.specification.e2-pins", "Specification_E2_pins", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E2", "e2-0/*", "engine/tests/native/specification_checks.cpp", "specification_checks", "E2-0"}, check_e2_pins},
        {{"native.specification.e2-reaction", "Specification_E2_reaction", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E2", "e2-2/*", "engine/tests/native/specification_checks.cpp", "specification_checks", "E2-2"}, check_e2_reaction},
        {{"native.specification.e2-support", "Specification_E2_support", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E2", "e2-3/*", "engine/tests/native/specification_checks.cpp", "specification_checks", "E2-3"}, check_e2_support},
        {{"native.specification.e2-charge", "Specification_E2_charge", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E2", "e2-4/*", "engine/tests/native/specification_checks.cpp", "specification_checks", "E2-4"}, check_e2_charge},
    });
}
