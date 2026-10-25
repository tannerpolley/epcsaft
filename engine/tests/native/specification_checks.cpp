#include "equilibrium/problem.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "core/physical_constants.hpp"
#include "equilibrium_support.hpp"
#include "recording_support.hpp"
#include "surface/parameter_builder.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
namespace {
// Two balance bases span the same row space when stacking them adds no rank (any basis, not only orthonormal).
bool same_subspace(const std::vector<std::vector<double>>& left, const std::vector<std::vector<double>>& right, double tolerance) {
    if (left.empty() || right.empty()) return left.empty() && right.empty();
    const auto rank = [&](std::size_t rows, std::size_t skip) {
        Eigen::MatrixXd matrix(static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(left.front().size()));
        for (std::size_t row = 0; row < rows; ++row)
            for (std::size_t column = 0; column < left.front().size(); ++column)
                matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column))
                    = (row + skip < left.size() ? left[row + skip] : right[row + skip - left.size()])[column];
        Eigen::FullPivLU<Eigen::MatrixXd> lu(matrix);
        return static_cast<std::size_t>(lu.setThreshold(tolerance).rank());
    };
    return left.size() == right.size() && rank(left.size(), 0) == left.size() && rank(right.size(), left.size()) == left.size()
        && rank(2 * left.size(), 0) == left.size();
}

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
                        accepted = true;
                        jacobian_rank = compiled.jacobian_rank;
                        generic_nullity = compiled.generic_nullity;
                        row_count = compiled.rows.size();
                        unknown_count = compiled.unknown_count();
                        charge_rows = compiled.charge_row_count;
                        for (const auto& item : compiled.problem.phases) support_count += item.support.size();
                        // Structural degeneracy is reported before any nonlinear work; every
                        // lattice all-pinned cell is rank deficient and returns the pre-solve infeasible record.
                        if (!compiled.feasible) {
                            const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
                            solver_status = result.solver_status;
                            CHECK(!result.success && result.residuals.empty() && !result.has_state()
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
                    // A vapor holds only zero-charge components (ALG-EQUILIBRIUM): ions in a lone vapor have no
                    // support, and beside a zero-pinned liquid a free vapor leaves the ion balances nullity 2.
                    const bool electrolyte = chemistry == equilibrium_support::Chemistry::Electrolyte;
                    const bool ions_without_phase = electrolyte && phase_set == equilibrium_support::PhaseSet::OneVapor;
                    const bool ions_on_zero_liquid = electrolyte && phase_set == equilibrium_support::PhaseSet::TwoVaporLiquid
                        && amount_pattern == equilibrium_support::AmountPattern::LastPinned;
                    const bool structurally_square = f == k && phases - k <= rank && !ions_without_phase;
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
                        "accepted iff structural rows are square and ions have a phase", "absolute", 0.0, 0.0,
                        accepted == structurally_square, "c1-structural", case_id.c_str(),
                        "accepted", "phase_set,amount,intensive,chemistry", "accepted=structural",
                        "case");
                    CHECK((f - k) == 0 || !accepted);
                    CHECK(accepted == structurally_square);
                    CHECK(!(accepted && k == phases) || generic_nullity > 0);
                    CHECK(!(accepted && k < phases) || generic_nullity == (ions_on_zero_liquid ? 2U : 0U));
                    CHECK(!(accepted && ions_on_zero_liquid)
                          || (row_count == 12 && unknown_count == 12 && jacobian_rank == 10));
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
    // That mass-proportional vector is not a fraction feed (D1): it is rejected, and the same
    // equimolar mixture is stated as true mass fractions.
    bool unnormalized_rejected = false;
    try { (void)epcsaft::equilibrium::compile_problem(model, mass); }
    catch (const std::invalid_argument& error) {
        unnormalized_rejected = std::string(error.what()).find("mass fractions must sum to one") != std::string::npos;
    }
    CHECK(unnormalized_rejected);
    const double feed_mass = mass.feed.values[0] + mass.feed.values[1];
    for (double& value : mass.feed.values) value /= feed_mass;
    auto molality = amounts;
    molality.feed.kind = epcsaft::equilibrium::FeedKind::Molalities;
    molality.feed.values = {0.0, 1.0 / model.parameters().molar_mass_kg_per_mol[0]};
    molality.feed.solvent = {0};
    const auto reference = epcsaft::equilibrium::compile_problem(model, amounts);
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
                    recording_support::field("feed_kind", equilibrium_support::feed_kind_name(item.feed.kind)),
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
    const bool equal_subspace = same_subspace(
        basis_coordinates, mixed_coordinates, 1.0e-12);
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
    CHECK(compiled.generic_nullity == 0);
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
    // The omitted vapor support is water alone, so the free liquid carries every ion without a charge row.
    recording_support::compare(
        "boundary-charge-row", "charge_row_count", "count", static_cast<double>(boundary.charge_row_count), 0.0,
        recording_support::fields({recording_support::field("T", second.temperature.value),
                                   recording_support::field("P", second.pressure.value),
                                   recording_support::field("x", second.feed.values),
                                   recording_support::field("phase_count", second.phases.size()),
                                   recording_support::field("vapor_amount", second.phases[1].amount.fraction),
                                   recording_support::field("primal_size", boundary.primal_size),
                                   recording_support::field("unknowns", boundary.unknown_count())}),
        "Relationship", "charge-balance", "boundary electrolyte topology with a water vapor adds no charge row",
        "reference_or_one", 0.0, 0.0, boundary.charge_row_count == 0,
        "c8-boundary-charge", "boundary", "N_q", "T,P,\\mathbf{x},topology", "N_q=0", "topology");
    CHECK(compiled.charge_row_count == 0);
    CHECK(boundary.charge_row_count == 0);
    CHECK(boundary.problem.phases[1].support == std::vector<std::size_t>{0});
    CHECK(boundary.primal_size == 8 && boundary.unknown_count() == 11);
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
        epcsaft::equilibrium::maximum_abs(left.evaluate(point).values) - epcsaft::equilibrium::maximum_abs(right.evaluate(permuted).values));
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
    }
    epcsaft::HelmholtzModel relabelled(equilibrium_support::relabel_neutral_binary());
    auto relabelled_problem = first;
    relabelled_problem.feed.values = {first.feed.values[1], first.feed.values[0]};
    const auto relabelled_system = epcsaft::equilibrium::compile_problem(relabelled, relabelled_problem);
    auto relabelled_point = point;
    std::swap(relabelled_point[0], relabelled_point[1]); std::swap(relabelled_point[4], relabelled_point[5]);
    const auto relabelled_values = relabelled_system.evaluate(relabelled_point).values;
    const double reference_norm = epcsaft::equilibrium::maximum_abs(left.evaluate(point).values);
    const double relabelled_error = std::abs(
        epcsaft::equilibrium::maximum_abs(relabelled_system.evaluate(relabelled_point).values)
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
            static_cast<double>(compiled.rows.size()) - static_cast<double>(compiled.unknown_count()), 0.0,
            recording_support::fields({recording_support::field("phase_set", phase_set_name(phase_set)),
                                       recording_support::field("T", problem.temperature.value),
                                       recording_support::field("P", problem.pressure.value),
                                       recording_support::field("x", problem.feed.values),
                                       recording_support::field("phase_count", problem.phases.size())}),
            "Relationship", "rank-count", "single-phase specification is square",
            "reference_or_one", 0.0, 0.0, compiled.rows.size() == compiled.unknown_count(),
            "c9-square", std::to_string(case_index).c_str(), "N_{rows}-N_{unknowns}", "phase,T,P,\\mathbf{x}",
            "N_{rows}=N_{unknowns}", "phase");
        ++case_index;
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
        record_rejection("e2-0/over-sum", "e2-0-pins", model, over,
            "pinned phase fractions exceed one", pin_inputs(over));
        bool generic = false;
        for (const auto& span : capture.trace().spans)
            generic = generic || std::string(span.name) == "equilibrium_generic_point"
                || std::string(span.name) == "boundary_start";
        record_bool("e2-0/over-sum-before-rank", "e2-0-pins",
            "rejection precedes generic-point and start work", !generic, pin_inputs(over));
    }
    auto under = two_phase({1.0, 1.0}, 0.4, 0.5);
    record_rejection("e2-0/all-pinned-under-sum", "e2-0-pins", model, under,
        "all-pinned phase fractions must sum to one", pin_inputs(under));
    auto bad = two_phase({1.0, 1.0}, 1.2, 0.0);
    bad.phases[1].amount = AmountSpec::free();
    bad.pressure = FreeValue::fixed(1.0e5);
    record_rejection("e2-0/above-one", "e2-0-pins", model, bad,
        "invalid pinned phase fraction in phase L", pin_inputs(bad));
    // Sum-consistent but rank-deficient: pin admission passes and the pre-solve record is returned.
    auto exact = two_phase({1.0, 1.0}, 0.5, 0.5);
    const auto compiled = epcsaft::equilibrium::compile_problem(model, exact);
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, exact);
    record_bool("e2-0/exact-sum-not-run", "e2-0-pins",
        "exact all-pinned sum passes admission and returns the pre-solve infeasible record",
        !compiled.feasible && result.convergence.stop == epcsaft::core::Stop::infeasible
            && result.message.find("invalid pinned") == std::string::npos,
        recording_support::fields({pin_inputs(exact),
            recording_support::field("generic_nullity", compiled.generic_nullity)}));
    // Reactive counterexamples: A <-> B conserves moles (rules apply, both orientations).
    for (const double sign : {1.0, -1.0}) {
        auto reactive_over = two_phase({1.0, 0.0}, 0.6, 0.5);
        reactive_over.reactions = {{"A=B", {-sign, sign}, sign > 0.0 ? 2.0 : 0.5, std::nullopt}};
        reactive_over.phases[0].support = {0, 1};
        reactive_over.phases[1].support = {0, 1};
        const auto id = std::string("e2-0/isomerisation-over-sum/") + (sign > 0.0 ? "forward" : "reversed");
        record_rejection(id.c_str(), "e2-0-pins", model, reactive_over,
            "pinned phase fractions exceed one", pin_inputs(reactive_over));
        auto reactive_under = two_phase({1.0, 0.0}, 0.4, 0.5);
        reactive_under.reactions = reactive_over.reactions;
        reactive_under.phases[0].support = {0, 1};
        reactive_under.phases[1].support = {0, 1};
        const auto under_id = std::string("e2-0/isomerisation-under-sum/") + (sign > 0.0 ? "forward" : "reversed");
        record_rejection(under_id.c_str(), "e2-0-pins", model, reactive_under,
            "all-pinned phase fractions must sum to one", pin_inputs(reactive_under));
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
    record_bool("e2-0/doubling-no-sum-rejection", "e2-0-pins",
        "nonconserving A=2B pins are not rejected by a total-mole sum rule",
        doubling_message.find("pinned phase fractions") == std::string::npos
            && doubling_message.find("all-pinned") == std::string::npos && doubling_balance <= 1.0e-12,
        recording_support::fields({pin_inputs(doubling), recording_support::field("balance_residual", doubling_balance),
            recording_support::field("message", doubling_message)}));
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
    record_bool("e2-0/halving-no-sum-rejection", "e2-0-pins",
        "nonconserving 2A=B pin passes the sum stage",
        halving_message.find("pinned phase fractions") == std::string::npos
            && halving_message.find("all-pinned") == std::string::npos && halving_balance <= 1.0e-12,
        recording_support::fields({pin_inputs(halving), recording_support::field("balance_residual", halving_balance),
            recording_support::field("message", halving_message)}));
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
    record_bool("e2-3/omitted-and-full", "e2-3-support", "omitted and explicit full support agree",
        omitted.problem.phases[0].support == std::vector<std::size_t>{0, 1}
            && explicit_full.problem.phases[0].support == std::vector<std::size_t>{0, 1}
            && omitted.removed_species.empty() && explicit_full.removed_species.empty(), inputs);
    auto empty = base;
    empty.phases[0].support_declared = true;
    record_rejection("e2-3/explicit-empty", "e2-3-support", model, empty, "empty support in phase L", inputs);
    auto partial = base;
    partial.phases[0].support = {0};
    record_rejection("e2-3/missing-species", "e2-3-support", model, partial, "infeasible_inventory", inputs);
    auto duplicate = base;
    duplicate.phases[0].support = {0, 0};
    record_rejection("e2-3/duplicate", "e2-3-support", model, duplicate, "duplicate support index", inputs);
    auto range = base;
    range.phases[0].support = {0, 5};
    record_rejection("e2-3/out-of-range", "e2-3-support", model, range, "support index out of range", inputs);
    auto pruned = full;
    pruned.feed.values = {0.5, 0.0};
    const auto pruned_compiled = epcsaft::equilibrium::compile_problem(model, pruned);
    auto pruned_omitted = base;
    pruned_omitted.feed.values = {0.5, 0.0};
    const auto pruned_omitted_compiled = epcsaft::equilibrium::compile_problem(model, pruned_omitted);
    record_bool("e2-3/zero-feed-pruning", "e2-3-support",
        "zero-feed nonparticipant is removed exactly once from explicit and omitted support",
        pruned_compiled.problem.phases[0].support == std::vector<std::size_t>{0}
            && pruned_compiled.removed_species == std::vector<std::size_t>{1}
            && pruned_omitted_compiled.problem.phases[0].support == std::vector<std::size_t>{0}
            && pruned_omitted_compiled.removed_species == std::vector<std::size_t>{1},
        recording_support::fields({recording_support::field("feed", pruned.feed.values)}));
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
    record_bool("e2-2/orientation", "e2-2-reaction",
        "support, rank, balance subspace and rows are independent of reaction orientation",
        left.problem.phases[0].support == std::vector<std::size_t>{0, 1}
            && right.problem.phases[0].support == std::vector<std::size_t>{0, 1}
            && left.removed_species == std::vector<std::size_t>{2} && right.removed_species == std::vector<std::size_t>{2}
            && left.rank == 1 && right.rank == 1 && same_subspace(left.B, right.B, 1.0e-12)
            && same_rows, inputs);
    auto dependent = forward;
    dependent.reactions.push_back({"twice", {-2.0, 2.0, 0.0}, 4.0, std::nullopt});
    record_rejection("e2-2/dependent", "e2-2-reaction", model, dependent, "dependent reaction twice", inputs);
    // Standard-state falsifier at T = 300 K, nu = (-1, +1, 0), K = 2.
    const double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const double g_b = -R * 300.0 * std::log(2.0);
    auto supplied = equilibrium_support::flash(3);
    supplied.phases[0].kind = epcsaft::equilibrium::PhaseKind::Vapor;  // vapor-side sole root at 300 K, 1e5 Pa (#423)
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
    CHECK(supplied_compiled.rank == gauge_compiled.rank && supplied_compiled.rows.size() == gauge_compiled.rows.size());
    auto zero = supplied;
    zero.standard_gibbs_j_per_mol = {0.0, 0.0, 0.0};
    record_rejection("e2-2/zero-vector-with-K", "e2-2-reaction", model, zero, "inconsistent standard state", inputs);
    auto both = supplied;
    both.reactions[0].affinity_j_per_mol = 12.0;
    record_rejection("e2-2/two-targets", "e2-2-reaction", model, both, "at most one target", inputs);
    auto free_t = supplied;
    free_t.temperature = FreeValue::variable(300.0);
    free_t.pressure = FreeValue::fixed(1.0e5);
    free_t.phases[0].amount = AmountSpec::pinned(1.0);
    record_rejection("e2-2/free-T-point-consistent", "e2-2-reaction", model, free_t,
        "temperature-dependent standard state required", inputs);
    auto affinity_free_t = free_t;
    affinity_free_t.reactions = {{"A=B", {-1.0, 1.0, 0.0}, std::nullopt, g_b}};
    record_bool("e2-2/free-T-affinity-admitted", "e2-2-reaction",
        "a constant affinity target admits a free temperature with a consistent supplied vector",
        rejection(model, affinity_free_t).empty(), inputs);
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
    record_bool("e2-2/charged-permutation", "e2-2-reaction",
        "charged reaction rank, balance subspace, participant support and rows follow species identity",
        charged_compiled.rank == 2 && permuted_compiled.rank == 2
            && charged_compiled.problem.phases[0].support == std::vector<std::size_t>{0, 1, 2}
            && permuted_compiled.problem.phases[0].support == std::vector<std::size_t>{0, 1, 2}
            && same_topology
            && charged_compiled.charge_row_count == permuted_compiled.charge_row_count
            && same_subspace(charged_compiled.B, mapped_basis, 1.0e-12)
            && rejection(salt, charge_feed).find("feed carries net z") != std::string::npos
            && rejection(salt, charge_reaction).find("violates z") != std::string::npos,
        recording_support::fields({recording_support::field("feed", charged.feed.values),
            recording_support::field("nu", charged.reactions[0].stoichiometry), recording_support::field("K", 4.0),
            recording_support::field("permutation", std::vector<double>{2, 0, 1})}));
    // The charge test scales sqrt(eps) by max(1, ||nu||_1): |z.nu| = 3 sqrt(eps) lies above the old signed-sum
    // scale (sum nu = -3 sqrt(eps), so 1) and below the L1 scale (4 + 3 sqrt(eps)), and is admitted.
    const double defect = 3.0 * std::sqrt(std::numeric_limits<double>::epsilon());
    auto near_neutral = charged;
    near_neutral.reactions[0].stoichiometry = {2.0, -1.0, -1.0 - defect};
    record_bool("e2-2/charge-scale-l1", "e2-2-reaction",
        "a reaction charge defect within sqrt(eps) of the L1 stoichiometric scale is admitted",
        rejection(salt, near_neutral).find("violates z") == std::string::npos,
        recording_support::fields({recording_support::field("nu", near_neutral.reactions[0].stoichiometry),
            recording_support::field("z_dot_nu", defect)}));
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
    };
    record_rows("e2-4/one-free", epcsaft::equilibrium::compile_problem(model, base), 0, {});
    auto free_plus_zero = base;
    free_plus_zero.phases.push_back(equilibrium_support::phase("L2", PhaseKind::Liquid, AmountSpec::pinned(0.0)));
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
    record_rejection("e2-4/one-sign", "e2-4-charge", model, one_sign, "infeasible_inventory", inputs);
    const auto solved = epcsaft::equilibrium::solve_equilibrium(model, base);
    double charge = 0.0;
    if (solved.success) for (std::size_t i = 0; i < 3; ++i) charge += model.parameters().charge_number[i] * solved.mole_fractions[i];
    recording_support::compare("e2-4/solved-neutrality", "phase_charge", "1", charge, 0.0, inputs,
        "Relationship", "charge-balance", "solved rowless ionic phase is electroneutral", "strict_absolute",
        epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value, 0.0,
        solved.success && std::abs(charge) <= epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value,
        "e2-4-charge", "solved", "q", "T,P,feed", "q=0", "case");
}

// S1 (#51): physical inputs through the shared admission owner. Manufactured verification only.
namespace s1 {

using Order = std::vector<std::size_t>;
using epcsaft::equilibrium::compile_problem;
constexpr std::array<double, 3> CHARGE = {0.0, 2.0, -1.0};
const Order IDENTITY = {0, 1, 2}, PERMUTED = {2, 0, 1};
const std::vector<double> X = {0.94, 0.02, 0.04}, FEED = {50.0, 1.0, 2.0};

// The water / cation / anion mixture of properties_support::aqueous_electrolyte() as a caller
// states it, with explicit units; position k holds species order[k]. The cation carries +2
// so that an exchanged ionic charge moves z^2 and the feed charge.
epcsaft::ParameterDocument document(const Order& order, const std::array<double, 3>& charge = CHARGE) {
    std::array<std::size_t, 3> pos{};
    for (std::size_t k = 0; k < 3; ++k) pos[order[k]] = k;
    const char* name[3] = {"water", "cation", "anion"};
    const double segment[3] = {1.2047, 1.0, 1.0}, sigma[3] = {2.7927, 2.8232, 2.756};
    const double epsilon[3] = {353.95, 230.0, 170.0}, solvation[3] = {1.5, 1.0, 1.0};
    const double mass[3] = {18.0153, 0.02298977, 0.0354527};
    const char* mass_unit[3] = {"gram / mole", "kilogram / mole", "kg / mol"};
    const double ion_diameter[3] = {0.0, 2.484416, 2.42528}, born[3] = {0.0, 3.445, 4.1};
    epcsaft::ParameterDocument doc;
    doc.component_ids.resize(3);
    doc.site_component = {static_cast<int>(pos[0]), static_cast<int>(pos[0])};
    doc.site_multiplicity = {1, 1};
    doc.base_choice = "pc-saft";
    doc.association_choice = "general-site";
    doc.electrolyte_choice = "born";
    doc.ionic_dispersion_choice = "exclude-same-sign-ion-pairs";
    doc.permittivity_choice = "ion-specific-suppression";
    const auto add = [&](const char* field, std::vector<std::size_t> index, double value, const char* unit) {
        doc.entries.push_back({field + std::string("#") + std::to_string(doc.entries.size()), field,
                               std::move(index), value, unit});
    };
    for (std::size_t c = 0; c < 3; ++c) {
        doc.component_ids[pos[c]] = name[c];
        add("molar_mass", {pos[c]}, mass[c], mass_unit[c]);
        add("charge_number", {pos[c]}, charge[c], "elementary-charge");
        add("segment_count", {pos[c]}, segment[c], "dimensionless");
        add("sigma_constant", {pos[c]}, sigma[c], "angstrom");
        add("dispersion_energy_over_k", {pos[c]}, epsilon[c], "kelvin");
        add("solvation_factor", {pos[c]}, solvation[c], "dimensionless");
        if (c == 0) continue;
        add("packing_diameter", {pos[c]}, ion_diameter[c], "angstrom");
        add("debye_huckel_diameter", {pos[c]}, ion_diameter[c], "angstrom");
        add("born_diameter", {pos[c]}, born[c], "angstrom");
        add("ion_suppression_coefficient", {pos[c]}, 7.01, "dimensionless");
    }
    add("relative_permittivity", {pos[0]}, 78.09, "dimensionless");
    add("sigma_term_amplitude", {0, pos[0]}, 10.11, "angstrom");
    add("sigma_term_exponent", {0, pos[0]}, -0.01775, "1 / kelvin");
    add("sigma_term_amplitude", {1, pos[0]}, -1.417, "angstrom");
    add("sigma_term_exponent", {1, pos[0]}, -0.01146, "1 / kelvin");
    add("k_ij", {pos[0], pos[1]}, -0.3, "dimensionless");
    add("k_ij", {pos[2], pos[0]}, 0.0, "dimensionless");
    add("k_ij", {pos[1], pos[2]}, 0.8, "dimensionless");
    add("ionic_region_relative_permittivity", {}, 8.0, "dimensionless");
    add("born_c_shell", {}, 0.0, "dimensionless");
    add("born_c_dielectric", {}, 0.0, "dimensionless");
    add("association_energy_over_k", {0, 1}, 2425.7, "kelvin");
    add("association_volume", {1, 0}, 0.04509, "dimensionless");
    return doc;
}

// The accepted relabelling criterion for associating and ionic mixtures (tests/eos/test_state.py).
bool agrees(const std::string& id, const char* quantity, const char* unit, double value,
           double reference, const char* claim, const std::string& inputs) {
    const bool verdict = recording_support::scalar_verdict(value, reference, "symmetric_close", 1.0e-12, 1.0e-12);
    recording_support::compare(id.c_str(), quantity, unit, value, reference, inputs,
        "Relationship", "input-mapping", claim, "symmetric_close", 1.0e-12, 1.0e-12, verdict, "s1");
    return verdict;
}

// Position i of `state` holds the species at position order[i] of `reference`.
void same_state(const std::string& id, const epcsaft::State& state, const epcsaft::State& reference,
                const Order& order, const char* claim) {
    const auto inputs = recording_support::fields({recording_support::field("T", 298.15),
        recording_support::field("rho", state.molar_density), recording_support::field("reference_x", X)});
    agrees(id + "/P", "pressure", "Pa", state.pressure, reference.pressure, claim, inputs);
    agrees(id + "/a_res", "residual_helmholtz", "1", state.residual_helmholtz, reference.residual_helmholtz, claim, inputs);
    agrees(id + "/Z", "compressibility_factor", "1", state.compressibility_factor, reference.compressibility_factor, claim, inputs);
    for (std::size_t i = 0; i < 3; ++i) {
        agrees(id + "/mu_res/" + std::to_string(i), "residual_chemical_potential_over_rt", "1",
              state.residual_chemical_potential_over_rt[i].value(),
              reference.residual_chemical_potential_over_rt[order[i]].value(), claim, inputs);
    }
}

// Two liquids that both carry the ions: only then does a charge row, and its potential, exist.
epcsaft::equilibrium::EquilibriumProblem two_liquids(const std::vector<double>& feed) {
    auto problem = equilibrium_support::flash(3);
    problem.temperature.value = 298.15;
    problem.phases.push_back(equilibrium_support::phase("L2", epcsaft::equilibrium::PhaseKind::Liquid));
    for (auto& phase : problem.phases) phase.support = {0, 1, 2};
    problem.feed.values = feed;
    return problem;
}

void same_rows(const std::string& id, const epcsaft::equilibrium::CompiledSystem& system,
               const std::vector<double>& point, const epcsaft::equilibrium::CompiledSystem& reference,
               const std::vector<double>& reference_point, const Order& order, const char* claim) {
    using epcsaft::equilibrium::RowKind;
    const auto values = system.evaluate(point).values, expected = reference.evaluate(reference_point).values;
    CHECK(system.charge_row_count == 1 && reference.charge_row_count == 1 && values.size() == expected.size());
    const auto inputs = recording_support::fields({recording_support::field("feed", system.problem.feed.values)});
    for (std::size_t row = 0; row < system.rows.size(); ++row) {
        const auto& item = system.rows[row];
        const bool by_species = item.kind == RowKind::Chemical || item.kind == RowKind::Moiety;
        const std::size_t index = by_species ? order[item.index] : item.index;
        std::size_t match = 0;
        while (match < reference.rows.size() && !(reference.rows[match].kind == item.kind
               && reference.rows[match].phase == item.phase && reference.rows[match].index == index)) ++match;
        CHECK(match < reference.rows.size());
        agrees(id + "/" + std::to_string(row), epcsaft::equilibrium::row_kind_name(item.kind), "1",
              values[row], expected[match], claim, inputs);
    }
}

void check_document() {
    const double rho = 5.7e4;  // compressed liquid: the 1 bar liquid root is 5.67e4 mol/m^3
    epcsaft::HelmholtzModel model(epcsaft::build_parameters(document(IDENTITY)));
    epcsaft::HelmholtzModel permuted(epcsaft::build_parameters(document(PERMUTED)));
    auto qualified = properties_support::aqueous_electrolyte();
    qualified.charge_number = {0, 2, -1};
    const auto state = model.state(298.15, rho, X);
    same_state("s1-f1", state, epcsaft::HelmholtzModel(qualified).state(298.15, rho, X), IDENTITY,
               "document input equals the qualified raw parameters");
    const auto x_permuted = equilibrium_support::permute_vector(X, PERMUTED);
    same_state("s1-f2/state", permuted.state(298.15, rho, x_permuted), state, PERMUTED,
               "species order does not change the physical state");

    const auto reference = compile_problem(model, two_liquids(FEED));
    const auto system = compile_problem(permuted, two_liquids(equilibrium_support::permute_vector(FEED, PERMUTED)));
    auto reference_point = reference.generic_point();
    for (std::size_t k = 0; k < reference.rank + reference.charge_row_count; ++k) {
        reference_point[reference.primal_size + k] = 0.3 + 0.2 * static_cast<double>(k);
    }
    const auto shape = reference.prepare(reference_point).shape, target = system.prepare(reference_point).shape;
    auto point = reference_point;
    for (std::size_t phase = 0; phase < shape.size(); ++phase) {
        for (std::size_t j = 0; j < 3; ++j) point[target[phase].x[j]] = reference_point[shape[phase].x[PERMUTED[j]]];
    }
    for (std::size_t j = 0; j < 3; ++j) point[system.primal_size + j] = reference_point[reference.primal_size + PERMUTED[j]];
    same_rows("s1-f2/row", system, point, reference, reference_point, PERMUTED,
              "species order does not change the equilibrium rows");

    // Only the phase amounts are extensive; ln x, ln eta and the multipliers are not.
    for (const double scale : {1024.0, 1.0e-3}) {
        auto feed = FEED;
        for (double& value : feed) value *= scale;
        auto scaled_point = reference_point;
        for (const auto& layout : shape) scaled_point[*layout.amount] *= scale;
        same_rows("s1-f3/" + recording_support::number(scale), compile_problem(model, two_liquids(feed)),
                  scaled_point, reference, reference_point, IDENTITY, "feed scale leaves the rows unchanged");
    }

    // Negative controls: each wrong mapping must move a value far beyond the criterion.
    const auto moved = [&](const char* id, const char* claim, double value, double reference) {
        const double change = std::abs(value - reference) / std::abs(reference);
        return e2::record_bool(id, "s1", claim, change > 1.0e-6,
            recording_support::fields({recording_support::field("x", X), recording_support::field("relative_change", change)}));
    };
    CHECK(moved("s1-f2n/misplaced-x", "composition permuted without its document moves P",
                model.state(298.15, rho, x_permuted).pressure, state.pressure));
    epcsaft::HelmholtzModel exchanged(epcsaft::build_parameters(document(IDENTITY, {0.0, -1.0, 2.0})));
    CHECK(moved("s1-f2n/exchanged-charge-state", "exchanged ionic charges move a_res",
                exchanged.state(298.15, rho, X).residual_helmholtz, state.residual_helmholtz));
    CHECK(e2::record_rejection("s1-f2n/exchanged-charge-feed", "s1", exchanged, two_liquids(FEED),
                               "feed carries net z", recording_support::fields({recording_support::field("feed", FEED)})));
}

void check_feed() {
    using epcsaft::equilibrium::FeedKind;
    epcsaft::HelmholtzModel model(epcsaft::build_parameters(document(IDENTITY)));
    const auto inputs = [](const epcsaft::equilibrium::EquilibriumProblem& problem) {
        return recording_support::fields({recording_support::field("feed", problem.feed.values),
            recording_support::field("feed_kind", equilibrium_support::feed_kind_name(problem.feed.kind))});
    };
    auto fractions = two_liquids(X);
    fractions.feed.kind = FeedKind::MoleFractions;
    CHECK(compile_problem(model, fractions).n_feed == X);
    fractions.feed.values[0] += 1.0e-13;
    CHECK(compile_problem(model, fractions).n_feed == fractions.feed.values);
    fractions.feed.values[0] = X[0] + 1.0e-9;
    CHECK(e2::record_rejection("s1-f5/mole-sum", "s1", model, fractions, "mole fractions must sum to one", inputs(fractions)));

    // n_i = w_i / M_i against literal kg/mol constants: the document states water in gram / mole.
    const double molar_mass[3] = {0.0180153, 0.02298977, 0.0354527};
    double feed_mass = 0.0;
    for (std::size_t i = 0; i < 3; ++i) feed_mass += FEED[i] * molar_mass[i];
    auto mass = two_liquids(FEED);
    mass.feed.kind = FeedKind::MassFractions;
    for (std::size_t i = 0; i < 3; ++i) mass.feed.values[i] = FEED[i] * molar_mass[i] / feed_mass;
    const auto resolved = compile_problem(model, mass).n_feed;
    for (std::size_t i = 0; i < 3; ++i) {
        agrees("s1-f5/mass/" + std::to_string(i), "feed_amount", "mol/kg", resolved[i],
              FEED[i] / feed_mass, "mass fractions resolve to w_i / M_i", inputs(mass));
    }
    mass.feed.values[0] += 1.0e-9;
    CHECK(e2::record_rejection("s1-f5/mass-sum", "s1", model, mass, "mass fractions must sum to one", inputs(mass)));

    // D2: relative solvent amounts are required for several solvents and keep their meaning.
    epcsaft::HelmholtzModel neutral(properties_support::neutral_ternary());
    const double m0 = 0.016043, m1 = 0.03007;
    auto molality = equilibrium_support::flash(3);
    molality.feed.kind = FeedKind::Molalities;
    molality.feed.solvent = {0, 1};
    molality.feed.values = {0.0, 0.0, 0.5};
    CHECK(e2::record_rejection("s1-f5/solvent-amounts", "s1", neutral, molality,
                               "molality feed omits the relative solvent amounts", inputs(molality)));
    molality.feed.values = {2.0, 0.0, 0.5};
    const auto one_supplied = compile_problem(neutral, molality).n_feed;
    agrees("s1-f5/molality/one-supplied", "feed_amount", "mol", one_supplied[0], 1.0 / m0,
          "a supplied solvent carries the solvent mass", inputs(molality));
    CHECK(one_supplied[1] == 0.0 && one_supplied[2] == 0.5);
    molality.feed.values = {1.0, 3.0, 0.5};
    const auto both_supplied = compile_problem(neutral, molality).n_feed;
    agrees("s1-f5/molality/ratio-0", "feed_amount", "mol", both_supplied[0], 1.0 / (m0 + 3.0 * m1),
          "supplied solvent ratios are rescaled to the solvent mass", inputs(molality));
    agrees("s1-f5/molality/ratio-1", "feed_amount", "mol", both_supplied[1], 3.0 / (m0 + 3.0 * m1),
          "supplied solvent ratios are rescaled to the solvent mass", inputs(molality));

    auto missing = document(IDENTITY), duplicate = missing, unconsumed = missing;
    missing.entries.erase(missing.entries.begin() + 2);
    duplicate.entries.push_back(duplicate.entries[2]);
    unconsumed.permittivity_choice = "solvent-only";
    const auto rejects = [&](const char* id, const epcsaft::ParameterDocument& doc, const char* phrase) {
        std::string message;
        try { (void)epcsaft::build_parameters(doc); }
        catch (const std::invalid_argument& error) { message = error.what(); }
        return e2::record_bool(id, "s1", phrase, message.find(phrase) != std::string::npos,
            recording_support::fields({recording_support::field("message", message)}));
    };
    CHECK(rejects("s1-f5/builder-missing", missing, "missing component/water/segment_count"));
    CHECK(rejects("s1-f5/builder-duplicate", duplicate, "duplicate"));
    CHECK(rejects("s1-f5/builder-unconsumed", unconsumed, "unconsumed"));
}

}  // namespace s1

// S2 (#52): source-to-native reaction reference mapping. Manufactured coefficients only: no
// application constant is stated here, and an EOS-basis polynomial appears only as the exact
// algebraic image of a declared source basis.
namespace s2 {

constexpr double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
const std::vector<double> NU = {-1.0, 2.0};  // A -> 2B, products positive, dnu = +1
const std::array<double, 3> GRID = {280.0, 350.0, 420.0};

void check_reference_mapping() {
    using epcsaft::equilibrium::EOS_STANDARD_STATE_ID;
    const double p0 = 1.0e5, tref = 300.0, dnu = NU[0] + NU[1];
    const epcsaft::ReactionLogPolynomial source{-4.0, 1500.0, 1.2, 0.003, tref, 250.0, 450.0, "manufactured-pressure-basis-P0"};
    // The EOS density-basis (rho0 = 1 mol/m^3) image of the pressure-basis source.
    const epcsaft::ReactionLogPolynomial density{source.a - dnu * std::log(R * tref * 1.0 / p0), source.b,
        source.c - dnu, source.d, tref, 250.0, 450.0, EOS_STANDARD_STATE_ID};
    for (const double t : GRID) {
        // G3: a correlated target carries the series of -R T ln K(T); a constant K equal to K(T) does not.
        const double ln_k = density.a + density.b / t + density.c * std::log(t / tref) + density.d * t;
        const double l1 = -density.b / (t * t) + density.c / t + density.d;
        const double l2 = 2.0 * density.b / (t * t * t) - density.c / (t * t);
        const double l3 = -6.0 * density.b / (t * t * t * t) + 2.0 * density.c / (t * t * t);
        const double analytic[4] = {-ln_k, -(ln_k + t * l1), -t * (2.0 * l1 + t * l2), -t * t * (3.0 * l2 + t * l3)};
        const double constant[4] = {-ln_k, -ln_k, 0.0, 0.0};
        epcsaft::equilibrium::ReactionSpec correlated{"A=2B", NU, std::nullopt, std::nullopt, density};
        epcsaft::equilibrium::ReactionSpec fixed{"A=2B", NU, std::exp(ln_k), std::nullopt, std::nullopt};
        const auto series = epcsaft::equilibrium::reaction_target_series(t, {correlated}, 3);
        const auto flat = epcsaft::equilibrium::reaction_target_series(t, {fixed}, 3);
        const auto inputs = recording_support::fields({recording_support::field("T", t), recording_support::field("ln_K", ln_k),
            recording_support::field("a", density.a), recording_support::field("b", density.b),
            recording_support::field("c", density.c), recording_support::field("d", density.d), recording_support::field("T_ref", tref)});
        double power = 1.0 / (R * t);  // t_k T^(k-1) / R is dimensionless
        for (std::size_t k = 0; k < 4; ++k, power *= t) {
            const std::string point = recording_support::number(t) + "/" + std::to_string(k);
            CHECK(s1::agrees("s2-g3/correlated/" + point, "target_series", "1", series[k][0] * power, analytic[k],
                             "a correlated target carries its temperature series", inputs));
            CHECK(s1::agrees("s2-g3/constant/" + point, "target_series", "1", flat[k][0] * power, constant[k],
                             "a scalar K is a temperature-independent target", inputs));
        }
        const double omitted = std::abs(series[1][0] - flat[1][0]);  // = |Dh / T|, J/(mol K)
        const double criterion = std::max(1.0e-12, 1.0e-12 * std::abs(series[1][0]));
        CHECK(e2::record_bool(("s2-g3n/" + recording_support::number(t)).c_str(), "s2",
            "a constant K omits the reaction enthalpy Dh/T from the first temperature order", omitted / criterion >= 1.0e6,
            recording_support::fields({inputs, recording_support::field("abs_dh_over_T", omitted),
                recording_support::field("ratio_to_criterion", omitted / criterion)})));
        CHECK(s1::agrees("s2-g3n/dh/" + recording_support::number(t), "omitted_first_order", "J/(mol K)",
                         flat[1][0] - series[1][0], R * t * l1, "the omitted term is the van 't Hoff Dh/T", inputs));
    }

    // G4: the unit vocabulary is closed; a unit is rejected, never reinterpreted.
    const auto rejected = [](std::size_t entry, const char* unit, const char* phrase) {
        auto document = s1::document(s1::IDENTITY);
        document.entries[entry].unit = unit;
        std::string message;
        try { (void)epcsaft::build_parameters(document); }
        catch (const std::invalid_argument& error) { message = error.what(); }
        return e2::record_bool((std::string("s2-g4/") + unit).c_str(), "s2", phrase, message.find(phrase) != std::string::npos,
            recording_support::fields({recording_support::field("field", document.entries[entry].field),
                recording_support::field("message", message)}));
    };
    CHECK(rejected(3, "meter", "unit meter for sigma_constant"));  // an unlisted unit
    CHECK(rejected(4, "angstrom", "unit angstrom for dispersion_energy_over_k"));  // a listed unit of the wrong dimension
}

}  // namespace s2

}  // namespace
namespace {
void check_polish_solve_admission() {
    // A zero computed correction and tiny backward residual can hide a 1e-6 forward correction.
    const auto bound = epcsaft::equilibrium::polish_forward_bound;
    const double ill_conditioned = bound({1.0, 0.0, 0.0, 1e-14}, {0.0, 1e-20}, {0.0, 0.0}, {1.0, 1.0});
    const double singular = bound({1.0, 0.0, 0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {1.0, 1.0});
    const double admitted = bound({1.0, 0.0, 0.0, 1.0}, {-1e-14, 0.0}, {1e-14, 0.0}, {1.0, 1.0});
    recording_support::compare("polish-solve/forward", "forward_correction_bound", "1", ill_conditioned, 1e-6,
        "J22:1e-14;r2:1e-20;computed_delta:0", "Identity", "diagonal-system",
        "a tiny backward residual does not bound the forward correction", "reference_scale", 0.0, 1e-14,
        std::abs(ill_conditioned / 1e-6 - 1.0) <= 1e-14, "polish-solve", "diagonal", "forward", "J,r,delta", "J^-1 r", "solve");
    CHECK(ill_conditioned > epcsaft::policy::EQUILIBRIUM_POLISH_CORRECTION.value);
    CHECK(std::isinf(singular));
    CHECK(admitted <= epcsaft::policy::EQUILIBRIUM_POLISH_CORRECTION.value);
}

void check_inventory_support() {
    epcsaft::equilibrium::CompiledSystem c;
    c.n_feed = {2.0, 1.0, 0.0};
    c.feed_total = 3.0;
    c.B = {{1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}};
    c.balance = {2.0, 1.0};
    c.rank = 2;
    epcsaft::equilibrium::PhaseSpec phase;
    phase.support = {0};
    c.problem.phases = {phase};
    for (const bool fraction : {false, true}) {
        const auto r = epcsaft::equilibrium::inventory_bounds(c, {0}, fraction, 0.0)[0];
        CHECK(r.status == "infeasible_inventory" && !r.maximum);
    }
    c.problem.phases[0].support = {0, 1};
    for (const bool fraction : {false, true}) {
        const auto r = epcsaft::equilibrium::inventory_bounds(c, {2}, fraction, 0.0)[0];
        CHECK(r.status == "structural_zero" && r.maximum && *r.maximum == 0.0 && !r.log_fraction);
        recording_support::compare(("inventory/structural/" + std::to_string(fraction)).c_str(),
            "maximum", fraction ? "1" : "mol", *r.maximum, 0.0, "{}", "Identity", "structural-zero",
            "a species absent from every support has zero attainable inventory", "absolute", 0.0, 0.0, true,
            "inventory", "structural", "maximum", "support", "n_2=0", "zero");
    }
}

void check_trace_inventory_witness() {
    using namespace epcsaft;
    using namespace epcsaft::equilibrium;
    // #717 S3's water/ammonia inventory subproblem, with inert trace Li and Na.
    auto parameters = properties_support::parameters(7);
    std::fill(parameters.segment_count.begin(), parameters.segment_count.end(), 1.6069);
    std::fill(parameters.sigma_constant.begin(), parameters.sigma_constant.end(), 3.5206);
    std::fill(parameters.dispersion_energy_over_k.begin(), parameters.dispersion_energy_over_k.end(), 191.42);
    std::fill(parameters.molar_mass_kg_per_mol.begin(), parameters.molar_mass_kg_per_mol.end(), 0.03007);
    const HelmholtzModel model(parameters);
    EquilibriumProblem problem;
    problem.temperature = FreeValue::fixed(298.15);
    problem.pressure = FreeValue::fixed(101325.0);
    problem.feed.kind = FeedKind::Amounts;
    problem.feed.values = {55.50791769892174, 2.089181436265685e-7, 2.773272490012367,
        6.627433285849546e-11, 8.230107589476124e-5, 5.704369091575686e-4, 8.399925272785547e-5};
    PhaseSpec phase;
    phase.kind = PhaseKind::Liquid;
    phase.amount = AmountSpec::free();
    phase.support = {0, 1, 2, 3, 4, 5, 6};
    phase.support_declared = true;
    problem.phases = {phase};
    problem.reactions = {{"water", {-1, 0, 0, 1, 1, 0, 0}, 1.0},
        {"ammonium", {0, 0, 0, 1, 0, 1, -1}, 1.0}};
    const auto compiled = compile_problem(model, problem);
    for (const auto& order : {std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6},
             std::vector<std::size_t>{1, 0, 2, 3, 4, 5, 6}}) {
        const auto bounds = inventory_bounds(compiled, order, false, 0.0);
        const auto& li = *std::find_if(bounds.begin(), bounds.end(), [](const auto& r) { return r.species == 1; });
        CHECK(li.status == "optimal" && li.maximum);
        const double value = li.maximum.value_or(NAN), expected = problem.feed.values[1];
        recording_support::compare_custom(("inventory/trace/" + std::to_string(order.front())).c_str(),
            "maximum_Li", "mol", value, expected,
            recording_support::field("feed_mol", problem.feed.values), "Relationship", "conserved-inventory",
            "inert lithium cannot exceed its feed inventory", "relative<=1e-10",
            std::abs(value / expected - 1.0) <= 1e-10, "inventory");
        CHECK(std::max({li.balance_residual, li.normalization_residual, li.lower_bound_residual}) <= 1e-10);
        CHECK(li.witness.size() == problem.feed.values.size());
    }
}
}
int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.specification.polish-solve", "Newton-solve"}, check_polish_solve_admission},
        {{"native.specification.inventory", "Inventory"}, check_inventory_support},
        {{"native.specification.trace-inventory-witness", "trace-inventory-witness"}, check_trace_inventory_witness},
        {{"native.specification.c1", "C1"}, check_c1},
        {{"native.specification.c2", "C2"}, check_c2},
        {{"native.specification.c3", "C3"}, check_c3},
        {{"native.specification.c4", "C4"}, check_c4},
        {{"native.specification.c5", "C5"}, check_c5},
        {{"native.specification.c6", "C6"}, check_c6},
        {{"native.specification.c7", "C7"}, check_c7},
        {{"native.specification.c8", "C8"}, check_c8},
        {{"native.specification.c9", "C9"}, check_c9},
        {{"native.specification.e2-pins", "E2-0"}, check_e2_pins},
        {{"native.specification.e2-reaction", "E2-2"}, check_e2_reaction},
        {{"native.specification.e2-support", "E2-3"}, check_e2_support},
        {{"native.specification.e2-charge", "E2-4"}, check_e2_charge},
        {{"native.specification.s1-document", "S1-document"}, s1::check_document},
        {{"native.specification.s1-feed", "S1-feed"}, s1::check_feed},
        {{"native.specification.s2-reference", "S2-reference"}, s2::check_reference_mapping},
    });
}
