#include "equilibrium/problem.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "equilibrium_support.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {
void check_c1() {
    for (const auto phase_set : equilibrium_support::phase_sets()) {
        for (const auto amount_pattern : equilibrium_support::amount_patterns()) {
            for (const auto intensive_pattern : equilibrium_support::intensive_patterns()) {
                for (const auto chemistry : equilibrium_support::chemistry_values()) {
                    auto problem = equilibrium_support::lattice_problem(
                        phase_set, amount_pattern, intensive_pattern, chemistry);
                    const auto parameters = chemistry == equilibrium_support::Chemistry::Neutral ? model_support::neutral_binary() : chemistry == equilibrium_support::Chemistry::Electrolyte ? model_support::aqueous_electrolyte() : equilibrium_support::ascani_reactive_ternary();
                    epcsaft::HelmholtzModel model(parameters);
                    const int k = equilibrium_support::pinned_phases(phase_set, amount_pattern);
                    const int f = equilibrium_support::free_intensives(intensive_pattern);
                    const int reaction_count = static_cast<int>(epcsaft::equilibrium::reaction_rank(problem.feed.values.size(), problem.reactions));
                    const int rank = static_cast<int>(problem.feed.values.size()) - reaction_count;
                    bool accepted = false;
                    try {
                        const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
                        accepted = compiled.rows_minus_unknowns == 0;
                    } catch (const std::invalid_argument&) {
                        accepted = false;
                    }
                    CHECK((f - k) == 0 || !accepted);
                    const int phases = equilibrium_support::pinned_phases(phase_set, equilibrium_support::AmountPattern::AllFree) == 0 ? (phase_set == equilibrium_support::PhaseSet::OneLiquid || phase_set == equilibrium_support::PhaseSet::OneVapor ? 1 : phase_set == equilibrium_support::PhaseSet::Two || phase_set == equilibrium_support::PhaseSet::TwoVaporLiquid ? 2 : 3) : 1;
                    const bool structurally_square = f == k && k < phases
                        && phases - k <= rank;
                    CHECK(accepted == structurally_square);
                }
            }
        }
    }
}
void check_c4() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
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
    for (const auto& item : {fractions, mass, molality}) {
        const auto compiled = epcsaft::equilibrium::compile_problem(model, item);
        CHECK(compiled.feed_composition.size() == reference.feed_composition.size());
        for (std::size_t i = 0; i < reference.feed_composition.size(); ++i)
            CHECK(std::abs(compiled.feed_composition[i] - reference.feed_composition[i]) < 1e-12);
        CHECK(compiled.rows.size() == reference.rows.size());
        for (std::size_t row = 0; row < reference.rows.size(); ++row) {
            CHECK(compiled.rows[row].kind == reference.rows[row].kind);
            CHECK(compiled.rows[row].phase == reference.rows[row].phase);
            CHECK(compiled.rows[row].index == reference.rows[row].index);
            CHECK(std::string(compiled.rows[row].eqid) == std::string(reference.rows[row].eqid));
        }
    }
}
void check_c5() {
    auto parameters = model_support::parameters(3);
    parameters.segment_count = {1.0, 1.0, 1.0};
    parameters.sigma_constant = {3.0, 3.1, 3.2};
    parameters.dispersion_energy_over_k = {100.0, 110.0, 120.0};
    parameters.molar_mass_kg_per_mol = {0.01, 0.02, 0.03};
    epcsaft::HelmholtzModel model(parameters);
    auto problem = equilibrium_support::flash(3);
    problem.feed.values = {1.0, 1.0, 1.0};
    problem.reactions = {{"exchange", {1.0, -1.0, 0.0}, 1.0, std::nullopt}};
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    CHECK(compiled.rank == 2);
    CHECK(compiled.generic_nullity == 0);
    const std::vector<epcsaft::equilibrium::ReactionSpec> basis = {
        {"r1", {1.0, -1.0, 0.0}, 1.0, std::nullopt},
        {"r2", {0.0, 1.0, -1.0}, 1.0, std::nullopt},
    };
    const std::vector<epcsaft::equilibrium::ReactionSpec> mixed = {
        {"sum", {1.0, 0.0, -1.0}, 1.0, std::nullopt},
        {"difference", {1.0, -2.0, 1.0}, 1.0, std::nullopt},
    };
    CHECK(epcsaft::equilibrium::same_subspace(epcsaft::equilibrium::balance_basis(3, basis), epcsaft::equilibrium::balance_basis(3, mixed), 1.0e-12));
    const auto standard = epcsaft::equilibrium::standard_state(300.0, 3, basis, {});
    CHECK(epcsaft::equilibrium::reaction_affinities(basis, standard) == epcsaft::equilibrium::reaction_affinities(mixed, standard));
    auto dependent = problem;
    dependent.reactions.push_back({"dependent", {2.0, -2.0, 0.0}, 1.0, std::nullopt});
    bool rejected = false;
    try { (void)epcsaft::equilibrium::compile_problem(model, dependent); }
    catch (const std::invalid_argument& error) { rejected = std::string(error.what()).find("dependent") != std::string::npos; }
    CHECK(rejected);
    auto three_phase = equilibrium_support::flash(2);
    three_phase.phases = {equilibrium_support::phase("L1", epcsaft::equilibrium::PhaseKind::Liquid), equilibrium_support::phase("L2", epcsaft::equilibrium::PhaseKind::Liquid), equilibrium_support::phase("L3", epcsaft::equilibrium::PhaseKind::Liquid)};
    rejected = false;
    epcsaft::HelmholtzModel binary(model_support::neutral_binary());
    try { (void)epcsaft::equilibrium::compile_problem(binary, three_phase); }
    catch (const std::invalid_argument& error) { rejected = std::string(error.what()).find("phase count") != std::string::npos; }
    CHECK(rejected);
}
void check_c8() {
    epcsaft::HelmholtzModel model(model_support::aqueous_electrolyte());
    auto problem = equilibrium_support::flash(3);
    problem.feed.values = {0.98, 0.01, 0.01};
    problem.phases.front().support = {0, 1, 2};
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    CHECK(compiled.charge_row_count == 0);
    auto second = problem;
    second.phases.push_back(equilibrium_support::phase(
        "V", epcsaft::equilibrium::PhaseKind::Vapor,
        epcsaft::equilibrium::AmountSpec::pinned(0.0)));
    second.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    const auto boundary = epcsaft::equilibrium::compile_problem(model, second);
    CHECK(boundary.charge_row_count == 1);
}
void check_c2() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
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
    CHECK(named);
    auto omitted = equilibrium_support::flash(2);
    omitted.topology_declared = false;
    try { (void)epcsaft::equilibrium::compile_problem(model, omitted); }
    catch (const std::invalid_argument& error) {
        CHECK(std::string(error.what()).find("phase count") != std::string::npos);
    }
}
void check_c3() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    const auto compiled = epcsaft::equilibrium::compile_problem(model, equilibrium_support::flash(2));
    const auto point = compiled.generic_point();
    const auto values = compiled.evaluate(point).values;
    const auto jacobian = compiled.jacobian(point);
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
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
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
    CHECK(flash_values.size() == boundary_values.size());
    for (std::size_t row = 0; row < flash_values.size(); ++row)
        CHECK(std::abs(flash_values[row] - boundary_values[row]) < 1.0e-12);
}
void check_c7() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
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
    CHECK(left.generic_nullity == right.generic_nullity && a.size() == b.size());
    CHECK(std::abs(left.evaluate(point).residual_norm() - right.evaluate(permuted).residual_norm()) <= 1.0e-14);
    std::sort(a.begin(), a.end()); std::sort(b.begin(), b.end());
    for (std::size_t i = 0; i < a.size(); ++i) CHECK(std::abs(std::abs(a[i]) - std::abs(b[i])) <= 1.0e-14);
    epcsaft::HelmholtzModel relabelled(equilibrium_support::relabel_neutral_binary());
    auto relabelled_problem = first;
    relabelled_problem.feed.values = {first.feed.values[1], first.feed.values[0]};
    const auto relabelled_system = epcsaft::equilibrium::compile_problem(relabelled, relabelled_problem);
    auto relabelled_point = point;
    std::swap(relabelled_point[0], relabelled_point[1]); std::swap(relabelled_point[4], relabelled_point[5]);
    CHECK(relabelled_system.evaluate(relabelled_point).values.size() == a.size());
    CHECK(std::abs(relabelled_system.evaluate(relabelled_point).residual_norm()
                   - left.evaluate(point).residual_norm()) <= 1.0e-14 * left.evaluate(point).residual_norm());
}
void check_c9() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    for (const auto phase_set : {equilibrium_support::PhaseSet::OneLiquid,
                                 equilibrium_support::PhaseSet::OneVapor}) {
        auto problem = equilibrium_support::lattice_problem(
            phase_set, equilibrium_support::AmountPattern::AllFree,
            equilibrium_support::IntensivePattern::BothFixed,
            equilibrium_support::Chemistry::Neutral
        );
        const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
        CHECK(compiled.rows_minus_unknowns == 0);
        CHECK(compiled.generic_nullity == 0);
    }
}

}  // namespace
int main() {
    equilibrium_support::run_check("C1", [] { check_c1(); });
    equilibrium_support::run_check("C2", [] { check_c2(); });
    equilibrium_support::run_check("C3", [] { check_c3(); });
    equilibrium_support::run_check("C4", [] { check_c4(); });
    equilibrium_support::run_check("C5", [] { check_c5(); });
    equilibrium_support::run_check("C6", [] { check_c6(); });
    equilibrium_support::run_check("C7", [] { check_c7(); });
    equilibrium_support::run_check("C8", [] { check_c8(); });
    equilibrium_support::run_check("C9", [] { check_c9(); });
    return equilibrium_support::failures == 0 ? 0 : 1;
}
