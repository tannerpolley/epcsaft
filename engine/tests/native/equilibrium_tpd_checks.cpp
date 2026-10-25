// #223: phase removal, tangent-plane and normalized reaction-affinity identities.
#include "equilibrium_checks_shared.hpp"
#include "surface/parameter_builder.hpp"

namespace equilibrium_checks {

// #223 acceptance 8/9: numerical conservation and tangent-plane identities in RT units.
void tpd_compare(const std::string& id, double value, double reference, double tolerance,
                 const std::string& inputs = {}, const char* unit = "1") {
    const bool ok = std::isfinite(value) && std::abs(value - reference) <= tolerance;
    recording_support::compare(id.c_str(), "difference", unit, value, reference, inputs,
        "Relationship", "native-thermodynamic-identity", "#223 numerical falsifier",
        "absolute", tolerance, 0.0, ok);
}

// #662: EqID system_gibbs_helmholtz_coordinates retains N*x_i ideal log-curvature at trace support.
void check_trace_ion_curvature() {
    using namespace epcsaft;
    using namespace epcsaft::equilibrium;
    auto parameters = properties_support::aqueous_electrolyte();
    parameters.permittivity_form = PermittivityForm::DensityTemperature;
    parameters.dielectric_volume_m3_per_mol = {3.777e-4, 0.0, 0.0};
    parameters.dielectric_temperature_k = {1403.0, 0.0, 0.0};
    const HelmholtzModel model(parameters);
    auto phase = equilibrium_support::phase("L", PhaseKind::Liquid); phase.support = {0, 1, 2};
    for (double q : {1e-170, 1e-230}) {
        const std::vector<double> x{1.0 - 2.0 * q, q, q};
        const auto root = solve_density(model, 298.15, 1000.0, x, DensityBranch::Vapor, std::nullopt);
        const auto derivative = model.evaluate(298.15, 1.0 / root.molar_density, x, DerivativeOrder::IsothermalPhaseHessian);
        std::vector<double> shift(3);
        for (std::size_t i = 0; i < 3; ++i)
            shift[i] = -std::log(x[i] * root.molar_density) - derivative.gradient[2 + i];
        std::vector<std::vector<double>> transform;
        const auto trial = compile_stability_trial(model, 298.15, 1000.0, phase, x, shift, transform);
        auto point = trial.generic_point();
        for (std::size_t i = 0; i < 3; ++i) point[i] = std::log(x[i]);
        point[3] = 1.0;
        point[4] = std::log(model.packing_fraction_coefficient(298.15, x) * root.molar_density);
        std::fill(point.begin() + trial.primal_size, point.end(), 0.0);
        core::Capture capture(true);
        core::Scope scope("trace-ion-curvature");
        const auto result = run_driver(model, trial, point, structural_result(trial, scope), scope);
        const auto id = "trace-ion-curvature/" + recording_support::number(q);
        const auto inputs = recording_support::fields({recording_support::field("q", q),
            recording_support::field("residual", result.convergence.residual),
            recording_support::field("iterations", result.iterations), recording_support::field("evaluations", result.evaluations),
            recording_support::field("hessian_calls", result.hessian_calls)});
        double minimum = std::numeric_limits<double>::quiet_NaN();
        for (const auto& event : capture.trace().events)
            if (std::string(event.name) == "reduced_hessian_equilibrated") minimum = event.values[0];
        tpd_compare(id + "/minimum", minimum, 1.0, 1e-8, inputs);
        tpd_compare(id + "/stationarity", result.convergence.residual, 0.0, 1e-10, inputs);
        CHECK(result.success);
    }
}

// #664: manufactured isomerisation; one-row closed form for ALG-STABILITY amount correction.
void check_stability_amount_correction() {
    using namespace epcsaft;
    using namespace epcsaft::equilibrium;
    const HelmholtzModel model(properties_support::neutral_binary());
    constexpr double temperature = 180.0, pressure = 1e5, chemical_scale = 7.0;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature;
    for (const auto& [fraction, scale] : {std::pair{1e-17, 1.0}, {0.37, 10.0}, {0.83, 1.0}}) {
        const std::vector<double> x{fraction, 1.0 - fraction};
        const auto root = solve_density(model, temperature, pressure, x, DensityBranch::Liquid, std::nullopt);
        const double volume = 1.0 / root.molar_density;
        const auto derivative = model.evaluate(temperature, volume, x, DerivativeOrder::IsothermalPhaseHessian);
        auto problem = equilibrium_support::flash(2);
        problem.temperature = FreeValue::fixed(temperature); problem.feed.values = {0.5 * scale, 0.5 * scale};
        for (std::size_t i = 0; i < 2; ++i)
            problem.standard_gibbs_j_per_mol.push_back(rt * (chemical_scale - std::log(x[i] / volume) - derivative.gradient[2 + i]));
        ReactionSpec reaction; reaction.name = "isomerisation"; reaction.stoichiometry = {1.0, -1.0};
        reaction.affinity_j_per_mol = problem.standard_gibbs_j_per_mol[0] - problem.standard_gibbs_j_per_mol[1];
        problem.reactions = {reaction};
        problem.phases[0].composition_guess = x;
        problem.phases[0].packing_guess = model.packing_fraction_coefficient(temperature, x) / volume;
        auto compiled = compile_problem(model, problem);
        const auto parent = solve_compiled(compiled); CHECK(parent.success);
        const auto at = compiled.prepare(parent.coordinates);
        const auto independent = model.evaluate(temperature, at.volume[0], at.composition[0], DerivativeOrder::IsothermalPhaseHessian);
        std::array<double, 2> amount{}, mu{}, intended{}, stored{};
        double rhs = compiled.balance[0];
        const std::string id = "stability-amount-correction/" + recording_support::number(fraction);
        for (std::size_t i = 0; i < 2; ++i) {
            amount[i] = at.amount[0] * at.composition[0][i];
            rhs -= compiled.B[0][i] * amount[i];
            mu[i] = std::log(at.composition[0][i] / at.volume[0]) + independent.gradient[2 + i]
                + problem.standard_gibbs_j_per_mol[i] / rt;
            tpd_compare(id + "/chemical-scale/" + std::to_string(i), mu[i], chemical_scale, 1e-10);
        }
        // Extended arithmetic and the scalar closed form are independent of Eigen's decomposition.
        const long double denominator = static_cast<long double>(amount[0]) * amount[0]
            + static_cast<long double>(amount[1]) * amount[1];
        for (std::size_t i = 0; i < 2; ++i) {
            intended[i] = static_cast<double>(static_cast<long double>(amount[i]) * amount[i]
                * (rhs / compiled.B[0][0]) / denominator);
            stored[i] = amount[i] + intended[i];
        }
        const double total = stored[0] + stored[1];
        double correction = 0.0, energy = 0.0, representation_energy = 0.0;
        for (std::size_t i = 0; i < 2; ++i) {
            const double physical = total * (stored[i] / total);
            const double logged = total * std::exp(std::log(stored[i] / total));
            const double applied = stored[i] - amount[i];
            const double representation = std::abs(physical - stored[i]) + std::abs(logged - physical);
            correction += std::abs(applied) + representation;
            energy += (std::abs(applied) + std::abs(intended[i] - applied) + representation) * std::abs(mu[i]);
            representation_energy += representation * std::abs(mu[i]);
        }
        StabilityOptions options; options.evaluate_supplied_trials = true;
        auto trial = problem.phases[0]; trial.composition_guess = {0.5, 0.5}; options.trials = {trial};
        const auto report = stability_test(compiled, parent, options);
        CHECK(report.completed && report.trials.size() == 1 && energy > 0.0);
        const auto inputs = recording_support::fields({recording_support::field("rhs_mol", rhs),
            recording_support::field("representation_energy", representation_energy)});
        const auto& quantities = report.present_phases[0];
        tpd_compare(id + "/correction", quantities.at("correction_l1_mol"), correction, 1e-12 * correction, inputs, "mol");
        tpd_compare(id + "/energy", quantities.at("feasibility_estimate_over_rt"), energy, 1e-12 * energy, inputs, "mol");
        tpd_compare(id + "/donor-energy", report.trials[0].quantities.at("feasibility_error_estimate_over_rt"), energy, 1e-12 * energy, inputs, "mol");
        auto coordinates = parent.coordinates;
        std::copy(report.candidate_coordinates.begin(), report.candidate_coordinates.end(), coordinates.begin());
        const auto final = compiled.prepare(coordinates);
        CHECK(final.amount[0] > 0.0 && final.composition[0][0] > 0.0 && final.composition[0][1] > 0.0);
        const auto rows = compiled.evaluate(coordinates);
        for (std::size_t k = 0; k < rows.values.size(); ++k)
            if (compiled.rows[k].kind != RowKind::Chemical && compiled.rows[k].kind != RowKind::Pressure)
                tpd_compare(id + "/feasibility/" + std::to_string(k), rows.values[k], 0.0, 1e-12);
        if (fraction == 1e-17) {
            // Positive supported amounts cannot sum to this deliberately impossible target.
            compiled.balance[0] = -1.0;
            auto impossible = parent; impossible.residuals = compiled.evaluate(parent.coordinates).values;
            bool refused = false;
            try { (void)stability_test(compiled, impossible, options); }
            catch (const std::invalid_argument&) { refused = true; }
            tpd_compare(id + "/impossible-target", refused ? 1.0 : 0.0, 1.0, 0.0,
                recording_support::fields({recording_support::field("target_mol", -1.0)}));
        }
    }
}

// #663: manufactured stationary vapor with two independent ionic composition directions.
void check_trial_refinement() {
    using namespace epcsaft;
    using namespace epcsaft::equilibrium;
    const auto source = properties_support::nine_species();
    auto p = properties_support::parameters(4);
    const std::vector<std::size_t> kept{0, 3, 4, 6};
    const auto take = [&](const auto& values) {
        auto result = values; result.clear();
        for (auto i : kept) result.push_back(values.at(i));
        return result;
    };
    p.segment_count = take(source.segment_count); p.sigma_constant = take(source.sigma_constant);
    p.sigma_term_component = {0, 0}; p.sigma_term_amplitude = source.sigma_term_amplitude;
    p.sigma_term_exponent = source.sigma_term_exponent;
    p.dispersion_energy_over_k = take(source.dispersion_energy_over_k);
    p.molar_mass_kg_per_mol = take(source.molar_mass_kg_per_mol); p.charge_number = take(source.charge_number);
    p.packing_diameter_angstrom = take(source.packing_diameter_angstrom);
    p.debye_huckel_diameter_angstrom = take(source.debye_huckel_diameter_angstrom);
    p.born_diameter_angstrom = take(source.born_diameter_angstrom); p.solvation_factor = take(source.solvation_factor);
    p.relative_permittivity = take(source.relative_permittivity);
    p.ion_suppression_coefficient = take(source.ion_suppression_coefficient); p.electrolyte_form = source.electrolyte_form;
    p.ionic_dispersion = source.ionic_dispersion; p.ionic_region_relative_permittivity = source.ionic_region_relative_permittivity;
    for (std::size_t i = 0; i < 4; ++i) for (std::size_t j = 0; j < 4; ++j)
        p.k_ij[i * 4 + j] = source.k_ij[kept[i] * 9 + kept[j]];
    p.site_component = {0, 0}; p.site_multiplicity = {1, 1}; p.association_active = {0, 1, 1, 0};
    p.association_energy_over_k = {0, 2425.7, 2425.7, 0}; p.association_volume = {0, .04509, .04509, 0};
    p.permittivity_form = PermittivityForm::DensityTemperature;
    p.dielectric_volume_m3_per_mol = {3.777e-4, 0.0, 0.0, 0.0};
    p.dielectric_temperature_k = {1403.0, 0.0, 0.0, 0.0};
    const HelmholtzModel model(p);
    const double temperature = 298.15, pressure = 1000.0, a = 1e-30, b = 1e-230;
    const std::vector<double> target{1.0 - 2.0 * a - 3.0 * b, a + 2.0 * b, a, b}, seed{.99, .006, .002, .002};
    auto phase = equilibrium_support::phase("V", PhaseKind::Vapor); phase.support = {0, 1, 2, 3};
    const auto root = solve_density(model, temperature, pressure, target, DensityBranch::Vapor, std::nullopt);
    const auto derivative = model.evaluate(temperature, 1.0 / root.molar_density, target, DerivativeOrder::IsothermalPhaseHessian);
    std::vector<double> shift(4);
    for (std::size_t i = 0; i < 4; ++i) shift[i] = -std::log(target[i] * root.molar_density) - derivative.gradient[2 + i];
    std::vector<std::vector<double>> transform;
    const auto trial = compile_stability_trial(model, temperature, pressure, phase, seed, shift, transform);
    const auto start = solve_density(model, temperature, pressure, seed, DensityBranch::Vapor, std::nullopt);
    auto point = trial.generic_point();
    for (std::size_t i = 0; i < 4; ++i) point[i] = std::log(seed[i]);
    point[4] = 1.0; point[5] = std::log(model.packing_fraction_coefficient(temperature, seed) * start.molar_density);
    std::fill(point.begin() + trial.primal_size, point.end(), 0.0);
    core::Capture capture(true); core::Scope scope("trial-refinement");
    std::vector<double> attempt;
    const auto result = run_driver(model, trial, point, structural_result(trial, scope), scope, &attempt, nullptr, nullptr, false, {}, false, true);
    const auto raw = trial.prepare(attempt, true);
    const auto composition = reconstruct_trial_composition(raw.composition[0], phase.support, p.charge_number, phase.kind);
    const auto inputs = recording_support::fields({recording_support::field("composition", composition),
        recording_support::field("target", target), recording_support::field("raw_coordinates", attempt),
        recording_support::field("raw_rows", trial.evaluate_at(raw).values),
        recording_support::field("numerical_rows", trial.numerical_rows(raw, trial.evaluate_at(raw).values)),
        recording_support::field("cause", result.message),
        recording_support::field("iterations", result.iterations), recording_support::field("evaluations", result.evaluations)});
    tpd_compare("trial-refinement/stationarity", result.convergence.residual, 0.0, 1e-10, inputs);
    tpd_compare("trial-refinement/charge", std::inner_product(composition.begin(), composition.end(), p.charge_number.begin(), 0.0), 0.0, 1e-12, inputs);
    tpd_compare("trial-refinement/normalization", std::accumulate(composition.begin(), composition.end(), 0.0), 1.0, 1e-12, inputs);
    for (std::size_t i = 1; i < target.size(); ++i)
        tpd_compare("trial-refinement/log-fraction/" + std::to_string(i), std::log(composition[i]), std::log(target[i]), 1e-7, inputs);
    for (double x : composition) CHECK(std::isfinite(x) && x > 0.0);
    CHECK(result.success);
}

void check_phase_removal_inventory() {
    using namespace epcsaft::equilibrium;
    // Gross/Sadowski Table 2 and k_ij=0, gross-2001-methane-ethane/3 packet.
    const epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto one = equilibrium_support::flash(2);
    one.phases[0].kind = PhaseKind::Vapor;
    one.phases[0].support = {0, 1};
    one.solver_options.values[epcsaft::policy::row] = 1e-10;
    const auto reference = solve_equilibrium(model, one);
    CHECK(reference.success);
    auto two = one;
    two.phases.push_back(two.phases[0]); two.phases[1].name = "V2";
    const auto compiled = compile_problem(model, two);
    for (const bool disappearing : {true, false}) {
        const std::string id = disappearing ? "removal/disappearance" : "removal/coalescence";
        auto point = compiled.generic_point();
        const auto offsets = compiled.phase_offsets();
        const double second = disappearing ? 5e-9 : 0.6;
        for (std::size_t phase = 0; phase < 2; ++phase) {
            point[offsets[phase]] = point[offsets[phase] + 1] = std::log(0.5);
            point[offsets[phase] + 2] = phase == 0 ? 1.0 - second : second;
            point[offsets[phase] + 3] = std::log(reference.packing_fractions[0]);
        }
        const auto event = check_topology(compiled, point);
        CHECK(event.event == (disappearing ? TopologyEvent::PhaseDisappeared : TopologyEvent::PhasesCoalesced));
        CHECK(disappearing ? event.phase == 1 : event.pair == std::make_pair(std::size_t{0}, std::size_t{1}));
        const auto before = compiled.prepare(point);
        std::map<std::string, double> work;
        auto [reduced, seed] = remove_discovery_phase(compiled, point, &work);
        const auto after_system = compile_problem(model, reduced);
        seed.resize(after_system.unknown_count(), 0.0); // The all-free driver recovers physical multipliers.
        const auto after = after_system.prepare(seed);
        CHECK(after.amount.size() == 1);
        tpd_compare(id + "/amount", after.amount[0], 1.0, 1e-12, {}, "mol");
        std::vector<double> inventory(2, 0.0);
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t phase = 0; phase < 2; ++phase)
                inventory[i] += before.amount[phase] * before.composition[phase][i];
            tpd_compare(id + "/before-component" + std::to_string(i), inventory[i], one.feed.values[i], 1e-12, {}, "mol");
            tpd_compare(id + "/after-component" + std::to_string(i), after.amount[0] * after.composition[0][i], inventory[i], 1e-12, {}, "mol");
        }
        for (std::size_t k = 0; k < compiled.rank; ++k) {
            double merged = 0.0;
            for (std::size_t i = 0; i < 2; ++i) merged += compiled.B[k][i] * after.amount[0] * after.composition[0][i];
            tpd_compare(id + "/moiety" + std::to_string(k), merged, compiled.balance[k], 1e-12, {}, "mol");
        }
        const auto final = solve_equilibrium_from_coordinates(model, reduced, seed);
        CHECK(final.success && final.phase_count == 1);
        for (std::size_t i = 0; i < 2; ++i)
            tpd_compare(id + "/final-x" + std::to_string(i), final.mole_fractions[i], reference.mole_fractions[i], 1e-9);
        tpd_compare(id + "/relative-g", (final.objective - reference.objective) / std::max(1.0, std::abs(reference.objective)), 0.0, 1e-10);
        std::vector<double> attempt;
        const auto rejected = solve_equilibrium_from_coordinates(model, two, point, &attempt);
        recording_support::observe((id + "/driver-event").c_str(), "topology_event", "1",
            static_cast<double>(rejected.topology_event), std::nullopt, recording_support::fields({
                recording_support::field("attempt", attempt), recording_support::field("success", rejected.success ? 1U : 0U)}),
            "Numerical", "native-driver", "actual optimizer topology event", "none", "", "", "");
        if (!disappearing) {
            CHECK(rejected.topology_event == TopologyEvent::PhasesCoalesced && !rejected.has_state());
            CHECK(!attempt.empty());
            CHECK(check_topology(compiled, attempt).event == TopologyEvent::PhasesCoalesced);
            const auto preserved = compiled.prepare(attempt);
            tpd_compare(id + "/attempt-inventory", preserved.amount[0] + preserved.amount[1], 1.0, 1e-12, {}, "mol");
        }
    }
    // #423: a failed attempt whose liquid recipient L1 holds the methane(1e-8)/ethane feed at 300 K and 100 Pa, where
    // the sole root 0.0401 mol/m3 lies below the 7087.88 mol/m3 negative-slope sample. V2 disappears; the recipient's
    // reclosure is refused instead of reseeding L1 on that vapor-side root.
    auto gas = one;
    gas.pressure = FreeValue::fixed(100.0); gas.feed.values = {1e-8, 1.0 - 1e-8};
    const auto vapor = solve_equilibrium(model, gas);
    CHECK(vapor.success);
    auto contradicted = gas;
    contradicted.phases = {gas.phases[0], gas.phases[0]};
    contradicted.phases[0].name = "L1"; contradicted.phases[0].kind = PhaseKind::Liquid; contradicted.phases[1].name = "V2";
    const auto system = compile_problem(model, contradicted);
    auto point = system.generic_point();
    for (std::size_t phase = 0; phase < 2; ++phase) {
        const auto offset = system.phase_offsets()[phase];
        point[offset] = std::log(1e-8); point[offset + 1] = std::log1p(-1e-8);
        point[offset + 2] = phase == 0 ? 1.0 - 5e-9 : 5e-9;
        point[offset + 3] = std::log(vapor.packing_fractions[0]);
    }
    const auto event = check_topology(system, point);
    CHECK(event.event == TopologyEvent::PhaseDisappeared && event.phase == 1 && system.problem.phases[1].name == "V2");
    std::string cause;
    try { (void)remove_discovery_phase(system, point); } catch (const std::exception& error) { cause = error.what(); }
    tpd_compare("removal/contradicted-recipient-refused",
        cause.find("liquid branch not admitted by the finite search: the sole stable root 0.04009108") != std::string::npos
            && cause.find("below the negative-slope sample 7087.88") != std::string::npos ? 1.0 : 0.0, 1.0, 0.0,
        recording_support::fields({recording_support::field("dropped", system.problem.phases[1].name),
            recording_support::field("cause", cause)}));
}

void check_tangent_plane_curvature() {
    using namespace epcsaft::equilibrium;
    const epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto phase = equilibrium_support::phase("V", PhaseKind::Vapor);
    phase.support = {0, 1};
    std::vector<std::vector<double>> transform;
    const auto shifted = compile_stability_trial(model, 300.0, 1e5, phase, {0.25, 0.75}, {2.0, -3.0}, transform);
    const auto unshifted = compile_stability_trial(model, 300.0, 1e5, phase, {0.25, 0.75}, {0.0, 0.0}, transform);
    auto full = shifted.generic_point();
    full[0] = std::log(0.25); full[1] = std::log(0.75); full[2] = 1.0; full[3] = std::log(0.1);
    const auto hs = shifted.hessian(full), h0 = unshifted.hessian(full);
    const std::size_t nfull = full.size();
    for (std::size_t i = 0; i < nfull; ++i) for (std::size_t j = 0; j < nfull; ++j) {
        const double term = i == j && i < 2 ? (i == 0 ? 0.5 : -2.25)
            : i == 2 && j < 2 ? (j == 0 ? 0.5 : -2.25)
            : j == 2 && i < 2 ? (i == 0 ? 0.5 : -2.25) : 0.0;
        tpd_compare("plane-curvature/" + std::to_string(i) + "/" + std::to_string(j), hs[i * nfull + j] - h0[i * nfull + j], term, 1e-12);
    }
    auto point = full; point.resize(shifted.primal_size);
    GibbsMinimizer nlp(shifted, full);
    Ipopt::Index variables = 0, constraints = 0, jacobian_entries = 0, hessian_entries = 0;
    Ipopt::TNLP::IndexStyleEnum index_style;
    CHECK(nlp.get_nlp_info(variables, constraints, jacobian_entries, hessian_entries, index_style));
    const std::size_t n = point.size(), m = static_cast<std::size_t>(constraints);
    std::vector<double> weights(m);
    for (std::size_t i = 0; i < m; ++i) weights[i] = 0.2 * (i + 1);
    std::vector<double> gradient(n), jac(m * n), h(hessian_entries);
    std::vector<Ipopt::Index> hr(h.size()), hc(h.size());
    CHECK(nlp.eval_grad_f(variables, point.data(), true, gradient.data()));
    CHECK(nlp.eval_jac_g(variables, point.data(), false, constraints, jacobian_entries, nullptr, nullptr, jac.data()));
    CHECK(nlp.eval_h(variables, point.data(), false, 0.7, constraints, weights.data(), true, hessian_entries, hr.data(), hc.data(), nullptr));
    CHECK(nlp.eval_h(variables, point.data(), false, 0.7, constraints, weights.data(), false, hessian_entries, nullptr, nullptr, h.data()));
    const auto rows = callback_rows(shifted);
    for (double step : {1e-4, 3e-5, 1e-5}) {
        double first_error = 0.0, second_error = 0.0;
        for (std::size_t column = 0; column < n; ++column) {
            auto plus = point, minus = point;
            plus[column] += step; minus[column] -= step;
            const auto up = shifted.evaluate(plus), down = shifted.evaluate(minus);
            const double df = (up.objective - down.objective) / (2.0 * step);
            first_error = std::max(first_error, std::abs(df - gradient[column]) / std::max(1.0, std::abs(gradient[column])));
            for (std::size_t row = 0; row < m; ++row) {
                const double dg = (up.values[rows[row]] - down.values[rows[row]]) / (2.0 * step);
                first_error = std::max(first_error, std::abs(dg - jac[row * n + column]) / std::max(1.0, std::abs(jac[row * n + column])));
            }
            const auto gp = weighted_gradient(shifted, plus, weights, 0.7), gm = weighted_gradient(shifted, minus, weights, 0.7);
            for (std::size_t entry = 0; entry < h.size(); ++entry) if (static_cast<std::size_t>(hc[entry]) == column) {
                const double dh = (gp[hr[entry]] - gm[hr[entry]]) / (2.0 * step);
                second_error = std::max(second_error, std::abs(dh - h[entry]) / std::max(1.0, std::abs(h[entry])));
            }
        }
        const auto inputs = recording_support::fields({recording_support::field("step", step), recording_support::field("point", point)});
        tpd_compare("plane-stencil/first/" + recording_support::number(step), first_error, 0.0, 1e-7, inputs);
        tpd_compare("plane-stencil/second/" + recording_support::number(step), second_error, 0.0, 1e-5, inputs);
    }
}

// EqID tangent_plane_distance and present_phase_tangent_plane_identity (#300).
void check_supplied_tangent_plane(const epcsaft::equilibrium::CompiledSystem& c,
    const epcsaft::equilibrium::EquilibriumResult& central,
    const epcsaft::equilibrium::StabilityReport& minima, const std::string& id, bool stationary = true) {
    using namespace epcsaft::equilibrium;
    const auto& model = *c.model;
    const auto at = c.prepare(central.coordinates);
    const double eps = std::numeric_limits<double>::epsilon();
    for (std::size_t i = 0; i < c.n_feed.size(); ++i) {
        double plane = 0.0;
        for (std::size_t k = 0; k < c.rank; ++k) plane += c.B[k][i] * central.multipliers[k];
        tpd_compare(id + "/chemical-plane/" + std::to_string(i), minima.plane_over_rt[i], plane,
            32.0 * eps * std::max(1.0, std::abs(plane)));
    }
    const auto evaluate = [&](PhaseSpec phase, const std::vector<double>& w, std::optional<double> eta) {
        phase.composition_guess = w; phase.packing_guess = eta;
        StabilityOptions options; options.evaluate_supplied_trials = true;
        options.trials = std::vector<PhaseSpec>{phase};
        return stability_test(c, central, options);
    };
    for (std::size_t i = 0; i < minima.trials.size(); ++i) {
        const auto& t = minima.trials[i];
        if (!t.quantities.count("tpd_over_rt")) continue;
        std::vector<std::vector<double>> transform;
        const auto trial = compile_stability_trial(model, at.temperature, at.pressure, t.phase,
            t.composition, minima.standard_over_rt, transform);
        const auto raw = trial.prepare(t.raw_coordinates, true);
        const auto replay = evaluate(t.phase, t.composition, raw.eta[0]);
        CHECK(replay.trials.size() == 1);
        const auto& e = replay.trials[0]; const auto& q = e.quantities;
        double p0 = 0.0, p1 = 0.0;
        for (auto j : t.phase.support) {
            p0 += t.composition[j] * minima.plane_over_rt[j];
            p1 += e.composition[j] * replay.plane_over_rt[j];
        }
        const double scale = std::max({1.0, std::abs(t.quantities.at("g_over_rt")), std::abs(q.at("g_over_rt")), std::abs(p0), std::abs(p1)});
        const auto inputs = recording_support::fields({recording_support::field("S", scale),
            recording_support::field("raw_eta", raw.eta[0]), recording_support::field("evaluated_eta", q.at("packing_fraction")),
            recording_support::field("input_l1", q.at("input_reconstruction_l1")), recording_support::field("minimum", t.quantities.at("tpd_over_rt"))});
        tpd_compare(id + "/minimum-replay/" + std::to_string(i), q.at("tpd_over_rt"), t.quantities.at("tpd_over_rt"), 32.0 * eps * scale, inputs);
        if (t.optimizer.success) tpd_compare(id + "/optimizer-eta/" + std::to_string(i), t.optimizer.packing_fractions[0], raw.eta[0], 32.0 * eps * raw.eta[0], inputs);
        if (e.composition == t.composition && q.at("packing_fraction") == raw.eta[0])
            tpd_compare(id + "/same-final-point/" + std::to_string(i), q.at("tpd_over_rt"), t.quantities.at("tpd_over_rt"), 0.0, inputs);
        CHECK(e.split_checks.empty() && e.split_coordinates.empty() && !e.resolved_descent);
    }
    if (!stationary) return;
    auto corrected = central.coordinates;
    std::copy(minima.candidate_coordinates.begin(), minima.candidate_coordinates.end(), corrected.begin());
    const auto nearby = c.prepare(corrected);
    for (std::size_t p = 0; stationary && p < at.amount.size(); ++p) if (at.amount[p] > 0.0) {
        for (bool repair : {false, true}) {
            const auto& point = repair ? nearby : at;
            const auto evaluated = evaluate(c.problem.phases[p], point.composition[p], point.eta[p]);
            const auto& t = evaluated.trials[0]; const auto& q = t.quantities;
            const std::string name = id + (repair ? "/corrected/" : "/public/") + std::to_string(p);
            const auto inputs = recording_support::fields({recording_support::field("tpd", q.at("tpd_over_rt")),
                recording_support::field("correction", q.at("donor_correction_over_rt")),
                recording_support::field("pressure_row", q.at("pressure_row")), recording_support::field("input_l1", q.at("input_reconstruction_l1"))});
            tpd_compare(name + "/zero", q.at("tpd_over_rt"), 0.0, q.at("descent_error_estimate_over_rt"), inputs);
            const double v = q.at("volume_m3_per_mol");
            const auto d = model.evaluate(at.temperature, v, t.composition, epcsaft::DerivativeOrder::IsothermalPhaseHessian);
            double chemical = 0.0, charge = 0.0, electric = 0.0;
            std::size_t offset = c.primal_size + c.rank;
            for (const auto& row : c.rows) if (row.kind == RowKind::Electro) {
                if (row.phase == p) electric = central.coordinates[offset];
                ++offset;
            }
            for (auto j : t.phase.support) {
                charge += model.parameters().charge_number[j] * t.composition[j];
                chemical += t.composition[j] * (std::log(t.composition[j] / v)
                    + d.gradient[2 + j] + evaluated.standard_over_rt[j] - evaluated.plane_over_rt[j]
                    - electric * model.parameters().charge_number[j]);
            }
            tpd_compare(name + "/full-identity", q.at("tpd_over_rt"), chemical + electric * charge - q.at("pressure_row"), 32.0 * eps * std::max(1.0, std::abs(q.at("g_over_rt"))), inputs);
        }
        const auto rooted = evaluate(c.problem.phases[p], at.composition[p], std::nullopt);
        const auto& t = rooted.trials[0]; const auto& q = t.quantities;
        const auto root = epcsaft::solve_density(model, at.temperature, at.pressure, t.composition,
            t.phase.kind == PhaseKind::Liquid ? epcsaft::DensityBranch::Liquid : epcsaft::DensityBranch::Vapor, std::nullopt);
        tpd_compare(id + "/root-density/" + std::to_string(p), q.at("molar_density_mol_m3") / root.molar_density, 1.0, 32.0 * eps,
            recording_support::fields({recording_support::field("composition", t.composition), recording_support::field("density", q.at("molar_density_mol_m3"))}));
        tpd_compare(id + "/root-pressure/" + std::to_string(p), q.at("pressure_row"), 0.0, 1e-10);
    }
    auto phase = c.problem.phases[0]; phase.composition_guess = at.composition[0]; phase.packing_guess = at.eta[0];
    auto invalid = phase; invalid.composition_guess[invalid.support.back()] = 0.0;
    StabilityOptions options; options.evaluate_supplied_trials = true; options.trials = std::vector<PhaseSpec>{phase, invalid, phase};
    const auto ordered = stability_test(c, central, options);
    CHECK(ordered.trials.size() == 3 && ordered.trials[0].quantities == ordered.trials[2].quantities);
    CHECK(ordered.trials[1].message.find("nonpositive/nonfinite trial composition") != std::string::npos && ordered.trials[1].quantities.empty());
    tpd_compare(id + "/ordered-repeat", ordered.trials[0].quantities.at("tpd_over_rt"), ordered.trials[2].quantities.at("tpd_over_rt"), 0.0);
    phase.packing_guess = -0.1; options.trials = std::vector<PhaseSpec>{phase};
    const auto unavailable = stability_test(c, central, options);
    CHECK(unavailable.trials[0].message.find("stability point outside packing domain") != std::string::npos && unavailable.trials[0].quantities.empty());
    phase.packing_guess = at.eta[0]; options.trials = std::vector<PhaseSpec>{phase};
    options.standard_uncertainty_over_rt.assign(c.n_feed.size(), 1e-7);
    const auto uncertainty = stability_test(c, central, options);
    tpd_compare(id + "/standard-uncertainty", uncertainty.trials[0].quantities.at("standard_error_estimate_over_rt"), 2e-7, 32.0 * eps * 2e-7);
}

void check_tangent_plane_admission() {
    using namespace epcsaft::equilibrium;
    tpd_compare("witness/pressure-identity", present_phase_residual_distance(0.0, 1e-10), -1e-10, 0.0);
    tpd_compare("witness/uncertain-descent", descent_resolved(-1e-10, -1e-10, 1e-9) ? 1.0 : 0.0, 0.0, 0.0);
    for (double error : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
        tpd_compare("witness/missing-error/" + recording_support::number(error), descent_resolved(-1e-10, -1e-10, error) ? 1.0 : 0.0, 0.0, 0.0);
    auto parameters = properties_support::aqueous_electrolyte();
    parameters.permittivity_form = epcsaft::PermittivityForm::DensityTemperature;
    parameters.dielectric_volume_m3_per_mol = {3.777e-4, 0.0, 0.0};
    parameters.dielectric_temperature_k = {1403.0, 0.0, 0.0};
    const epcsaft::HelmholtzModel model(parameters);
    // #675: the admitted aqueous near-pure trial retains its unresolved slope and error.
    auto aqueous = equilibrium_support::flash(3); aqueous.feed.values = {0.98, 0.01, 0.01};
    aqueous.temperature = FreeValue::fixed(298.15); aqueous.pressure = FreeValue::fixed(1000.0);
    aqueous.solver_options.values[epcsaft::policy::row] = 1e-10;
    const auto aqueous_system = compile_problem(model, aqueous);
    const auto aqueous_state = solve_compiled(aqueous_system); CHECK(aqueous_state.success);
    const auto aqueous_report = stability_test(aqueous_system, aqueous_state);
    const auto near_pure = std::find_if(aqueous_report.trials.begin(), aqueous_report.trials.end(), [](const auto& trial) {
        return trial.seed == "near-pure-0" && trial.phase.support.size() == 3;
    });
    CHECK(near_pure != aqueous_report.trials.end());
    CHECK(near_pure->evidence == "no_resolved_descent" && !near_pure->resolved_descent);
    const auto& aqueous_q = near_pure->quantities;
    CHECK(std::isfinite(aqueous_q.at("tpd_over_rt")) && aqueous_q.at("tpd_over_rt") < 0.0);
    CHECK(std::isfinite(aqueous_q.at("donor_slope_over_rt")) && std::isfinite(aqueous_q.at("descent_error_estimate_over_rt")));
    CHECK(aqueous_q.at("descent_error_estimate_over_rt") >= 0.0
        && aqueous_q.at("donor_slope_over_rt") >= -aqueous_q.at("descent_error_estimate_over_rt"));
    const std::vector<double> raw{0.98, 0.0100000001, 0.0099999999};
    const auto rebuilt = reconstruct_trial_composition(raw, {0, 1, 2}, parameters.charge_number);
    auto off_normalized = raw; for (double& amount : off_normalized) amount *= 1.001;
    const auto normalized = reconstruct_trial_composition(off_normalized, {0, 1, 2}, parameters.charge_number);
    tpd_compare("witness/off-normalized-reconstruction", normalized[0] + normalized[1] + normalized[2], 1.0, 1e-12);
    const auto inputs = recording_support::fields({recording_support::field("raw", raw), recording_support::field("reconstructed", rebuilt)});
    tpd_compare("witness/reconstructed-charge", rebuilt[1] - rebuilt[2], 0.0, 1e-12, inputs);
    tpd_compare("witness/reconstructed-normalization", rebuilt[0] + rebuilt[1] + rebuilt[2], 1.0, 1e-12, inputs);
    // Reevaluate the same physical packing, with the spurious charge-plane term removed.
    auto phase = equilibrium_support::phase("L", PhaseKind::Liquid); phase.support = {0, 1, 2};
    std::vector<std::vector<double>> transform;
    const auto trial = compile_stability_trial(model, 298.15, 1e5, phase, rebuilt, {0.0, -100.0, 100.0}, transform);
    auto point = trial.generic_point();
    for (std::size_t i = 0; i < 3; ++i) point[i] = std::log(raw[i]);
    point[3] = 1.0; point[4] = std::log(0.3);
    const auto raw_energy = trial.evaluate(point);
    const auto no_charge_plane = compile_stability_trial(model, 298.15, 1e5, phase, rebuilt, {0.0, 0.0, 0.0}, transform);
    tpd_compare("witness/raw-apparent-charge-descent", raw_energy.objective - no_charge_plane.evaluate(point).objective, -2e-8, 1e-12);
    for (std::size_t i = 0; i < 3; ++i) point[i] = std::log(rebuilt[i]);
    const auto admitted_energy = trial.evaluate(point);
    tpd_compare("witness/reconstructed-charge-descent", admitted_energy.objective - no_charge_plane.evaluate(point).objective, 0.0, 1e-12);
    CHECK(std::isfinite(raw_energy.objective) && std::isfinite(admitted_energy.objective));
    const auto at = trial.prepare(point);
    tpd_compare("witness/reevaluated-charge-plane", -100.0 * at.composition[0][1] + 100.0 * at.composition[0][2], 0.0, 1e-12,
        recording_support::fields({inputs, recording_support::field("raw_energy_over_rt", raw_energy.objective),
            recording_support::field("reevaluated_energy_over_rt", admitted_energy.objective), recording_support::field("volume_m3_per_mol", at.volume[0])}));
    tpd_compare("support/singleton", reconstruct_trial_composition({7.0, 0.0, 0.0}, {0}, parameters.charge_number)[0], 1.0, 0.0);
    for (const auto& support : std::vector<std::vector<std::size_t>>{{}, {0, 1}}) {
        bool rejected = false;
        try { (void)reconstruct_trial_composition(raw, support, parameters.charge_number); }
        catch (const std::invalid_argument&) { rejected = true; }
        tpd_compare("support/reconstruction-rejected/" + std::to_string(support.size()), rejected ? 1.0 : 0.0, 1.0, 0.0);
    }
    auto problem = equilibrium_support::flash(3); problem.feed.values = {0.98, 0.01, 0.01};
    problem.temperature = FreeValue::fixed(298.15); problem.solver_options.values[epcsaft::policy::row] = 1e-10;
    const auto candidate = solve_equilibrium(model, problem); CHECK(candidate.success);
    const auto compiled = compile_problem(model, problem);
    const auto report = stability_test(compiled, candidate);
    check_supplied_tangent_plane(compiled, candidate, report, "supplied/charged");
    // One declared-state owner (#341): a changed row tolerance and phase name leave the state and its rows
    // unchanged, so stability admits them and reports the same trials; a changed feed is another state.
    auto relabeled = problem; relabeled.solver_options.values[epcsaft::policy::row] = 1e-9; relabeled.phases[0].name = "L2";
    const auto relabeled_report = stability_test(compile_problem(model, relabeled), candidate);
    double report_difference = relabeled_report.trials.size() == report.trials.size() ? 0.0 : 1.0;
    const auto widen = [&](const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() != y.size()) { report_difference = 1.0; return; }
        for (std::size_t i = 0; i < x.size(); ++i) report_difference = std::max(report_difference, std::abs(x[i] - y[i]));
    };
    widen(relabeled_report.plane_over_rt, report.plane_over_rt);
    for (std::size_t i = 0; i < std::min(report.trials.size(), relabeled_report.trials.size()); ++i) {
        widen(relabeled_report.trials[i].composition, report.trials[i].composition);
        const auto& quantities = relabeled_report.trials[i].quantities;
        for (const auto& [key, value] : report.trials[i].quantities)
            if (quantities.count(key)) widen({quantities.at(key)}, {value}); else report_difference = 1.0;
    }
    tpd_compare("declared-state/relabeled-report", report_difference, 0.0, 0.0);
    auto other_feed = problem; other_feed.feed.values = {0.97, 0.015, 0.015};
    std::string refusal;
    try { (void)stability_test(compile_problem(model, other_feed), candidate); }
    catch (const std::invalid_argument& error) { refusal = error.what(); }
    tpd_compare("declared-state/feed-value", refusal == "stability result/declaration mismatch: feed values differ" ? 1.0 : 0.0, 1.0, 0.0,
        recording_support::fields({recording_support::field("refusal", refusal)}));
    CHECK(!report.trials.empty());
    for (std::size_t p = 0; p < report.present_phases.size(); ++p)
        tpd_compare("witness/present-phase/" + std::to_string(p), report.present_phases[p].at("present_identity_defect_over_rt"), 0.0, 1e-10);
    std::size_t physical_trials = 0;
    for (std::size_t i = 0; i < report.trials.size(); ++i) {
        const auto& item = report.trials[i];
        const auto definition = recording_support::fields({recording_support::field("label", equilibrium_support::phase_kind_name(item.phase.kind)),
            recording_support::field("support", std::vector<int>(item.phase.support.begin(), item.phase.support.end())), recording_support::field("seed", item.seed),
            recording_support::field("evidence", item.evidence), recording_support::field("composition", item.composition),
            recording_support::field("association_floor", item.quantities.count("gradient_floor_max_si_columns")
                ? item.quantities.at("gradient_floor_max_si_columns") : std::numeric_limits<double>::quiet_NaN())});
        if (item.phase.kind == PhaseKind::Vapor)
            for (auto species : item.phase.support) CHECK(parameters.charge_number[species] == 0);
        if (item.composition.empty()) continue;
        ++physical_trials;
        double charge = 0.0, sum = 0.0;
        for (auto species : item.phase.support) { CHECK(item.composition[species] > 0.0); charge += parameters.charge_number[species] * item.composition[species]; sum += item.composition[species]; }
        tpd_compare("support/trial-charge/" + std::to_string(i), charge, 0.0, 1e-12, definition);
        tpd_compare("support/trial-normalization/" + std::to_string(i), sum, 1.0, 1e-12, definition);
        if (item.optimizer.success) tpd_compare("witness/trial-dual/" + std::to_string(i), item.quantities.at("trial_identity_defect_over_rt"), 0.0, 1e-10, definition);
    }
    CHECK(physical_trials > 0);
    for (const auto& support : std::vector<std::vector<std::size_t>>{{}, {0, 1}}) {
        auto invalid = phase; invalid.support = support; invalid.support_declared = true;
        StabilityOptions options; options.trials = std::vector<PhaseSpec>{invalid};
        bool rejected = false;
        try { (void)stability_test(compiled, candidate, options); } catch (const std::invalid_argument&) { rejected = true; }
        tpd_compare("support/physics-rejected/" + std::to_string(support.size()), rejected ? 1.0 : 0.0, 1.0, 0.0);
    }
    // #771: supplied charged Vapor points retain their kind and match their own root, never its sibling.
    const std::vector<double> trace{1.0 - 2e-8, 1e-8, 1e-8};
    const auto vapor_root = epcsaft::solve_density(model, 298.15, 1000.0, trace, epcsaft::DensityBranch::Vapor, std::nullopt);
    const auto liquid_root = epcsaft::solve_density(model, 298.15, 1000.0, trace, epcsaft::DensityBranch::Liquid, std::nullopt);
    CHECK(liquid_root.molar_density > vapor_root.molar_density);
    StabilityOptions supplied; supplied.evaluate_supplied_trials = true;
    std::vector<PhaseSpec> points;
    for (double rho : {vapor_root.molar_density, 1.1 * vapor_root.molar_density, liquid_root.molar_density}) {
        auto point = phase; point.kind = PhaseKind::Vapor; point.support_declared = true;
        point.composition_guess = trace; point.packing_guess = model.packing_fraction_coefficient(298.15, trace) * rho;
        points.push_back(point);
    }
    supplied.trials = points;
    const auto evaluated = stability_test(aqueous_system, aqueous_state, supplied);
    CHECK(evaluated.trials.size() == points.size());
    for (std::size_t j = 0; j < points.size(); ++j) {
        const auto& trial = evaluated.trials[j];
        CHECK(trial.phase.kind == PhaseKind::Vapor && trial.phase.support == phase.support);
        const auto id = "charged-vapor/supplied/" + std::to_string(j);
        const auto inputs = recording_support::fields({recording_support::field("composition", trial.composition),
            recording_support::field("raw_rows", trial.raw_residuals), recording_support::field("cause", trial.message)});
        for (const auto& [name, value] : trial.quantities)
            recording_support::observe((id + "/" + name).c_str(), name.c_str(), "native-field-units", value, std::nullopt, inputs,
                "Relationship", "native-supplied-point", "#771 pressure-root admission");
        if (j == 0) {
            CHECK(trial.evidence == "evaluated" && std::isfinite(trial.quantities.at("tpd_over_rt")));
            tpd_compare(id + "/pressure", trial.quantities.at("pressure_pa"), 1000.0, 1e-4, inputs, "Pa");
            CHECK(trial.quantities.at("pressure_density_derivative") > 0.0);
        } else CHECK(trial.evidence == "unresolved" && trial.message.find("pressure-root window") != std::string::npos);
    }
    auto dipolar = parameters; dipolar.permittivity_form = epcsaft::PermittivityForm::DipolarPerturbation;
    dipolar.dipole_scaling_debye2 = {5.199, 5.407, 2.563}; dipolar.correlation_integral_parameter = {.1276, .215, .138};
    const epcsaft::HelmholtzModel dipolar_model(dipolar);
    auto charged = aqueous; charged.phases = {points[0]};
    epcsaft::core::Capture refusal_capture(true); std::string cause;
    try { (void)compile_problem(dipolar_model, charged); }
    catch (const epcsaft::core::DeclarationError& error) { cause = error.what(); }
    CHECK(cause.find("not numerically qualified") != std::string::npos);
    for (const auto& span : refusal_capture.trace().spans) CHECK(std::string(span.name) != "evaluate");
    for (bool referenced : {false, true}) {
        auto omitted = charged; omitted.pressure = FreeValue::fixed(1e5); omitted.feed.values = {1.0, 1e-6, 1e-6};
        if (referenced) {
            ReactionSpec reaction; reaction.stoichiometry = {-1, 1, 1};
            reaction.correlation = epcsaft::ReactionLogPolynomial{std::log(1e-8), 0, 0, 0, 298.15, 250, 350, "analytic-source"};
            reaction.reference = epcsaft::ReactionReference{std::nullopt, epcsaft::ReferenceSourceBasis::RawMoleFractionInfiniteDilution};
            omitted.reactions = {reaction};
            epcsaft::NeutralReferenceDeclaration reference; reference.solvent_index = 0;
            reference.solute_molality_ratios = {0, 1, 1}; reference.standard_molality_mol_per_kg = 1;
            omitted.source_reference = reference;
        }
        epcsaft::core::Capture omitted_capture(true); std::string cause;
        try { (void)compile_problem(model, omitted); }
        catch (const std::invalid_argument& error) { cause = error.what(); }
        CHECK(cause == "charged vapor cannot be the omitted-charge phase");
        std::size_t eos_evaluations = 0, density_solves = 0;
        for (const auto& span : omitted_capture.trace().spans) {
            eos_evaluations += std::string(span.name) == "evaluate";
            density_solves += std::string(span.name) == "solve_density";
        }
        const auto id = std::string("charged-vapor/omitted-charge-phase/") + (referenced ? "referenced/" : "unreferenced/");
        tpd_compare(id + "eos-evaluations", eos_evaluations, 0.0, 0.0);
        tpd_compare(id + "density-solves", density_solves, 0.0, 0.0);
    }
    // A declared electrolyte bubble has ions absent from its first, neutral vapor phase.
    // Explicit trials must still seed their ionic support with strictly positive amounts.
    auto bubble = problem;
    bubble.phases = {equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.0)), phase};
    bubble.phases[0].support = {0}; bubble.phases[1].name = "L";
    // Zero vapor and component balances fix these declared input compositions.
    bubble.phases[0].composition_guess = {1.0};
    bubble.phases[1].composition_guess = bubble.feed.values;
    bubble.pressure = FreeValue::variable(1e5);
    const auto bubble_state = solve_equilibrium(model, bubble); CHECK(bubble_state.success);
    const auto bubble_system = compile_problem(model, bubble);
    StabilityOptions restricted; restricted.trials = std::vector<PhaseSpec>{phase};
    const auto restricted_report = stability_test(bubble_system, bubble_state, restricted);
    check_supplied_tangent_plane(bubble_system, bubble_state, restricted_report, "supplied/bubble", false);
    CHECK(!restricted_report.trials.empty());
    std::size_t positive_trials = 0;
    for (std::size_t i = 0; i < restricted_report.trials.size(); ++i) {
        const auto& item = restricted_report.trials[i];
        if (item.composition.empty()) continue;
        ++positive_trials;
        for (auto species : phase.support) CHECK(item.composition[species] > 0.0);
        tpd_compare("support/absent-first-phase-charge/" + std::to_string(i), item.composition[1] - item.composition[2], 0.0, 1e-12,
            recording_support::fields({recording_support::field("composition", item.composition), recording_support::field("seed", item.seed),
                recording_support::field("label", equilibrium_support::phase_kind_name(item.phase.kind)), recording_support::field("support", std::vector<int>{0, 1, 2})}));
    }
    CHECK(positive_trials > 0);
    const auto charged_parameters = equilibrium_support::permute_parameters(parameters, {1, 2, 0});
    const epcsaft::HelmholtzModel charged_model(charged_parameters);
    const auto charged_trial = compile_stability_trial(charged_model, 298.15, 1e5, phase, {0.01, 0.01, 0.98}, {0.0, 0.0, 0.0}, transform);
    double dual = 0.0;
    for (std::size_t k = 0; k < 2; ++k) dual += (k == 0 ? 2.0 : 4.0) * charged_trial.balance[k];
    tpd_compare("witness/charged-first-dual", dual, 3.0, 1e-12);
    // The numerical trial uses the candidate's physical plane, not the arbitrary
    // multipliers in the arithmetic example or a zero shift that empties solvent.
    auto charged_problem = problem; charged_problem.feed.values = {0.01, 0.01, 0.98};
    const auto charged_system = compile_problem(charged_model, charged_problem);
    const auto charged_candidate = solve_compiled(charged_system); CHECK(charged_candidate.success);
    const auto charged_point = charged_system.prepare(charged_candidate.coordinates);
    auto charged_shift = charged_point.standard.derivative[0];
    for (std::size_t j = 0; j < charged_shift.size(); ++j) for (std::size_t k = 0; k < charged_system.rank; ++k)
        charged_shift[j] -= charged_system.B[k][j] * charged_candidate.multipliers[k];
    const auto actual_trial = compile_stability_trial(charged_model, 298.15, 1e5, phase,
        {0.01, 0.01, 0.98}, charged_shift, transform);
    const auto stationary = solve_compiled(actual_trial); CHECK(stationary.success);
    double stationary_dual = 0.0;
    for (std::size_t k = 0; k < actual_trial.rank; ++k) stationary_dual += stationary.multipliers[k] * actual_trial.balance[k];
    tpd_compare("witness/charged-first-stationary-dual", stationary.objective, stationary_dual, 1e-10);
}

void check_normalized_reaction_affinity() {
    using namespace epcsaft::equilibrium;
    auto parameters = properties_support::parameters(2);
    parameters.segment_count = {1.0, 1.0}; parameters.sigma_constant = {3.0, 3.0};
    parameters.dispersion_energy_over_k = {200.0, 200.0}; parameters.molar_mass_kg_per_mol = {0.030, 0.015};
    const epcsaft::HelmholtzModel model(parameters);
    auto problem = equilibrium_support::flash(2); problem.feed.values = {1.0, 0.0}; problem.solver_options.values[epcsaft::policy::row] = 1e-10;
    problem.reactions = {{"A=2B", {-1.0, 2.0}, 1.0, std::nullopt}};
    const auto candidate = solve_equilibrium(model, problem); CHECK(candidate.success);
    const auto compiled = compile_problem(model, problem);
    StabilityOptions options; auto phase = problem.phases[0]; phase.support = {0, 1}; options.trials = std::vector<PhaseSpec>{phase};
    const auto report = stability_test(compiled, candidate, options);
    check_supplied_tangent_plane(compiled, candidate, report, "supplied/reacting");
    std::size_t stationary_trials = 0;
    for (std::size_t i = 0; i < report.trials.size(); ++i) {
        const auto& item = report.trials[i];
        if (!item.optimizer.success || item.composition.empty()) continue;
        ++stationary_trials;
        const auto derivative = model.evaluate(300.0, item.quantities.at("volume_m3_per_mol"), item.composition, epcsaft::DerivativeOrder::Hessian);
        double affinity = 0.0;
        for (std::size_t species = 0; species < 2; ++species)
            affinity += (species == 0 ? -1.0 : 2.0) * (std::log(item.composition[species] / item.quantities.at("volume_m3_per_mol"))
                + derivative.gradient[2 + species] + report.standard_over_rt[species]);
        tpd_compare("reaction/normalized-affinity/" + std::to_string(i), affinity, item.quantities.at("tpd_over_rt"), 1e-10,
            recording_support::fields({recording_support::field("composition", item.composition), recording_support::field("seed", item.seed)}));
    }
    CHECK(stationary_trials > 0);
    const auto physical = compiled.prepare(candidate.coordinates);
    double affinity = 0.0;
    for (std::size_t species = 0; species < 2; ++species)
        affinity += (species == 0 ? -1.0 : 2.0) * (std::log(physical.composition[0][species] / physical.volume[0])
            + physical.derivative[0].gradient[2 + species] + report.standard_over_rt[species]);
    tpd_compare("reaction/candidate-affinity", affinity, 0.0, 1e-10);
    // Algebraic k=-0.1, pi=lambda*(2,1), mu+s=pi+k*e gives nu.(mu+s)=k*sum(nu).
    const std::array<double, 2> algebraic{2.0, 1.0};
    tpd_compare("reaction/negative-normalized-affinity", -(2.0 * algebraic[0] - 0.1) + 2.0 * (2.0 * algebraic[1] - 0.1), -0.1, 1e-15);
}

// #390: Ascani 2023 esterification, acid/alcohol/ester/water, printed parameters of #387
// (system 1 1-pentanol at 318.15 K, system 2 1-hexanol at 353.6 K). Discovery returns the
// declared lower-G distinct two-liquid state; references are the retained declared solves.
epcsaft::Parameters ascani_2023(int system) {
    const bool s1 = system == 1;
    epcsaft::ParameterDocument d;
    d.component_ids = {"acetic-acid", "alcohol", "ester", "water"};
    d.site_component = {0, 0, 1, 1, 2, 2, 3, 3}; d.site_multiplicity.assign(8, 1);
    d.base_choice = "pc-saft"; d.association_choice = "general-site";
    d.electrolyte_choice = d.permittivity_choice = d.polar_choice = "none";
    const std::array<double, 4> M{60.05196, s1 ? 88.14818 : 102.17476, s1 ? 130.18486 : 144.21144, 18.01528},
        m{1.3402, s1 ? 3.626 : 3.5146, s1 ? 4.7077 : 4.8847, 1.2047}, eps{311.59, s1 ? 247.28 : 262.32, s1 ? 234.57 : 241.42, 353.95},
        sigma{3.8582, s1 ? 3.4508 : 3.6735, s1 ? 3.4729 : 3.5834, 2.7927}, self_e{3044.4, s1 ? 2252.1 : 2538.9, 0.0, 2425.7},
        self_k{0.07555, s1 ? 0.01033 : 0.00575, 0.04509, 0.04509};
    for (std::size_t i = 0; i < 4; ++i) {
        d.entries.push_back({"molar mass", "molar_mass", {i}, M[i], "gram / mole"});
        d.entries.push_back({"charge", "charge_number", {i}, 0.0, "elementary-charge"});
        d.entries.push_back({"segments", "segment_count", {i}, m[i], "dimensionless"});
        d.entries.push_back({"dispersion", "dispersion_energy_over_k", {i}, eps[i], "kelvin"});
        d.entries.push_back({"diameter", "sigma_constant", {i}, sigma[i], "angstrom"});
        d.entries.push_back({"association energy", "association_energy_over_k", {2 * i, 2 * i + 1}, self_e[i], "kelvin"});
        d.entries.push_back({"association volume", "association_volume", {2 * i, 2 * i + 1}, self_k[i], "dimensionless"});
    }
    d.entries.push_back({"water diameter", "sigma_term_amplitude", {0, 3}, 10.11, "angstrom"});
    d.entries.push_back({"water diameter", "sigma_term_exponent", {0, 3}, -0.01775, "1 / kelvin"});
    d.entries.push_back({"water diameter", "sigma_term_amplitude", {1, 3}, -1.417, "angstrom"});
    d.entries.push_back({"water diameter", "sigma_term_exponent", {1, 3}, -0.01146, "1 / kelvin"});
    const std::array<std::size_t, 6> a{0, 0, 0, 1, 1, 2}, b{1, 2, 3, 2, 3, 3};
    const std::array<double, 6> kij = s1 ? std::array<double, 6>{-0.1, -0.1, -0.1247, -0.0095, 0.004804, -0.0228}
        : std::array<double, 6>{-0.033, -0.10218000000000002, -0.1247, -0.0042, 0.032506800000000016, 0.07317500000000007};
    const std::array<double, 6> cross_e = s1 ? std::array<double, 6>{2648.25, 1522.2, 2735.05, 1126.05, 2338.8999999999996, 1212.85}
        : std::array<double, 6>{2791.65, 1522.2, 2735.05, 1269.45, 2482.3, 1212.85};
    const std::array<double, 6> cross_k = s1 ? std::array<double, 6>{0.02784940934583932, 0.05820443410811945,
        0.05686339209809261, 0.02158171005511797, 0.021341087457223958, 0.04455631309050417}
        : std::array<double, 6>{0.020830030355047026, 0.05828606348969585, 0.056849265487925495,
        0.016099303503434192, 0.015798666029097522, 0.044385444116459734};
    for (std::size_t k = 0; k < 6; ++k) {
        d.entries.push_back({"k_ij", "k_ij", {a[k], b[k]}, kij[k], "dimensionless"});
        for (const auto& sites : {std::vector<std::size_t>{2 * a[k], 2 * b[k] + 1}, std::vector<std::size_t>{2 * a[k] + 1, 2 * b[k]}}) {
            d.entries.push_back({"association energy", "association_energy_over_k", sites, cross_e[k], "kelvin"});
            d.entries.push_back({"association volume", "association_volume", sites, cross_k[k], "dimensionless"});
        }
    }
    return epcsaft::build_parameters(d);
}

// Feed [z, 1-z, 0, 0] mol at 1e5 Pa with the Table 4 constant K on the raw mole-fraction basis.
epcsaft::equilibrium::EquilibriumProblem ascani_2023_problem(int system, double z) {
    using namespace epcsaft; using namespace epcsaft::equilibrium;
    const double T = system == 1 ? 318.15 : 353.6;
    auto problem = equilibrium_support::flash(4);
    problem.temperature = FreeValue::fixed(T); problem.pressure = FreeValue::fixed(1e5);
    problem.feed.values = {z, 1.0 - z, 0.0, 0.0}; problem.solver_options.values[epcsaft::policy::row] = 1e-10;
    ReactionSpec reaction; reaction.name = "esterification"; reaction.stoichiometry = {-1, -1, 1, 1};
    reaction.correlation = ReactionLogPolynomial{std::log(system == 1 ? 43.99 : 22.92), 0, 0, 0, T, T, T, "ascani-2023"};
    reaction.reference = ReactionReference{std::nullopt, ReferenceSourceBasis::RawMoleFractionInfiniteDilution};
    problem.reactions = {reaction};
    NeutralReferenceDeclaration reference; reference.solvent_index = 3; reference.solute_molality_ratios = {1, 1, 1, 0};
    reference.standard_molality_mol_per_kg = 1; reference.species_reference.assign(4, SpeciesReference::PureLiquid);
    problem.source_reference = reference;
    return problem;
}

// Reaction-invariant balances B n with B nu = 0: acid+water, alcohol+water, ester-water.
constexpr std::array<std::array<double, 4>, 3> ascani_2023_balances{{{1, 0, 0, 1}, {0, 1, 0, 1}, {0, 0, 1, -1}}};

void check_reacting_liquid_split() {
    using namespace epcsaft; using namespace epcsaft::equilibrium;
    struct Case { int system; double z, g_two, minor; };
    for (const Case& item : {Case{1, 0.66, -1.7399460167689516, 0.007122036772532215},
             Case{2, 0.75, -0.065662761908543296, 0.16528120690812875}, Case{2, 0.76, -0.04111209796253866, 0.13472846849100448},
             Case{2, 0.78, 0.012659963549595034, 0.05694860636520922}}) {
        const std::string id = "reacting-liquid-split/" + std::to_string(item.system) + "/" + std::to_string(item.z).substr(0, 4);
        const HelmholtzModel model(ascani_2023(item.system));
        auto problem = ascani_2023_problem(item.system, item.z);
        problem.phases.clear(); problem.topology_declared = false; problem.discovery = StabilityOptions{};
        const auto result = discover_phases(model, problem);
        CHECK(result.success); CHECK(result.phase_count == 2); CHECK(result.convergence.residual <= 1e-10);
        tpd_compare(id + "/g", result.objective, item.g_two, 1e-10 * std::max(1.0, std::abs(item.g_two)), {}, "mol");
        tpd_compare(id + "/minor-amount", std::min(result.amounts[0], result.amounts[1]), item.minor, 1e-10, {}, "mol");
        double distance = 0.0;
        for (std::size_t i = 0; i < 4; ++i) distance = std::max(distance, std::abs(result.mole_fractions[i] - result.mole_fractions[4 + i]));
        for (std::size_t k = 0; k < 3; ++k) {
            double balance = 0.0;
            for (std::size_t i = 0; i < 4; ++i)
                balance += ascani_2023_balances[k][i] * (result.amounts[0] * result.mole_fractions[i] + result.amounts[1] * result.mole_fractions[4 + i]);
            tpd_compare(id + "/balance/" + std::to_string(k), balance,
                ascani_2023_balances[k][0] * item.z + ascani_2023_balances[k][1] * (1.0 - item.z), 1e-12, {}, "mol");
        }
        CHECK(distance > 0.3);
    }
}

// #407: the declared aq/org starts of the Ascani reproduction, at 1, 100 and 10,000 mol from the cold
// default start and from the accepted one-mol state, return the same lower-G distinct liquids as
// the discovery references above (pre-#243 declared solves agree), with aq water-rich.
void check_declared_reacting_liquid_split() {
    using namespace epcsaft; using namespace epcsaft::equilibrium;
    struct Case { int system; double z, g_two, minor; std::vector<double> aq, org; };
    for (const Case& item : {Case{1, 0.66, -1.7399460167689516, 0.007122036772532215, {.3, .003, .003, .69}, {.4, .05, .27, .28}},
             Case{2, 0.78, 0.012659963549595034, 0.05694860636520922, {.48, .0013, .0055, .51}, {.58, .012, .22, .19}}}) {
        const HelmholtzModel model(ascani_2023(item.system));
        std::vector<double> accepted;  // primal coordinates of the one-mol cold solve
        for (const double N : {1.0, 100.0, 10000.0}) {
            auto problem = ascani_2023_problem(item.system, item.z);
            for (double& n : problem.feed.values) n *= N;
            auto aq = equilibrium_support::phase("aq", PhaseKind::Liquid), org = equilibrium_support::phase("org", PhaseKind::Liquid);
            aq.composition_guess = item.aq; org.composition_guess = item.org;
            problem.phases = {aq, org};
            const auto compiled = compile_problem(model, problem);
            for (const bool cold : {true, false}) {
                const std::string id = "declared-reacting-liquid-split/" + std::to_string(item.system) + "/"
                    + std::to_string(static_cast<int>(N)) + (cold ? "/cold" : "/accepted");
                auto start = accepted;
                if (!cold) for (const std::size_t i : {std::size_t{4}, std::size_t{10}}) start[i] *= N;  // aq, org amounts
                const auto result = cold ? solve_equilibrium(model, problem)
                                         : solve_equilibrium_from_coordinates(model, problem, add_multipliers(compiled, start));
                CHECK(result.success); CHECK(result.convergence.residual <= 1e-10);
                if (accepted.empty()) accepted.assign(result.coordinates.begin(), result.coordinates.begin() + compiled.primal_size);
                tpd_compare(id + "/g", result.objective / N, item.g_two, 1e-10 * std::max(1.0, std::abs(item.g_two)));
                tpd_compare(id + "/aq-amount", result.amounts[0] / N, item.minor, 1e-12);
                CHECK(result.mole_fractions[3] - result.mole_fractions[7] > 0.3);  // aq water-rich, distinct liquids
                for (std::size_t k = 0; k < 3; ++k) {
                    double balance = 0.0;
                    for (std::size_t i = 0; i < 4; ++i)
                        balance += ascani_2023_balances[k][i] * (result.amounts[0] * result.mole_fractions[i] + result.amounts[1] * result.mole_fractions[4 + i]);
                    tpd_compare(id + "/balance/" + std::to_string(k), balance / N,
                        ascani_2023_balances[k][0] * item.z + ascani_2023_balances[k][1] * (1.0 - item.z), 1e-12);
                }
            }
        }
    }
}

// #391/#410: Gross 2002 Fig. 6 butane(0.42)/1-butanol bubble at 473.15 K, k_ij 0.015, from the
// generated start. The inverse-K finite liquid failed here; the feed liquid on the resolved
// liquid side reaches the declared-start state.
void check_gross_2002_bubble() {
    using namespace epcsaft; using namespace epcsaft::equilibrium;
    ParameterDocument d;
    d.component_ids = {"butane", "1-butanol"};
    d.site_component = {1, 1}; d.site_multiplicity = {1, 1};
    d.base_choice = "pc-saft"; d.association_choice = "general-site";
    d.electrolyte_choice = d.permittivity_choice = d.polar_choice = "none";
    const std::array<double, 2> M{58.123, 74.123}, m{2.3316, 2.7515}, eps{222.88, 259.59}, sigma{3.7086, 3.6139};
    for (std::size_t i = 0; i < 2; ++i) {
        d.entries.push_back({"molar mass", "molar_mass", {i}, M[i], "gram / mole"});
        d.entries.push_back({"charge", "charge_number", {i}, 0.0, "elementary-charge"});
        d.entries.push_back({"segments", "segment_count", {i}, m[i], "dimensionless"});
        d.entries.push_back({"dispersion", "dispersion_energy_over_k", {i}, eps[i], "kelvin"});
        d.entries.push_back({"diameter", "sigma_constant", {i}, sigma[i], "angstrom"});
    }
    d.entries.push_back({"association energy", "association_energy_over_k", {0, 1}, 2544.6, "kelvin"});
    d.entries.push_back({"association volume", "association_volume", {0, 1}, 0.006692, "dimensionless"});
    d.entries.push_back({"k_ij", "k_ij", {0, 1}, 0.015, "dimensionless"});
    const HelmholtzModel model(build_parameters(d));
    auto problem = equilibrium_support::boundary();
    problem.temperature = FreeValue::fixed(473.15); problem.pressure = FreeValue::variable(3.9e6);
    problem.feed.values = {0.42, 0.58}; problem.solver_options.values[epcsaft::policy::row] = 1e-8;
    const auto result = solve_equilibrium(model, problem);
    CHECK(result.success); CHECK(result.convergence.residual <= 1e-8);
    tpd_compare("gross-2002-bubble/P", *result.pressure, 3842628.8934729425, 1e-10 * 3842628.8934729425, {}, "Pa");
    tpd_compare("gross-2002-bubble/rho-L", result.molar_densities.at(0), 7161.966436021788, 1e-10 * 7161.966436021788, {}, "mol/m3");
    tpd_compare("gross-2002-bubble/rho-V", result.molar_densities.at(1), 1639.905384348973, 1e-10 * 1639.905384348973, {}, "mol/m3");
    tpd_compare("gross-2002-bubble/y-butane", result.mole_fractions.at(2), 0.7115607362253475, 1e-10, {}, "1");
}

void gross_2002_component(epcsaft::ParameterDocument& d, std::size_t i, double M, double m, double eps, double sigma) {
    d.entries.push_back({"molar mass", "molar_mass", {i}, M, "kg / mol"});
    d.entries.push_back({"charge", "charge_number", {i}, 0.0, "elementary-charge"});
    d.entries.push_back({"segments", "segment_count", {i}, m, "dimensionless"});
    d.entries.push_back({"dispersion", "dispersion_energy_over_k", {i}, eps, "kelvin"});
    d.entries.push_back({"diameter", "sigma_constant", {i}, sigma, "angstrom"});
}

epcsaft::ParameterDocument gross_2002_associating() {
    epcsaft::ParameterDocument d;
    d.base_choice = "pc-saft"; d.association_choice = "general-site";
    d.electrolyte_choice = d.permittivity_choice = d.polar_choice = "none";
    return d;
}

// Gross 2002 Fig. 10 water/1-pentanol, k_ij 0.016.
epcsaft::Parameters gross_2002_water_pentanol() {
    auto d = gross_2002_associating();
    d.component_ids = {"water", "1-pentanol"};
    d.site_component = {0, 0, 1, 1}; d.site_multiplicity = {1, 1, 1, 1};
    const std::array<double, 2> M{0.018015, 0.08815}, m{1.0656, 3.6260}, eps{366.51, 247.28}, sigma{3.0007, 3.4508},
        eps_ab{2500.7, 2252.1}, kappa{0.034868, 0.010319};
    for (std::size_t i = 0; i < 2; ++i) {
        gross_2002_component(d, i, M[i], m[i], eps[i], sigma[i]);
        d.entries.push_back({"association energy", "association_energy_over_k", {2 * i, 2 * i + 1}, eps_ab[i], "kelvin"});
        d.entries.push_back({"association volume", "association_volume", {2 * i, 2 * i + 1}, kappa[i], "dimensionless"});
    }
    d.entries.push_back({"k_ij", "k_ij", {0, 1}, 0.016, "dimensionless"});
    d.entries.push_back({"wolbach-sandler", "association_combining_rule", {0, 3, 0, 1, 2, 3}, 0.0, "dimensionless"});
    d.entries.push_back({"wolbach-sandler", "association_combining_rule", {1, 2, 0, 1, 2, 3}, 0.0, "dimensionless"});
    return epcsaft::build_parameters(d);
}

// #391: Gross 2002 Fig. 10 water(0.87)/1-pentanol dew at 101,300 Pa, k_ij 0.016, temperature
// started at 375 K. The ordinary solve fails; the cooling onset reaches the first liquid on
// cooling, the pentanol-rich declared-start dew, not the lower 369.425 K root. Water 0.8875
// started at 475 K: the ordinary solve succeeds on the metastable 368.061 K root, whose endpoint
// liquid trials descend, so the onset search returns the 369.966 K declared-start dew.
// Gross 2002 row 149 at 100 times its pressure, 1-propanol(0.413)/ethylbenzene, k_ij 0.023, 10,132,000 Pa
// from 392.05 K: at the ordinary dew's 227.63 K the feed's sole stable density root lies above the scan's
// negative-slope sample, so its vapor is refused (#423) and the onset returns that refusal with no state.
void check_gross_2002_cooling_onset() {
    using namespace epcsaft; using namespace epcsaft::equilibrium;
    const HelmholtzModel model(gross_2002_water_pentanol());
    auto e = gross_2002_associating();
    e.component_ids = {"1-propanol", "ethylbenzene"};
    e.site_component = {0, 0}; e.site_multiplicity = {1, 1};
    gross_2002_component(e, 0, 0.060096, 2.9997, 233.40, 3.2522);
    gross_2002_component(e, 1, 0.106167, 3.0799, 287.35, 3.7974);
    e.entries.push_back({"association energy", "association_energy_over_k", {0, 1}, 2276.8, "kelvin"});
    e.entries.push_back({"association volume", "association_volume", {0, 1}, 0.015268, "dimensionless"});
    e.entries.push_back({"k_ij", "k_ij", {0, 1}, 0.023, "dimensionless"});
    const HelmholtzModel propanol(build_parameters(e));
    const auto run = [&](const HelmholtzModel& system, double first, double start, double pressure) {
        EquilibriumProblem problem;
        problem.phases = {equilibrium_support::phase("L", PhaseKind::Liquid, AmountSpec::pinned(0.0)),
                          equilibrium_support::phase("V", PhaseKind::Vapor)};
        problem.temperature = FreeValue::variable(start); problem.pressure = FreeValue::fixed(pressure);
        problem.feed.kind = FeedKind::Amounts; problem.feed.values = {first, 1.0 - first}; problem.solver_options.values[epcsaft::policy::row] = 1e-8;
        return cooling_onset(compile_problem(system, problem));
    };
    const auto onset = [&](double water, double start) {
        auto onset = run(model, water, start, 1.013 * 1e5);
        CHECK(onset.first.success); CHECK(!onset.second.empty()); CHECK(onset.first.convergence.residual <= 1e-8);
        return onset;
    };
    const auto [result, reports] = onset(0.87, 375.0);
    tpd_compare("gross-2002-cooling-onset/T", *result.temperature, 369.6782514610692, 1e-10 * 369.6782514610692, {}, "K");
    tpd_compare("gross-2002-cooling-onset/x-water", result.mole_fractions.at(0), 0.46929206305581567, 1e-10, {}, "1");
    tpd_compare("gross-2002-cooling-onset/rho-L", result.molar_densities.at(0), 13988.729693124995, 1e-10 * 13988.729693124995, {}, "mol/m3");
    const auto [metastable, checks] = onset(0.8875, 475.0);
    CHECK(std::count_if(checks.front().trials.begin(), checks.front().trials.end(), [](const auto& t) { return t.resolved_descent; }) == 1);
    tpd_compare("gross-2002-cooling-onset/T-475", *metastable.temperature, 369.9663336068358, 1e-10 * 369.9663336068358, {}, "K");
    const auto refused = run(propanol, 0.413, 392.05, 1.0132e7).first;
    CHECK(!refused.success); CHECK(!refused.temperature && refused.molar_densities.empty());
    CHECK(refused.message.find("vapor branch not admitted by the finite search: the sole stable root 10325.2153") != std::string::npos);
}

// #423: a sole cold root lies on the other side of its negative-slope sample from the declared phase.
// Declared phases are refused; discovery refuses that initial candidate and relabels an added phase to
// the established side with its Gibbs decrease kept; a finite-split donor on the other side is refused.
// Water(0.87)/1-pentanol at 150 K, 101,300 Pa with a vapor trial family; Gross 2001 methane/ethane z=1e-8
// at 300 K, 100 Pa and z=0.5 at 240 K, 3 MPa. Values from main 4363d103 and #423 slice 2.
void check_density_branch_admission() {
    using namespace epcsaft; using namespace epcsaft::equilibrium;
    const HelmholtzModel water_pentanol(gross_2002_water_pentanol()), methane_ethane(properties_support::neutral_binary());
    const auto problem = [](double T, double P, std::vector<double> z) {
        auto p = equilibrium_support::flash(2);
        p.temperature = FreeValue::fixed(T); p.pressure = FreeValue::fixed(P); p.feed.values = std::move(z);
        return p;
    };
    const auto liquid = problem(150.0, 101300.0, {0.87, 0.13}), gas = problem(300.0, 100.0, {1e-8, 1.0 - 1e-8});
    const auto refused = [](const HelmholtzModel& model, EquilibriumProblem p, PhaseKind kind, const char* cause) {
        p.phases[0].kind = kind;
        const auto result = solve_equilibrium(model, p);
        CHECK(!result.success); CHECK(result.molar_densities.empty()); CHECK(result.message.find(cause) != std::string::npos);
    };
    refused(water_pentanol, liquid, PhaseKind::Vapor, "vapor branch not admitted by the finite search: the sole stable root 35911.87");
    refused(methane_ethane, gas, PhaseKind::Liquid, "liquid branch not admitted by the finite search: the sole stable root 0.04009108");
    const auto discover = [](const HelmholtzModel& model, EquilibriumProblem p) {
        auto trial = equilibrium_support::phase("trial-V", PhaseKind::Vapor); trial.support = {0, 1};
        StabilityOptions options; options.trials = std::vector<PhaseSpec>{trial}; options.declaration_budget = 4;
        p.phases.clear(); p.topology_declared = false; p.discovery = options;
        auto result = solve_equilibrium(model, p);
        CHECK(result.success && result.declaration);
        return result;
    };
    const auto split = discover(water_pentanol, liquid);
    CHECK(split.phase_count == 2);
    for (const auto& phase : split.declaration->phases) CHECK(phase.kind == PhaseKind::Liquid);
    tpd_compare("branch-admission/split-g", split.objective, -18.07210327606597, 1e-10 * 18.07210327606597, {}, "mol");
    tpd_compare("branch-admission/split-rho", split.molar_densities.at(1), 11193.38560166406, 1e-10 * 11193.38560166406, {}, "mol/m3");
    tpd_compare("branch-admission/split-amount", split.amounts.at(1), 0.13612504455029806, 1e-10, {}, "mol");
    const auto vapor = discover(methane_ethane, gas);
    CHECK(vapor.phase_count == 1); CHECK(vapor.declaration->phases[0].kind == PhaseKind::Vapor);
    tpd_compare("branch-admission/vapor-rho", vapor.molar_densities.at(0), 0.040091082039588319, 1e-10 * 0.040091082039588319, {}, "mol/m3");
    // Donor: the declared liquid at 240 K has two stable feed roots. Withdrawing its x_methane=0.4623 liquid trial at
    // half and quarter capacity leaves donors whose sole roots lie below their samples; an eighth is admitted.
    const auto compiled = compile_problem(methane_ethane, problem(240.0, 3e6, {0.5, 0.5}));
    const auto declared = solve_compiled(compiled);
    CHECK(declared.success);
    const auto report = stability_test(compiled, declared, StabilityOptions{});
    const auto donor = std::find_if(report.trials.begin(), report.trials.end(), [](const auto& t) {
        return !t.composition.empty() && std::abs(t.composition[0] - 0.46231861363169502) < 1e-8; });
    CHECK(donor != report.trials.end());
    const auto beta = donor->quantities.find("split_beta_mol");
    tpd_compare("branch-admission/donor-beta", beta == donor->quantities.end() ? NAN : beta->second, 0.11623984312000767,
                1e-10 * 0.11623984312000767, {}, "mol");
    CHECK(donor->resolved_descent && donor->split_checks.size() == 3);
    CHECK(donor->split_checks[0].at("admitted") == 0.0 && donor->split_checks[0].at("component_error_mol") <= 1e-12);
    CHECK(donor->message.find("finite split: liquid branch not admitted by the finite search: the sole stable root 2124.87") != std::string::npos
          && donor->message.find("below the negative-slope sample 8026.59") != std::string::npos);
}

// #647 falsifier 4, EqID michelsen_modified_tangent_plane: J_ij = n dln(phi_i)/dn_j at fixed T, P from one
// IsothermalPhaseHessian evaluation against exact pressure-chart ln(phi_i) actions along e_j - x, and
// Gibbs-Duhem sum_i x_i J_ij = 0. Gross 2001 methane/ethane x = 0.7 on its two stable roots; the associating binary liquid.
void check_tangent_plane_jacobian() {
    using namespace epcsaft;
    struct Point { const char* id; Parameters parameters; double T, P; std::vector<double> x; DensityBranch branch; };
    for (const auto& [id, parameters, T, P, x, branch] : std::vector<Point>{
             {"methane-ethane-220K-4MPa/liquid", properties_support::neutral_binary(), 220.0, 4e6, {0.7, 0.3}, DensityBranch::Liquid},
             {"methane-ethane-220K-4MPa/vapor", properties_support::neutral_binary(), 220.0, 4e6, {0.7, 0.3}, DensityBranch::Vapor},
             {"associating-binary-298K-0.1MPa/liquid", properties_support::associating_binary(), 298.15, 1e5, {0.3, 0.7},
              DensityBranch::Liquid}}) {
        const HelmholtzModel model(parameters);
        const double v = 1.0 / solve_density(model, T, P, x, branch, std::nullopt).molar_density;
        const auto J = equilibrium::tangent_plane_jacobian(model.evaluate(T, v, x, DerivativeOrder::IsothermalPhaseHessian), v, {0, 1});
        const auto inputs = recording_support::fields({recording_support::field("T_K", T), recording_support::field("P_Pa", P),
            recording_support::field("x", x), recording_support::field("v_m3_per_mol", v)});
        for (std::size_t j = 0; j < 2; ++j) {
            PressureActionDirection direction{0.0, 0.0, {-x[0], -x[1]}, {}};
            direction.composition[j] += 1.0;
            const auto batch = model.pressure_actions(T, P, x, branch, {PropertyObservable::LogFugacityCoefficient}, {direction});
            double duhem = 0.0, scale = 1.0;
            for (std::size_t i = 0; i < 2; ++i) {
                const auto& action = batch.results.at(i);
                CHECK(action.component == i && action.action);
                const double value = J[i * 2 + j];
                tpd_compare(std::string("tangent-plane-jacobian/") + id + "/J" + std::to_string(i) + std::to_string(j), value,
                    action.action.value_or(NAN), 1e-10 * std::max(1.0, std::abs(value)), inputs);
                duhem += x[i] * value; scale = std::max(scale, std::abs(value));
            }
            tpd_compare(std::string("tangent-plane-jacobian/") + id + "/gibbs-duhem/" + std::to_string(j), duhem, 0.0, 1e-10 * scale, inputs);
        }
    }
}

}  // namespace equilibrium_checks
