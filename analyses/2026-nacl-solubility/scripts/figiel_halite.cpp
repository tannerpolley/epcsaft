// Fixed Clarke–Glew Table 34 law and packet /2: #902 binary halite applicability.
#include "equilibrium/problem.hpp"
#include "surface/parameter_builder.hpp"
#include "figiel_inputs.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace epcsaft;
using namespace epcsaft::equilibrium;

double maximum(const std::vector<double>& rows) {
    double v = rows.empty() ? INFINITY : 0;
    for (double r : rows) v = std::isfinite(r) ? std::max(v, std::abs(r)) : INFINITY;
    return v;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const HelmholtzModel model(build_parameters(input_document()));
    const double mw = model.parameters().molar_mass_kg_per_mol[0], R = GAS_CONSTANT_J_PER_MOL_K;
    NeutralReferenceDeclaration neutral;
    neutral.solvent_index = 0;
    neutral.solute_molality_ratios = {0, 1, 1};
    neutral.standard_molality_mol_per_kg = 1;
    std::ofstream out(std::string(argv[1]) + "/attempts.tsv");
    out << std::setprecision(17)
        << "P_Pa\tT_K\tdG_J_mol\tfeed_mol\tstart_mol_kg\tsuccess\tmessage\tm_mol_kg\tw_salt\tsolid_mol\tSI_decade\tmaterial_scaled\tcharge_scaled\tclassified_max\tactivity_contraction_ln\taffinity_difference_ln\trho_mol_m3\tdp_drho\terror_m\trel_error_m\terror_w\trel_error_w\tdoubled_feed_relative\tnumerical_pass\tphysical_pass\tsource_interval_pass\tdelta_m_mol_kg\n";
    std::ofstream diagnostics(std::string(argv[1]) + "/diagnostics.tsv");
    diagnostics << std::setprecision(17) << "attempt\tquantity\tindex\tvalue\n";
    int failures = 0, attempt = 0;
    for (double pressure : {reference_pressure, pinho_pressure}) {
        double central = NAN;
        for (int realization = 0; realization < 5; ++realization) {
            const double t = reference_temperature + (realization == 3 ? -0.1 : realization == 4 ? 0.1 : 0);
            const double shift = realization == 1 ? -g_se : realization == 2 ? g_se : 0;
            const double b = -(dH - reference_temperature * dCp) / R;
            ReactionLogPolynomial law{-(dG + shift) / (R * reference_temperature) - b / reference_temperature,
                b, dCp / R, 0, reference_temperature, 288.15, 308.15, "molality-infinite-dilution"};
            for (double feed : {12., 24.}) {
                if (realization && feed == 24) continue;
                for (double start : {1., 12.}) {
                    out << pressure << '\t' << t << '\t' << shift << '\t' << feed << '\t' << start << '\t';
                    ++attempt;
                    try {
                        EquilibriumProblem p;
                        p.temperature = FreeValue::fixed(t);
                        p.pressure = FreeValue::fixed(pressure);
                        p.feed = {FeedKind::Molalities, {0, feed, feed}, {0}, 1};
                        PhaseSpec liquid;
                        liquid.name = "aqueous";
                        liquid.composition_guess = {1 / mw, start, start};
                        const double total = 1 / mw + 2 * start;
                        for (double& x : liquid.composition_guess) x /= total;
                        p.phases = {liquid};
                        p.source_reference = neutral;
                        SolidSpec solid;
                        solid.name = "halite";
                        solid.stoichiometry = {0, 1, 1};
                        solid.correlation = law;
                        solid.reference = ReactionReference{reference_pressure, ReferenceSourceBasis::CommonMolalityInfiniteDilution};
                        solid.molar_volume_m3_per_mol = solid_volume;
                        p.solids = {solid};
                        const auto result = solve_equilibrium(model, p);
                        for (const auto& item : std::vector<std::pair<std::string, std::vector<double>>>{
                            {"amounts", result.amounts}, {"mole_fractions", result.mole_fractions},
                            {"residuals", result.residuals}, {"classified_residuals", result.classified_residuals},
                            {"solid_amounts", result.solid_amounts}, {"solid_affinity_decades", result.solid_affinity_decades}})
                            for (std::size_t k = 0; k < item.second.size(); ++k)
                                diagnostics << attempt << '\t' << item.first << '\t' << k << '\t' << item.second[k] << '\n';
                        std::string message = result.message;
                        std::replace(message.begin(), message.end(), '\n', ' ');
                        std::replace(message.begin(), message.end(), '\t', ' ');
                        out << result.success << '\t' << message;
                        if (!result.success) { ++failures; out << '\n'; continue; }
                        const auto& x = result.mole_fractions;
                        const double n = result.amounts.at(0), ns = result.solid_amounts.at(0);
                        const double water = n * x[0] * mw, m = n * x[1] / water;
                        const double w = m * salt_mass / (1 + m * salt_mass);
                        const double material = std::max({std::abs(water - 1), std::abs(n*x[1]+ns-feed)/feed,
                            std::abs(n*x[2]+ns-feed)/feed});
                        const double charge = std::abs(n * (x[1]-x[2])) / (2 * feed);
                        const double si = result.solid_affinity_decades.at(0), classified = maximum(result.classified_residuals);
                        const auto ref = infinite_dilution_reference(model, t, pressure, {1,0,0}, {1,2}, {1,1});
                        const auto activity = electrolyte_activity(model, t, pressure, {1,0,0}, {1,2}, {1,1}, m, ref);
                        const auto source = reaction_reference(model, t, pressure, neutral, {0,1,1}, *solid.reference);
                        const auto system = reaction_reference(model, t, pressure, neutral, {0,1,1},
                            ReactionReference{pressure, ReferenceSourceBasis::CommonMolalityInfiniteDilution});
                        if (!source.terminal_contraction || !system.terminal_contraction) throw std::runtime_error("reference unavailable");
                        const double contraction = 2 * (std::log(m) + activity.mean_ionic_log_molality_scale)
                            - reaction_log_value(law, t) + *system.terminal_contraction - *source.terminal_contraction
                            - solid_volume * (pressure-reference_pressure) / (R*t);
                        const double difference = contraction - std::log(10.) * si;
                        const auto state = model.state(t, result.molar_densities.at(0), x);
                        if (realization == 0 && feed == 12 && start == 1) central = m;
                        const double feed_change = realization == 0 ? std::abs(m/central-1) : NAN;
                        const double observed_m = pressure == reference_pressure ? measured_m : measured_w / (salt_mass * (1-measured_w));
                        const double em = m-observed_m, ew = w-measured_w;
                        const bool numerical = std::isfinite(m) && ns >= 0 && material <= 1e-10 && charge <= 1e-10
                            && classified <= 1e-10 && (ns > 0 ? std::abs(si) <= 1e-9 : si <= 1e-9)
                            && std::abs(difference) <= 1e-8 && state.molar_density > 0 && state.pressure_density_derivative > 0
                            && (realization != 0 || feed_change < 1e-8);
                        out << '\t' << m << '\t' << w << '\t' << ns << '\t' << si << '\t' << material << '\t' << charge
                            << '\t' << classified << '\t' << contraction << '\t' << difference << '\t' << state.molar_density
                            << '\t' << state.pressure_density_derivative << '\t' << em << '\t' << em/observed_m
                            << '\t' << ew << '\t' << ew/measured_w << '\t' << feed_change << '\t' << numerical
                            << '\t' << (std::abs(em/observed_m) <= .03)
                            << '\t' << (pressure == reference_pressure ? std::abs(em) <= measured_m_se : std::abs(ew) <= measured_w_u) << '\t' << m-central << '\n';
                        out.flush();
                    } catch (const std::exception& error) {
                        ++failures;
                        out << "\tEXCEPTION: " << error.what() << '\n';
                    }
                }
            }
        }
    }
    return failures ? 1 : 0;
}
