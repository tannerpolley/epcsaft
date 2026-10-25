// Native finite-concentration states for Bülow 2021 Figs. 4, 7 and 9: solvent activity coefficients from
// the Engine's electrolyte_activity osmotic coefficient, and sum_i nu_i ln(x_i phi_i) of a 1:1 salt whose
// ions are the last two components, at an explicit composition.
// run.py writes the parameter documents and molalities, compiles this file and reads its CSV:
//   tools/native-driver regression_checks analyses/reproduction/2021-bulow/scripts/bulow2021_states.cpp
//   build/native-drivers/bulow2021_states/driver <job file>
#include "properties/helmholtz_model.hpp"
#include "surface/parameter_builder.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>
using namespace epcsaft;

ParameterDocument read_document(std::istream& in) {
    ParameterDocument doc;
    std::size_t n, sites, entries;
    in >> n >> sites >> entries;
    doc.component_ids.resize(n); doc.site_component.resize(sites); doc.site_multiplicity.resize(sites);
    for (auto& id : doc.component_ids) in >> std::quoted(id);
    for (std::size_t i = 0; i < sites; ++i) in >> doc.site_component[i] >> doc.site_multiplicity[i];
    in >> std::quoted(doc.base_choice) >> std::quoted(doc.association_choice) >> std::quoted(doc.electrolyte_choice)
       >> std::quoted(doc.ionic_dispersion_choice) >> std::quoted(doc.permittivity_choice) >> std::quoted(doc.polar_choice);
    for (std::size_t i = 0; i < entries; ++i) {
        ParameterEntry entry;
        std::size_t width;
        in >> std::quoted(entry.identity) >> std::quoted(entry.field) >> width;
        entry.index.resize(width);
        for (auto& index : entry.index) in >> index;
        in >> entry.magnitude >> std::quoted(entry.unit);
        doc.entries.push_back(entry);
    }
    if (!in) throw std::runtime_error("incomplete parameter document");
    return doc;
}

double log_phi(const State& state, std::size_t i) { return state.log_fugacity_coefficient.value()[i].value(); }

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    std::ifstream in(argv[1]);
    double T, P;
    std::size_t blocks;
    in >> T >> P >> blocks;
    // activity: components (solvent, cation, anion); molality m sets x_ion = mM/(1 + 2mM), x_solvent = 1/(1 + 2mM).
    // sum: any component list with the salt ions last, at the listed mole fractions.
    const std::vector<double> solvent{1.0, 0.0, 0.0};
    std::cout << std::setprecision(17)
              << "label,kind,molality,osmotic_coefficient,log_solvent_activity,log_solvent_activity_fugacity_ratio,"
                 "gamma_solvent_species_basis,gamma_solvent_formula_basis,log_phi_cation,log_phi_anion,sum_log_x_phi,"
                 "bulk_relative_permittivity,molar_density,pressure_residual_relative\n";
    for (std::size_t b = 0; b < blocks; ++b) {
        std::string label;
        in >> std::quoted(label);
        const HelmholtzModel model(build_parameters(read_document(in)));
        const std::size_t n = model.parameters().molar_mass_kg_per_mol.size();
        const double M = model.parameters().molar_mass_kg_per_mol[0];
        // only activity states need the neutral-solvent reference; all-ion sets refuse it (#409)
        std::optional<decltype(infinite_dilution_reference(model, T, P, solvent, {1, 2}, {1, 1}))> reference;
        std::size_t count;
        in >> count;
        for (std::size_t k = 0; k < count; ++k) {
            std::string kind;
            in >> kind;
            if (kind == "activity") {
                double m;
                in >> m;
                if (!in) throw std::runtime_error("incomplete state list");
                const double s = 2.0 * m * M;
                std::cout << std::quoted(label) << ',' << kind << ',' << m << ',';
                if (!reference) reference = infinite_dilution_reference(model, T, P, solvent, {1, 2}, {1, 1});
                const auto activity = electrolyte_activity(model, T, P, solvent, {1, 2}, {1, 1}, m, *reference);
                const double phi = std::get<double>(activity.osmotic_coefficient.value());
                const double log_a = -s * phi;  // practical osmotic coefficient definition
                // identity check: ln(x_s phi_s / phi_s,pure) from the same states
                const double log_a_fugacity = -std::log1p(s) + log_phi(activity.state, 0) - log_phi(*reference->pure_solvent_state, 0);
                std::cout << phi << ',' << log_a << ',' << log_a_fugacity << ','
                          << std::exp(log_a) * (1.0 + s) << ',' << std::exp(log_a) * (1.0 + 0.5 * s) << ",,,,,,\n";
            } else if (kind == "sum") {
                std::vector<double> x(n);
                for (auto& v : x) in >> v;
                if (!in) throw std::runtime_error("incomplete state list");
                const auto state = model.state_at_pressure(T, P, x, DensityBranch::Liquid);
                const double residual = (state.pressure - P) / P;
                if (std::abs(residual) > 1e-6 || state.pressure_density_derivative <= 0.0)
                    throw std::runtime_error("liquid pressure closure failed");
                const double cation = log_phi(state, n - 2), anion = log_phi(state, n - 1);
                std::cout << std::quoted(label) << ",sum,,,,,,," << cation << ',' << anion << ','
                          << std::log(x[n - 2]) + cation + std::log(x[n - 1]) + anion << ',' << state.bulk_relative_permittivity
                          << ',' << state.molar_density << ',' << residual << '\n';
            } else {
                throw std::runtime_error("unknown state kind " + kind);
            }
        }
    }
    return 0;
}
