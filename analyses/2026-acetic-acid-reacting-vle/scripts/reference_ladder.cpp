// Independent verification of EqID reference_limit_terminal_pair; no production reference replacement.
#include "properties/helmholtz_model.hpp"
#include "surface/parameter_builder.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    std::ifstream in(argv[1]);
    epcsaft::ParameterDocument doc;
    std::size_t n, sites, entries;
    in >> n >> sites >> entries;
    doc.component_ids.resize(n);
    doc.site_component.resize(sites);
    doc.site_multiplicity.resize(sites);
    for (auto& id : doc.component_ids) in >> std::quoted(id);
    for (std::size_t i = 0; i < sites; ++i) in >> doc.site_component[i] >> doc.site_multiplicity[i];
    in >> std::quoted(doc.base_choice) >> std::quoted(doc.association_choice) >> std::quoted(doc.electrolyte_choice)
       >> std::quoted(doc.ionic_dispersion_choice) >> std::quoted(doc.permittivity_choice) >> std::quoted(doc.polar_choice);
    for (std::size_t i = 0; i < entries; ++i) {
        epcsaft::ParameterEntry e;
        std::size_t width;
        in >> std::quoted(e.identity) >> std::quoted(e.field) >> width;
        e.index.resize(width);
        for (auto& index : e.index) in >> index;
        in >> e.magnitude >> std::quoted(e.unit);
        doc.entries.push_back(e);
    }
    const auto parameters = epcsaft::build_parameters(doc);
    const epcsaft::HelmholtzModel model(parameters);
    double temperature, solvent_mass;
    std::size_t reactions, rungs;
    in >> temperature >> solvent_mass >> reactions >> rungs;
    std::vector<double> ratios(n), molalities(rungs);
    for (auto& v : ratios) in >> v;
    for (auto& v : molalities) in >> v;
    std::cout << "reaction_index,source_pressure_Pa,analytic_contraction,molality_mol_kg,finite_contraction,terminal_pair,log_shift\n" << std::setprecision(17);
    for (std::size_t reaction = 0; reaction < reactions; ++reaction) {
        double pressure, analytic, previous = 0.0;
        in >> pressure >> analytic;
        std::vector<double> nu(n);
        for (auto& v : nu) in >> v;
        if (!in) throw std::runtime_error("Incomplete reference declarations");
        for (std::size_t rung = 0; rung < rungs; ++rung) {
            std::vector<double> x(n);
            for (std::size_t i = 0; i < n; ++i) x[i] = i == 0 ? 1 / solvent_mass : molalities[rung] * ratios[i];
            const auto total = std::accumulate(x.begin(), x.end(), 0.0);
            for (auto& v : x) v /= total;
            const auto state = model.state_at_pressure(temperature, pressure, x, epcsaft::DensityBranch::Liquid);
            long double contraction = 0.0L;
            for (std::size_t i = 0; i < n; ++i)
                contraction += nu[i] * (std::log(state.molar_density) + epcsaft::require_component_property(state.residual_chemical_potential_over_rt[i]));
            const double current = static_cast<double>(contraction);
            if (!std::isfinite(current)) throw std::runtime_error("Nonfinite finite-dilution contraction");
            std::cout << reaction << ',' << pressure << ',' << analytic << ',' << molalities[rung] << ',' << current << ',';
            if (rung) {
                const double ratio = std::sqrt(molalities[rung - 1] / molalities[rung]);
                const double limit = (ratio * current - previous) / (ratio - 1);
                std::cout << limit << ',' << limit - analytic;
            } else std::cout << ',';
            std::cout << '\n';
            previous = current;
        }
    }
}
