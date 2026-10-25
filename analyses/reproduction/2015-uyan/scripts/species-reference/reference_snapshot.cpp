#include "properties_support.hpp"
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
using namespace epcsaft;
void value(const std::optional<double>& x) {
    if (x) std::cout << std::hexfloat << *x; else std::cout << "unavailable";
    std::cout << ',';
}
void dump(const NeutralReferenceReactionDiagnostic& d) {
    std::cout << int(d.status) << ',' << d.message << ',';
    for (auto x : {d.terminal_contraction, d.observed_terminal_change,
        d.terminal_contraction_pressure_derivative, d.terminal_contraction_pressure_second_derivative,
        d.terminal_contraction_temperature_derivative, d.observed_temperature_derivative_change,
        d.terminal_contraction_temperature_second_derivative, d.observed_temperature_second_derivative_change,
        d.terminal_contraction_temperature_pressure_derivative}) value(x);
    for (const auto* v : {&d.terminal_contraction_parameter_derivative,
        &d.terminal_contraction_temperature_parameter_derivative, &d.observed_temperature_parameter_change})
        for (auto x : *v) value(x);
    for (const auto* columns : {&d.lower_derivatives, &d.derivatives}) for (const auto& c : *columns) {
        std::cout << c.temperature_order << ',' << c.pressure_order << ',';
        if (c.parameter) std::cout << *c.parameter;
        std::cout << ','; value(c.value);
        for (auto x : c.ladder) value(x);
        for (auto x : c.scaled_changes) value(x);
        value(c.criterion);
        if (c.action_status) std::cout << int(*c.action_status);
        std::cout << ',' << c.message << ',';
    }
    std::cout << '\n';
}
int main() {
    for (bool ionic : {false, true}) {
        auto p = ionic ? properties_support::aqueous_electrolyte() : properties_support::parameters(2);
        if (!ionic) {
            p.segment_count = {1.6069,1.6069}; p.sigma_constant = {3.5206,3.5206};
            p.dispersion_energy_over_k = {191.42,191.42}; p.molar_mass_kg_per_mol = {0.03007,0.03007};
        }
        HelmholtzModel model(p, {{ParameterFamily::DispersionEnergyOverK,{0}}});
        NeutralReferenceDeclaration decl{0, ionic ? std::vector<double>{0,1,1} : std::vector<double>{0,1},1};
        const std::vector<double> nu = ionic ? std::vector<double>{-1,1,1} : std::vector<double>{-2,1};
        for (auto tp : {std::pair{260.0,3e6},std::pair{280.0,5e6},std::pair{298.15,1e5}}) {
            if (ionic && tp.first != 298.15) continue;
            for (auto basis : {ReferenceSourceBasis::CommonMolalityInfiniteDilution,ReferenceSourceBasis::RawMoleFractionInfiniteDilution}) {
                for (bool fixed : {false,true}) {
                    ReactionReference ref{fixed ? std::optional<double>(tp.second*2) : std::nullopt,basis};
                    const auto base = reaction_references(model,tp.first,tp.second,decl,{nu},{ref},true,std::nullopt,2,true,4)[0];
                    dump(base);
#ifdef MIXED_REFERENCE
                    auto full = decl;
                    full.species_reference.assign(p.component_count(), SpeciesReference::InfiniteDilution);
                    full.species_reference[0] = SpeciesReference::PureLiquid;
                    const auto explicit_default = reaction_references(model,tp.first,tp.second,full,{nu},{ref},true,std::nullopt,2,true,4)[0];
                    std::ostringstream a,b;
                    auto* stream = std::cout.rdbuf(a.rdbuf()); dump(base);
                    std::cout.rdbuf(b.rdbuf()); dump(explicit_default); std::cout.rdbuf(stream);
                    if (a.str() != b.str()) throw std::runtime_error("explicit default differs");
#endif
                }
            }
        }
    }
}
