#pragma once
#include "properties_support.hpp"
#include "caloric_specification_inputs.hpp"

// Gross--Vrabec (2006), Table 2; Gross (2002), Table 1, 2B methanol.
namespace pure_saturation_checks {
inline epcsaft::Parameters acetone() {
    auto p = properties_support::parameters(1);
    p.segment_count = {2.7447}; p.sigma_constant = {3.2742}; p.dispersion_energy_over_k = {232.99};
    p.molar_mass_kg_per_mol = {0.05808}; p.dipole_moment_debye = {2.88}; p.dipole_count = {1};
    return p;
}
inline epcsaft::Parameters methanol() {
    auto p = properties_support::parameters(1);
    p.segment_count = {1.5255}; p.sigma_constant = {3.2300}; p.dispersion_energy_over_k = {188.90};
    p.molar_mass_kg_per_mol = {0.032042}; p.site_component = {0, 0}; p.site_multiplicity = {1, 1};
    p.association_active = {0, 1, 1, 0}; p.association_energy_over_k = {0, 2899.5, 2899.5, 0};
    p.association_volume = {0, 0.035176, 0.035176, 0}; return p;
}
inline epcsaft::equilibrium::EquilibriumProblem saturation(double t) {
    using namespace epcsaft::equilibrium;
    EquilibriumProblem p;
    PhaseSpec l; l.name = "L"; l.kind = PhaseKind::Liquid; l.support = {0};
    auto v = l; v.name = "V"; v.kind = PhaseKind::Vapor; v.amount = AmountSpec::pinned(0);
    p.phases = {l, v}; p.feed.values = {1}; p.temperature = FreeValue::fixed(t); p.pressure = FreeValue::variable();
    p.solver_options.values[epcsaft::policy::row] = 1e-12; return p;
}
}
