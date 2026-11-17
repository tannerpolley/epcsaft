// Published MEA Part I census state 19 (record be66f6cf), Amine-Thermodynamics #339, issue #21.
ParameterDocument parameters;
parameters.component_ids = {"carbon-dioxide", "monoethanolamine", "water", "protonated-monoethanolamine", "carbamate-anion", "bicarbonate-anion", "carbonate-anion", "hydronium-cation", "hydroxide-anion"};
parameters.site_ids = {"monoethanolamine/a", "monoethanolamine/b", "water/a", "water/b", "carbon-dioxide/a", "carbon-dioxide/b"};
parameters.site_component = {1, 1, 2, 2, 0, 0};
parameters.site_multiplicity = {1, 1, 1, 1, 1, 1};
parameters.base_choice = "pc-saft";
parameters.association_choice = "general-site";
parameters.electrolyte_choice = "born";
parameters.ionic_dispersion_choice = "exclude-same-sign-ion-pairs";
parameters.permittivity_choice = "solvent-only";
parameters.polar_choice = "none";
parameters.entries = {
    {"component/carbon-dioxide/molar_mass", "molar_mass", {0}, 0.04401, "kilogram / mole"},
    {"component/carbon-dioxide/charge_number", "charge_number", {0}, 0.0, "elementary-charge"},
    {"component/carbon-dioxide/dispersion_energy_over_k", "dispersion_energy_over_k", {0}, 173.44025, "kelvin"},
    {"component/carbon-dioxide/segment_count", "segment_count", {0}, 2.0729, "dimensionless"},
    {"component/carbon-dioxide/segment_diameter", "sigma_constant", {0}, 2.7852, "angstrom"},
    {"component/carbon-dioxide/solvation_factor", "solvation_factor", {0}, 1.0, "dimensionless"},
    {"component/monoethanolamine/molar_mass", "molar_mass", {1}, 0.06108, "kilogram / mole"},
    {"component/monoethanolamine/charge_number", "charge_number", {1}, 0.0, "elementary-charge"},
    {"component/monoethanolamine/dispersion_energy_over_k", "dispersion_energy_over_k", {1}, 277.174, "kelvin"},
    {"component/monoethanolamine/relative_permittivity", "relative_permittivity", {1}, 32.0, "dimensionless"},
    {"component/monoethanolamine/segment_count", "segment_count", {1}, 3.0353, "dimensionless"},
    {"component/monoethanolamine/segment_diameter", "sigma_constant", {1}, 3.0435, "angstrom"},
    {"component/monoethanolamine/solvation_factor", "solvation_factor", {1}, 1.0, "dimensionless"},
    {"component/water/molar_mass", "molar_mass", {2}, 0.01801528, "kilogram / mole"},
    {"component/water/charge_number", "charge_number", {2}, 0.0, "elementary-charge"},
    {"component/water/dispersion_energy_over_k", "dispersion_energy_over_k", {2}, 353.95, "kelvin"},
    {"component/water/segment_count", "segment_count", {2}, 1.2047, "dimensionless"},
    {"component/water/solvation_factor", "solvation_factor", {2}, 1.5, "dimensionless"},
    {"component/protonated-monoethanolamine/molar_mass", "molar_mass", {3}, 0.06209, "kilogram / mole"},
    {"component/protonated-monoethanolamine/charge_number", "charge_number", {3}, 1.0, "elementary-charge"},
    {"component/protonated-monoethanolamine/born_diameter", "born_diameter", {3}, 3.53322927146, "angstrom"},
    {"component/protonated-monoethanolamine/debye_huckel_diameter", "debye_huckel_diameter", {3}, 3.0668752979568, "angstrom"},
    {"component/protonated-monoethanolamine/dispersion_energy_over_k", "dispersion_energy_over_k", {3}, 232.687201645, "kelvin"},
    {"component/protonated-monoethanolamine/packing_diameter", "packing_diameter", {3}, 3.0668752979568, "angstrom"},
    {"component/protonated-monoethanolamine/segment_count", "segment_count", {3}, 1.0, "dimensionless"},
    {"component/protonated-monoethanolamine/segment_diameter", "sigma_constant", {3}, 3.48508556586, "angstrom"},
    {"component/protonated-monoethanolamine/solvation_factor", "solvation_factor", {3}, 1.0, "dimensionless"},
    {"component/carbamate-anion/molar_mass", "molar_mass", {4}, 0.10408, "kilogram / mole"},
    {"component/carbamate-anion/charge_number", "charge_number", {4}, -1.0, "elementary-charge"},
    {"component/carbamate-anion/born_diameter", "born_diameter", {4}, 3.54107030822, "angstrom"},
    {"component/carbamate-anion/debye_huckel_diameter", "debye_huckel_diameter", {4}, 3.1111830263448, "angstrom"},
    {"component/carbamate-anion/dispersion_energy_over_k", "dispersion_energy_over_k", {4}, 453.265244384, "kelvin"},
    {"component/carbamate-anion/packing_diameter", "packing_diameter", {4}, 3.1111830263448, "angstrom"},
    {"component/carbamate-anion/segment_count", "segment_count", {4}, 1.0, "dimensionless"},
    {"component/carbamate-anion/segment_diameter", "sigma_constant", {4}, 3.53543525721, "angstrom"},
    {"component/carbamate-anion/solvation_factor", "solvation_factor", {4}, 1.0, "dimensionless"},
    {"component/bicarbonate-anion/molar_mass", "molar_mass", {5}, 0.0610168, "kilogram / mole"},
    {"component/bicarbonate-anion/charge_number", "charge_number", {5}, -1.0, "elementary-charge"},
    {"component/bicarbonate-anion/born_diameter", "born_diameter", {5}, 3.0, "angstrom"},
    {"component/bicarbonate-anion/debye_huckel_diameter", "debye_huckel_diameter", {5}, 2.5780480000000003, "angstrom"},
    {"component/bicarbonate-anion/dispersion_energy_over_k", "dispersion_energy_over_k", {5}, 70.0, "kelvin"},
    {"component/bicarbonate-anion/packing_diameter", "packing_diameter", {5}, 2.5780480000000003, "angstrom"},
    {"component/bicarbonate-anion/segment_count", "segment_count", {5}, 1.0, "dimensionless"},
    {"component/bicarbonate-anion/segment_diameter", "sigma_constant", {5}, 2.9296, "angstrom"},
    {"component/bicarbonate-anion/solvation_factor", "solvation_factor", {5}, 1.0, "dimensionless"},
    {"component/carbonate-anion/molar_mass", "molar_mass", {6}, 0.06001, "kilogram / mole"},
    {"component/carbonate-anion/charge_number", "charge_number", {6}, -2.0, "elementary-charge"},
    {"component/carbonate-anion/born_diameter", "born_diameter", {6}, 3.0, "angstrom"},
    {"component/carbonate-anion/debye_huckel_diameter", "debye_huckel_diameter", {6}, 2.149136, "angstrom"},
    {"component/carbonate-anion/dispersion_energy_over_k", "dispersion_energy_over_k", {6}, 249.26, "kelvin"},
    {"component/carbonate-anion/packing_diameter", "packing_diameter", {6}, 2.149136, "angstrom"},
    {"component/carbonate-anion/segment_count", "segment_count", {6}, 1.0, "dimensionless"},
    {"component/carbonate-anion/segment_diameter", "sigma_constant", {6}, 2.4422, "angstrom"},
    {"component/carbonate-anion/solvation_factor", "solvation_factor", {6}, 1.0, "dimensionless"},
    {"component/hydronium-cation/molar_mass", "molar_mass", {7}, 0.01902, "kilogram / mole"},
    {"component/hydronium-cation/charge_number", "charge_number", {7}, 1.0, "elementary-charge"},
    {"component/hydronium-cation/born_diameter", "born_diameter", {7}, 1.218, "angstrom"},
    {"component/hydronium-cation/debye_huckel_diameter", "debye_huckel_diameter", {7}, 3.049552, "angstrom"},
    {"component/hydronium-cation/dispersion_energy_over_k", "dispersion_energy_over_k", {7}, 500.0, "kelvin"},
    {"component/hydronium-cation/packing_diameter", "packing_diameter", {7}, 3.049552, "angstrom"},
    {"component/hydronium-cation/segment_count", "segment_count", {7}, 1.0, "dimensionless"},
    {"component/hydronium-cation/segment_diameter", "sigma_constant", {7}, 3.4654, "angstrom"},
    {"component/hydronium-cation/solvation_factor", "solvation_factor", {7}, 1.0, "dimensionless"},
    {"component/hydroxide-anion/molar_mass", "molar_mass", {8}, 0.01701, "kilogram / mole"},
    {"component/hydroxide-anion/charge_number", "charge_number", {8}, -1.0, "elementary-charge"},
    {"component/hydroxide-anion/born_diameter", "born_diameter", {8}, 3.08107689400404, "angstrom"},
    {"component/hydroxide-anion/debye_huckel_diameter", "debye_huckel_diameter", {8}, 1.775576, "angstrom"},
    {"component/hydroxide-anion/dispersion_energy_over_k", "dispersion_energy_over_k", {8}, 650.0, "kelvin"},
    {"component/hydroxide-anion/packing_diameter", "packing_diameter", {8}, 1.775576, "angstrom"},
    {"component/hydroxide-anion/segment_count", "segment_count", {8}, 1.0, "dimensionless"},
    {"component/hydroxide-anion/segment_diameter", "sigma_constant", {8}, 2.0177, "angstrom"},
    {"component/hydroxide-anion/solvation_factor", "solvation_factor", {8}, 1.0, "dimensionless"},
    {"component/water/segment_diameter/constant-plus-sum-of-exponentials/constant", "sigma_constant", {2}, 2.7927, "angstrom"},
    {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/amplitude", "sigma_term_amplitude", {0, 2}, 10.11, "angstrom"},
    {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/exponent_coefficient", "sigma_term_exponent", {0, 2}, -0.01775, "1 / kelvin"},
    {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/amplitude", "sigma_term_amplitude", {1, 2}, -1.417, "angstrom"},
    {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/exponent_coefficient", "sigma_term_exponent", {1, 2}, -0.01146, "1 / kelvin"},
    {"component/water/relative_permittivity/constant-plus-polynomial-in-temperature/constant", "relative_permittivity", {2}, 78.41043478578715, "dimensionless"},
    {"component/water/relative_permittivity/constant-plus-polynomial-in-temperature/reference_temperature", "permittivity_reference_temperature", {2}, 298.15, "kelvin"},
    {"component/water/relative_permittivity/constant-plus-polynomial-in-temperature/term-0/coefficient", "permittivity_term_coefficient", {0, 2, 1}, -0.36133766233691506, "1 / kelvin"},
    {"component/water/relative_permittivity/constant-plus-polynomial-in-temperature/term-1/coefficient", "permittivity_term_coefficient", {1, 2, 2}, 0.00076555618295, "1 / kelvin ** 2"},
    {"pair/bicarbonate-anion/hydronium-cation/k_ij", "k_ij", {5, 7}, 0.0, "dimensionless"},
    {"pair/carbamate-anion/hydronium-cation/k_ij", "k_ij", {4, 7}, 0.0, "dimensionless"},
    {"pair/bicarbonate-anion/carbon-dioxide/k_ij", "k_ij", {5, 0}, 0.0, "dimensionless"},
    {"pair/carbamate-anion/carbon-dioxide/k_ij", "k_ij", {4, 0}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/carbonate-anion/k_ij", "k_ij", {0, 6}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/hydronium-cation/k_ij", "k_ij", {0, 7}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/hydroxide-anion/k_ij", "k_ij", {0, 8}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/monoethanolamine/k_ij", "k_ij", {0, 1}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/monoethanolamine/k_ij/reciprocal_temperature_slope", "k_ij_reciprocal_temperature_slope", {0, 1}, 0.0, "kelvin"},
    {"pair/carbon-dioxide/monoethanolamine/k_ij/reciprocal_temperature_slope/reference_temperature", "k_ij_reference_temperature", {0, 1}, 313.15, "kelvin"},
    {"pair/carbon-dioxide/protonated-monoethanolamine/k_ij", "k_ij", {0, 3}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/water/k_ij", "k_ij", {0, 2}, 0.013262879176919628, "dimensionless"},
    {"pair/carbonate-anion/hydronium-cation/k_ij", "k_ij", {6, 7}, 0.0, "dimensionless"},
    {"pair/hydronium-cation/hydroxide-anion/k_ij", "k_ij", {7, 8}, 0.0, "dimensionless"},
    {"pair/bicarbonate-anion/monoethanolamine/k_ij", "k_ij", {5, 1}, 0.0, "dimensionless"},
    {"pair/carbamate-anion/monoethanolamine/k_ij", "k_ij", {4, 1}, 0.0, "dimensionless"},
    {"pair/carbonate-anion/monoethanolamine/k_ij", "k_ij", {6, 1}, 0.0, "dimensionless"},
    {"pair/hydronium-cation/monoethanolamine/k_ij", "k_ij", {7, 1}, 0.0, "dimensionless"},
    {"pair/hydroxide-anion/monoethanolamine/k_ij", "k_ij", {8, 1}, 0.0, "dimensionless"},
    {"pair/monoethanolamine/protonated-monoethanolamine/k_ij", "k_ij", {1, 3}, 0.0, "dimensionless"},
    {"pair/monoethanolamine/water/k_ij", "k_ij", {1, 2}, -0.07352749874985018, "dimensionless"},
    {"pair/bicarbonate-anion/protonated-monoethanolamine/k_ij", "k_ij", {5, 3}, 0.0, "dimensionless"},
    {"pair/carbamate-anion/protonated-monoethanolamine/k_ij", "k_ij", {4, 3}, -0.00201813457644, "dimensionless"},
    {"pair/carbonate-anion/protonated-monoethanolamine/k_ij", "k_ij", {6, 3}, 0.0, "dimensionless"},
    {"pair/hydroxide-anion/protonated-monoethanolamine/k_ij", "k_ij", {8, 3}, 0.0, "dimensionless"},
    {"pair/bicarbonate-anion/water/k_ij", "k_ij", {5, 2}, 0.0, "dimensionless"},
    {"pair/carbamate-anion/water/k_ij", "k_ij", {4, 2}, 0.0, "dimensionless"},
    {"pair/carbonate-anion/water/k_ij", "k_ij", {6, 2}, -0.25, "dimensionless"},
    {"pair/hydronium-cation/water/k_ij", "k_ij", {7, 2}, 0.25, "dimensionless"},
    {"pair/hydroxide-anion/water/k_ij", "k_ij", {8, 2}, -0.25, "dimensionless"},
    {"pair/protonated-monoethanolamine/water/k_ij", "k_ij", {3, 2}, 0.0, "dimensionless"},
    {"model/electrolyte/c_shell", "born_c_shell", {}, 1.0, "dimensionless"},
    {"model/electrolyte/c_dielectric", "born_c_dielectric", {}, 1.0, "dimensionless"},
    {"model/relative_permittivity/ionic_region_relative_permittivity", "ionic_region_relative_permittivity", {}, 8.0, "dimensionless"},
    {"association/monoethanolamine/a/monoethanolamine/b/energy_over_k", "association_energy_over_k", {0, 1}, 2586.3, "kelvin"},
    {"association/monoethanolamine/a/monoethanolamine/b/volume", "association_volume", {0, 1}, 0.03747, "dimensionless"},
    {"association/water/a/water/b/energy_over_k", "association_energy_over_k", {2, 3}, 2425.7, "kelvin"},
    {"association/water/a/water/b/volume", "association_volume", {2, 3}, 0.04509, "dimensionless"},
    {"association/carbon-dioxide/a/water/b/energy_over_k", "association_energy_over_k", {4, 3}, 1212.85, "kelvin"},
    {"association/carbon-dioxide/a/water/b/volume", "association_volume", {4, 3}, 0.04509, "dimensionless"},
    {"association/carbon-dioxide/b/water/a/energy_over_k", "association_energy_over_k", {5, 2}, 1212.85, "kelvin"},
    {"association/carbon-dioxide/b/water/a/volume", "association_volume", {5, 2}, 0.04509, "dimensionless"},
    {"arithmetic-energy-geometric-volume", "association_combining_rule", {0, 3, 0, 1, 2, 3}, 0.0, "dimensionless"},
    {"arithmetic-energy-geometric-volume", "association_combining_rule", {1, 2, 0, 1, 2, 3}, 0.0, "dimensionless"},
    {"no_cross_association/carbon-dioxide/monoethanolamine", "no_cross_association", {0, 1}, 0.0, "dimensionless"},
};
const HelmholtzModel model(build_parameters(parameters));
EquilibriumProblem problem;
problem.temperature = FreeValue::fixed(313.15);
problem.pressure = FreeValue::fixed(5137.128272820414);
problem.feed.values = {0.48646410500000004, 9.7292821, 76.9751925, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
PhaseSpec liquid;
liquid.name = "L";
liquid.support = {0, 1, 2, 3, 4, 5, 6, 7, 8};
liquid.support_declared = true;
liquid.composition_guess = {};
liquid.packing_guess = std::nullopt;
PhaseSpec vapor;
vapor.kind = PhaseKind::Vapor;
vapor.name = "V";
vapor.support = {0, 1, 2};
vapor.support_declared = true;
vapor.composition_guess = {};
vapor.packing_guess = std::nullopt;
problem.phases = {liquid, vapor};
NeutralReferenceDeclaration reference;
reference.solvent_index = 2;
reference.solute_molality_ratios = {0.125, 0.125, 0.0, 0.625, 0.25, 0.125, 0.125, 0.125, 0.125};
reference.standard_molality_mol_per_kg = 1.0;
reference.species_reference = {};
problem.source_reference = reference;
{
    ReactionSpec reaction;
    reaction.name = "R1";
    reaction.stoichiometry = {0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0};
    ReactionLogPolynomial law;
    law.a = 3.7290977365566675;
    law.b = -13445.9;
    law.c = -22.4773;
    law.d = 0.0;
    law.reference_temperature = 313.15;
    law.temperature_min = 273.15;
    law.temperature_max = 498.15;
    law.standard_state_id = "aqueous-mole-fraction-infinite-dilution-water-v1";
    law.e = 0.0;
    law.f = 0.0;
    reaction.correlation = law;
    reaction.reference = ReactionReference{std::nullopt, ReferenceSourceBasis::RawMoleFractionInfiniteDilution};
    problem.reactions.push_back(std::move(reaction));
}
{
    ReactionSpec reaction;
    reaction.name = "R2";
    reaction.stoichiometry = {-1.0, 0.0, -2.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0};
    ReactionLogPolynomial law;
    law.a = 20.959245323180852;
    law.b = -11105.640030520277;
    law.c = -36.7816;
    law.d = 0.0;
    law.reference_temperature = 313.15;
    law.temperature_min = 293.15;
    law.temperature_max = 393.15;
    law.standard_state_id = "aqueous-molality-infinite-dilution-water-v1";
    law.e = 0.0;
    law.f = 0.0;
    reaction.correlation = law;
    reaction.reference = ReactionReference{std::nullopt, ReferenceSourceBasis::CommonMolalityInfiniteDilution};
    problem.reactions.push_back(std::move(reaction));
}
{
    ReactionSpec reaction;
    reaction.name = "R3";
    reaction.stoichiometry = {0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 1.0, 1.0, 0.0};
    ReactionLogPolynomial law;
    law.a = 12.145792981306897;
    law.b = -12431.7;
    law.c = -35.4819;
    law.d = 0.0;
    law.reference_temperature = 313.15;
    law.temperature_min = 273.15;
    law.temperature_max = 498.15;
    law.standard_state_id = "aqueous-mole-fraction-infinite-dilution-water-v1";
    law.e = 0.0;
    law.f = 0.0;
    reaction.correlation = law;
    reaction.reference = ReactionReference{std::nullopt, ReferenceSourceBasis::RawMoleFractionInfiniteDilution};
    problem.reactions.push_back(std::move(reaction));
}
{
    ReactionSpec reaction;
    reaction.name = "R4";
    reaction.stoichiometry = {0.0, 1.0, -1.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0};
    ReactionLogPolynomial law;
    law.a = 1.505015374192114;
    law.b = -1317.0489707842564;
    law.c = 0.0;
    law.d = 0.0;
    law.reference_temperature = 313.15;
    law.temperature_min = 293.15;
    law.temperature_max = 393.15;
    law.standard_state_id = "aqueous-molality-infinite-dilution-water-v1";
    law.e = 0.0;
    law.f = 0.0;
    reaction.correlation = law;
    reaction.reference = ReactionReference{100000.0, ReferenceSourceBasis::CommonMolalityInfiniteDilution};
    problem.reactions.push_back(std::move(reaction));
}
{
    ReactionSpec reaction;
    reaction.name = "R5";
    reaction.stoichiometry = {0.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
    ReactionLogPolynomial law;
    law.a = 2.3424545466714632;
    law.b = -6994.424474742253;
    law.c = 0.0;
    law.d = -0.0009848156442735533;
    law.reference_temperature = 313.15;
    law.temperature_min = 293.15;
    law.temperature_max = 393.15;
    law.standard_state_id = "aqueous-molality-infinite-dilution-water-v1";
    law.e = 0.0;
    law.f = 0.0;
    reaction.correlation = law;
    reaction.reference = ReactionReference{100000.0, ReferenceSourceBasis::CommonMolalityInfiniteDilution};
    problem.reactions.push_back(std::move(reaction));
}
problem.solver_options.set(policy::row, 1e-10);
