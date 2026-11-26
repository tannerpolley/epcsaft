// Amine-Thermodynamics 046e7462f, MDEA Endpoint A, 100 kPa reaction reference.
// Endpoint SHA-256 3c75b7c2475112d64d5845e04cda9a0a0f83b0ce6baa3ab32a36bcadc14ec950.
ParameterDocument document;
document.component_ids = {"carbon-dioxide", "water", "bicarbonate-anion", "carbonate-anion", "hydronium-cation", "hydroxide-anion", "methyldiethanolamine", "protonated-methyldiethanolamine"};
document.site_component = {6, 6, 7, 7, 1, 1, 0, 0};
document.site_multiplicity = {2, 2, 2, 2, 1, 1, 1, 1};
document.base_choice = "pc-saft";
document.association_choice = "general-site";
document.electrolyte_choice = "born";
document.ionic_dispersion_choice = "exclude-same-sign-ion-pairs";
document.permittivity_choice = "solvent-density-temperature";
document.polar_choice = "none";
document.entries.push_back({"component/carbon-dioxide/molar_mass", "molar_mass", {0}, 0.04401, "kilogram / mole"});
document.entries.push_back({"component/carbon-dioxide/charge_number", "charge_number", {0}, 0.0, "elementary-charge"});
document.entries.push_back({"component/carbon-dioxide/dispersion_energy_over_k", "dispersion_energy_over_k", {0}, 169.21, "kelvin"});
document.entries.push_back({"component/carbon-dioxide/segment_count", "segment_count", {0}, 2.0729, "dimensionless"});
document.entries.push_back({"component/carbon-dioxide/segment_diameter", "sigma_constant", {0}, 2.7852, "angstrom"});
document.entries.push_back({"component/carbon-dioxide/solvation_factor", "solvation_factor", {0}, 1.0, "dimensionless"});
document.entries.push_back({"component/carbon-dioxide/dielectric_volume", "dielectric_volume", {0}, 1e-08, "liter / mole"});
document.entries.push_back({"component/carbon-dioxide/dielectric_temperature", "dielectric_temperature", {0}, 1403.0, "kelvin"});
document.entries.push_back({"component/water/molar_mass", "molar_mass", {1}, 0.01801528, "kilogram / mole"});
document.entries.push_back({"component/water/charge_number", "charge_number", {1}, 0.0, "elementary-charge"});
document.entries.push_back({"component/water/dispersion_energy_over_k", "dispersion_energy_over_k", {1}, 432.3788130932931, "kelvin"});
document.entries.push_back({"component/water/segment_count", "segment_count", {1}, 1.252753837934407, "dimensionless"});
document.entries.push_back({"component/water/solvation_factor", "solvation_factor", {1}, 1.5, "dimensionless"});
document.entries.push_back({"component/water/dielectric_volume", "dielectric_volume", {1}, 0.3777, "liter / mole"});
document.entries.push_back({"component/water/dielectric_temperature", "dielectric_temperature", {1}, 1403.0, "kelvin"});
document.entries.push_back({"component/bicarbonate-anion/molar_mass", "molar_mass", {2}, 0.0610168, "kilogram / mole"});
document.entries.push_back({"component/bicarbonate-anion/charge_number", "charge_number", {2}, -1.0, "elementary-charge"});
document.entries.push_back({"component/bicarbonate-anion/born_diameter", "born_diameter", {2}, 3.0, "angstrom"});
document.entries.push_back({"component/bicarbonate-anion/debye_huckel_diameter", "debye_huckel_diameter", {2}, 2.9296, "angstrom"});
document.entries.push_back({"component/bicarbonate-anion/dispersion_energy_over_k", "dispersion_energy_over_k", {2}, 70.0, "kelvin"});
document.entries.push_back({"component/bicarbonate-anion/packing_diameter", "packing_diameter", {2}, 2.5780480000000003, "angstrom"});
document.entries.push_back({"component/bicarbonate-anion/segment_count", "segment_count", {2}, 1.0, "dimensionless"});
document.entries.push_back({"component/bicarbonate-anion/segment_diameter", "sigma_constant", {2}, 2.9296, "angstrom"});
document.entries.push_back({"component/bicarbonate-anion/solvation_factor", "solvation_factor", {2}, 1.0, "dimensionless"});
document.entries.push_back({"component/carbonate-anion/molar_mass", "molar_mass", {3}, 0.06001, "kilogram / mole"});
document.entries.push_back({"component/carbonate-anion/charge_number", "charge_number", {3}, -2.0, "elementary-charge"});
document.entries.push_back({"component/carbonate-anion/born_diameter", "born_diameter", {3}, 3.0, "angstrom"});
document.entries.push_back({"component/carbonate-anion/debye_huckel_diameter", "debye_huckel_diameter", {3}, 2.4422, "angstrom"});
document.entries.push_back({"component/carbonate-anion/dispersion_energy_over_k", "dispersion_energy_over_k", {3}, 249.26, "kelvin"});
document.entries.push_back({"component/carbonate-anion/packing_diameter", "packing_diameter", {3}, 2.149136, "angstrom"});
document.entries.push_back({"component/carbonate-anion/segment_count", "segment_count", {3}, 1.0, "dimensionless"});
document.entries.push_back({"component/carbonate-anion/segment_diameter", "sigma_constant", {3}, 2.4422, "angstrom"});
document.entries.push_back({"component/carbonate-anion/solvation_factor", "solvation_factor", {3}, 1.0, "dimensionless"});
document.entries.push_back({"component/hydronium-cation/molar_mass", "molar_mass", {4}, 0.01902, "kilogram / mole"});
document.entries.push_back({"component/hydronium-cation/charge_number", "charge_number", {4}, 1.0, "elementary-charge"});
document.entries.push_back({"component/hydronium-cation/born_diameter", "born_diameter", {4}, 1.218, "angstrom"});
document.entries.push_back({"component/hydronium-cation/debye_huckel_diameter", "debye_huckel_diameter", {4}, 3.4654, "angstrom"});
document.entries.push_back({"component/hydronium-cation/dispersion_energy_over_k", "dispersion_energy_over_k", {4}, 500.0, "kelvin"});
document.entries.push_back({"component/hydronium-cation/packing_diameter", "packing_diameter", {4}, 3.049552, "angstrom"});
document.entries.push_back({"component/hydronium-cation/segment_count", "segment_count", {4}, 1.0, "dimensionless"});
document.entries.push_back({"component/hydronium-cation/segment_diameter", "sigma_constant", {4}, 3.4654, "angstrom"});
document.entries.push_back({"component/hydronium-cation/solvation_factor", "solvation_factor", {4}, 1.0, "dimensionless"});
document.entries.push_back({"component/hydroxide-anion/molar_mass", "molar_mass", {5}, 0.01701, "kilogram / mole"});
document.entries.push_back({"component/hydroxide-anion/charge_number", "charge_number", {5}, -1.0, "elementary-charge"});
document.entries.push_back({"component/hydroxide-anion/born_diameter", "born_diameter", {5}, 3.08107689400404, "angstrom"});
document.entries.push_back({"component/hydroxide-anion/debye_huckel_diameter", "debye_huckel_diameter", {5}, 2.0177, "angstrom"});
document.entries.push_back({"component/hydroxide-anion/dispersion_energy_over_k", "dispersion_energy_over_k", {5}, 650.0, "kelvin"});
document.entries.push_back({"component/hydroxide-anion/packing_diameter", "packing_diameter", {5}, 1.775576, "angstrom"});
document.entries.push_back({"component/hydroxide-anion/segment_count", "segment_count", {5}, 1.0, "dimensionless"});
document.entries.push_back({"component/hydroxide-anion/segment_diameter", "sigma_constant", {5}, 2.0177, "angstrom"});
document.entries.push_back({"component/hydroxide-anion/solvation_factor", "solvation_factor", {5}, 1.0, "dimensionless"});
document.entries.push_back({"component/methyldiethanolamine/molar_mass", "molar_mass", {6}, 0.1191622, "kilogram / mole"});
document.entries.push_back({"component/methyldiethanolamine/charge_number", "charge_number", {6}, 0.0, "elementary-charge"});
document.entries.push_back({"component/methyldiethanolamine/segment_count", "segment_count", {6}, 3.675, "dimensionless"});
document.entries.push_back({"component/methyldiethanolamine/segment_diameter", "sigma_constant", {6}, 3.563, "angstrom"});
document.entries.push_back({"component/methyldiethanolamine/dispersion_energy_over_k", "dispersion_energy_over_k", {6}, 228.711, "kelvin"});
document.entries.push_back({"component/methyldiethanolamine/solvation_factor", "solvation_factor", {6}, 1.0, "dimensionless"});
document.entries.push_back({"component/methyldiethanolamine/dielectric_volume", "dielectric_volume", {6}, 0.8630669728921845, "liter / mole"});
document.entries.push_back({"component/methyldiethanolamine/dielectric_temperature", "dielectric_temperature", {6}, 1121.7734731422595, "kelvin"});
document.entries.push_back({"component/protonated-methyldiethanolamine/molar_mass", "molar_mass", {7}, 0.12017014, "kilogram / mole"});
document.entries.push_back({"component/protonated-methyldiethanolamine/charge_number", "charge_number", {7}, 1.0, "elementary-charge"});
document.entries.push_back({"component/protonated-methyldiethanolamine/charged_molecule_ion", "charged_molecule_ion", {7}, 1.0, "dimensionless"});
document.entries.push_back({"component/protonated-methyldiethanolamine/segment_count", "segment_count", {7}, 3.675, "dimensionless"});
document.entries.push_back({"component/protonated-methyldiethanolamine/segment_diameter", "sigma_constant", {7}, 3.5959113120691675, "angstrom"});
document.entries.push_back({"component/protonated-methyldiethanolamine/dispersion_energy_over_k", "dispersion_energy_over_k", {7}, 243.52151216647658, "kelvin"});
document.entries.push_back({"component/protonated-methyldiethanolamine/solvation_factor", "solvation_factor", {7}, 1.0, "dimensionless"});
document.entries.push_back({"component/protonated-methyldiethanolamine/debye_huckel_diameter", "debye_huckel_diameter", {7}, 3.5959113120691675, "angstrom"});
document.entries.push_back({"component/protonated-methyldiethanolamine/born_diameter", "born_diameter", {7}, 3.5959113120691675, "angstrom"});
document.entries.push_back({"component/water/segment_diameter/constant-plus-sum-of-exponentials/constant", "sigma_constant", {1}, 2.780432614942962, "angstrom"});
document.entries.push_back({"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/amplitude", "sigma_term_amplitude", {0, 1}, 10.11, "angstrom"});
document.entries.push_back({"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/exponent_coefficient", "sigma_term_exponent", {0, 1}, -0.01775, "1 / kelvin"});
document.entries.push_back({"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/amplitude", "sigma_term_amplitude", {1, 1}, -1.417, "angstrom"});
document.entries.push_back({"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/exponent_coefficient", "sigma_term_exponent", {1, 1}, -0.01146, "1 / kelvin"});
document.entries.push_back({"pair/bicarbonate-anion/hydronium-cation/k_ij", "k_ij", {2, 4}, 0.0, "dimensionless"});
document.entries.push_back({"pair/bicarbonate-anion/carbon-dioxide/k_ij", "k_ij", {2, 0}, 0.0, "dimensionless"});
document.entries.push_back({"pair/carbon-dioxide/carbonate-anion/k_ij", "k_ij", {0, 3}, 0.0, "dimensionless"});
document.entries.push_back({"pair/carbon-dioxide/hydronium-cation/k_ij", "k_ij", {0, 4}, 0.0, "dimensionless"});
document.entries.push_back({"pair/carbon-dioxide/hydroxide-anion/k_ij", "k_ij", {0, 5}, 0.0, "dimensionless"});
document.entries.push_back({"pair/carbon-dioxide/water/k_ij", "k_ij", {0, 1}, -0.046416901120577676, "dimensionless"});
document.entries.push_back({"pair/carbon-dioxide/water/k_ij/reciprocal_temperature_slope", "k_ij_reciprocal_temperature_slope", {0, 1}, -35.66505220645838, "kelvin"});
document.entries.push_back({"pair/carbon-dioxide/water/k_ij/reciprocal_temperature_slope/reference_temperature", "k_ij_reference_temperature", {0, 1}, 313.15, "kelvin"});
document.entries.push_back({"pair/carbonate-anion/hydronium-cation/k_ij", "k_ij", {3, 4}, 0.0, "dimensionless"});
document.entries.push_back({"pair/hydronium-cation/hydroxide-anion/k_ij", "k_ij", {4, 5}, 0.0, "dimensionless"});
document.entries.push_back({"pair/bicarbonate-anion/water/k_ij", "k_ij_temperature_law", {0, 2, 1}, 313.15, "kelvin"});
document.entries.push_back({"pair/bicarbonate-anion/water/k_ij/a0", "temperature_law_coefficient", {0, 0}, 0.1706931974036258, "dimensionless"});
document.entries.push_back({"pair/bicarbonate-anion/water/k_ij/a1", "temperature_law_coefficient", {0, 1}, 0.0039311876328716865, "1 / kelvin"});
document.entries.push_back({"pair/carbonate-anion/water/k_ij", "k_ij_temperature_law", {1, 3, 1}, 313.15, "kelvin"});
document.entries.push_back({"pair/carbonate-anion/water/k_ij/a0", "temperature_law_coefficient", {1, 0}, -0.9705531093245405, "dimensionless"});
document.entries.push_back({"pair/carbonate-anion/water/k_ij/a1", "temperature_law_coefficient", {1, 1}, 0.0019631260450306285, "1 / kelvin"});
document.entries.push_back({"pair/hydronium-cation/water/k_ij", "k_ij", {4, 1}, 0.25, "dimensionless"});
document.entries.push_back({"pair/hydroxide-anion/water/k_ij", "k_ij", {5, 1}, -0.25, "dimensionless"});
document.entries.push_back({"pair/carbon-dioxide/methyldiethanolamine/k_ij", "k_ij", {0, 6}, 0.0, "dimensionless"});
document.entries.push_back({"pair/carbon-dioxide/protonated-methyldiethanolamine/k_ij", "k_ij", {0, 7}, 0.3784408781206014, "dimensionless"});
document.entries.push_back({"pair/mdea-water/k_ij", "k_ij_temperature_law", {2, 1, 6}, 423.0, "kelvin"});
document.entries.push_back({"pair/mdea-water/k_ij/a0", "temperature_law_coefficient", {2, 0}, -0.06711345443725586, "dimensionless"});
document.entries.push_back({"pair/mdea-water/k_ij/a1", "temperature_law_coefficient", {2, 1}, 0.0011694693134037676, "1 / kelvin"});
document.entries.push_back({"pair/water/protonated-methyldiethanolamine/k_ij", "k_ij", {1, 7}, -0.2518751165455549, "dimensionless"});
document.entries.push_back({"pair/water/protonated-methyldiethanolamine/k_ij/reciprocal_temperature_slope", "k_ij_reciprocal_temperature_slope", {1, 7}, -31.02545057904153, "kelvin"});
document.entries.push_back({"pair/water/protonated-methyldiethanolamine/k_ij/reciprocal_temperature_slope/reference_temperature", "k_ij_reference_temperature", {1, 7}, 313.15, "kelvin"});
document.entries.push_back({"pair/bicarbonate-anion/methyldiethanolamine/k_ij", "k_ij", {2, 6}, 0.0, "dimensionless"});
document.entries.push_back({"pair/bicarbonate-anion/protonated-methyldiethanolamine/k_ij", "k_ij", {2, 7}, -0.03126710546865702, "dimensionless"});
document.entries.push_back({"pair/carbonate-anion/methyldiethanolamine/k_ij", "k_ij", {3, 6}, 0.0, "dimensionless"});
document.entries.push_back({"pair/carbonate-anion/protonated-methyldiethanolamine/k_ij", "k_ij", {3, 7}, 0.0, "dimensionless"});
document.entries.push_back({"pair/hydronium-cation/methyldiethanolamine/k_ij", "k_ij", {4, 6}, 0.0, "dimensionless"});
document.entries.push_back({"pair/hydroxide-anion/methyldiethanolamine/k_ij", "k_ij", {5, 6}, 0.0, "dimensionless"});
document.entries.push_back({"pair/hydroxide-anion/protonated-methyldiethanolamine/k_ij", "k_ij", {5, 7}, 0.0, "dimensionless"});
document.entries.push_back({"pair/methyldiethanolamine/protonated-methyldiethanolamine/k_ij", "k_ij", {6, 7}, 0.0, "dimensionless"});
document.entries.push_back({"model/electrolyte/c_shell", "born_c_shell", {}, 1.0, "dimensionless"});
document.entries.push_back({"model/electrolyte/c_dielectric", "born_c_dielectric", {}, 1.0, "dimensionless"});
document.entries.push_back({"model/relative_permittivity/ionic_region_relative_permittivity", "ionic_region_relative_permittivity", {}, 8.0, "dimensionless"});
document.entries.push_back({"association/methyldiethanolamine/energy_over_k", "association_energy_over_k", {0, 1}, 2046.624, "kelvin"});
document.entries.push_back({"association/methyldiethanolamine/volume", "association_volume", {0, 1}, 0.11800686, "dimensionless"});
document.entries.push_back({"wolbach-sandler", "association_combining_rule", {4, 1, 4, 5, 0, 1}, 0.0, "dimensionless"});
document.entries.push_back({"wolbach-sandler", "association_combining_rule", {0, 5, 4, 5, 0, 1}, 0.0, "dimensionless"});
document.entries.push_back({"association/protonated-methyldiethanolamine/energy_over_k", "association_energy_over_k", {2, 3}, 0.0, "kelvin"});
document.entries.push_back({"association/protonated-methyldiethanolamine/volume", "association_volume", {2, 3}, 0.11800686, "dimensionless"});
document.entries.push_back({"wolbach-sandler", "association_combining_rule", {2, 5, 2, 3, 4, 5}, 0.0, "dimensionless"});
document.entries.push_back({"wolbach-sandler", "association_combining_rule", {4, 3, 2, 3, 4, 5}, 0.0, "dimensionless"});
document.entries.push_back({"wolbach-sandler", "association_combining_rule", {2, 1, 2, 3, 0, 1}, 0.0, "dimensionless"});
document.entries.push_back({"wolbach-sandler", "association_combining_rule", {0, 3, 2, 3, 0, 1}, 0.0, "dimensionless"});
document.entries.push_back({"association/water/a/water/b/energy_over_k", "association_energy_over_k", {4, 5}, 1756.6497964508999, "kelvin"});
document.entries.push_back({"association/water/a/water/b/volume", "association_volume", {4, 5}, 0.0376197279502706, "dimensionless"});
document.entries.push_back({"association/carbon-dioxide/a/water/b/energy_over_k", "association_energy_over_k", {6, 5}, 878.3248982254499, "kelvin"});
document.entries.push_back({"association/carbon-dioxide/a/water/b/volume", "association_volume", {6, 5}, 0.0376197279502706, "dimensionless"});
document.entries.push_back({"association/carbon-dioxide/b/water/a/energy_over_k", "association_energy_over_k", {7, 4}, 878.3248982254499, "kelvin"});
document.entries.push_back({"association/carbon-dioxide/b/water/a/volume", "association_volume", {7, 4}, 0.0376197279502706, "dimensionless"});
document.entries.push_back({"no_cross_association/methyldiethanolamine/carbon-dioxide", "no_cross_association", {6, 0}, 0.0, "dimensionless"});
document.entries.push_back({"no_cross_association/protonated-methyldiethanolamine/carbon-dioxide", "no_cross_association", {7, 0}, 0.0, "dimensionless"});

#include "mdea-thermochemistry.inc"
const auto parameters = build_parameters(document);
const HelmholtzModel model(parameters, {}, thermal);
EquilibriumProblem cap;
cap.temperature = FreeValue::fixed(333.15);
cap.pressure = FreeValue::fixed(1121000.0);
cap.feed.values = {0.9112803667763014, 9.9217608607804, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
PhaseSpec liquid;
liquid.name = "L";
liquid.support = {0, 1, 2, 3, 4, 5, 6, 7};
liquid.support_declared = true;
cap.phases = {liquid};
cap.reactions = {
    {"R1", {0.0, -2.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, std::nullopt, std::nullopt,
        ReactionLogPolynomial{3.7290977365566675, -13445.9, -22.4773, 0.0, 313.15, 288.71, 422.04, "aqueous-mole-fraction-infinite-dilution-water-v1"},
        ReactionReference{100000.0, ReferenceSourceBasis::RawMoleFractionInfiniteDilution}},
    {"R2", {-1.0, -2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0}, std::nullopt, std::nullopt,
        ReactionLogPolynomial{24.109329984336796, -12092.1, -36.7816, 0.0, 313.15, 288.71, 422.04, "aqueous-molality-infinite-dilution-water-v1"},
        ReactionReference{100000.0, ReferenceSourceBasis::CommonMolalityInfiniteDilution}},
    {"R3", {0.0, -1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0}, std::nullopt, std::nullopt,
        ReactionLogPolynomial{12.145792981306897, -12431.7, -35.4819, 0.0, 313.15, 288.71, 422.04, "aqueous-mole-fraction-infinite-dilution-water-v1"},
        ReactionReference{100000.0, ReferenceSourceBasis::RawMoleFractionInfiniteDilution}},
    {"MDEA-dissociation", {0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 1.0, -1.0}, std::nullopt, std::nullopt,
        ReactionLogPolynomial{-16.398914976729678, -819.7, 10.976, 0.0, 313.15, 288.71, 422.04, "aqueous-molality-infinite-dilution-water-v1"},
        ReactionReference{100000.0, ReferenceSourceBasis::CommonMolalityInfiniteDilution}},
};
cap.source_reference = NeutralReferenceDeclaration{1, {0.125, 0.0, 0.125, 0.03125, 0.125, 0.0625, 0.125, 0.125},
    1.0, {}};
cap.solver_options.set(policy::row, policy::EQUILIBRIUM_ROW_RESIDUAL.value);
