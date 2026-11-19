// Amine-Thermodynamics f69a7d8, reactive-A12s-stage1-endpoint.json and Hilliard D.5.
// The eight reaction declarations match the retained row 55 bubble replay.
ParameterDocument document;
document.component_ids = {"carbon-dioxide", "water", "bicarbonate-anion", "carbonate-anion", "hydronium-cation", "hydroxide-anion", "piperazine", "piperazinium", "piperazinediium", "piperazine-carbamate", "piperazine-dicarbamate", "protonated-piperazine-carbamate"};
document.site_component = {1, 1, 0, 0, 6, 6};
document.site_multiplicity = {1, 1, 1, 1, 1, 1};
document.base_choice = "pc-saft";
document.association_choice = "general-site";
document.electrolyte_choice = "born";
document.permittivity_choice = "solvent-density-temperature";
document.ionic_dispersion_choice = "exclude-same-sign-ion-pairs";
document.polar_choice = "none";
document.entries = {
    {"component/carbon-dioxide/charge_number", "charge_number", {0}, 0, "elementary-charge"},
    {"component/carbon-dioxide/molar_mass", "molar_mass", {0}, 0.04401, "kilogram / mole"},
    {"component/carbon-dioxide/dispersion_energy_over_k", "dispersion_energy_over_k", {0}, 169.21, "kelvin"},
    {"component/carbon-dioxide/segment_count", "segment_count", {0}, 2.0729, "dimensionless"},
    {"component/carbon-dioxide/segment_diameter", "sigma_constant", {0}, 2.7852, "angstrom"},
    {"component/carbon-dioxide/solvation_factor", "solvation_factor", {0}, 1.0, "dimensionless"},
    {"component/carbon-dioxide/dielectric_volume", "dielectric_volume", {0}, 1e-08, "liter / mole"},
    {"component/carbon-dioxide/dielectric_temperature", "dielectric_temperature", {0}, 1403.0, "kelvin"},
    {"component/water/charge_number", "charge_number", {1}, 0, "elementary-charge"},
    {"component/water/molar_mass", "molar_mass", {1}, 0.01801528, "kilogram / mole"},
    {"component/water/dispersion_energy_over_k", "dispersion_energy_over_k", {1}, 432.3788130932931, "kelvin"},
    {"component/water/segment_count", "segment_count", {1}, 1.252753837934407, "dimensionless"},
    {"component/water/solvation_factor", "solvation_factor", {1}, 1.5, "dimensionless"},
    {"component/water/dielectric_volume", "dielectric_volume", {1}, 0.3777, "liter / mole"},
    {"component/water/dielectric_temperature", "dielectric_temperature", {1}, 1403.0, "kelvin"},
    {"component/bicarbonate-anion/charge_number", "charge_number", {2}, -1, "elementary-charge"},
    {"component/bicarbonate-anion/molar_mass", "molar_mass", {2}, 0.0610168, "kilogram / mole"},
    {"component/bicarbonate-anion/born_diameter", "born_diameter", {2}, 3.0, "angstrom"},
    {"component/bicarbonate-anion/debye_huckel_diameter", "debye_huckel_diameter", {2}, 2.9296, "angstrom"},
    {"component/bicarbonate-anion/dispersion_energy_over_k", "dispersion_energy_over_k", {2}, 70.0, "kelvin"},
    {"component/bicarbonate-anion/packing_diameter", "packing_diameter", {2}, 2.5780480000000003, "angstrom"},
    {"component/bicarbonate-anion/segment_count", "segment_count", {2}, 1.0, "dimensionless"},
    {"component/bicarbonate-anion/segment_diameter", "sigma_constant", {2}, 2.9296, "angstrom"},
    {"component/bicarbonate-anion/solvation_factor", "solvation_factor", {2}, 1.0, "dimensionless"},
    {"component/carbonate-anion/charge_number", "charge_number", {3}, -2, "elementary-charge"},
    {"component/carbonate-anion/molar_mass", "molar_mass", {3}, 0.06001, "kilogram / mole"},
    {"component/carbonate-anion/born_diameter", "born_diameter", {3}, 3.0, "angstrom"},
    {"component/carbonate-anion/debye_huckel_diameter", "debye_huckel_diameter", {3}, 2.4422, "angstrom"},
    {"component/carbonate-anion/dispersion_energy_over_k", "dispersion_energy_over_k", {3}, 249.26, "kelvin"},
    {"component/carbonate-anion/packing_diameter", "packing_diameter", {3}, 2.149136, "angstrom"},
    {"component/carbonate-anion/segment_count", "segment_count", {3}, 1.0, "dimensionless"},
    {"component/carbonate-anion/segment_diameter", "sigma_constant", {3}, 2.4422, "angstrom"},
    {"component/carbonate-anion/solvation_factor", "solvation_factor", {3}, 1.0, "dimensionless"},
    {"component/hydronium-cation/charge_number", "charge_number", {4}, 1, "elementary-charge"},
    {"component/hydronium-cation/molar_mass", "molar_mass", {4}, 0.01902, "kilogram / mole"},
    {"component/hydronium-cation/born_diameter", "born_diameter", {4}, 1.218, "angstrom"},
    {"component/hydronium-cation/debye_huckel_diameter", "debye_huckel_diameter", {4}, 3.4654, "angstrom"},
    {"component/hydronium-cation/dispersion_energy_over_k", "dispersion_energy_over_k", {4}, 500.0, "kelvin"},
    {"component/hydronium-cation/packing_diameter", "packing_diameter", {4}, 3.049552, "angstrom"},
    {"component/hydronium-cation/segment_count", "segment_count", {4}, 1.0, "dimensionless"},
    {"component/hydronium-cation/segment_diameter", "sigma_constant", {4}, 3.4654, "angstrom"},
    {"component/hydronium-cation/solvation_factor", "solvation_factor", {4}, 1.0, "dimensionless"},
    {"component/hydroxide-anion/charge_number", "charge_number", {5}, -1, "elementary-charge"},
    {"component/hydroxide-anion/molar_mass", "molar_mass", {5}, 0.01701, "kilogram / mole"},
    {"component/hydroxide-anion/born_diameter", "born_diameter", {5}, 3.08107689400404, "angstrom"},
    {"component/hydroxide-anion/debye_huckel_diameter", "debye_huckel_diameter", {5}, 2.0177, "angstrom"},
    {"component/hydroxide-anion/dispersion_energy_over_k", "dispersion_energy_over_k", {5}, 650.0, "kelvin"},
    {"component/hydroxide-anion/packing_diameter", "packing_diameter", {5}, 1.775576, "angstrom"},
    {"component/hydroxide-anion/segment_count", "segment_count", {5}, 1.0, "dimensionless"},
    {"component/hydroxide-anion/segment_diameter", "sigma_constant", {5}, 2.0177, "angstrom"},
    {"component/hydroxide-anion/solvation_factor", "solvation_factor", {5}, 1.0, "dimensionless"},
    {"component/piperazine/charge_number", "charge_number", {6}, 0, "elementary-charge"},
    {"component/piperazine/molar_mass", "molar_mass", {6}, 0.08613, "kilogram / mole"},
    {"component/piperazine/segment_count", "segment_count", {6}, 2.6519, "dimensionless"},
    {"component/piperazine/segment_diameter", "sigma_constant", {6}, 3.7136, "angstrom"},
    {"component/piperazine/dispersion_energy_over_k", "dispersion_energy_over_k", {6}, 313.20597765826267, "kelvin"},
    {"component/piperazine/solvation_factor", "solvation_factor", {6}, 1, "dimensionless"},
    {"component/piperazine/dielectric_volume", "dielectric_volume", {6}, 0.7309865858310542, "liter / mole"},
    {"component/piperazine/dielectric_temperature", "dielectric_temperature", {6}, 1019.0802597119892, "kelvin"},
    {"component/piperazinium/charge_number", "charge_number", {7}, 1, "elementary-charge"},
    {"component/piperazinium/molar_mass", "molar_mass", {7}, 0.087146, "kilogram / mole"},
    {"component/piperazinium/segment_count", "segment_count", {7}, 1, "dimensionless"},
    {"component/piperazinium/segment_diameter", "sigma_constant", {7}, 4.290267650574509, "angstrom"},
    {"component/piperazinium/dispersion_energy_over_k", "dispersion_energy_over_k", {7}, 100.71654266353482, "kelvin"},
    {"component/piperazinium/born_diameter", "born_diameter", {7}, 5, "angstrom"},
    {"component/piperazinium/packing_diameter", "packing_diameter", {7}, 3.775435532505568, "angstrom"},
    {"component/piperazinium/debye_huckel_diameter", "debye_huckel_diameter", {7}, 4.290267650574509, "angstrom"},
    {"component/piperazinium/solvation_factor", "solvation_factor", {7}, 1, "dimensionless"},
    {"component/piperazinediium/charge_number", "charge_number", {8}, 2, "elementary-charge"},
    {"component/piperazinediium/molar_mass", "molar_mass", {8}, 0.088154, "kilogram / mole"},
    {"component/piperazinediium/segment_count", "segment_count", {8}, 1, "dimensionless"},
    {"component/piperazinediium/segment_diameter", "sigma_constant", {8}, 4.290267650574509, "angstrom"},
    {"component/piperazinediium/dispersion_energy_over_k", "dispersion_energy_over_k", {8}, 100.71654266353482, "kelvin"},
    {"component/piperazinediium/born_diameter", "born_diameter", {8}, 5, "angstrom"},
    {"component/piperazinediium/packing_diameter", "packing_diameter", {8}, 3.775435532505568, "angstrom"},
    {"component/piperazinediium/debye_huckel_diameter", "debye_huckel_diameter", {8}, 4.290267650574509, "angstrom"},
    {"component/piperazinediium/solvation_factor", "solvation_factor", {8}, 1, "dimensionless"},
    {"component/piperazine-carbamate/charge_number", "charge_number", {9}, -1, "elementary-charge"},
    {"component/piperazine-carbamate/molar_mass", "molar_mass", {9}, 0.129132, "kilogram / mole"},
    {"component/piperazine-carbamate/segment_count", "segment_count", {9}, 1, "dimensionless"},
    {"component/piperazine-carbamate/segment_diameter", "sigma_constant", {9}, 4.907681021003549, "angstrom"},
    {"component/piperazine-carbamate/dispersion_energy_over_k", "dispersion_energy_over_k", {9}, 691.9841751513, "kelvin"},
    {"component/piperazine-carbamate/born_diameter", "born_diameter", {9}, 5, "angstrom"},
    {"component/piperazine-carbamate/packing_diameter", "packing_diameter", {9}, 4.318759298483124, "angstrom"},
    {"component/piperazine-carbamate/debye_huckel_diameter", "debye_huckel_diameter", {9}, 4.907681021003549, "angstrom"},
    {"component/piperazine-carbamate/solvation_factor", "solvation_factor", {9}, 1, "dimensionless"},
    {"component/piperazine-dicarbamate/charge_number", "charge_number", {10}, -2, "elementary-charge"},
    {"component/piperazine-dicarbamate/molar_mass", "molar_mass", {10}, 0.172134, "kilogram / mole"},
    {"component/piperazine-dicarbamate/segment_count", "segment_count", {10}, 1, "dimensionless"},
    {"component/piperazine-dicarbamate/segment_diameter", "sigma_constant", {10}, 4.907681021003549, "angstrom"},
    {"component/piperazine-dicarbamate/dispersion_energy_over_k", "dispersion_energy_over_k", {10}, 691.9841751513, "kelvin"},
    {"component/piperazine-dicarbamate/born_diameter", "born_diameter", {10}, 5, "angstrom"},
    {"component/piperazine-dicarbamate/packing_diameter", "packing_diameter", {10}, 4.318759298483124, "angstrom"},
    {"component/piperazine-dicarbamate/debye_huckel_diameter", "debye_huckel_diameter", {10}, 4.907681021003549, "angstrom"},
    {"component/piperazine-dicarbamate/solvation_factor", "solvation_factor", {10}, 1, "dimensionless"},
    {"component/protonated-piperazine-carbamate/charge_number", "charge_number", {11}, 0, "elementary-charge"},
    {"component/protonated-piperazine-carbamate/molar_mass", "molar_mass", {11}, 0.13014, "kilogram / mole"},
    {"component/protonated-piperazine-carbamate/segment_count", "segment_count", {11}, 2.6519, "dimensionless"},
    {"component/protonated-piperazine-carbamate/segment_diameter", "sigma_constant", {11}, 3.838681667122176, "angstrom"},
    {"component/protonated-piperazine-carbamate/dispersion_energy_over_k", "dispersion_energy_over_k", {11}, 948.1409230864688, "kelvin"},
    {"component/protonated-piperazine-carbamate/solvation_factor", "solvation_factor", {11}, 1, "dimensionless"},
    {"component/protonated-piperazine-carbamate/dielectric_volume", "dielectric_volume", {11}, 0.7309865858310542, "liter / mole"},
    {"component/protonated-piperazine-carbamate/dielectric_temperature", "dielectric_temperature", {11}, 1019.0802597119892, "kelvin"},
    {"component/water/segment_diameter/constant-plus-sum-of-exponentials/constant", "sigma_constant", {1}, 2.780432614942962, "angstrom"},
    {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/amplitude", "sigma_term_amplitude", {0, 1}, 10.11, "angstrom"},
    {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/exponent_coefficient", "sigma_term_exponent", {0, 1}, -0.01775, "1 / kelvin"},
    {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/amplitude", "sigma_term_amplitude", {1, 1}, -1.417, "angstrom"},
    {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/exponent_coefficient", "sigma_term_exponent", {1, 1}, -0.01146, "1 / kelvin"},
    {"pair/bicarbonate-anion/hydronium-cation/k_ij", "k_ij", {2, 4}, 0.0, "dimensionless"},
    {"pair/bicarbonate-anion/carbon-dioxide/k_ij", "k_ij", {2, 0}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/carbonate-anion/k_ij", "k_ij", {0, 3}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/hydronium-cation/k_ij", "k_ij", {0, 4}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/hydroxide-anion/k_ij", "k_ij", {0, 5}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/water/k_ij", "k_ij", {0, 1}, -0.046416901120577676, "dimensionless"},
    {"pair/carbon-dioxide/water/k_ij/reciprocal_temperature_slope", "k_ij_reciprocal_temperature_slope", {0, 1}, -35.66505220645838, "kelvin"},
    {"pair/carbon-dioxide/water/k_ij/reciprocal_temperature_slope/reference_temperature", "k_ij_reference_temperature", {0, 1}, 313.15, "kelvin"},
    {"pair/carbonate-anion/hydronium-cation/k_ij", "k_ij", {3, 4}, 0.0, "dimensionless"},
    {"pair/hydronium-cation/hydroxide-anion/k_ij", "k_ij", {4, 5}, 0.0, "dimensionless"},
    {"pair/bicarbonate-anion/water/k_ij", "k_ij", {2, 1}, 1.0, "dimensionless"},
    {"pair/carbonate-anion/water/k_ij", "k_ij", {3, 1}, -0.25, "dimensionless"},
    {"pair/hydronium-cation/water/k_ij", "k_ij", {4, 1}, 0.25, "dimensionless"},
    {"pair/hydroxide-anion/water/k_ij", "k_ij", {5, 1}, -0.25, "dimensionless"},
    {"pair/pz-water/k", "k_ij_temperature_law", {0, 6, 1}, 353.15, "kelvin"},
    {"pair/pz-water/k0", "temperature_law_coefficient", {0, 0}, -0.19764540162608946, "dimensionless"},
    {"pair/pz-water/k1", "temperature_law_coefficient", {0, 1}, 0.0011250147848721607, "1 / kelvin"},
    {"pair/carbon-dioxide/piperazine/k_ij", "k_ij", {0, 6}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/piperazinium/k_ij", "k_ij", {0, 7}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/piperazinediium/k_ij", "k_ij", {0, 8}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/piperazine-carbamate/k_ij", "k_ij", {0, 9}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/piperazine-dicarbamate/k_ij", "k_ij", {0, 10}, 0.0, "dimensionless"},
    {"pair/carbon-dioxide/protonated-piperazine-carbamate/k_ij", "k_ij", {0, 11}, 0.0, "dimensionless"},
    {"pair/water/piperazinium/k_ij", "k_ij", {1, 7}, -1.0, "dimensionless"},
    {"pair/water/piperazinediium/k_ij", "k_ij", {1, 8}, -1.0, "dimensionless"},
    {"pair/water/piperazine-carbamate/k_ij", "k_ij", {1, 9}, 0.06561976300270762, "dimensionless"},
    {"pair/water/piperazine-dicarbamate/k_ij", "k_ij", {1, 10}, 0.06561976300270762, "dimensionless"},
    {"pair/water/protonated-piperazine-carbamate/k_ij", "k_ij", {1, 11}, 0.0, "dimensionless"},
    {"pair/bicarbonate-anion/piperazine/k_ij", "k_ij", {2, 6}, 0.0, "dimensionless"},
    {"pair/bicarbonate-anion/piperazinium/k_ij", "k_ij", {2, 7}, 0.33836633155927126, "dimensionless"},
    {"pair/bicarbonate-anion/piperazinediium/k_ij", "k_ij", {2, 8}, 0.33836633155927126, "dimensionless"},
    {"pair/bicarbonate-anion/protonated-piperazine-carbamate/k_ij", "k_ij", {2, 11}, 0.0, "dimensionless"},
    {"pair/carbonate-anion/piperazine/k_ij", "k_ij", {3, 6}, 0.0, "dimensionless"},
    {"pair/carbonate-anion/piperazinium/k_ij", "k_ij", {3, 7}, 0.0, "dimensionless"},
    {"pair/carbonate-anion/piperazinediium/k_ij", "k_ij", {3, 8}, 0.0, "dimensionless"},
    {"pair/carbonate-anion/protonated-piperazine-carbamate/k_ij", "k_ij", {3, 11}, 0.0, "dimensionless"},
    {"pair/hydronium-cation/piperazine/k_ij", "k_ij", {4, 6}, 0.0, "dimensionless"},
    {"pair/hydronium-cation/piperazine-carbamate/k_ij", "k_ij", {4, 9}, 0.0, "dimensionless"},
    {"pair/hydronium-cation/piperazine-dicarbamate/k_ij", "k_ij", {4, 10}, 0.0, "dimensionless"},
    {"pair/hydronium-cation/protonated-piperazine-carbamate/k_ij", "k_ij", {4, 11}, 0.0, "dimensionless"},
    {"pair/hydroxide-anion/piperazine/k_ij", "k_ij", {5, 6}, 0.0, "dimensionless"},
    {"pair/hydroxide-anion/piperazinium/k_ij", "k_ij", {5, 7}, 0.0, "dimensionless"},
    {"pair/hydroxide-anion/piperazinediium/k_ij", "k_ij", {5, 8}, 0.0, "dimensionless"},
    {"pair/hydroxide-anion/protonated-piperazine-carbamate/k_ij", "k_ij", {5, 11}, 0.0, "dimensionless"},
    {"pair/piperazine/piperazinium/k_ij", "k_ij", {6, 7}, 0.0, "dimensionless"},
    {"pair/piperazine/piperazinediium/k_ij", "k_ij", {6, 8}, 0.0, "dimensionless"},
    {"pair/piperazine/piperazine-carbamate/k_ij", "k_ij", {6, 9}, 0.0, "dimensionless"},
    {"pair/piperazine/piperazine-dicarbamate/k_ij", "k_ij", {6, 10}, 0.0, "dimensionless"},
    {"pair/piperazine/protonated-piperazine-carbamate/k_ij", "k_ij", {6, 11}, 0.0, "dimensionless"},
    {"pair/piperazinium/piperazine-carbamate/k_ij", "k_ij", {7, 9}, 0.0, "dimensionless"},
    {"pair/piperazinium/piperazine-dicarbamate/k_ij", "k_ij", {7, 10}, 0.0, "dimensionless"},
    {"pair/piperazinium/protonated-piperazine-carbamate/k_ij", "k_ij", {7, 11}, 0.0, "dimensionless"},
    {"pair/piperazinediium/piperazine-carbamate/k_ij", "k_ij", {8, 9}, 0.0, "dimensionless"},
    {"pair/piperazinediium/piperazine-dicarbamate/k_ij", "k_ij", {8, 10}, 0.0, "dimensionless"},
    {"pair/piperazinediium/protonated-piperazine-carbamate/k_ij", "k_ij", {8, 11}, 0.0, "dimensionless"},
    {"pair/piperazine-carbamate/protonated-piperazine-carbamate/k_ij", "k_ij", {9, 11}, 0.0, "dimensionless"},
    {"pair/piperazine-dicarbamate/protonated-piperazine-carbamate/k_ij", "k_ij", {10, 11}, 0.0, "dimensionless"},
    {"model/electrolyte/c_shell", "born_c_shell", {}, 1.0, "dimensionless"},
    {"model/electrolyte/c_dielectric", "born_c_dielectric", {}, 1.0, "dimensionless"},
    {"model/relative_permittivity/ionic_region_relative_permittivity", "ionic_region_relative_permittivity", {}, 8.0, "dimensionless"},
    {"association/water/a/water/b/energy_over_k", "association_energy_over_k", {0, 1}, 1756.6497964508999, "kelvin"},
    {"association/water/a/water/b/volume", "association_volume", {0, 1}, 0.0376197279502706, "dimensionless"},
    {"association/carbon-dioxide/a/water/b/energy_over_k", "association_energy_over_k", {2, 1}, 878.3248982254499, "kelvin"},
    {"association/carbon-dioxide/a/water/b/volume", "association_volume", {2, 1}, 0.0376197279502706, "dimensionless"},
    {"association/carbon-dioxide/b/water/a/energy_over_k", "association_energy_over_k", {3, 0}, 878.3248982254499, "kelvin"},
    {"association/carbon-dioxide/b/water/a/volume", "association_volume", {3, 0}, 0.0376197279502706, "dimensionless"},
    {"association/piperazine/a/b/energy", "association_energy_over_k", {4, 5}, 1469.820092254, "kelvin"},
    {"association/piperazine/a/b/volume", "association_volume", {4, 5}, 0.027229855340530382, "dimensionless"},
    {"wolbach-sandler", "association_combining_rule", {4, 1, 4, 5, 0, 1}, 0.0, "dimensionless"},
    {"wolbach-sandler", "association_combining_rule", {5, 0, 4, 5, 0, 1}, 0.0, "dimensionless"},
    {"no_cross_association/carbon-dioxide/piperazine", "no_cross_association", {0, 6}, 0.0, "dimensionless"},
};
const HelmholtzModel model(build_parameters(document));
EquilibriumProblem problem;
problem.topology_declared = true;
problem.solver_options.values[policy::row] = 1e-10;
problem.source_reference = NeutralReferenceDeclaration{1,
    {1.0, 0, 1.0, 1.0, 4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, 1.0, {}};
PhaseSpec liquid;
liquid.name = "L";
liquid.support = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
liquid.support_declared = true;
PhaseSpec vapor;
vapor.name = "V";
vapor.support = {1, 0, 6, 11};
vapor.support_declared = true;
vapor.kind = PhaseKind::Vapor;
vapor.amount.kind = AmountKind::Pinned;
vapor.amount.fraction = 0.0;
vapor.composition_guess = {0.4999, 0.5, 0, 0, 0, 0, 5e-05, 0, 0, 0, 0, 5e-05};
problem.phases = {liquid, vapor};
ReactionSpec r1;
r1.name = "R1";
r1.stoichiometry = {0, -2.0, 0, 0, 1.0, 1.0, 0, 0, 0, 0, 0, 0};
r1.correlation = ReactionLogPolynomial{3.7290977365566675, -13445.9, -22.4773, 0.0,
    313.15, 283.15, 498.15, "aqueous-mole-fraction-infinite-dilution-water-v1"};
r1.reference = ReactionReference{std::nullopt, ReferenceSourceBasis::RawMoleFractionInfiniteDilution};
ReactionSpec r2;
r2.name = "R2";
r2.stoichiometry = {-1.0, -2.0, 1.0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0};
r2.correlation = ReactionLogPolynomial{23.87092013089665, -12017.430996556475, -36.7816, 0.0,
    313.15, 283.15, 393.15, "aqueous-molality-infinite-dilution-water-v1"};
r2.reference = ReactionReference{std::nullopt, ReferenceSourceBasis::CommonMolalityInfiniteDilution};
ReactionSpec r3;
r3.name = "R3";
r3.stoichiometry = {0, -1.0, -1.0, 1.0, 1.0, 0, 0, 0, 0, 0, 0, 0};
r3.correlation = ReactionLogPolynomial{12.145792981306897, -12431.7, -35.4819, 0.0,
    313.15, 283.15, 498.15, "aqueous-mole-fraction-infinite-dilution-water-v1"};
r3.reference = ReactionReference{std::nullopt, ReferenceSourceBasis::RawMoleFractionInfiniteDilution};
ReactionSpec p1;
p1.name = "P1";
p1.stoichiometry = {0, -1, 0, 0, 1, 0, 1, -1, 0, 0, 0, 0};
p1.correlation = ReactionLogPolynomial{-11.559187724019623, -3123.5331542628023, 6.786961056344094, 0.0,
    313.15, 283.15, 353.0, "aqueous-molality-infinite-dilution-water-v1"};
p1.reference = ReactionReference{100000.0, ReferenceSourceBasis::CommonMolalityInfiniteDilution};
ReactionSpec p2;
p2.name = "P2";
p2.stoichiometry = {0, -1, 0, 0, 1, 0, 0, 1, -1, 0, 0, 0};
p2.correlation = ReactionLogPolynomial{-2.804785141571699, -2824.70341953464, 3.563321538783252, 0.0,
    313.15, 283.15, 353.0, "aqueous-molality-infinite-dilution-water-v1"};
p2.reference = ReactionReference{100000.0, ReferenceSourceBasis::CommonMolalityInfiniteDilution};
ReactionSpec p3;
p3.name = "P3";
p3.stoichiometry = {0, -1, 1, 0, 0, 0, 1, 0, 0, -1, 0, 0};
p3.correlation = ReactionLogPolynomial{8.778505217580486, -3615.97963003785, 0.0, 0.0,
    313.15, 283.15, 333.3, "aqueous-molality-infinite-dilution-water-v1"};
p3.reference = ReactionReference{100000.0, ReferenceSourceBasis::CommonMolalityInfiniteDilution};
ReactionSpec p4;
p4.name = "P4";
p4.stoichiometry = {0, -1, 0, 0, 1, 0, 0, 0, 0, 1, 0, -1};
p4.correlation = ReactionLogPolynomial{-10.042964504571048, -3493.0215860719677, 0.0, 0.0,
    313.15, 283.15, 333.3, "aqueous-molality-infinite-dilution-water-v1"};
p4.reference = ReactionReference{100000.0, ReferenceSourceBasis::CommonMolalityInfiniteDilution};
ReactionSpec p5;
p5.name = "P5";
p5.stoichiometry = {0, -1, 1, 0, 0, 0, 0, 0, 0, 1, -1, 0};
p5.correlation = ReactionLogPolynomial{3.3900402036432573, -1322.1443603971813, 0.0, 0.0,
    313.15, 283.15, 333.3, "aqueous-molality-infinite-dilution-water-v1"};
p5.reference = ReactionReference{100000.0, ReferenceSourceBasis::CommonMolalityInfiniteDilution};
problem.reactions = {r1, r2, r3, p1, p2, p3, p4, p5};
