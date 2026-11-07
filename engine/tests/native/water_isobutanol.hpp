// Khudaida et al. (2026), Table 5/7 and the retained parameter packet.
// Packet: analyses/2026-khudaida-idaes-settler/inputs/data/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json
// SHA-256: 8babed899e7347e58c403c7aff24a7fdd6d60411af3b4d81eb6aa90c2431dfa4
epcsaft::Parameters water_isobutanol() {
    epcsaft::ParameterDocument d;
    d.component_ids = {"water", "isobutanol"};
    d.site_component = {1, 1, 0, 0};
    d.site_multiplicity = {1, 1, 1, 1};
    d.base_choice = "pc-saft";
    d.association_choice = "general-site";
    d.electrolyte_choice = "born";
    d.ionic_dispersion_choice = "exclude-same-sign-ion-pairs";
    d.permittivity_choice = "ion-fraction-suppression";
    d.polar_choice = "none";
    d.entries = {
        {"component/water/molar_mass", "molar_mass", {0}, 18.01528, "gram / mole"},
        {"component/water/charge_number", "charge_number", {0}, 0.0, "elementary-charge"},
        {"component/water/dispersion_energy_over_k", "dispersion_energy_over_k", {0}, 353.95, "kelvin"},
        {"component/water/segment_count", "segment_count", {0}, 1.2047, "dimensionless"},
        {"component/water/solvation_factor", "solvation_factor", {0}, 1.5, "dimensionless"},
        {"component/isobutanol/molar_mass", "molar_mass", {1}, 74.1216, "gram / mole"},
        {"component/isobutanol/charge_number", "charge_number", {1}, 0.0, "elementary-charge"},
        {"component/isobutanol/dispersion_energy_over_k", "dispersion_energy_over_k", {1}, 257.13, "kelvin"},
        {"component/isobutanol/segment_count", "segment_count", {1}, 2.844, "dimensionless"},
        {"component/isobutanol/segment_diameter", "sigma_constant", {1}, 3.5561, "angstrom"},
        {"component/isobutanol/solvation_factor", "solvation_factor", {1}, 1.5, "dimensionless"},
        {"component/isobutanol/relative_permittivity/constant/constant", "relative_permittivity", {1}, 17.2, "dimensionless"},
        {"component/water/relative_permittivity/constant/constant", "relative_permittivity", {0}, 78.09, "dimensionless"},
        {"component/water/segment_diameter/constant-plus-sum-of-exponentials/constant", "sigma_constant", {0}, 2.7927, "angstrom"},
        {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/amplitude", "sigma_term_amplitude", {0, 0}, 10.11, "angstrom"},
        {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/exponent_coefficient", "sigma_term_exponent", {0, 0}, -0.01775, "1 / kelvin"},
        {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/amplitude", "sigma_term_amplitude", {1, 0}, -1.417, "angstrom"},
        {"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/exponent_coefficient", "sigma_term_exponent", {1, 0}, -0.01146, "1 / kelvin"},
        {"pair/isobutanol/water/k_ij", "k_ij", {1, 0}, -0.0071, "dimensionless"},
        {"model/electrolyte/c_shell", "born_c_shell", {}, 1.0, "dimensionless"},
        {"model/electrolyte/c_dielectric", "born_c_dielectric", {}, 1.0, "dimensionless"},
        {"model/relative_permittivity/ion_fraction_suppression_coefficient", "ion_suppression_coefficient", {}, 7.01, "dimensionless"},
        {"model/relative_permittivity/ionic_region_relative_permittivity", "ionic_region_relative_permittivity", {}, 8.0, "dimensionless"},
        {"association/isobutanol/a/isobutanol/b/energy_over_k", "association_energy_over_k", {0, 1}, 2514.3, "kelvin"},
        {"association/isobutanol/a/isobutanol/b/volume", "association_volume", {0, 1}, 0.00391, "dimensionless"},
        {"wolbach-sandler", "association_combining_rule", {0, 3, 0, 1, 2, 3}, 0.0, "dimensionless"},
        {"wolbach-sandler", "association_combining_rule", {1, 2, 0, 1, 2, 3}, 0.0, "dimensionless"},
        {"association/water/a/water/b/energy_over_k", "association_energy_over_k", {2, 3}, 2425.7, "kelvin"},
        {"association/water/a/water/b/volume", "association_volume", {2, 3}, 0.04509, "dimensionless"},
    };
    return epcsaft::build_parameters(d);
}
