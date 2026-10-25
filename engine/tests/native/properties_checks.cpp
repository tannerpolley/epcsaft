#include "recording_support.hpp"

// The case registry. Each case function is defined in one properties_*_checks.cpp unit.
namespace properties_checks {
void residual_cv_temperature_extraction();
void born_parameter_columns();
void reference_f8();
void reference_f9();
void reference_limits();
void residual_properties_and_caloric_ladders();
void residual_internal_energy_identity_ladder();
void derivatives_and_limits();
void finite_difference_gradient();
void hessian_reference();
void hessian_fd_and_pressure_slope();
void same_sign_pair_is_inert();
void parameter_derivative_and_helpers();
void core_c2_parameter_coordinates();
void density_and_reference();
void properties_p2();
void properties_thermochemistry_p3();
void properties_p4();
void cache_invariance();
void modular_permittivity_actions();
void electrolyte_parameter_actions();
void energy_e5_ladder();
void energy_action_ladders();
void neutral_property_pressure_hessian();
void log_fugacity_actions();
void liquid_reference_branch_evidence();
void joule_thomson_identities();
void temperature_laws();
void cation_mixing();
}  // namespace properties_checks

int main(int argc, char** argv) {
    using namespace properties_checks;
    return recording_support::run(argc, argv, {
        {{"native.properties.residual-cv-temperature-extraction", "ResidualCvTemperature"}, residual_cv_temperature_extraction},
        {{"native.properties.born-parameter", "BornParameter"}, born_parameter_columns},
        {{"native.properties.reference-f8", "F8"}, reference_f8, "F8"},
        {{"native.properties.reference-f9", "F9"}, reference_f9, "F9"},
        {{"native.properties.reference-limits", "ReferenceLimits"}, reference_limits, "ReferenceLimits"},
        {{"native.properties.p1", "P1"}, residual_properties_and_caloric_ladders, "P1"},
        {{"native.properties.f1", "F1"}, residual_internal_energy_identity_ladder, "F1"},
        {{"native.properties.m2", "M2"}, derivatives_and_limits},
        {{"native.properties.m3", "M3"}, finite_difference_gradient},
        {{"native.properties.m4", "M4"}, hessian_reference},
        {{"native.properties.m5", "M5"}, hessian_fd_and_pressure_slope},
        {{"native.properties.m6", "M6"}, same_sign_pair_is_inert},
        {{"native.properties.m7", "M7"}, parameter_derivative_and_helpers},
        {{"native.properties.c2", "C2"}, core_c2_parameter_coordinates},
        {{"native.properties.m8", "M8"}, density_and_reference},
        {{"native.properties.p2", "P2"}, properties_p2, "P2"},
        {{"native.properties.p3", "P3"}, properties_thermochemistry_p3, "P3"},
        {{"native.properties.p4", "P4"}, properties_p4, "P4"},
        {{"native.properties.m9", "M9"}, cache_invariance},
        {{"native.properties.issue26", "Issue26"}, modular_permittivity_actions},
        {{"native.properties.electrolyte-coordinates", "ElectrolyteCoordinates"}, electrolyte_parameter_actions},
        {{"native.properties.e5-ladder", "C4"}, energy_e5_ladder},
        {{"native.properties.issue26e", "Issue26E"}, energy_action_ladders},
        {{"native.properties.issue10", "Issue10"}, neutral_property_pressure_hessian},
        {{"native.properties.log-fugacity", "LogFugacity"}, log_fugacity_actions},
        {{"native.properties.cation-mixing", "CationMixing"}, cation_mixing},
        {{"native.properties.temperature-laws", "TemperatureLaws"}, temperature_laws},
        {{"native.properties.issue423", "Issue423"}, liquid_reference_branch_evidence},
        {{"native.properties.joule-thomson", "JouleThomson"}, joule_thomson_identities},
    });
}
