// Amine-Thermodynamics c58e337ea, MDEA_slope_bound declaration and warm-start-reproduction (#275).
#include "kim_pure_water.inc"
#include "mdea-bubble.inc"
const auto coefficient = [&](const char* identity, double value) {
    const auto entry = std::find_if(document.entries.begin(), document.entries.end(),
        [&](const ParameterEntry& e) { return e.identity == identity; });
    if (entry == document.entries.end()) throw std::runtime_error(identity);
    entry->magnitude = value;
};
coefficient("component/protonated-methyldiethanolamine/segment_diameter", 3.5959113120691675);
coefficient("component/protonated-methyldiethanolamine/dispersion_energy_over_k", 243.52151216647658);
coefficient("component/protonated-methyldiethanolamine/debye_huckel_diameter", 3.5959113120691675);
coefficient("component/protonated-methyldiethanolamine/born_diameter", 3.5959113120691675);
coefficient("pair/carbon-dioxide/protonated-methyldiethanolamine/k_ij", 0.3784408781206014);
coefficient("pair/mdea-water/k_ij/a1", 0.0011694693134037676);
coefficient("pair/water/protonated-methyldiethanolamine/k_ij", -0.2518751165455549);
coefficient("pair/water/protonated-methyldiethanolamine/k_ij/reciprocal_temperature_slope", -31.02545057904153);
coefficient("pair/bicarbonate-anion/protonated-methyldiethanolamine/k_ij", -0.03126710546865702);
document.entries.erase(std::remove_if(document.entries.begin(), document.entries.end(),
    [](const ParameterEntry& e) { return e.identity == "pair/bicarbonate-anion/water/k_ij"
        || e.identity == "pair/carbonate-anion/water/k_ij"; }), document.entries.end());
document.entries.insert(document.entries.end(), {
    {"pair/bicarbonate-anion/water/k_ij", "k_ij_temperature_law", {1, 2, 1}, 313.15, "kelvin"},
    {"pair/bicarbonate-anion/water/k_ij/a0", "temperature_law_coefficient", {1, 0}, 0.1706931974036258, "dimensionless"},
    {"pair/bicarbonate-anion/water/k_ij/a1", "temperature_law_coefficient", {1, 1}, 0.0039311876328716865, "1 / kelvin"},
    {"pair/carbonate-anion/water/k_ij", "k_ij_temperature_law", {2, 3, 1}, 313.15, "kelvin"},
    {"pair/carbonate-anion/water/k_ij/a0", "temperature_law_coefficient", {2, 0}, -0.9705531093245405, "dimensionless"},
    {"pair/carbonate-anion/water/k_ij/a1", "temperature_law_coefficient", {2, 1}, 0.0019631260450306285, "1 / kelvin"},
});
const HelmholtzModel model(build_parameters(document));
for (auto& reaction : bubble.reactions) reaction.reference->reference_pressure_pa = 100000.0;
bubble.temperature = FreeValue::fixed(313.15);
bubble.feed.values = {1e-5, 211.7659574468085, 0, 0, 0, 0, 1, 0};
