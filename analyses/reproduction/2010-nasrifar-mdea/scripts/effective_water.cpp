// Nasrifar and Tafazzol (2010), section 6, p.7627: effective water = true water plus every declared ion,
// free molecular species preserved. Solves the R1 water case and the two Kamps (2001) MDEA feeds through the
// public native ideal mole-fraction solve, applies the conversion and prints one JSON document on stdout.
#include "core/numerical_policy.hpp"
#include "equilibrium/problem.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace epcsaft;
using namespace epcsaft::equilibrium;

std::map<std::string, double> effective_water_amounts(const IdealMoleFraction& model, const EquilibriumResult& result,
    const std::string& water, const std::vector<std::string>& ions, const std::vector<std::string>& molecular_species) {
    if (result.ideal_identity != &model || !result.has_state() || !result.requested_tolerance_met.value_or(false)
        || result.phase_count != 1 || result.mole_fractions.size() != model.component_ids.size())
        throw std::invalid_argument("effective water requires this ideal model's criterion-meeting liquid result");
    std::vector<bool> seen(model.component_ids.size(), false);
    std::map<std::string, double> converted;
    const auto take = [&](const std::string& name, bool ion) {
        const auto found = std::find(model.component_ids.begin(), model.component_ids.end(), name);
        if (found == model.component_ids.end()) throw std::invalid_argument("unknown conversion species: " + name);
        const auto i = static_cast<std::size_t>(found - model.component_ids.begin());
        if (seen[i] || (model.charge_number[i] != 0) != ion) throw std::invalid_argument("conversion species duplicate or charge role: " + name);
        seen[i] = true;
        const double amount = result.amounts[0] * result.mole_fractions[i];
        if (!std::isfinite(amount) || amount < 0.0) throw std::invalid_argument("invalid conversion amount: " + name);
        return amount;
    };
    double effective = take(water, false);
    for (const auto& ion : ions) effective += take(ion, true);
    if (!std::isfinite(effective)) throw std::invalid_argument("nonfinite effective water");
    converted[water] = effective;
    for (const auto& name : molecular_species) converted[name] = take(name, false);
    if (std::find(seen.begin(), seen.end(), false) != seen.end()) throw std::invalid_argument("incomplete conversion species partition");
    return converted;
}

// Nasrifar Table 1 ln K = a + b/T + c ln T + d T with its fit interval; the reactions the Engine check declares.
const std::array<std::array<double, 6>, 6> coefficients = {{
    {132.899, -13445.9, -22.4773, 0, 273, 498}, {231.456, -12092.1, -36.7816, 0, 273, 498},
    {216.049, -12431.7, -35.4819, 0, 273, 498}, {214.582, -12995.4, -33.5471, 0, 273, 423},
    {-32, -3338, 0, 0, 287, 343}, {-9.4165, -4234.98, 0, 0, 298, 333}}};

ReactionSpec reaction(std::size_t r, std::vector<double> nu) {
    const auto& a = coefficients[r];
    ReactionSpec spec; spec.name = "R" + std::to_string(r + 1); spec.stoichiometry = std::move(nu);
    spec.correlation = ReactionLogPolynomial{a[0], a[1], a[2], a[3], 1.0, a[4], a[5], IDEAL_MOLE_FRACTION_STANDARD_STATE_ID};
    return spec;
}

EquilibriumProblem liquid(double temperature, std::vector<double> feed) {
    EquilibriumProblem p;
    PhaseSpec phase; phase.name = "L"; phase.kind = PhaseKind::Liquid;
    p.phases = {phase};
    p.temperature = FreeValue::fixed(temperature);
    p.pressure = FreeValue::fixed(1.0e5);
    p.feed.kind = FeedKind::Amounts;
    p.feed.values = std::move(feed);
    p.solver_options.values[policy::row] = 1e-10;
    return p;
}

std::string records;
bool passed = true;
void record(const std::string& name, double model, double reference, double tolerance, const char* unit = "mol") {
    const bool ok = std::abs(model - reference) <= tolerance;
    passed &= ok;
    char line[512];
    std::snprintf(line, sizeof(line), "%s\n  {\"comparison_id\": \"%s\", \"model\": %.17g, \"reference\": %.17g, \"atol\": %.17g, \"unit\": \"%s\", \"verdict\": %s}",
        records.empty() ? "" : ",", name.c_str(), model, reference, tolerance, unit, ok ? "true" : "false");
    records += line;
}

int main() {
    IdealMoleFraction water({"water", "OH-", "H3O+"}, {0, -1, 1});
    auto p = liquid(298.15, {1.0, 0.0, 0.0});
    p.reactions = {reaction(0, {-2, 1, 1})};
    const auto r1 = solve_equilibrium(water, p);
    record("R1-effective-water", effective_water_amounts(water, r1, "water", {"OH-", "H3O+"}, {}).at("water"), 1.0, 1e-10);
    std::string message;
    try { effective_water_amounts(water, r1, "water", {"OH-"}, {}); } catch (const std::invalid_argument& error) { message = error.what(); }
    record("refusal-incomplete conversion", message.find("incomplete conversion") != std::string::npos, 1.0, 0.0, "1");
    std::string effective;
    for (bool sulfide : {false, true}) {
        IdealMoleFraction model({"water", "MDEA", sulfide ? "H2S" : "CO2", "OH-", "H3O+",
            sulfide ? "HS-" : "HCO3-", sulfide ? "S--" : "CO3--", "MDEAH+"}, {0,0,0,-1,1,-1,-2,1});
        const double W = 1.0 / 0.01801528, M = sulfide ? 8.001 : 7.994, G = sulfide ? 6.254 : 6.363;
        auto problem = liquid(313.11, {W,M,G,0,0,0,0,0});
        problem.reactions = {reaction(0, {-2,0,0,1,1,0,0,0}),
            reaction(sulfide ? 3 : 1, {sulfide ? -1.0 : -2.0,0,-1,0,1,1,0,0}),
            reaction(sulfide ? 4 : 2, {-1,0,0,0,1,-1,1,0}), reaction(5, {-1,1,0,0,1,0,0,-1})};
        const auto result = solve_equilibrium(model, problem);
        std::vector<double> n = result.mole_fractions;
        for (double& value : n) value *= result.amounts[0];
        const auto out = effective_water_amounts(model, result, "water",
            {"OH-", "H3O+", model.component_ids[5], model.component_ids[6], "MDEAH+"}, {"MDEA", model.component_ids[2]});
        const std::string gas = sulfide ? "H2S" : "CO2";
        record(gas + "-converted-total", out.at("water")+out.at("MDEA")+out.at(gas), result.amounts[0], 1e-10);
        record(gas + "-free-MDEA", out.at("MDEA"), n[1], 1e-10);
        record(gas + "-free-gas", out.at(gas), n[2], 1e-10);
        record(gas + "-effective-water", out.at("water"), W+n[7]+(sulfide ? n[5]+n[6] : 0.0), 1e-10);
        record(gas + "-effective-inventory", out.at("water")+n[1]+(sulfide ? n[2] : 0.0), W+M+(sulfide ? G : 0.0), 1e-10);
        char value[96];
        std::snprintf(value, sizeof(value), "%s\"%s\": %.17g", effective.empty() ? "" : ", ", gas.c_str(), out.at("water"));
        effective += value;
    }
    std::printf("{\"effective_water_mol\": {%s},\n \"comparisons\": [%s\n]}\n", effective.c_str(), records.c_str());
    return passed ? 0 : 1;
}
