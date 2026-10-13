// #73: pressure-chart ln(phi_i) actions against exact identities, the neutral infinite-dilution
// limit and insertion statuses; numerical verification only. Refinement rows for ln(phi_i) and
// the density action through third order are in P4.
#include "properties_checks_shared.hpp"

namespace properties_checks {

namespace {

using namespace epcsaft;
constexpr auto LOG_PHI = PropertyObservable::LogFugacityCoefficient;
constexpr double NOT_AVAILABLE = std::numeric_limits<double>::quiet_NaN();

void verdict_row(const std::string& id, const std::string& quantity, const char* unit, double model, double reference,
                 double atol, double rtol, const std::string& inputs, const char* method, const char* claim) {
    recording_support::compare(("log-fugacity/" + id).c_str(), quantity.c_str(), unit, model, reference, inputs,
        "Relationship", method, claim, "symmetric_close", atol, rtol,
        recording_support::scalar_verdict(model, reference, "symmetric_close", atol, rtol), "log-fugacity");
}

// Sum_i x_i D ln(phi_i) and its term scale, from the rows after `first`.
std::pair<double, double> weighted(const PropertyActionBatch& batch, const std::vector<double>& x, std::size_t first) {
    double sum = 0.0, scale = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const auto& item = batch.results[first + i];
        const double term = item.status == PropertyActionStatus::Available && item.component == i && item.action
            ? x[i] * *item.action : NOT_AVAILABLE;
        sum += term;
        scale += std::abs(term);
    }
    return {sum, scale};
}

// Fixed T and x: sum_i x_i ln(phi_i) = g_res/RT, so its T action is -h_res/(R T^2) and its P actions
// are 1/(rho R T) - 1/P and their P derivatives; at fixed T and P, sum_i x_i dln(phi_i) = 0 (Gibbs-Duhem).
void exact_identities() {
    const HelmholtzModel gross(properties_support::neutral_binary());  // Gross 2001 methane/ethane, k_ij = 0
    const double R = GAS_CONSTANT_J_PER_MOL_K;
    const PressureActionDirection t{1.0, 0.0, {0.0, 0.0}, {}}, p{0.0, 1.0, {0.0, 0.0}, {}};
    const PressureActionDirection c{0.0, 0.0, {0.02, -0.02}, {}};
    // 180 K, 3 MPa: a one-pascal direction lifts dV/V ~ 3e-9, which the power-of-two direction scaling
    // keeps exact. 150 K, 0.1 MPa: a dense liquid with Z = 4.1e-3, where ln Z = ln(P/(R T rho)).
    using Point = std::tuple<const char*, double, double, std::vector<double>>;
    for (const auto& [name, T, P, x] : std::array<Point, 2>{{
             {"compressed-liquid", 180.0, 3.0e6, {0.4, 0.6}}, {"dense-liquid", 150.0, 1.0e5, {0.05, 0.95}}}}) {
        const auto state = gross.state_at_pressure(T, P, x, DensityBranch::Liquid);
        const auto batch = [&, T = T, P = P, x = x](const std::vector<PressureActionDirection>& directions) {
            return gross.pressure_actions(T, P, x, DensityBranch::Liquid,
                {PropertyObservable::ResidualEnthalpy, LOG_PHI}, directions);
        };
        const double rho = state.molar_density;
        const auto inputs = recording_support::fields({recording_support::field("T_K", T),
            recording_support::field("P_Pa", P), recording_support::field("x", x),
            recording_support::field("rho", rho), recording_support::field("Z", state.compressibility_factor)});
        const auto first_t = batch({t}), first_p = batch({p}), second_p = batch({p, p});
        for (std::size_t i = 0; i < x.size(); ++i) {
            const auto& item = first_t.results[1 + i];
            verdict_row(std::string(name) + "/value/" + std::to_string(i), "ln_phi_" + std::to_string(i), "1",
                item.component == i ? item.value.value_or(NOT_AVAILABLE) : NOT_AVAILABLE,
                (*state.log_fugacity_coefficient)[i].value(), 1.0e-12, 1.0e-10, inputs, "state-at-pressure",
                "the component-tagged action value is the fixed-pressure ln(phi_i)");
        }
        const double rho_p = first_p.density_action.at(0), rho_pp = second_p.density_mixed_action.value_or(NOT_AVAILABLE);
        const double h_t = first_t.results[0].value.value_or(NOT_AVAILABLE);
        const double h_p = first_p.results[0].action.value_or(NOT_AVAILABLE);
        const double h_pp = second_p.results[0].action.value_or(NOT_AVAILABLE);
        // Unit directions dT = 1 K and dP = 1 Pa, so each sum carries 1/K and 1/Pa per direction.
        const std::array<std::tuple<const char*, PropertyActionBatch, double, double, const char*>, 6> sums = {{
            {"T", first_t, -h_t / (R * T * T), 1.0e-10, "1/K"},
            {"P", first_p, 1.0 / (rho * R * T) - 1.0 / P, 1.0e-10, "1/Pa"},
            {"PP", second_p, 1.0 / (P * P) - rho_p / (rho * rho * R * T), 1.0e-10, "1/Pa^2"},
            {"TP", batch({t, p}), -h_p / (R * T * T), 1.0e-10, "1/(K Pa)"},
            {"PPP", batch({p, p, p}), -rho_pp / (rho * rho * R * T) + 2.0 * rho_p * rho_p / (rho * rho * rho * R * T)
                - 2.0 / (P * P * P), 1.0e-10, "1/Pa^3"},
            {"TPP", batch({t, p, p}), -h_pp / (R * T * T), 1.0e-9, "1/(K Pa^2)"}}};
        for (const auto& [id, actions, reference, rtol, unit] : sums)
            verdict_row(std::string(name) + "/sum/" + id, std::string("sum_x_D") + id + "_ln_phi", unit,
                weighted(actions, x, 1).first, reference, 0.0, rtol, inputs, "exact-identity",
                "sum_i x_i D ln(phi_i) equals the D action of g_res/RT");
        using Directions = std::vector<PressureActionDirection>;
        for (const auto& [id, directions, unit] : std::array<std::tuple<const char*, Directions, const char*>, 4>{{
                 {"c", {c}, "1"}, {"Tc", {t, c}, "1/K"}, {"Pc", {p, c}, "1/Pa"}, {"TPc", {t, p, c}, "1/(K Pa)"}}}) {
            const auto [sum, scale] = weighted(batch(directions), x, 1);
            verdict_row(std::string(name) + "/gibbs-duhem/" + id, std::string("sum_x_D") + id + "_ln_phi", unit,
                sum, 0.0, 1.0e-10 * scale, 0.0,
                recording_support::fields({inputs, recording_support::field("term_scale", scale)}), "Gibbs-Duhem",
                "at fixed T and P, sum_i x_i D ln(phi_i) vanishes to the rounding of its terms");
        }
    }
}

// Owner decision 2026-09-26: a direction entering a zero-fraction neutral component follows the exact
// infinite-dilution limit, f(0) = 2 f(eps) - f(2 eps) + O(eps^2).
void neutral_infinite_dilution_limit() {
    const HelmholtzModel model(properties_support::neutral_ternary());  // methane, ethane, propane
    const PressureActionDirection d{0.0, 0.0, {0.02, -0.01, -0.01}, {}}, t{1.0, 0.0, {0.0, 0.0, 0.0}, {}};
    constexpr double T = 250.0, P = 2.0e6, eps = 1.0e-7;
    // Outputs: density action, h_res, ln(phi_i) values, then the h_res and ln(phi_i) actions.
    const auto outputs = [&](double trace, const std::vector<PressureActionDirection>& directions) {
        const auto batch = model.pressure_actions(T, P, {trace, 0.5 - trace / 2.0, 0.5 - trace / 2.0},
            DensityBranch::Liquid, {PropertyObservable::ResidualEnthalpy, LOG_PHI}, directions);
        std::vector<double> result = {directions.size() == 1 ? batch.density_action.at(0)
            : batch.density_mixed_action.value_or(NOT_AVAILABLE)};
        for (const auto& item : batch.results) result.push_back(item.value.value_or(NOT_AVAILABLE));
        for (const auto& item : batch.results)
            result.push_back(item.status == PropertyActionStatus::Available ? item.action.value_or(NOT_AVAILABLE)
                : NOT_AVAILABLE);
        return result;
    };
    using Directions = std::vector<PressureActionDirection>;
    // d is a composition direction (dimensionless) and t is dT = 1 K, so Td actions carry 1/K.
    for (const auto& [id, directions, per] : std::array<std::tuple<const char*, Directions, std::string>, 3>{{
             {"d", {d}, ""}, {"dd", {d, d}, ""}, {"Td", {t, d}, "/K"}}}) {
        const auto exact = outputs(0.0, directions), one = outputs(eps, directions), two = outputs(2.0 * eps, directions);
        const std::array<std::pair<std::string, std::string>, 9> labels = {{
            {"density_action", "mol/m^3" + per}, {"residual_enthalpy", "J/mol"}, {"ln_phi_0", "1"}, {"ln_phi_1", "1"},
            {"ln_phi_2", "1"}, {"residual_enthalpy_action", "J/mol" + per}, {"ln_phi_0_action", per.empty() ? "1" : "1/K"},
            {"ln_phi_1_action", per.empty() ? "1" : "1/K"}, {"ln_phi_2_action", per.empty() ? "1" : "1/K"}}};
        for (std::size_t k = 0; k < exact.size(); ++k)
            verdict_row(std::string("neutral-limit/") + id + "/" + std::to_string(k), labels.at(k).first,
                labels.at(k).second.c_str(), exact[k], 2.0 * one[k] - two[k], 1.0e-12, 1.0e-9,
                recording_support::fields({recording_support::field("T_K", T), recording_support::field("P_Pa", P),
                    recording_support::field("eps", eps), recording_support::field("output", k)}),
                "infinite-dilution-extrapolation", "x_CH4 = 0 outputs equal their trace-methane limit");
    }
}

// Only the x ln x ideal entropy is singular at a neutral insertion; an ionic insertion is singular for
// every action, so a simultaneous neutral insertion must not hide it (support_status precedence).
void insertion_statuses() {
    const HelmholtzModel neutral(properties_support::neutral_binary(), {}, p3_binary_record());
    const std::vector<PropertyObservable> calorics = {PropertyObservable::IdealEntropy,
        PropertyObservable::TotalEntropy, PropertyObservable::ResidualEntropy,
        PropertyObservable::IdealEnthalpy, PropertyObservable::TotalEnthalpy, PropertyObservable::ResidualEnthalpy};
    const DensityActionDirection insert{0.0, 0.0, {-0.01, 0.01}, {}};
    const double pressure = neutral.state(300.0, 1000.0, {1.0, 0.0}).pressure;
    auto with_log_phi = calorics;
    with_log_phi.push_back(LOG_PHI);
    const auto density = neutral.density_actions(300.0, 1000.0, {1.0, 0.0}, calorics, {insert});
    const auto fixed_p = neutral.pressure_actions(300.0, pressure, {1.0, 0.0}, DensityBranch::Unspecified,
        with_log_phi, {{0.0, 0.0, insert.composition, {}}});
    bool neutral_ok = fixed_p.density_status == PropertyActionStatus::Available && fixed_p.density_action.size() == 1
        && std::isfinite(fixed_p.density_action[0]) && fixed_p.results.size() == calorics.size() + 2;
    for (const auto* batch : {&density, &fixed_p})
        for (std::size_t k = 0; k < batch->results.size(); ++k) {
            const auto& item = batch->results[k];
            neutral_ok = neutral_ok && item.value && std::isfinite(*item.value)
                && (k < 2 ? item.status == PropertyActionStatus::SupportChange && !item.action
                          : item.status == PropertyActionStatus::Available && item.action && std::isfinite(*item.action));
        }
    recording_support::compare_custom("log-fugacity/insertion/neutral", "status", "bool", neutral_ok ? 1.0 : 0.0, 1.0,
        "methane/ethane x=(1,0), dx=(-.01,.01), T=300 K, rho=1000 mol/m3 and its pressure", "Availability",
        "insertion-status", "only the ideal and total entropy are singular at a neutral insertion", "absolute",
        neutral_ok, "log-fugacity", "insertion/neutral", "status", "T,rho|P,x", "support_change");
    // MEA system: water/MEA with CO2 (neutral) and the MEAH+/MEACOO- pair entering together; the
    // entropy rows carry the precedence (ionic before neutral).
    ThermochemistryRecord record;
    record.standard_pressure = 1.0e5;
    record.component_intervals.assign(9, {IdealInterval{250.0, 350.0, true, true, properties_support::p3_forms()[0].second}});
    const HelmholtzModel mea(properties_support::nine_species(), {}, record);
    const std::vector<double> x = {0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const std::vector<double> dx = {-0.03, 0.0, 0.01, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0};
    const auto liquid = mea.state_at_pressure(313.15, 1.0e5, x, DensityBranch::Liquid);
    const std::vector<PropertyObservable> mixed = {PropertyObservable::ResidualEnthalpy,
        PropertyObservable::ResidualEntropy, PropertyObservable::IdealEntropy, PropertyObservable::TotalEntropy};
    const auto ionic_rho = mea.density_actions(313.15, liquid.molar_density, x, mixed, {{0.0, 0.0, dx, {}}});
    auto with_ln_phi = mixed;
    with_ln_phi.push_back(LOG_PHI);
    const auto ionic_p = mea.pressure_actions(313.15, 1.0e5, x, DensityBranch::Liquid, with_ln_phi, {{0.0, 0.0, dx, {}}});
    bool ionic_ok = ionic_p.density_status == PropertyActionStatus::IonicInsertion && ionic_p.density_action.empty()
        && ionic_p.results.size() == mixed.size() + x.size();
    for (const auto* batch : {&ionic_rho, &ionic_p})
        for (const auto& item : batch->results)
            ionic_ok = ionic_ok && item.status == PropertyActionStatus::IonicInsertion && !item.action;
    recording_support::compare_custom("log-fugacity/insertion/neutral-and-ionic", "status", "bool",
        ionic_ok ? 1.0 : 0.0, 1.0, "nine-species x=(0.9,0.1,0...), dx=(-.03,0,.01,.01,.01,0...), 313.15 K, 0.1 MPa",
        "Availability", "insertion-status", "a charge-balanced ionic insertion stays singular beside a neutral one",
        "absolute", ionic_ok, "log-fugacity", "insertion/neutral-and-ionic", "status", "T,rho|P,x", "ionic_insertion");
    // Without charged species the state has no ion ln(phi_i); the neutral rows keep the state's values.
    const auto neutral_state = mea.pressure_actions(313.15, 1.0e5, x, DensityBranch::Liquid, {LOG_PHI},
        {{1.0, 0.0, std::vector<double>(x.size(), 0.0), {}}});
    bool absent_ions = neutral_state.results.size() == x.size();
    for (std::size_t i = 0; absent_ions && i < x.size(); ++i) {
        const auto& item = neutral_state.results[i];
        const auto& reference = (*liquid.log_fugacity_coefficient)[i];
        absent_ions = item.component == i && (reference
            ? item.status == PropertyActionStatus::Available && item.action
                && std::abs(*item.value - *reference) <= 1.0e-10 * std::max(1.0, std::abs(*reference))
            : item.status == PropertyActionStatus::IonicInsertion && !item.value && !item.action);
    }
    recording_support::compare_custom("log-fugacity/insertion/absent-ions", "status", "bool",
        absent_ions ? 1.0 : 0.0, 1.0, "nine-species x=(0.9,0.1,0...), dT=1, 313.15 K, 0.1 MPa", "Availability",
        "state-at-pressure", "ion rows without charged species are ionic_insertion with no value, as in the state",
        "absolute", absent_ions, "log-fugacity", "insertion/absent-ions", "status", "T,P,x", "ionic_insertion");
    // Gross-Vrabec dipolar term: with every dipole absent its A2-only branch omits A3 beyond second
    // composition order, so a neutral insertion in a model with the term stays a support change.
    auto dipolar = properties_support::neutral_binary();
    dipolar.dipole_moment_debye = {1.5, 0.0};
    dipolar.dipole_count = {1.0, 0.0};
    const HelmholtzModel polar(dipolar);
    const auto polar_insert = polar.pressure_actions(180.0, 3.0e6, {0.0, 1.0}, DensityBranch::Liquid,
        {PropertyObservable::ResidualEnthalpy, LOG_PHI}, {{0.0, 0.0, {0.01, -0.01}, {}}});
    bool dipolar_ok = polar_insert.density_status == PropertyActionStatus::SupportChange
        && polar_insert.density_action.empty() && polar_insert.results.size() == 3;
    for (const auto& item : polar_insert.results)
        dipolar_ok = dipolar_ok && item.status == PropertyActionStatus::SupportChange && item.value && !item.action;
    recording_support::compare_custom("log-fugacity/insertion/absent-dipolar-family", "status", "bool",
        dipolar_ok ? 1.0 : 0.0, 1.0, "methane(1.5 D)/ethane x=(0,1), dx=(.01,-.01), 180 K, 3 MPa", "Availability",
        "insertion-status", "a neutral insertion in a model with the dipolar term stays a support change", "absolute",
        dipolar_ok, "log-fugacity", "insertion/absent-dipolar-family", "status", "T,P,x", "support_change");
}

}  // namespace

void log_fugacity_actions() {
    exact_identities();
    neutral_infinite_dilution_limit();
    insertion_statuses();
}

}  // namespace properties_checks
