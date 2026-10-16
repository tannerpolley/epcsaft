// Temperature laws q(T) = sum_p a_p (T - T_r)^p + s (1/T - 1/T_r) for k_ij, l_ij, association energy and
// volume, and the Wolbach--Sandler volume correction (#567): published values, the pinned four-site case,
// exact derivatives against step-refined central differences (verification only) and the domain failure.
#include "properties_checks_shared.hpp"
#include "regression/observations.hpp"
#include "surface/parameter_builder.hpp"

namespace properties_checks {
namespace {

using epcsaft::ParameterEntry;
using epcsaft::ParameterFamily;
using recording_support::number;

constexpr std::array<double, 3> TEMPERATURES = {300.0, 320.0, 340.0};
constexpr double RHO = 100.0;  // mol/m3
const std::vector<double> X = {0.4, 0.6};

// |model - reference| <= tolerance max(scale, |reference|), recorded with its inputs.
void normalized(const std::string& id, const char* quantity, const char* unit, double model, double reference,
                double scale, double tolerance, const std::string& inputs, const char* claim) {
    const double allowance = tolerance * std::max(scale, std::abs(reference));
    recording_support::compare(id.c_str(), quantity, unit, model, reference, inputs, "Identity", "temperature_law",
        claim, "absolute", allowance, 0.0, std::abs(model - reference) <= allowance, "temperature-laws", id.c_str());
}

ParameterEntry entry(std::string identity, std::string field, std::vector<std::size_t> index, double magnitude,
                     std::string unit) {
    return {std::move(identity), std::move(field), std::move(index), magnitude, std::move(unit)};
}

// Units of a_0..a_3 and s for a dimensionless law and for an energy law (K).
const std::array<const char*, 5> DIMENSIONLESS = {"dimensionless", "1 / kelvin", "1 / kelvin ** 2", "1 / kelvin ** 3", "kelvin"};
const std::array<const char*, 5> ENERGY = {"kelvin", "dimensionless", "1 / kelvin", "1 / kelvin ** 2", "kelvin ** 2"};

// One law [law, a, b] of `field` with the listed coefficients (a_0..a_3, then s; NaN omits a term).
void add_law(std::vector<ParameterEntry>& entries, const std::string& name, const std::string& field, std::size_t law,
             std::size_t a, std::size_t b, double reference, const std::array<double, 5>& coefficients) {
    entries.push_back(entry(name, field + "_temperature_law", {law, a, b}, reference, "kelvin"));
    const auto& units = field == "association_energy_over_k" ? ENERGY : DIMENSIONLESS;
    for (std::size_t p = 0; p < 5; ++p) {
        if (std::isnan(coefficients[p])) continue;
        entries.push_back(p < 4 ? entry(name + "/a" + std::to_string(p), "temperature_law_coefficient", {law, p}, coefficients[p], units[p])
                                : entry(name + "/s", "temperature_law_reciprocal_slope", {law}, coefficients[p], units[p]));
    }
}

epcsaft::ParameterDocument binary_document(bool sites) {
    epcsaft::ParameterDocument document;
    document.component_ids = {"0", "1"};
    document.base_choice = "pc-saft";
    document.association_choice = sites ? "general-site" : "none";
    document.electrolyte_choice = document.permittivity_choice = "none";
    const double m[] = {4.5354, 2.1945}, sigma[] = {2.6019, 2.229}, eps[] = {204.0438, 141.66}, mass[] = {0.06108, 0.01801528};
    for (std::size_t i = 0; i < 2; ++i) {
        auto& e = document.entries;
        e.push_back(entry("", "molar_mass", {i}, mass[i], "kilogram / mole"));
        e.push_back(entry("", "charge_number", {i}, 0.0, "elementary-charge"));
        e.push_back(entry("", "segment_count", {i}, m[i], "dimensionless"));
        e.push_back(entry("", "sigma_constant", {i}, sigma[i], "angstrom"));
        e.push_back(entry("", "dispersion_energy_over_k", {i}, eps[i], "kelvin"));
    }
    if (sites) {
        document.site_component = {0, 0, 1, 1};
        document.site_multiplicity = {2, 1, 2, 2};
    }
    return document;
}

const double NAN_TERM = std::numeric_limits<double>::quiet_NaN();
constexpr double T_R = 300.0;
// The design's six law vectors [a0, a1, a2, a3, s]; slot of coefficient p of law n is 5 n + p.
const std::array<double, 5> K01 = {-0.052, 1e-4, 2e-6, 1e-8, 10.0}, L01 = {0.01, 2e-5, 1e-7, -1e-9, 1.0},
    E00 = {2383.4744, 0.8, 0.002, 1e-5, 3000.0}, V00 = {0.118488, 2e-4, -5e-7, 2e-9, 0.05},
    E11 = {1804.17, -0.4, 0.001, -3e-6, 1500.0}, C01 = {0.15, 1e-3, -2e-6, 1e-8, 3.0};
const std::array<double, 5> DIMENSIONLESS_SCALE = {0.01, 1e-3, 1e-5, 1e-7, 10.0}, ENERGY_SCALE = {100.0, 1.0, 0.01, 1e-4, 3000.0};
const std::array<const char*, 6> LAW_NAMES = {"k01", "l01", "E00", "v00", "E11", "c01"};

// The pinned four-site mixture: explicit self edges (0,1), (2,3); thermal Wolbach--Sandler cross edges
// (0,3), (1,2) sharing c01. `legacy` replaces every law by the given constants and drops c01.
epcsaft::Parameters four_site(const std::array<double, 5>& correction, const std::array<double, 6>* legacy = nullptr) {
    auto document = binary_document(true);
    auto& e = document.entries;
    if (legacy == nullptr) {
        add_law(e, "manufactured/k01", "k_ij", 0, 0, 1, T_R, K01);
        add_law(e, "manufactured/l01", "l_ij", 1, 0, 1, T_R, L01);
        add_law(e, "manufactured/E00", "association_energy_over_k", 2, 0, 1, T_R, E00);
        add_law(e, "manufactured/v00", "association_volume", 3, 0, 1, T_R, V00);
        add_law(e, "manufactured/E11", "association_energy_over_k", 4, 2, 3, T_R, E11);
        add_law(e, "manufactured/c01", "association_volume_correction", 5, 0, 3, T_R, correction);
        e.push_back(entry("manufactured/c01", "association_volume_correction_temperature_law", {5, 1, 2}, T_R, "kelvin"));
    } else {
        const auto& v = *legacy;  // k, l, E00, v00, E11 constants
        e.push_back(entry("k01", "k_ij", {0, 1}, v[0], "dimensionless"));
        e.push_back(entry("l01", "l_ij", {0, 1}, v[1], "dimensionless"));
        e.push_back(entry("E00", "association_energy_over_k", {0, 1}, v[2], "kelvin"));
        e.push_back(entry("v00", "association_volume", {0, 1}, v[3], "dimensionless"));
        e.push_back(entry("E11", "association_energy_over_k", {2, 3}, v[4], "kelvin"));
    }
    e.push_back(entry("v11", "association_volume", {2, 3}, 0.2039, "dimensionless"));
    e.push_back(entry("wolbach-sandler", "association_combining_rule", {0, 3, 0, 1, 2, 3}, 0.0, "dimensionless"));
    e.push_back(entry("wolbach-sandler", "association_combining_rule", {1, 2, 0, 1, 2, 3}, 0.0, "dimensionless"));
    return epcsaft::build_parameters(document);
}

double law_scale(std::size_t slot) { return (slot / 5 == 2 || slot / 5 == 4 ? ENERGY_SCALE : DIMENSIONLESS_SCALE)[slot % 5]; }

std::vector<epcsaft::ActiveParameter> coefficients() {
    std::vector<epcsaft::ActiveParameter> result;
    for (std::size_t slot = 0; slot < 30; ++slot) result.push_back({ParameterFamily::TemperatureLawCoefficient, {slot}});
    return result;
}

// One output (A^res/(nRT), P/Pa, h^res/(J/mol)) and its exact actions d_T^n and d_u, d_T d_u per coefficient
// coordinate u = (a - origin)/scale, from the model's own derivative owners.
struct Actions final {
    std::array<double, 3> value{};
    std::array<std::array<double, 3>, 3> thermal{};  // [output][order - 1]
    std::array<std::array<double, 30>, 3> coefficient{}, mixed{};
};

Actions actions(const epcsaft::HelmholtzModel& model, double T) {
    using epcsaft::PropertyObservable;
    constexpr double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const double V = 1.0 / RHO;
    Actions out;
    const auto state = model.state(T, RHO, X);
    const auto d = model.evaluate(T, V, X, epcsaft::DerivativeOrder::Hessian);
    const std::size_t width = d.gradient.size();
    const auto unit = [&](std::size_t i) { std::vector<double> e(width, 0.0); e[i] = 1.0; return e; };
    const auto energy = [&](std::vector<std::vector<double>> directions) {
        return model.residual_energy_action(T, V, X, directions).action.value();
    };
    const auto h = [&](std::vector<epcsaft::DensityActionDirection> directions) {
        return model.density_actions(T, RHO, X, {PropertyObservable::ResidualEnthalpy}, directions).results.front().action.value();
    };
    const epcsaft::DensityActionDirection along_t{1.0, 0.0, {0.0, 0.0}, {}};
    const double phi_v = d.gradient[1], phi_tv = d.hessian[1], phi_ttv = energy({unit(0), unit(0), unit(1)});
    out.value = {d.value, state.pressure, state.residual_enthalpy};
    out.thermal[0] = {d.gradient[0], d.hessian[0], energy({unit(0), unit(0), unit(0)})};
    out.thermal[1] = {R / V - R * phi_v - R * T * phi_tv, -2.0 * R * phi_tv - R * T * phi_ttv,
                      -3.0 * R * phi_ttv - R * T * energy({unit(0), unit(0), unit(0), unit(1)})};
    out.thermal[2] = {h({along_t}), h({along_t, along_t}), h({along_t, along_t, along_t})};
    for (std::size_t slot = 0; slot < 30; ++slot) {
        const std::size_t u = 4 + slot;
        const double scale = law_scale(slot), phi_vu = d.hessian[1 * width + u];
        const epcsaft::DensityActionDirection along_u{0.0, 0.0, {0.0, 0.0},
            {{ParameterFamily::TemperatureLawCoefficient, {slot}, scale}}};
        out.coefficient[0][slot] = scale * d.gradient[u];
        out.mixed[0][slot] = scale * d.hessian[u];
        out.coefficient[1][slot] = -R * T * scale * phi_vu;
        out.mixed[1][slot] = scale * (-R * phi_vu - R * T * energy({unit(0), unit(1), unit(u)}));
        out.coefficient[2][slot] = h({along_u});
        out.mixed[2][slot] = h({along_t, along_u});
    }
    return out;
}

// Central differences of `f` over the step ladder; the best normalized error must be <= 1e-6 and the second
// rung must improve on the first as second-order differencing does, unless the first is already at that floor.
void ladder(const std::string& id, const char* unit, double exact, const std::vector<double>& steps,
            const std::function<double(double)>& f, const std::string& inputs) {
    std::vector<double> differences, errors;
    for (const double step : steps) {
        differences.push_back((f(step) - f(-step)) / (2.0 * step));
        errors.push_back(std::abs(differences.back() - exact) / std::max(1.0, std::abs(exact)));
    }
    const auto best = static_cast<std::size_t>(std::min_element(errors.begin(), errors.end()) - errors.begin());
    const bool second_order = errors[0] <= 1.0e-6 || errors[1] <= 0.35 * errors[0];
    CHECK(second_order);
    normalized(id, "derivative", unit, differences[best], exact, 1.0, 1.0e-6,
        recording_support::fields({inputs, recording_support::field("steps", steps),
            recording_support::field("differences", differences), recording_support::field("normalized_errors", errors)}),
        "exact action equals step-refined central differences (verification only)");
}

epcsaft::Parameters shifted(epcsaft::Parameters parameters, std::size_t slot, double u) {
    parameters.temperature_law_coefficient[slot] += law_scale(slot) * u;
    return parameters;
}

}  // namespace

void temperature_laws() {
    using recording_support::field;
    using recording_support::fields;
    // 1. Published laws. Held 2014 Table 2 water/Na+ 2.37999 - 0.007981 T/K is recentred exactly at 298.15 K:
    // a0 = 2.37999 - 0.007981(298.15) = 0.00045485. Derived values are equation values, not observations.
    struct Source final { const char* id; double reference; std::array<double, 5> law; double T, expected, tolerance; };
    const std::vector<Source> sources = {
        {"held2014/water-benzene/298.15", 298.15, {0.0255, 6.07e-4, NAN_TERM, NAN_TERM, NAN_TERM}, 298.15, 0.0255, 1e-14},
        {"held2014/water-benzene/273.15", 298.15, {0.0255, 6.07e-4, NAN_TERM, NAN_TERM, NAN_TERM}, 273.15, 0.010325, 1e-12},
        {"held2014/water-na/273.15", 298.15, {0.00045485, -0.007981, NAN_TERM, NAN_TERM, NAN_TERM}, 273.15, 0.19997985, 1e-12},
        {"held2014/water-na/298.15", 298.15, {0.00045485, -0.007981, NAN_TERM, NAN_TERM, NAN_TERM}, 298.15, 0.00045485, 1e-12},
        {"uyan2015/co2-water/304", 298.0, {-0.022, 0.00042, -0.0000017, NAN_TERM, NAN_TERM}, 304.0, -0.0195412, 1e-12},
        {"uyan2015/co2-water/348", 298.0, {-0.022, 0.00042, -0.0000017, NAN_TERM, NAN_TERM}, 348.0, -0.00525, 1e-12},
        {"uyan2015/co2-water/422", 298.0, {-0.022, 0.00042, -0.0000017, NAN_TERM, NAN_TERM}, 422.0, 0.0039408, 1e-12},
        {"uyan2015/mdea-water/313", 298.0, {-0.015275, -0.0000524, 0.0000079, NAN_TERM, NAN_TERM}, 313.0, -0.0142835, 1e-12},
        {"uyan2015/mdea-water/353", 298.0, {-0.015275, -0.0000524, 0.0000079, NAN_TERM, NAN_TERM}, 353.0, 0.0057405, 1e-12},
        {"bulow2021/water-mdea/298.15", 298.15, {-0.118, -7.59e-4, -1.32e-5, 1.25e-7, NAN_TERM}, 298.15, -0.118, 1e-14},
        {"bulow2021/water-mdea/333.15", 298.15, {-0.118, -7.59e-4, -1.32e-5, 1.25e-7, NAN_TERM}, 333.15, -0.155375625, 1e-12},
    };
    for (const auto& source : sources) {
        auto document = binary_document(false);
        add_law(document.entries, source.id, "k_ij", 0, 0, 1, source.reference, source.law);
        const auto parameters = epcsaft::build_parameters(document);
        const double k = epcsaft::temperature_law_values(source.T, parameters).at(0);
        recording_support::compare(source.id, "k_ij", "1", k, source.expected, field("T_K", source.T), "Identity",
            "temperature_law", "a published k_ij(T) law gives its tabulated or derived value", "absolute",
            source.tolerance, 0.0, std::abs(k - source.expected) <= source.tolerance, "temperature-laws", source.id);
    }
    // Uyan 2015 Table 5: constant water/MDEA volume correction c = 0.015 on a Wolbach--Sandler edge.
    {
        auto document = binary_document(true);
        auto& e = document.entries;
        e.push_back(entry("E00", "association_energy_over_k", {0, 1}, 2383.4744, "kelvin"));
        e.push_back(entry("v00", "association_volume", {0, 1}, 0.118488, "dimensionless"));
        e.push_back(entry("E11", "association_energy_over_k", {2, 3}, 1804.17, "kelvin"));
        e.push_back(entry("v11", "association_volume", {2, 3}, 0.2039, "dimensionless"));
        e.push_back(entry("k01", "k_ij", {0, 1}, 0.0, "dimensionless"));
        e.push_back(entry("uyan2015/water-mdea/c", "association_volume_correction", {0, 0, 3}, 0.015, "dimensionless"));
        e.push_back(entry("wolbach-sandler", "association_combining_rule", {0, 3, 0, 1, 2, 3}, 0.0, "dimensionless"));
        const auto parameters = epcsaft::build_parameters(document);
        const double c = epcsaft::temperature_law_values(298.0, parameters).at(0);
        recording_support::compare("uyan2015/water-mdea/c", "c", "1", c, 0.015, field("T_K", 298.0), "Identity",
            "temperature_law", "a constant source volume correction is admitted unchanged", "absolute", 1e-14, 0.0,
            std::abs(c - 0.015) <= 1e-14, "temperature-laws", "uyan2015/water-mdea/c");
    }

    // 2. Pinned four-site case at 320 K: law values, absolute Delta_03 and the corrected/legacy ratio.
    const auto parameters = four_site(C01);
    const double T = 320.0;
    const auto laws = epcsaft::temperature_law_values(T, parameters);
    const std::array<double, 6> at_320 = {-0.0512033333333333333333333333333333, 0.0102236666666666666666666666666667,
                                          2399.7294, 1467523.0 / 12000000.0, 1796.2335, 0.168655};
    for (std::size_t law = 0; law < 6; ++law)
        normalized(std::string("four-site/law/") + LAW_NAMES[law], LAW_NAMES[law], law == 2 || law == 4 ? "K" : "1",
                   laws[law], at_320[law], 1.0, 1e-12, field("T_K", T), "the law equals its exact decimal value");
    const auto sigma = epcsaft::segment_sigma(T, parameters), diameter = epcsaft::segment_diameter(T, parameters, sigma);
    const double delta_03 = epcsaft::association_delta(T, RHO, X, parameters, sigma, diameter)[3];
    const double d0 = 2.6019 * (1.0 - 0.12 * std::exp(-3.0 * 204.0438 / 320.0));
    const double d1 = 2.229 * (1.0 - 0.12 * std::exp(-3.0 * 141.66 / 320.0));
    const double rho_hat = 100.0 * 6.02214076e-7, pi = 3.14159265358979323846;
    const double zeta2 = pi * rho_hat / 6.0 * (0.4 * 4.5354 * d0 * d0 + 0.6 * 2.1945 * d1 * d1);
    const double zeta3 = pi * rho_hat / 6.0 * (0.4 * 4.5354 * d0 * d0 * d0 + 0.6 * 2.1945 * d1 * d1 * d1);
    const double u = d0 * d1 / (d0 + d1), b = 1.0 - zeta3;
    const double g01 = 1.0 / b + 3.0 * u * zeta2 / (b * b) + 2.0 * u * u * zeta2 * zeta2 / (b * b * b);
    const double size = std::sqrt(2.6019 * 2.229) * 2969329.0 / 3000000.0;
    const double reference = g01 * (std::exp(2097.98145 / 320.0) - 1.0) * size * size * size
        * std::sqrt(1467523.0 / 12000000.0 * 0.2039) * 0.831345;
    normalized("four-site/delta03/absolute", "Delta_03", "angstrom^3", delta_03, reference, 1.0, 1e-12,
        fields({field("T_K", T), field("rho_mol_m3", RHO), field("reference_angstrom3", reference)}),
        "Delta_03 equals the closed Uyan Eq. 9 expression with one size factor and (1 - c) outside the root");
    const std::array<double, 6> legacy = {-0.052, 30671.0 / 3000000.0, 2399.7294, 1467523.0 / 12000000.0, 1796.2335};
    const auto legacy_parameters = four_site(C01, &legacy);
    const double legacy_delta = epcsaft::association_delta(T, RHO, X, legacy_parameters, sigma, diameter)[3];
    normalized("four-site/delta03/ratio", "ratio", "1", delta_03 / legacy_delta, 0.831345, 1.0, 1e-12,
        field("legacy_delta_angstrom3", legacy_delta), "the corrected strength is (1 - c) times the legacy rule");

    // Scalar law derivatives on the tape against q^(n) = sum_{p>=n} p!/(p-n)! a_p d^(p-n) + (-1)^n n! s/T^(n+1),
    // dq/da_p = d^p, dq/ds = 1/T - 1/T_r and their T derivatives, d = T - T_r.
    for (const double temperature : TEMPERATURES) {
        using AD = CppAD::AD<double>;
        std::vector<double> point = {temperature};
        point.insert(point.end(), parameters.temperature_law_coefficient.begin(), parameters.temperature_law_coefficient.end());
        std::vector<AD> independent(point.begin(), point.end());
        CppAD::Independent(independent);
        const auto record = epcsaft::detail::make_tape_parameters(parameters, coefficients(), {independent.begin() + 1, independent.end()});
        CppAD::ADFun<double> tape(independent, epcsaft::temperature_law_values(independent[0], record));
        std::vector<double> along_t(point.size(), 0.0), zero(point.size(), 0.0);
        along_t[0] = 1.0;
        (void)tape.Forward(0, point);
        const std::array<std::vector<double>, 3> taylor = {tape.Forward(1, along_t), tape.Forward(2, zero), tape.Forward(3, zero)};
        const auto jacobian = tape.Jacobian(point);
        const double d = temperature - T_R;
        for (std::size_t law = 0; law < 6; ++law) {
            const auto tag = "law-derivative/T=" + number(temperature) + "/" + LAW_NAMES[law];
            std::vector<double> weight(6, 0.0);
            weight[law] = 1.0;
            const auto hessian = tape.Hessian(point, weight);
            const std::array<const std::array<double, 5>*, 6> vectors = {&K01, &L01, &E00, &V00, &E11, &C01};
            const auto* a = vectors[law];
            const char* unit = law == 2 || law == 4 ? "K/K^n" : "1/K^n";
            for (int n = 1, factorial = 1; n <= 3; factorial *= ++n) {
                double analytic = (n % 2 ? -1.0 : 1.0) * factorial * (*a)[4] / std::pow(temperature, n + 1);
                for (int q = n; q <= 3; ++q) {
                    double falling = 1.0;
                    for (int k = 0; k < n; ++k) falling *= q - k;
                    analytic += falling * (*a)[static_cast<std::size_t>(q)] * std::pow(d, q - n);
                }
                normalized(tag + "/dT" + std::to_string(n), "d^n q/dT^n", unit, factorial * taylor[n - 1][law], analytic, 1.0,
                           1e-12, field("T_K", temperature), "tape thermal derivative equals the analytic law");
            }
            for (std::size_t p = 0; p < 5; ++p) {
                const std::size_t column = 1 + 5 * law + p;
                const double basis = p < 4 ? std::pow(d, static_cast<double>(p)) : 1.0 / temperature - 1.0 / T_R;
                const double thermal = p < 4 ? (p == 0 ? 0.0 : p * std::pow(d, static_cast<double>(p) - 1.0))
                                             : -1.0 / (temperature * temperature);
                const auto name = tag + "/" + (p < 4 ? "a" + std::to_string(p) : std::string("s"));
                normalized(name + "/dq", "dq/dc", "1", jacobian[law * point.size() + column], basis, 1.0, 1e-12,
                           field("T_K", temperature), "tape coefficient derivative equals its basis function");
                normalized(name + "/dTdq", "d2q/dTdc", "1/K", hessian[column], thermal, 1.0, 1e-12,
                           field("T_K", temperature), "tape mixed derivative equals the basis derivative");
            }
        }
    }

    // 3. Sites, exact actions and central differences at every temperature; pressure witnesses at 320 K.
    const auto active = coefficients();
    const epcsaft::HelmholtzModel model(parameters, active);
    const std::array<const char*, 3> outputs = {"a_res", "P", "h_res"}, units = {"1", "Pa", "J/mol"};
    const std::vector<double> t_steps = {0.4, 0.2, 0.1, 0.05, 0.025}, u_steps = {1e-3, 5e-4, 2.5e-4, 1.25e-4};
    // At T = T_r the a2/a3 mixed actions vanish and their T differences carry pure h^2 truncation, so this
    // ladder continues halving to 0.00625 K (lead decision 2026-10-16); the tolerance is unchanged.
    const std::vector<double> mixed_t_steps = {0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625};
    for (const double temperature : TEMPERATURES) {
        const auto tag = "four-site/T=" + number(temperature);
        const auto state = model.state(temperature, RHO, X);
        const auto s = epcsaft::segment_sigma(temperature, parameters), dd = epcsaft::segment_diameter(temperature, parameters, s);
        const auto residuals = epcsaft::association_residuals(RHO, X, parameters, state.site_fractions,
            epcsaft::association_delta(temperature, RHO, X, parameters, s, dd));
        double worst = 0.0;
        for (const double r : residuals) worst = std::max(worst, std::abs(r));
        recording_support::compare((tag + "/site-residual").c_str(), "max|F_a|", "1", worst, 0.0,
            field("site_fractions", state.site_fractions), "Identity", "temperature_law", "solved sites satisfy mass action",
            "absolute", 1e-12, 0.0, worst <= 1e-12, "temperature-laws");
        for (const double fraction : state.site_fractions) CHECK(fraction > 0.0 && fraction <= 1.0);
        const auto inputs = fields({field("T_K", temperature), field("rho_mol_m3", RHO), field("site_fractions", state.site_fractions)});

        const auto exact = actions(model, temperature);
        std::map<double, Actions> thermal_shift;
        const auto at_t = [&](double step) -> const Actions& {
            auto found = thermal_shift.find(step);
            if (found == thermal_shift.end()) found = thermal_shift.emplace(step, actions(model, temperature + step)).first;
            return found->second;
        };
        if (temperature == 320.0) {
            // Nonzero witnesses through the implicit site chain: a lost dependency cannot pass by returning zero.
            const double floor = 1e-6 * std::max(1.0, std::abs(exact.value[1]));
            for (const std::size_t slot : {25u, 10u, 15u, 5u}) {
                const double first = exact.coefficient[1][slot], mixed = temperature * exact.mixed[1][slot];
                const bool ok = std::abs(first) > floor && std::abs(mixed) > floor;
                const auto id = tag + "/witness/" + LAW_NAMES[slot / 5] + "/a0";
                recording_support::compare(id.c_str(), "verdict", "bool", ok ? 1.0 : 0.0, 1.0,
                    fields({inputs, field("dP_du_Pa", first), field("T_dTdu_P_Pa", mixed), field("floor_Pa", floor)}),
                    "Relationship", "temperature_law", "|dP/du| and |T d_T d_u P| exceed 1e-6 max(1 Pa, |P|)",
                    "absolute", 0.0, 0.0, ok, "temperature-laws");
            }
        }
        for (std::size_t o = 0; o < 3; ++o) {
            const auto value = [&](double t) {
                return o == 0 ? model.evaluate(t, 1.0 / RHO, X, epcsaft::DerivativeOrder::Value).value
                    : o == 1 ? model.state(t, RHO, X).pressure : model.state(t, RHO, X).residual_enthalpy;
            };
            ladder(tag + "/" + outputs[o] + "/dT", units[o], exact.thermal[o][0], t_steps,
                [&](double step) { return value(temperature + step); }, inputs);
            for (std::size_t n = 1; n < 3; ++n)
                ladder(tag + "/" + outputs[o] + "/dT" + std::to_string(n + 1), units[o], exact.thermal[o][n], t_steps,
                    [&](double step) { return at_t(step).thermal[o][n - 1]; }, inputs);
        }
        for (std::size_t slot = 0; slot < 30; ++slot) {
            std::map<double, std::array<std::array<double, 2>, 3>> rebuilt;  // [output] {value, d_T}
            const auto at_u = [&](double step) -> const auto& {
                auto found = rebuilt.find(step);
                if (found == rebuilt.end()) {
                    const epcsaft::HelmholtzModel other(shifted(parameters, slot, step));
                    const auto other_state = other.state(temperature, RHO, X);
                    const auto d = other.evaluate(temperature, 1.0 / RHO, X, epcsaft::DerivativeOrder::Hessian);
                    constexpr double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
                    const epcsaft::DensityActionDirection along_t{1.0, 0.0, {0.0, 0.0}, {}};
                    const double h_t = other.density_actions(temperature, RHO, X, {epcsaft::PropertyObservable::ResidualEnthalpy},
                        {along_t}).results.front().action.value();
                    found = rebuilt.emplace(step, std::array<std::array<double, 2>, 3>{{
                        {d.value, d.gradient[0]},
                        {other_state.pressure, R * RHO - R * d.gradient[1] - R * temperature * d.hessian[1]},
                        {other_state.residual_enthalpy, h_t}}}).first;
                }
                return found->second;
            };
            const auto name = tag + "/" + LAW_NAMES[slot / 5] + "/" + (slot % 5 == 4 ? "s" : "a" + std::to_string(slot % 5));
            for (std::size_t o = 0; o < 3; ++o) {
                ladder(name + "/" + outputs[o] + "/du", units[o], exact.coefficient[o][slot], u_steps,
                    [&](double step) { return at_u(step)[o][0]; }, inputs);
                ladder(name + "/" + outputs[o] + "/dTdu/over-u", units[o], exact.mixed[o][slot], u_steps,
                    [&](double step) { return at_u(step)[o][1]; }, inputs);
                ladder(name + "/" + outputs[o] + "/dTdu/over-T", units[o], exact.mixed[o][slot], mixed_t_steps,
                    [&](double step) { return at_t(step).coefficient[o][slot]; }, inputs);
            }
        }
    }

    // The fit's density columns: P* = P(T, 100 mol/m3, x), vapor branch, observed 101, Difference scale 100.
    // r(u=0) = -0.01 when the closure reselects the same root; dr/du against rebuilt-model differences.
    for (const double temperature : TEMPERATURES) {
        using namespace epcsaft::regression;
        Observation observation;
        observation.kind = ObservationKind::Density;
        observation.temperature = temperature;
        observation.pressure = model.state(temperature, RHO, X).pressure;
        observation.composition = X;
        observation.branch = epcsaft::DensityBranch::Vapor;
        observation.observed = 101.0;
        observation.scale = 100.0;
        const auto prediction = predict(model, observation);
        CHECK(prediction.status == ObservationStatus::Available && prediction.state->pressure_density_derivative > 0.0);
        const auto tag = "density-column/T=" + number(temperature);
        const auto inputs = fields({field("P_star_Pa", observation.pressure), field("rho_mol_m3", *prediction.value),
                                    field("dP_drho", prediction.state->pressure_density_derivative)});
        recording_support::compare((tag + "/residual").c_str(), "r", "1", *prediction.residual, -0.01, inputs, "Identity",
            "temperature_law", "the pressure closure reselects the vapor root at 100 mol/m3", "absolute", 1e-10, 0.0,
            std::abs(*prediction.residual + 0.01) <= 1e-10, "temperature-laws");
        const auto gradient = differentiate(model, observation, prediction);
        CHECK(gradient.status == ObservationStatus::Available && gradient.value.size() == 30);
        for (std::size_t slot = 0; slot < 30; ++slot)
            ladder(tag + "/" + LAW_NAMES[slot / 5] + "/" + std::to_string(slot % 5) + "/dr_du", "1",
                   law_scale(slot) * gradient.value[slot] / observation.scale, u_steps, [&](double step) {
                       const auto other = predict(epcsaft::HelmholtzModel(shifted(parameters, slot, step)), observation);
                       if (other.status != ObservationStatus::Available) throw std::runtime_error("rebuilt density prediction failed");
                       return *other.residual;
                   }, inputs);
    }

    // 4. A law crossing to negative volume on a replayed model is refused with its target and temperature.
    const auto refused = [&](const std::string& id, const std::function<void()>& call, const std::string& cause) {
        std::string message;
        try { call(); } catch (const std::invalid_argument& error) { message = error.what(); }
        const bool ok = message.find(cause) != std::string::npos;
        recording_support::compare(id.c_str(), "verdict", "bool", ok ? 1.0 : 0.0, 1.0, field("message", message),
            "Relationship", "temperature_law", "the named refusal is preserved", "absolute", 0.0, 0.0, ok, "temperature-laws");
    };
    {
        const epcsaft::HelmholtzModel crossing(four_site({0.9, 0.01, NAN_TERM, NAN_TERM, NAN_TERM}));
        (void)crossing.state(300.0, RHO, X);
        refused("domain/negative-volume", [&] { (void)crossing.state(320.0, RHO, X); },
                "invalid evaluated association_volume for pair 0/3 at T = 320");
    }
    {   // l01 = 0.5 + 0.01 (T - 300 K) reaches 1.1 at 360 K on a pair without sites: no double site solve precedes the replay.
        auto document = binary_document(false);
        document.entries.push_back(entry("k01", "k_ij", {0, 1}, 0.0, "dimensionless"));
        add_law(document.entries, "l01", "l_ij", 0, 0, 1, T_R, {0.5, 0.01, NAN_TERM, NAN_TERM, NAN_TERM});
        const epcsaft::HelmholtzModel crossing(epcsaft::build_parameters(document));
        const auto h = epcsaft::PropertyObservable::ResidualEnthalpy;
        (void)crossing.density_actions(300.0, RHO, X, {h}, {{1.0, 0.0, {0.0, 0.0}, {}}});
        refused("domain/l-crossing/density-actions", [&] { (void)crossing.density_actions(360.0, RHO, X, {h}, {{1.0, 0.0, {0.0, 0.0}, {}}}); },
                "invalid evaluated l_ij for pair 0/1 at T = 360");
        refused("domain/l-crossing/pressure-actions", [&] { (void)crossing.pressure_actions(360.0, 1.0e5, X,
                    epcsaft::DensityBranch::Vapor, {h}, {{1.0, 0.0, {0.0, 0.0}, {}}}, std::nullopt); },
                "invalid evaluated l_ij for pair 0/1 at T = 360");
    }
    // 5. Law-owned constants are unconsumed; one of value or law per quantity; declared units.
    refused("refuse/kij-under-law", [&] { (void)epcsaft::HelmholtzModel(parameters, {{ParameterFamily::Kij, {1, 2}}}); },
            "unconsumed law-owned coordinate");
    refused("refuse/rule-volume", [&] { (void)epcsaft::HelmholtzModel(parameters, {{ParameterFamily::AssociationVolume, {3, 12}}}); },
            "unconsumed law-owned coordinate");
    const auto document_with = [](const std::function<void(std::vector<ParameterEntry>&)>& edit) {
        auto document = binary_document(false);
        edit(document.entries);
        return [document] { (void)epcsaft::build_parameters(document); };
    };
    refused("refuse/value-and-law", document_with([](auto& e) {
        add_law(e, "k", "k_ij", 0, 0, 1, T_R, {0.01, NAN_TERM, NAN_TERM, NAN_TERM, NAN_TERM});
        e.push_back(entry("k-value", "k_ij", {1, 0}, 0.01, "dimensionless")); }), "duplicate");
    refused("refuse/coefficient-unit", document_with([](auto& e) {
        add_law(e, "k", "k_ij", 0, 0, 1, T_R, {0.01, NAN_TERM, NAN_TERM, NAN_TERM, NAN_TERM});
        e.push_back(entry("k/a1", "temperature_law_coefficient", {0, 1}, 1e-4, "dimensionless")); }), "expected 1 / kelvin");
    refused("refuse/shared-coefficient-identity", document_with([](auto& e) {
        add_law(e, "shared", "k_ij", 0, 0, 1, T_R, {0.03, NAN_TERM, NAN_TERM, NAN_TERM, NAN_TERM});
        e.push_back(entry("l", "l_ij_temperature_law", {1, 0, 1}, T_R, "kelvin"));
        e.push_back(entry("shared/a0", "temperature_law_coefficient", {1, 0}, 0.01, "dimensionless")); }),
        "repeated temperature law coefficient shared/a0");
}

}  // namespace properties_checks
