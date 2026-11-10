// Hong et al. 2027 (Sep. Purif. Technol. 422, 140522) pathway-weighted reactive cross-association on the non-ionic
// path, falsifiers F1-F4 of its design at states L, V and Z: pathway_weighted_association_helmholtz,
// symmetric_reactive_site_correction, pathway_weight and pathway_weighted_cross_dispersion. Hong Tables 1-3 values;
// every linear law has T_r = 298.15 K. Test topology: each reactive CO2 site bonds MEA site b only. That is a
// falsifier declaration, not a chemistry choice. Differences verify exact derivatives only.
#include "equilibrium_support.hpp"
#include "properties_checks_shared.hpp"
#include "surface/parameter_builder.hpp"

#include <functional>

namespace properties_checks {
namespace {

using epcsaft::ActiveParameter;
using epcsaft::DensityBranch;
using epcsaft::ParameterEntry;
using epcsaft::ParameterFamily;
using epcsaft::PropertyObservable;
using recording_support::field;
using recording_support::fields;

constexpr double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
constexpr double T_R = 298.15;
constexpr double IDENTITY = 1e-10;    // F1, F2, F2b, F3c, F4(i): scaled identity tolerance
constexpr double DIFFERENCE = 1e-7;   // F3: best of three scaled difference errors

// Single, Symmetric (nu = 2, ordinary association) and Corrected (nu = 2 with omega_2) are one block each;
// Blend weighs a corrected symmetric block by w and a single-site block by 1 - w; Uncorrected drops the correction.
enum class Form { Single, Symmetric, Corrected, Blend, Uncorrected };

struct Pathway final {
    double k_slope, k_intercept, kappa, c, d;  // k = slope T + intercept (Eq. 12); eps/k = c T + d (Eq. 13)
    std::size_t nu;
};
constexpr Pathway SYMMETRIC{3.675e-3, -2.3981, 1.555e-3, 4.1160, 3546.55, 2};
constexpr Pathway SINGLE{5.939e-3, -1.9778, 6.5e-5, 12.5398, 1197.42, 1};

struct Weight final {
    double intercept, loading, temperature;  // b_0, b_alpha, b_T
};
constexpr Weight LOGISTIC{1.7, -10.0, -1.5};  // illustrative, shaped on Table S1, not fitted

ParameterEntry entry(std::string identity, std::string field, std::vector<std::size_t> index, double magnitude,
                     std::string unit) {
    return {std::move(identity), std::move(field), std::move(index), magnitude, std::move(unit)};
}

// q = a_0 + a_1 (T - T_r) for q = slope T + intercept.
void linear_law(std::vector<ParameterEntry>& e, const std::string& field, std::size_t law, std::size_t a, std::size_t b,
                double slope, double intercept) {
    const bool energy = field == "association_energy_over_k";
    const std::string name = "law/" + std::to_string(law);
    e.push_back(entry(name, field + "_temperature_law", {law, a, b}, T_R, "kelvin"));
    e.push_back(entry(name + "/a0", "temperature_law_coefficient", {law, 0}, intercept + slope * T_R, energy ? "kelvin" : "dimensionless"));
    e.push_back(entry(name + "/a1", "temperature_law_coefficient", {law, 1}, slope, energy ? "dimensionless" : "1 / kelvin"));
}

// Components MEA (2B), water (4C), CO2; per block sites MEA a, b, water a, b, CO2 r. order[p] is the component
// declared at position p (the same model in another component order).
using Order = std::array<std::size_t, 3>;
epcsaft::Parameters hong(Form form, Weight weight = LOGISTIC, Order order = {0, 1, 2}) {
    Order at{};
    for (std::size_t p = 0; p < 3; ++p) at[order[p]] = p;
    epcsaft::ParameterDocument document;
    for (const auto component : order) document.component_ids.push_back(std::array{"mea", "water", "co2"}[component]);
    document.base_choice = "pc-saft";
    document.association_choice = "general-site";
    document.electrolyte_choice = "none";
    document.permittivity_choice = "none";
    auto& e = document.entries;
    const std::array<double, 3> m = {3.1200, 2.1945, 2.0729};
    const std::array<double, 3> sigma = {3.0229, 2.2290, 2.7852};
    const std::array<double, 3> epsilon = {287.94, 141.66, 169.21};
    const std::array<double, 3> mass = {0.06108352, 0.01801528, 0.0440095};
    for (std::size_t i = 0; i < 3; ++i) {
        e.push_back(entry("", "molar_mass", {at[i]}, mass[i], "kilogram / mole"));
        e.push_back(entry("", "charge_number", {at[i]}, 0.0, "elementary-charge"));
        e.push_back(entry("", "segment_count", {at[i]}, m[i], "dimensionless"));
        e.push_back(entry("", "sigma_constant", {at[i]}, sigma[i], "angstrom"));
        e.push_back(entry("", "dispersion_energy_over_k", {at[i]}, epsilon[i], "kelvin"));
    }
    const bool blend = form == Form::Blend || form == Form::Uncorrected;
    const std::vector<Pathway> blocks = blend ? std::vector<Pathway>{SYMMETRIC, SINGLE}
        : std::vector<Pathway>{form == Form::Single ? SINGLE : SYMMETRIC};
    linear_law(e, "k_ij", 0, at[0], at[1], 9.0e-5, -0.0603);  // Table 2, x_MEA <= 0.20 set
    linear_law(e, "k_ij", 1, at[1], at[2], 2.46e-4, -0.0633);
    linear_law(e, "k_ij", 2, at[0], at[2], blocks[0].k_slope, blocks[0].k_intercept);
    std::size_t law = 3;
    if (blend) linear_law(e, "k_ij_complement", law++, at[0], at[2], SINGLE.k_slope, SINGLE.k_intercept);
    // Eq. 6: eps_MW = alpha_HB sqrt(eps_MEA eps_W), one power-0 law on every MEA-water edge; Eq. 5 kappa is fixed.
    const std::size_t cross = law++;
    e.push_back(entry("law/cross/a0", "temperature_law_coefficient", {cross, 0}, 1.0374 * std::sqrt(2355.03 * 1804.17), "kelvin"));
    const double kappa_cross = std::pow((std::cbrt(0.0313) + std::cbrt(0.2039)) / 2.0, 3.0);
    for (std::size_t block = 0; block < blocks.size(); ++block) {
        const std::size_t s = 5 * block;
        for (const std::size_t component : {0, 0, 1, 1, 2}) document.site_component.push_back(static_cast<int>(at[component]));
        for (const std::size_t nu : {std::size_t{1}, std::size_t{1}, std::size_t{2}, std::size_t{2}, blocks[block].nu})
            document.site_multiplicity.push_back(static_cast<int>(nu));
        e.push_back(entry("", "association_energy_over_k", {s, s + 1}, 2355.03, "kelvin"));
        e.push_back(entry("", "association_volume", {s, s + 1}, 0.0313, "dimensionless"));
        e.push_back(entry("", "association_energy_over_k", {s + 2, s + 3}, 1804.17, "kelvin"));
        e.push_back(entry("", "association_volume", {s + 2, s + 3}, 0.2039, "dimensionless"));
        for (const auto& [a, b] : {std::pair{s, s + 3}, std::pair{s + 1, s + 2}}) {
            e.push_back(entry("law/cross", "association_energy_over_k_temperature_law", {cross, a, b}, T_R, "kelvin"));
            e.push_back(entry("", "association_volume", {a, b}, kappa_cross, "dimensionless"));
        }
        linear_law(e, "association_energy_over_k", law++, s + 4, s + 1, blocks[block].c, blocks[block].d);
        e.push_back(entry("", "association_volume", {s + 4, s + 1}, blocks[block].kappa, "dimensionless"));
        for (std::size_t site = s; blend && site < s + 5; ++site)
            e.push_back(entry("", "site_pathway", {site}, static_cast<double>(block + 1), "dimensionless"));
    }
    if (form == Form::Corrected || form == Form::Blend) e.push_back(entry("", "site_correction", {4}, 1.0, "dimensionless"));
    e.push_back(entry("", "no_cross_association", {at[1], at[2]}, 0.0, "dimensionless"));
    if (blend) {
        e.push_back(entry("", "pathway_loading", {at[2], at[0]}, 0.0, "dimensionless"));
        e.push_back(entry("b0", "pathway_weight", {0}, weight.intercept, "dimensionless"));
        e.push_back(entry("b_alpha", "pathway_weight", {1}, weight.loading, "dimensionless"));
        e.push_back(entry("b_T", "pathway_weight", {2}, weight.temperature, "dimensionless"));
    }
    return epcsaft::build_parameters(document);
}

// The 16 blend coordinates: k_MW a0, a1; k_sym a0, a1; k_single a0, a1; eps_MW; eps_sym a0, a1; eps_single a0, a1;
// kappa_sym; kappa_single; b_0, b_alpha, b_T. A single-block model has k_MW, k_MEA-CO2, eps_MW, eps_r and kappa_r.
std::vector<ActiveParameter> coordinates(Form form, bool weights = true) {
    std::vector<ActiveParameter> result;
    const bool blend = form == Form::Blend || form == Form::Uncorrected;
    for (const std::size_t slot : blend ? std::vector<std::size_t>{0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12}
                                        : std::vector<std::size_t>{0, 1, 4, 5, 6, 7, 8})
        result.push_back({ParameterFamily::TemperatureLawCoefficient, {slot}});
    if (!blend) {
        result.push_back({ParameterFamily::AssociationVolume, {9, 21}});
        return result;
    }
    result.push_back({ParameterFamily::AssociationVolume, {14, 41}});
    result.push_back({ParameterFamily::AssociationVolume, {69, 96}});
    for (std::size_t k = 0; weights && k < 3; ++k) result.push_back({ParameterFamily::PathwayWeight, {k}});
    return result;
}
// Blend coordinate of each single-block coordinate, for the w = 1 (symmetric) and w = 0 (single) limits.
constexpr std::array<std::size_t, 8> SYMMETRIC_SHARED = {0, 1, 2, 3, 6, 7, 8, 11};
constexpr std::array<std::size_t, 8> SINGLE_SHARED = {0, 1, 4, 5, 6, 9, 10, 12};

struct StatePoint final {
    const char* name;
    double temperature, density;
    std::vector<double> x;
    DensityBranch branch;
};
// The six component orders of the same model, for the double rounding estimate e_A.
constexpr std::array<Order, 6> ORDERS{{{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}};
const std::vector<StatePoint> STATES = {
    {"L", 333.15, 42424.79, {0.10404183733215595, 0.8231288765353348, 0.07282928613250916}, DensityBranch::Liquid},
    {"V", 353.15, 34.49684, {1.0e-5, 0.6, 0.39999}, DensityBranch::Vapor},
    {"Z", 333.15, 43629.45, {0.11221432663480933, 0.8877856733651908, 0.0}, DensityBranch::Liquid}};

std::string inputs(const StatePoint& s, const char* model) {
    return fields({field("state", s.name), field("model", model), field("T_K", s.temperature),
        field("rho_mol_m3", s.density), field("x", s.x)});
}

// |model - reference| <= tolerance max(|reference|, scale).
void identity(const std::string& id, const char* quantity, double model, double reference, double scale,
              double tolerance, const std::string& context, const char* claim) {
    const double allowance = tolerance * std::max(std::abs(reference), scale);
    recording_support::compare(id.c_str(), quantity, "scaled", model, reference, context, "Identity", "pathway_weight",
        claim, "absolute", allowance, 0.0, std::abs(model - reference) <= allowance, "pathway-weight", id.c_str());
}

// Named state quantities and their scales (S = 1 reduced, RT for J/mol, R for J/(mol K)).
std::vector<std::tuple<std::string, double, double>> quantities(const epcsaft::State& state, const StatePoint& s,
    std::size_t block, bool solvent_only = false) {
    const double rt = R * s.temperature;
    std::vector<std::tuple<std::string, double, double>> out = {
        {"a_res", state.residual_helmholtz, 1.0}, {"Z", state.compressibility_factor, 1.0},
        {"u_res", state.residual_internal_energy, rt}, {"h_res", state.residual_enthalpy, rt},
        {"cv_res", state.residual_isochoric_heat_capacity, R}, {"cp_res", state.residual_isobaric_heat_capacity.value(), R}};
    for (std::size_t i = 0; i < (solvent_only ? 2 : 3); ++i)
        out.emplace_back("mu_res/RT[" + std::to_string(i) + "]", state.residual_chemical_potential_over_rt[i].value(), 1.0);
    for (std::size_t site = 0; !solvent_only && site < 5; ++site)
        out.emplace_back("X[" + std::to_string(site) + "]", state.site_fractions[5 * block + site], 1.0);
    return out;
}

// F1 and F2: a saturated-weight blend equals the matching single-block model at each state, including B and the
// fixed T, P, x log-fugacity sensitivities of the shared coordinates (scaled by |theta|).
void limit(const char* name, Form standalone, Form blend, double intercept, std::size_t block,
           const std::array<std::size_t, 8>& shared, const char* claim) {
    const epcsaft::HelmholtzModel one(hong(standalone), coordinates(standalone));
    const epcsaft::HelmholtzModel two(hong(blend, {intercept, 0.0, 0.0}), coordinates(blend, false));
    for (const auto& s : STATES) {
        const auto a = two.state(s.temperature, s.density, s.x, true), b = one.state(s.temperature, s.density, s.x, true);
        const auto qa = quantities(a, s, block), qb = quantities(b, s, 0);
        for (std::size_t q = 0; q < qa.size(); ++q)
            identity(std::string(name) + "/" + s.name + "/" + std::get<0>(qa[q]), std::get<0>(qa[q]).c_str(),
                std::get<1>(qa[q]), std::get<1>(qb[q]), std::get<2>(qb[q]), IDENTITY, inputs(s, name), claim);
        const double second = two.virial_coefficients(s.temperature, s.x).second;
        const double reference = one.virial_coefficients(s.temperature, s.x).second;
        identity(std::string(name) + "/" + s.name + "/B", "B", second, reference, 0.0, IDENTITY, inputs(s, name), claim);
        for (std::size_t c = 0; c < shared.size(); ++c) {
            const auto coordinate = coordinates(blend, false)[shared[c]];
            const double theta = std::abs(epcsaft::detail::family_values(two.parameters(), coordinate.family).first[coordinate.indices[0]]);
            for (std::size_t i = 0; i < 3; ++i)
                identity(std::string(name) + "/" + s.name + "/dlnphi" + std::to_string(i) + "/theta" + std::to_string(shared[c]),
                    "theta dlnphi/dtheta", theta * a.log_fugacity_parameter_derivative[i * 13 + shared[c]].value(),
                    theta * b.log_fugacity_parameter_derivative[i * 8 + c].value(), 1.0, IDENTITY, inputs(s, name), claim);
        }
    }
}

// Best of three scaled errors |D_k - A| / max(|A|, scale) at most 1e-7, with every error recorded. A step that
// would take an amount negative uses the forward second-order difference.
double ladder(const std::string& id, double exact, double scale, const std::array<double, 3>& steps,
              const std::function<double(double)>& f, const std::function<bool(double)>& forward, const std::string& context,
              double limit = DIFFERENCE) {
    const double centre = f(0.0);
    std::vector<double> errors, differences;
    for (const double h : steps) {
        const double d = forward(h) ? (-3.0 * centre + 4.0 * f(h) - f(2.0 * h)) / (2.0 * h) : (f(h) - f(-h)) / (2.0 * h);
        differences.push_back(d);
        errors.push_back(std::abs(d - exact) / std::max(std::abs(exact), scale));
    }
    const double best = *std::min_element(errors.begin(), errors.end());
    if (limit > 0.0)
        recording_support::compare(id.c_str(), "scaled difference error", "1", best, 0.0,
            fields({context, field("exact", exact), field("steps", std::vector<double>(steps.begin(), steps.end())),
                field("differences", differences), field("errors", errors)}),
            "Relationship", "three-step differences", "exact derivative equals refined differences (verification only)",
            "absolute", limit, 0.0, best <= limit, "pathway-weight", id.c_str());
    return best;
}

epcsaft::Parameters shifted(epcsaft::Parameters p, const ActiveParameter& coordinate, double step) {
    auto values = epcsaft::detail::family_values(p, coordinate.family);
    for (const auto index : coordinate.indices) values.first[index] += step;
    return p;
}

// d(N a)/dn_i with every site fraction frozen at its root: the error an omitted implicit site chain would leave.
double frozen_chemical_potential(const epcsaft::Parameters& p, const StatePoint& s, const std::vector<double>& sites, std::size_t i) {
    using Series = epcsaft::Taylor<1>;
    std::vector<Series> amounts(s.x.begin(), s.x.end()), x;
    amounts[i].c[1] = 1.0;
    Series total(0.0);
    for (const auto& n : amounts) total += n;
    for (const auto& n : amounts) x.push_back(n / total);
    const Series density = total * s.density, temperature(s.temperature);
    const auto sigma = epcsaft::segment_sigma(temperature, p);
    const auto diameter = epcsaft::segment_diameter(temperature, p, sigma);
    const std::vector<Series> site(sites.begin(), sites.end());
    const std::vector<Series> energy(p.association_energy_over_k.begin(), p.association_energy_over_k.end());
    const std::vector<Series> volume(p.association_volume.begin(), p.association_volume.end());
    const auto delta = epcsaft::association_delta_values(temperature, density, x, p, sigma, diameter, energy, volume);
    const auto terms = epcsaft::residual_helmholtz(temperature, density, x, p, sigma, diameter, site, Series(0.0), Series(0.0),
        epcsaft::association_products(density, x, p, site, delta));
    return (total * (terms.hard_chain + terms.dispersion + terms.association)).c[1];
}

// The z of pathway_weight in its own operation order, so the fixed-weight model has the logistic's w bit for bit.
double logistic_argument(const StatePoint& s) {
    double z = LOGISTIC.intercept + LOGISTIC.temperature * (s.temperature - 353.15) / 40.0;
    return z + LOGISTIC.loading * (s.x[2] / s.x[0] - 0.5);
}

void f3_state(const epcsaft::HelmholtzModel& model, const StatePoint& s) {
    const double rt = R * s.temperature, volume = 1.0 / s.density;
    const auto base = model.state(s.temperature, s.density, s.x, true);
    const std::string context = inputs(s, "blend");
    const std::array<double, 3> relative = {1e-3, 1e-4, 1e-5}, kelvin = {0.1, 0.01, 0.001};
    const std::array<double, 3> amount = {5e-5, 5e-6, 5e-7};
    const auto never = [](double) { return false; };
    const auto value = [&](double t, double v, const std::vector<double>& n) {
        return model.evaluate(t, v, n, epcsaft::DerivativeOrder::Value).value;
    };
    const std::string id = std::string("F3/") + s.name;
    ladder(id + "/Z-1", base.compressibility_factor - 1.0, 1.0, relative,
        [&](double h) { return value(s.temperature, volume / (1.0 + h), s.x); }, never, context);
    ladder(id + "/da_dT", -base.residual_internal_energy / (rt * s.temperature), 0.01, kelvin,
        [&](double h) { return value(s.temperature + h, volume, s.x); }, never, context);
    ladder(id + "/cv_res", base.residual_isochoric_heat_capacity, rt / 100.0, kelvin,
        [&](double h) { return model.state(s.temperature + h, s.density, s.x).residual_internal_energy; }, never, context);
    const auto hessian = model.evaluate(s.temperature, volume, s.x, epcsaft::DerivativeOrder::Hessian);
    const std::size_t width = hessian.gradient.size();
    for (std::size_t i = 0; i < 3; ++i) {
        const auto along = [&, i](double h) { auto n = s.x; n[i] += h; return n; };
        const auto forward = [&, i](double h) { return s.x[i] - h < 0.0; };
        ladder(id + "/mu" + std::to_string(i), base.residual_chemical_potential_over_rt[i].value(), 1.0, amount,
            [&](double h) { return value(s.temperature, volume, along(h)); }, forward, context);
        for (std::size_t j = 0; j < 3; ++j)
            ladder(id + "/dmu" + std::to_string(i) + "_dn" + std::to_string(j), hessian.hessian[(2 + i) * width + 2 + j], 1.0, amount,
                [&, i, j](double h) {
                    auto n = s.x; n[j] += h;
                    return model.evaluate(s.temperature, volume, n, epcsaft::DerivativeOrder::Gradient).gradient[2 + i];
                }, [&, j](double h) { return s.x[j] - h < 0.0; }, context);
    }
    // c_p,res from h_res at fixed P through the native difference owner.
    std::vector<double> errors;
    for (const double step : kelvin) {
        const auto check = model.pressure_action_difference(s.temperature, base.pressure, s.x, s.branch,
            {PropertyObservable::ResidualEnthalpy}, {1.0, 0.0, {0.0, 0.0, 0.0}, {}}, step, 1e-8, 10.0).front();
        errors.push_back(std::abs(check.difference - check.exact) / std::max(std::abs(check.exact), rt / 100.0));
    }
    const double best = *std::min_element(errors.begin(), errors.end());
    recording_support::compare((id + "/cp_res").c_str(), "scaled difference error", "1", best, 0.0,
        fields({context, field("errors", errors)}), "Relationship", "pressure_action_difference",
        "c_p,res equals differences of h_res at fixed P (verification only)", "absolute", DIFFERENCE, 0.0,
        best <= DIFFERENCE, "pathway-weight", (id + "/cp_res").c_str());
    // Each coordinate: d ln phi_i/d theta at fixed T, P, x against shifted models solved at the same P.
    const auto active = coordinates(Form::Blend);
    for (std::size_t c = 0; c < active.size(); ++c) {
        const double theta = std::abs(epcsaft::detail::family_values(model.parameters(), active[c].family).first[active[c].indices[0]]);
        std::map<double, std::vector<double>> cache;
        const auto lnphi = [&](double h) -> const std::vector<double>& {
            auto found = cache.find(h);
            if (found != cache.end()) return found->second;
            const epcsaft::HelmholtzModel moved(shifted(model.parameters(), active[c], h));
            const auto state = moved.state_at_pressure(s.temperature, base.pressure, s.x, s.branch, s.density);
            std::vector<double> values;
            for (const auto& v : state.log_fugacity_coefficient.value()) values.push_back(v.value_or(0.0));
            return cache.emplace(h, std::move(values)).first->second;
        };
        const std::array<double, 3> steps = {1e-3 * theta, 1e-4 * theta, 1e-5 * theta};
        for (std::size_t i = 0; i < 3; ++i) {
            if (s.x[i] == 0.0) continue;  // ln phi of an absent component is not compared
            ladder(id + "/dlnphi" + std::to_string(i) + "/theta" + std::to_string(c),
                base.log_fugacity_parameter_derivative[i * active.size() + c].value(), 1.0 / theta, steps,
                [&, i](double h) { return lnphi(h)[i]; }, never, context);
        }
    }
}

}  // namespace

void pathway_weight() {
    // F1: w = 0 against the single-site model; w = 1 against the standalone corrected symmetric model.
    limit("F1/w0", Form::Single, Form::Blend, -800.0, 1, SINGLE_SHARED, "the w = 0 blend equals the single-site model");
    limit("F1/w1", Form::Corrected, Form::Blend, 40.0, 0, SYMMETRIC_SHARED, "the w = 1 blend equals the corrected symmetric model");
    // F2: correction off, w = 1, against ordinary two-site association.
    limit("F2", Form::Symmetric, Form::Uncorrected, 40.0, 0, SYMMETRIC_SHARED, "the uncorrected w = 1 blend equals ordinary association");

    const epcsaft::HelmholtzModel blend(hong(Form::Blend), coordinates(Form::Blend));
    const epcsaft::HelmholtzModel uncorrected(hong(Form::Uncorrected));
    const epcsaft::HelmholtzModel single(hong(Form::Single));
    const auto& L = STATES[0];
    const auto& Z = STATES[2];

    // F2b: on the carbon-free face dA/dw = 0, so the logistic blend equals the standalone values (sym = single there);
    // the correction moves the infinite-dilution CO2 potential by -w X nu (ln X - X/2 + 1/2) = 0.1188144166.
    {
        const auto a = blend.state(Z.temperature, Z.density, Z.x, true), b = single.state(Z.temperature, Z.density, Z.x, true);
        const auto qa = quantities(a, Z, 0, true), qb = quantities(b, Z, 0, true);
        for (std::size_t q = 0; q < qa.size(); ++q)
            identity("F2b/Z/" + std::get<0>(qa[q]), std::get<0>(qa[q]).c_str(), std::get<1>(qa[q]), std::get<1>(qb[q]),
                std::get<2>(qb[q]), IDENTITY, inputs(Z, "blend"), "the carbon-free blend equals the standalone model");
        const auto off = uncorrected.state(Z.temperature, Z.density, Z.x);
        identity("F2b/Z/correction", "mu_res,CO2/RT on - off", a.residual_chemical_potential_over_rt[2].value()
            - off.residual_chemical_potential_over_rt[2].value(), 0.1188144166, 1.0, IDENTITY, inputs(Z, "blend"),
            "the correction shifts the infinite-dilution CO2 potential by -w X nu (ln X - X/2 + 1/2)");
    }

    // F3 at L, V and Z, then the failures it must catch at L.
    for (const auto& s : STATES) f3_state(blend, s);
    {
        const double volume = 1.0 / L.density;
        const auto base = blend.state(L.temperature, L.density, L.x);
        const double z = logistic_argument(L);
        const epcsaft::HelmholtzModel fixed(hong(Form::Blend, {z, 0.0, 0.0}), {{ParameterFamily::PathwayWeight, {0}}});
        const auto held = fixed.state(L.temperature, L.density, L.x);
        const auto value = [&](double t, const std::vector<double>& n) {
            return blend.evaluate(t, volume, n, epcsaft::DerivativeOrder::Value).value;
        };
        const auto never = [](double) { return false; };
        for (std::size_t i = 0; i < 3; ++i) {
            const auto difference = [&, i](double h) { auto n = L.x; n[i] += h; return value(L.temperature, n); };
            const double frozen = frozen_chemical_potential(blend.parameters(), L, base.site_fractions, i);
            const double omitted = held.residual_chemical_potential_over_rt[i].value();
            const std::string id = "F3/L/caught/mu" + std::to_string(i);
            const double frozen_error = ladder(id, frozen, 1.0, {5e-5, 5e-6, 5e-7}, difference, never, "", 0.0);
            const double omitted_error = ladder(id, omitted, 1.0, {5e-5, 5e-6, 5e-7}, difference, never, "", 0.0);
            const bool caught = frozen_error > 1e3 * DIFFERENCE && (i == 1 || omitted_error > 1e3 * DIFFERENCE);
            recording_support::compare_custom(id.c_str(), "scaled difference error", "1", std::min(frozen_error, omitted_error), 0.0,
                fields({inputs(L, "blend"), field("frozen_site_mu", frozen), field("frozen_site_error", frozen_error),
                    field("exact_mu", base.residual_chemical_potential_over_rt[i].value()),
                    field("omitted_dw_dn_mu", omitted), field("omitted_dw_dn_error", omitted_error)}),
                "Relationship", "three-step differences", "frozen sites and an omitted dw/dn each fail the F3 amount check",
                "each error > 1e3 x 1e-7 (water has dw/dn = 0)", caught);
        }
        const double omitted = -held.residual_internal_energy / (R * L.temperature * L.temperature);
        const double error = ladder("F3/L/caught/da_dT", omitted, 0.01, {0.1, 0.01, 0.001},
            [&](double h) { return value(L.temperature + h, L.x); }, never, "", 0.0);
        recording_support::compare_custom("F3/L/caught/da_dT", "scaled difference error", "1", error, 0.0,
            fields({inputs(L, "blend"), field("omitted_dw_dT_u_res", held.residual_internal_energy),
                field("exact_u_res", base.residual_internal_energy)}), "Relationship", "three-step differences",
            "an omitted dw/dT fails the F3 temperature check", "error > 1e3 x 1e-7", error > 1e3 * DIFFERENCE);

        // F3c: the logistic and the fixed weight differ only by the weight chains, both Engine outputs.
        const auto gradient = fixed.evaluate(L.temperature, volume, L.x, epcsaft::DerivativeOrder::Gradient).gradient;
        const double db0 = gradient[2 + 3];
        const std::array<double, 3> dalpha = {-L.x[2] / (L.x[0] * L.x[0]), 0.0, 1.0 / L.x[0]};
        for (std::size_t i = 0; i < 3; ++i)
            identity("F3c/L/mu" + std::to_string(i), "mu_b - mu_a", base.residual_chemical_potential_over_rt[i].value()
                - held.residual_chemical_potential_over_rt[i].value(), db0 * LOGISTIC.loading * dalpha[i], 1.0, IDENTITY,
                inputs(L, "blend"), "the chemical potentials differ by dA/db0 b_alpha dalpha/dn_i");
        identity("F3c/L/u_res", "u_b - u_a", base.residual_internal_energy - held.residual_internal_energy,
            -R * L.temperature * L.temperature * db0 * LOGISTIC.temperature / 40.0, R * L.temperature, IDENTITY,
            inputs(L, "blend"), "the internal energies differ by -R T^2 dA/db0 b_T / 40 K");

        // Third order, mixed: d c_p,res / d b_T at fixed T, P, x.
        const std::size_t bt = 15;
        const epcsaft::PressureActionDirection along_t{1.0, 0.0, {0.0, 0.0, 0.0}, {}};
        const epcsaft::PressureActionDirection along_b{0.0, 0.0, {0.0, 0.0, 0.0}, {{ParameterFamily::PathwayWeight, {2}, 1.0}}};
        const double exact = blend.pressure_actions(L.temperature, base.pressure, L.x, L.branch,
            {PropertyObservable::ResidualEnthalpy}, {along_t, along_b}, L.density).results.front().action.value();
        const double theta = std::abs(LOGISTIC.temperature);
        ladder("F3/L/dcp_db_T", exact, R / theta, {1e-3 * theta, 1e-4 * theta, 1e-5 * theta}, [&](double h) {
            const epcsaft::HelmholtzModel moved(shifted(blend.parameters(), coordinates(Form::Blend)[bt], h));
            return moved.state_at_pressure(L.temperature, base.pressure, L.x, L.branch, L.density, true).residual_isobaric_heat_capacity.value();
        }, never, inputs(L, "blend"));
    }

    // F3b: density series against the tape, the round trip, the double floor against wider replay, and vacuum B.
    for (const auto& s : STATES) {
        const auto state = blend.state(s.temperature, s.density, s.x);
        const auto back = blend.state_at_pressure(s.temperature, state.pressure, s.x, s.branch);
        identity(std::string("F3b/") + s.name + "/round-trip", "rho_back/rho", back.molar_density / s.density, 1.0, 1.0, 1e-12,
            inputs(s, "blend"), "the density series returns the tape state's density at its pressure");
        double floor = 0.0;
        const double volume = 1.0 / s.density;
        const auto series = epcsaft::properties::density_volume_taylor(blend, s.temperature, s.density, s.x, floor);
        const auto tape = blend.evaluate(s.temperature, volume, s.x, epcsaft::DerivativeOrder::Hessian);
        std::vector<double> direction(5 + 16, 0.0);
        direction[1] = 1.0;
        const double third = blend.residual_energy_action(s.temperature, volume, s.x, {direction, direction, direction}).action.value();
        const std::size_t width = tape.gradient.size();
        const std::array<double, 3> exact = {volume * tape.gradient[1], volume * volume * tape.hessian[width + 1], volume * volume * volume * third};
        const std::array<double, 3> computed = {volume * series.c[1], 2.0 * volume * volume * series.c[2], 6.0 * volume * volume * volume * series.c[3]};
        for (std::size_t k = 0; k < 3; ++k)
            identity(std::string("F3b/") + s.name + "/series" + std::to_string(k + 1), "V^k d^kA/dV^k", computed[k], exact[k], 1.0,
                IDENTITY, fields({inputs(s, "blend"), field("double_floor", floor)}), "volume series coefficients equal the tape's");
        // |p_double - p_wider| <= RT floor + 4 eps |p| + RT e_A with p = RT (rho - A_V) (lead design revision, 2026-11-09).
        // The difference is double rounding of the whole A_V (about 50 ulp at L); the block site floor under exact
        // propagation is at most 4.0e-9 Pa (blend at L). e_A is the arithmetic-noise spread of the double A_V over the six component orders
        // x +-8 one-ulp density steps, its exact linear change (from the exact A_VV of the canonical order) removed.
        const double rt = R * s.temperature, eps = std::numeric_limits<double>::epsilon();
        double floor_correction_off = NAN;
        for (const auto& [name, form] : {std::pair{"correction-off", Form::Uncorrected}, std::pair{"blend", Form::Blend},
                                         std::pair{"single", Form::Single}}) {
            long double low = 0.0L, high = 0.0L;
            double own_floor = 0.0, a_v = NAN, a_vv = NAN;
            for (const auto& order : ORDERS) {
                std::vector<double> x(3);
                for (std::size_t p = 0; p < 3; ++p) x[p] = s.x[order[p]];
                const epcsaft::HelmholtzModel ordered(hong(form, LOGISTIC, order));
                for (int j = 0; j <= 16; ++j) {
                    const int k = j % 2 ? (j + 1) / 2 : -j / 2;  // 0, 1, -1, ..., 8, -8: the canonical center first
                    double density = s.density, unused = 0.0;
                    for (int step = 0; step < std::abs(k); ++step) density = std::nextafter(density, k * INFINITY);
                    const auto series = epcsaft::properties::density_volume_taylor(ordered, s.temperature, density, x,
                        &order == &ORDERS[0] && k == 0 ? own_floor : unused);
                    if (&order == &ORDERS[0] && k == 0) a_v = series.c[1], a_vv = 2.0 * series.c[2];
                    const long double change = a_vv * (1.0L / density - 1.0L / s.density);
                    const long double noise = static_cast<long double>(series.c[1]) - a_v - change;
                    low = std::min(low, noise);
                    high = std::max(high, noise);
                }
            }
            const double e_a = static_cast<double>(high - low);
            if (form == Form::Uncorrected) floor_correction_off = own_floor;
            const epcsaft::HelmholtzModel model(hong(form));
            long double p_wider = NAN;
            std::string stop;
            try {
                const auto wide = model.isothermal<epcsaft::WiderScalar>(s.temperature, volume, {s.x.begin(), s.x.end()});
                p_wider = rt * (static_cast<long double>(s.density) - static_cast<long double>(wide.gradient[1]));
                stop = wide.association_refinement->cause;
            } catch (const std::exception& error) {
                stop = error.what();
            }
            const double p = rt * (s.density - a_v);
            const double difference = static_cast<double>(p_wider - static_cast<long double>(p));
            const double bound = rt * own_floor + 4.0 * eps * std::abs(p) + rt * e_a;
            const auto id = std::string("F3b/") + s.name + "/wider/" + name;
            const auto context = fields({inputs(s, name), field("p_double", p), field("p_wider_minus_p_double", difference),
                field("A_V", a_v), field("floor", own_floor), field("floor_correction_off", floor_correction_off), field("e_A", e_a),
                field("bound", bound), field("wider_stop", stop)});
            if (std::abs(difference) <= bound)
                recording_support::compare_custom(id.c_str(), "|p_double - p_wider|", "Pa", std::abs(difference), bound, context,
                    "Relationship", "wider replay", "the double pressure agrees with the wider replay within its rounding estimates",
                    "error <= RT floor + 4 eps |p| + RT e_A", true);
            else recording_support::observe(id.c_str(), "|p_double - p_wider|", "Pa", std::abs(difference), bound, context,
                "Observation", "wider replay", "exceeds RT floor + 4 eps |p| + RT e_A: reported (lead, 2026-11-09)");
            if (s.name[0] != 'L' || form != Form::Blend) continue;
            // The bound still separates a site defect: a 1e-4 association tolerance moves p by more than it allows.
            epcsaft::policy::SolverOptions options;
            options.set(epcsaft::policy::association, 1e-4);
            const epcsaft::HelmholtzModel loose(hong(form), {}, std::nullopt, options);
            double loose_floor = 0.0;
            const double moved = rt * (s.density - epcsaft::properties::density_volume_taylor(loose, s.temperature, s.density, s.x,
                loose_floor).c[1]);
            const double defect = static_cast<double>(std::fabs(p_wider - static_cast<long double>(moved)));
            recording_support::compare_custom("F3b/L/wider/site-defect", "|p_loose - p_wider|", "Pa", defect, bound,
                fields({context, field("association_tolerance", 1e-4), field("loose_floor", loose_floor)}), "Relationship",
                "controlled site defect", "a site defect exceeds the tight solve's bound", "defect > RT floor + 4 eps |p| + RT e_A",
                defect > bound);
        }
        if (s.name[0] == 'L') {
            // Rounding allowance on the pressure scale rho RT (p/RT is far smaller than A_V in a liquid). The correction-off
            // blend is the floor without h, reported beside it.
            for (const auto& [name, form] : {std::pair{"blend", Form::Blend}, std::pair{"correction-off", Form::Uncorrected}})
                for (const double tolerance : {1e-4, 1e-7}) {
                    epcsaft::policy::SolverOptions options;
                    options.set(epcsaft::policy::association, tolerance);
                    const epcsaft::HelmholtzModel tight(hong(form));
                    const epcsaft::HelmholtzModel loose(hong(form), {}, std::nullopt, options);
                    double tight_floor = 0.0, loose_floor = 0.0;
                    const double exact = epcsaft::properties::density_volume_taylor(tight, s.temperature, s.density, s.x, tight_floor).c[1];
                    const double moved = epcsaft::properties::density_volume_taylor(loose, s.temperature, s.density, s.x, loose_floor).c[1];
                    const double bound = loose_floor + tight_floor + 4.0 * std::numeric_limits<double>::epsilon() * s.density;
                    const std::string id = std::string("F3b/L/floor/") + name + "/" + recording_support::number(tolerance);
                    const auto context = fields({inputs(s, name), field("association_tolerance", tolerance), field("loose_floor", loose_floor),
                        field("tight_floor", tight_floor), field("bound", bound)});
                    if (form == Form::Uncorrected)
                        recording_support::observe(id.c_str(), "|A_V,loose - A_V,tight|", "mol/m3", std::abs(moved - exact), bound, context,
                            "Observation", "controlled site defect", "the floor without the correction's h terms, reported");
                    else recording_support::compare_custom(id.c_str(), "|A_V,loose - A_V,tight|", "mol/m3", std::abs(moved - exact), bound,
                        context, "Relationship", "controlled site defect", "the double floor bounds the A_V error a site defect leaves",
                        "error <= floor_loose + floor_tight + 4 eps rho", std::abs(moved - exact) <= bound);
                }
            const auto virial = blend.virial_coefficients(s.temperature, s.x);
            const double rho = 0.01, z = blend.state(s.temperature, rho, s.x).compressibility_factor;
            identity("F3b/L/B", "B", virial.second, (z - 1.0) / rho - virial.third * rho, 0.0, 1e-6,
                fields({inputs(s, "blend"), field("B_correction_off", uncorrected.virial_coefficients(s.temperature, s.x).second)}),
                "B equals the low-density slope (the remainder is O(rho^2))");
        }
        // A site defect against the exact-gradient floor (PressureSlope: A_V; Hessian: T, V, n) and against the
        // gradient coefficient error estimate of each column, both formed from the defective sites. Both are local
        // first-order estimates with a worst-sign allowance for the achieved residual, not rigorous enclosures: at Z with a
        // requested association tolerance of 0.1, the n_CO2 error 1.3086e-3 is about 20% above the first-order estimate
        // without the allowance (1.0920e-3); the complete estimate over the error is 8.53 (floor) and 8.58 (coefficient).
        if (s.name[0] != 'V') {
            const double tolerance = s.name[0] == 'L' ? 1e-4 : 0.1;
            epcsaft::policy::SolverOptions options;
            options.set(epcsaft::policy::association, tolerance);
            const epcsaft::HelmholtzModel tight(hong(Form::Blend));
            const epcsaft::HelmholtzModel loose(hong(Form::Blend), {}, std::nullopt, options);
            const double rounding = 4.0 * std::numeric_limits<double>::epsilon() * s.density;
            const auto check = [&](const std::string& id, std::size_t column, double error, double bound) {
                recording_support::compare_custom(id.c_str(), "|dA/dq_j,loose - dA/dq_j,tight|", "1", error, bound,
                    fields({inputs(s, "blend"), field("association_tolerance", tolerance), field("column", static_cast<double>(column)),
                        field("bound", bound)}),
                    "Relationship", "controlled site defect", "the error estimate bounds the gradient error a site defect leaves",
                    "error <= estimate_loose + estimate_tight + 4 eps rho", error <= bound);
            };
            for (const auto order : {epcsaft::DerivativeOrder::PressureSlope, epcsaft::DerivativeOrder::Hessian}) {
                const auto exact = tight.evaluate(s.temperature, volume, s.x, order);
                const auto moved = loose.evaluate(s.temperature, volume, s.x, order);
                const bool slope = order == epcsaft::DerivativeOrder::PressureSlope;
                for (std::size_t j = slope ? 1 : 0; j < (slope ? 2 : 5); ++j)
                    check(std::string("F3b/") + s.name + "/gradient-floor/" + (slope ? "slope/" : "hessian/") + std::to_string(j), j,
                        std::abs(moved.gradient[j] - exact.gradient[j]), moved.gradient_floor + exact.gradient_floor + rounding);
            }
            const auto exact = tight.evaluate(s.temperature, volume, s.x, epcsaft::DerivativeOrder::Gradient);
            const auto moved = loose.evaluate(s.temperature, volume, s.x, epcsaft::DerivativeOrder::Gradient);
            const auto defective = loose.state(s.temperature, s.density, s.x).site_fractions;
            const std::vector<double> point{s.temperature, volume, s.x[0], s.x[1], s.x[2]};
            for (std::size_t j = 0; j < 5; ++j) {
                std::vector<double> unit(5, 0.0);
                unit[j] = 1.0;
                const double bound = tight.residual_coefficient_errors({point, unit}, defective)[1]
                    + tight.residual_coefficient_errors({point, unit}, exact.site_fractions)[1] + rounding;
                check(std::string("F3b/") + s.name + "/coefficient-error/" + std::to_string(j), j, std::abs(moved.gradient[j] - exact.gradient[j]), bound);
            }
        }
    }

    // F4: one consistent energy. (i) partial molar residual enthalpy at 333.15 K, 2e5 Pa, x_L, liquid; (ii) Euler.
    {
        const double pressure = 2e5;
        const epcsaft::PressureActionDirection along_t{1.0, 0.0, {0.0, 0.0, 0.0}, {}};
        const auto lnphi = blend.pressure_actions(L.temperature, pressure, L.x, L.branch,
            {PropertyObservable::LogFugacityCoefficient}, {along_t});
        for (std::size_t i = 0; i < 3; ++i) {
            std::vector<double> d(3);
            for (std::size_t j = 0; j < 3; ++j) d[j] = (i == j ? 1.0 : 0.0) - L.x[j];
            const auto h = blend.pressure_actions(L.temperature, pressure, L.x, L.branch,
                {PropertyObservable::ResidualEnthalpy}, {{0.0, 0.0, d, {}}}).results.front();
            identity("F4/partial-enthalpy/" + std::to_string(i), "(h_i + R T^2 dlnphi_i/dT)/RT", (h.value.value() + h.action.value())
                / (R * L.temperature) + L.temperature * lnphi.results[i].action.value(), 0.0, 1.0, IDENTITY, inputs(L, "blend"),
                "the caloric partial molar residual enthalpy equals -R T^2 (d ln phi_i/dT)_P,x");
        }
    }
    for (const auto& s : STATES) {
        const auto state = blend.state(s.temperature, s.density, s.x);
        double sum = 0.0;
        for (std::size_t i = 0; i < 3; ++i) sum += s.x[i] * state.residual_chemical_potential_over_rt[i].value();
        identity(std::string("F4/euler/") + s.name, "a - sum x mu + Z - 1", state.residual_helmholtz - sum
            + state.compressibility_factor - 1.0, 0.0, std::max(1.0, std::abs(state.residual_helmholtz)), 1e-12, inputs(s, "blend"),
            "Euler homogeneity at fixed T");
        const auto terms = blend.log_fugacity_terms(s.temperature, s.density, s.x);
        for (std::size_t i = 0; i < 3 && s.x[i] > 0.0; ++i) {
            const double total = terms.hard_chain[i].value() + terms.dispersion[i].value() + terms.association[i].value()
                - terms.log_compressibility_factor;
            identity(std::string("F4/terms/") + s.name + "/" + std::to_string(i), "sum of ln phi terms", total,
                state.log_fugacity_coefficient.value()[i].value(), 1.0, IDENTITY, inputs(s, "blend"),
                "the per-term log-fugacity decomposition sums to ln phi, the implicit site chain included");
        }
    }

    // F4(iii): tangent-plane stability of the declared w(alpha, T) liquid on the Table S1 grid (30 wt% MEA on a CO2-free
    // basis, 313.15-393.15 K, alpha 0.1-1.0) at the F4 pressure; reported, not asserted. A negative liquid-trial minimum
    // is a predicted split of this fluid, a negative vapor-trial minimum a liquid above its bubble pressure at 2e5 Pa.
    // The minima are over the Engine's local trial set: no descent found is not a proof of stability.
    for (std::size_t t = 0; t < 9; ++t)
        for (std::size_t a = 1; a <= 10; ++a) {
            const double temperature = 313.15 + 10.0 * static_cast<double>(t), loading = 0.1 * static_cast<double>(a);
            const double amine = Z.x[0], total = 1.0 + amine * loading;
            const std::vector<double> x{amine / total, (1.0 - amine) / total, amine * loading / total};
            auto problem = equilibrium_support::flash(3);
            problem.temperature = epcsaft::equilibrium::FreeValue::fixed(temperature);
            problem.pressure = epcsaft::equilibrium::FreeValue::fixed(2e5);
            problem.feed.values = x;
            const auto compiled = epcsaft::equilibrium::compile_problem(blend, problem);
            const auto parent = epcsaft::equilibrium::solve_compiled(compiled);
            std::map<std::string, double> least{{"liquid", NAN}, {"vapor", NAN}};
            std::string cause = parent.success ? "" : "liquid solve: " + parent.message;
            if (parent.success) {
                const auto measure = epcsaft::equilibrium::stability_measure(compiled, parent);
                cause = measure.cause;
                for (const auto& trial : measure.report.trials) {
                    const auto tpd = trial.quantities.find("tpd_over_rt");
                    auto& slot = least[equilibrium_support::phase_kind_name(trial.phase.kind)];
                    if (tpd != trial.quantities.end() && !(tpd->second >= slot)) slot = tpd->second;
                }
            }
            const auto id = "F4/stability/" + recording_support::number(temperature) + "/" + recording_support::number(loading);
            recording_support::observe(id.c_str(), "min tpd/RT", "1", std::fmin(least["liquid"], least["vapor"]), std::nullopt,
                fields({field("model", "blend"), field("T_K", temperature), field("p_Pa", 2e5), field("alpha", loading), field("x", x),
                    field("liquid_trial_min", least["liquid"]), field("vapor_trial_min", least["vapor"]), field("cause", cause)}),
                "Observation", "tangent-plane stability", "reported, not asserted: instability is a model prediction");
        }

    // Domain: no MEA with CO2 present needs the loading only when b_alpha is nonzero or active.
    const std::vector<double> amine_free = {0.0, 0.6, 0.4};
    const epcsaft::HelmholtzModel zero_slope(hong(Form::Blend, {1.7, 0.0, -1.5}));
    CHECK(std::isfinite(zero_slope.evaluate(353.15, 1e-3, amine_free, epcsaft::DerivativeOrder::Hessian).value));
    for (const auto& model : {&blend, &uncorrected}) {
        bool refused = false;
        try { (void)model->evaluate(353.15, 1e-3, amine_free, epcsaft::DerivativeOrder::Gradient); }
        catch (const std::invalid_argument& error) { refused = std::string(error.what()).find("loading undefined") != std::string::npos; }
        CHECK(refused);
    }
    bool active_zero_refused = false;
    try {
        const epcsaft::HelmholtzModel active(hong(Form::Blend, {1.7, 0.0, -1.5}), {{ParameterFamily::PathwayWeight, {1}}});
        (void)active.evaluate(353.15, 1e-3, amine_free, epcsaft::DerivativeOrder::Gradient);
    } catch (const std::invalid_argument& error) {
        active_zero_refused = std::string(error.what()).find("loading undefined") != std::string::npos;
    }
    CHECK(active_zero_refused);
}

}  // namespace properties_checks
