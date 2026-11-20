#include "properties_checks_shared.hpp"

#include <cstring>
#include <limits>

namespace properties_checks {
namespace {

using namespace epcsaft;
using Wide = WiderScalar;
using recording_support::Rule;

void compare(const std::string& id, const char* quantity, const char* unit, double actual, double expected,
             double absolute, double relative, const char* method) {
    recording_support::compare("msa/" + id, quantity, unit, actual, expected, Rule::ReferenceScale,
        absolute, relative, method);
}

struct Ion { Wide n, d, z; };
using Ions = std::vector<Ion>;

struct PrimitiveReference { Wide gamma, energy, internal, pressure; };

// Blum and Høye (2.13–2.19), (3.1), (3.4), (3.7), evaluated independently with bisection.
PrimitiveReference primitive_reference(Wide ell, const Ions& ions) {
    using std::acos;
    const Wide pi = acos(Wide(-1));
    Wide packing{}, s0{};
    for (const auto& [n, d, z] : ions) {
        packing += n * d * d * d;
        s0 += n * z * z;
    }
    const Wide b = pi / (2 * (1 - pi * packing / 6));
    const auto moments = [&](Wide gamma) {
        Wide num = 0, omega = 1, energy = 0, q = 0;
        for (const auto& [n, d, z] : ions) {
            const Wide w = n / (1 + gamma * d);
            num += w * z * d;
            omega += b * w * d * d * d;
            energy += w * z * z;
        }
        const Wide pn = num / omega;
        for (const auto& [n, d, z] : ions) {
            const Wide e = (z - b * pn * d * d) / (1 + gamma * d);
            q += n * e * e;
        }
        return std::array<Wide, 3>{q, -ell * (gamma * energy + b * omega * pn * pn), pn};
    };
    Wide lo = 0, hi = wide_math::sqrt(pi * ell * s0);
    while (hi * hi < pi * ell * moments(hi)[0]) hi *= 2;
    for (int i = 0; i < std::numeric_limits<Wide>::digits; ++i) {
        const Wide mid = (lo + hi) / 2;
        if (mid * mid < pi * ell * moments(mid)[0]) lo = mid;
        else hi = mid;
    }
    const Wide gamma = (lo + hi) / 2;
    const auto [q, internal, pn] = moments(gamma);
    const Wide entropy = gamma * gamma * gamma / (3 * pi);
    const Wide y = b * pn;
    return {gamma, internal + entropy, internal, -entropy - 2 * ell * y * y / pi};
}

constexpr double root_rtol = 2.0 * policy::DENSITY_ROOT_RESIDUAL.value;
constexpr double roundoff = 32.0 * std::numeric_limits<double>::epsilon();

void primitive_values() {
    for (const double n : {1e-11, 1e-6, 0.004}) {
        const double ell = 7.14, d = 4.25;
        const std::vector<double> q{ell, n, d, 2.0, 2.0 * n, d, -1.0};
        const double kappa = std::sqrt(4 * PI * ell * 6 * n);
        const double gamma = msa_gamma(q);
        const double exact = kappa / (1 + std::sqrt(1 + 2 * kappa * d));
        const double energy = -exact * exact * exact * (2 + 3 * exact * d) / (3 * PI);
        compare("equal/gamma", "Gamma", "1/angstrom", gamma, exact, 0, root_rtol, "equal-diameter-closed-form");
        compare("equal/energy", "beta A/V", "1/angstrom^3", msa_helmholtz_density(gamma, q),
            energy, 0, roundoff, "equal-diameter-closed-form");
    }
    const std::vector<double> counterexample{1e-9, 0.004, 1.0, 1.0, 0.001, 5.0, 1.0, 0.005, 2.0, -1.0};
    const double gamma = msa_gamma(counterexample);
    CHECK(gamma > std::sqrt(PI * 1e-9 * 0.01));
    compare("bracket/gamma", "Gamma L", "1", gamma, 5.6120130e-6, 0.5e-13, 0, "review-reference");

    constexpr double concentration = 0.0256, ell = 7.14;
    const auto cation_density = [](double c) { return 1000 * c * NUMBER_DENSITY_PER_ANGSTROM3; };
    const double n = cation_density(concentration);
    const auto source = [&](double c) {
        const Wide nplus = cation_density(c);
        return primitive_reference(ell, {{nplus, 20, 10}, {10 * nplus, 4.25, -1}});
    };
    const auto nominal = source(concentration);
    const auto q = std::vector<Taylor<1>>{ell, n, 20, 10, 10 * n, 4.25, -1};
    auto seeded = q;
    seeded[0].c[1] = ell;
    const auto actual_gamma = msa_gamma(seeded);
    const auto actual_energy = msa_helmholtz_density(actual_gamma, seeded);
    const double interval = std::max(
        std::abs(static_cast<double>(source(concentration - 0.0002).internal / (11 * cation_density(concentration - 0.0002))
            - nominal.internal / (11 * n))),
        std::abs(static_cast<double>(source(concentration + 0.0002).internal / (11 * cation_density(concentration + 0.0002))
            - nominal.internal / (11 * n))));
    compare("mohoric/internal", "beta U/N_ions", "1", actual_energy.c[1] / (11 * n), -2.431,
        0.0005 + interval, 0, "Mohoric-Table-6");
    compare("mohoric/gamma", "Gamma", "1/angstrom", actual_gamma.c[0], static_cast<double>(nominal.gamma),
        0, root_rtol, "wider-precision-bisection");
    compare("mohoric/energy", "beta A/N_ions", "1", actual_energy.c[0] / (11 * n),
        static_cast<double>(nominal.energy / (11 * n)), 0, roundoff, "wider-precision-bisection");
}

constexpr double bjerrum_angstrom_kelvin = ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C
    / (4 * PI * VACUUM_PERMITTIVITY_F_PER_M * BOLTZMANN_J_PER_K * ANGSTROM_M);

void dilute_limit() {
    auto p = properties_support::aqueous_electrolyte();
    p.ion_suppression_coefficient = {0, 0, 0};
    p.debye_huckel_diameter_angstrom = p.sigma_constant;
    const double T = 298.15, rho = 24000;
    const double ell = bjerrum_angstrom_kelvin / (T * p.relative_permittivity[0]);
    const double density = rho * NUMBER_DENSITY_PER_ANGSTROM3;
    for (const double xion : {1e-4, 1e-6, 1e-8, 1e-10}) {
        const std::vector<double> x{1 - 2 * xion, xion, xion};
        p.electrolyte_form = ElectrolyteForm::Msa;
        const double msa = electrolyte_contribution(T, rho, x, x, p).screening * density;
        p.electrolyte_form = ElectrolyteForm::DebyeHuckel;
        const double dh = electrolyte_contribution(T, rho, x, x, p).screening * density;
        const double n = density * xion, da = p.sigma_constant[1], db = p.sigma_constant[2];
        const double kappa = std::sqrt(8 * PI * ell * n);
        const double leading = -kappa * kappa * kappa / (12 * PI);
        const double b = PI / (2 * (1 - PI * n * (da * da * da + db * db * db) / 6));
        const double size_bound = (kappa * kappa / 4) * n * (da + db);
        const double coupling_bound = b * (n * (da + db)) * (n * (da + db));
        compare("dilute/finite-size-dh", "MSA minus finite-size DH / |DH limit|", "1",
            (msa - dh) / std::abs(leading), 0, ell * (2 * size_bound + coupling_bound) / std::abs(leading) + roundoff, 0,
            "energy-minimum-and-DH-series-bound");
        compare("dilute/leading", "MSA / DH limit", "1", msa / leading, 1,
            ell * (size_bound + coupling_bound) / std::abs(leading) + roundoff, 0, "point-ion-DH-leading-coefficient");
    }
}

using Point = std::array<Wide, 4>; // dimensionless T, V, salt amount and cation diameter
using Directions = std::vector<std::size_t>;

struct Difference { Wide value, floor; };

// Tensor product of centered first differences. Repetition gives a centered higher difference.
template <class Function>
Difference centered(const Function& f, Point point, const Directions& directions, Wide h, std::size_t i = 0) {
    if (i == directions.size()) return f(point);
    const auto axis = directions[i];
    point[axis] += h;
    const auto plus = centered(f, point, directions, h, i + 1);
    point[axis] -= 2 * h;
    const auto minus = centered(f, point, directions, h, i + 1);
    return {(plus.value - minus.value) / (2 * h), (plus.floor + minus.floor) / (2 * h)};
}

template <class Function>
void verify_derivative(const char* id, double actual, const Function& f, const Directions& axes,
                       double subtraction_scale = 0.0, const char* unit = "1") {
    const auto a = centered(f, {1, 1, 1, 1}, axes, Wide(0.004));
    const auto b = centered(f, {1, 1, 1, 1}, axes, Wide(0.002));
    const auto c = centered(f, {1, 1, 1, 1}, axes, Wide(0.001));
    const Wide reference = (4 * c.value - b.value) / 3;
    const Wide previous = (4 * b.value - a.value) / 3;
    const double allowance = static_cast<double>(wide_math::abs(reference - previous) + (4 * c.floor + b.floor) / 3)
        + roundoff * subtraction_scale;
    compare(std::string("derivative/") + id, "scaled derivative", unit, actual,
        static_cast<double>(reference), allowance, roundoff, "centered-difference-refinement");
}

void derivatives_and_neutral_support() {
    auto dh = properties_support::aqueous_electrolyte();
    dh.sigma_constant[1] = 4.4;
    dh.debye_huckel_diameter_angstrom = dh.sigma_constant;
    dh.permittivity_form = PermittivityForm::DensityTemperatureArithmetic;
    dh.dielectric_volume_m3_per_mol[0] = 0.03;
    dh.dielectric_temperature_k[0] = 1000;
    auto msa = dh;
    msa.electrolyte_form = ElectrolyteForm::MsaBorn;
    const std::vector<ActiveParameter> active{{ParameterFamily::SigmaConstant, {1}}};
    const HelmholtzModel m(msa, active), d(dh, active);
    const std::vector<double> x{0.88, 0.06, 0.06};
    const double d0 = dh.sigma_constant[1], d1 = dh.sigma_constant[2];
    for (const auto& [T, rho] : {std::pair{390.0, 600.0}, std::pair{430.0, 900.0}}) {
        const double V = 1 / rho;
        const auto physical = [&](const Point& q) {
            const Wide ell = Wide(bjerrum_angstrom_kelvin) / (T * q[0] * (1 + Wide(x[0]) / (V * q[1])
                * dh.dielectric_volume_m3_per_mol[0] * (dh.dielectric_temperature_k[0] / (T * q[0]) - 1)));
            const Wide n = Wide(x[1]) * q[2] * NUMBER_DENSITY_PER_ANGSTROM3 / (V * q[1]);
            return std::pair<Wide, Ions>{ell, {{n, d0 * q[3], 1}, {n, d1, -1}}};
        };
        const auto energy_difference = [&](const Point& q) {
            const auto [ell, ions] = physical(q);
            const Wide k = wide_math::sqrt(8 * Wide(PI) * ell * ions.front().n);
            const auto chi = [](Wide y) { return 3 * (wide_math::log1p(y) - y + y * y / 2) / (y * y * y); };
            const Wide dh_density = -k * ell * ions.front().n * (chi(k * d0) + chi(k * d1)) / 3;
            const Wide msa_density = primitive_reference(ell, ions).energy;
            const Wide conversion = V * q[1] / NUMBER_DENSITY_PER_ANGSTROM3;
            // The independently evaluated terms cancel, so the reference floor uses their magnitudes.
            const Wide floor = 32 * std::numeric_limits<Wide>::epsilon() * conversion
                * (wide_math::abs(msa_density) + wide_math::abs(dh_density));
            return Difference{(msa_density - dh_density) * conversion, floor};
        };
        const auto a = m.evaluate(T, V, x, DerivativeOrder::Hessian);
        const auto b = d.evaluate(T, V, x, DerivativeOrder::Hessian);
        const State ma = m.state(T, rho, x), db = d.state(T, rho, x);
        CHECK(ma.born == db.born && ma.msa != 0 && ma.debye_huckel == 0);
        const auto terms = m.log_fugacity_terms(T, rho, x);
        CHECK(terms.msa[1].value() != 0 && terms.debye_huckel[1] == 0.0);
        verify_derivative("temperature", T * (a.gradient[0] - b.gradient[0]), energy_difference, {0},
            T * std::max(std::abs(a.gradient[0]), std::abs(b.gradient[0])));
        verify_derivative("pressure", -(ma.pressure - db.pressure) * V / (GAS_CONSTANT_J_PER_MOL_K * T),
            energy_difference, {1}, V * std::max(std::abs(ma.pressure), std::abs(db.pressure)) / (GAS_CONSTANT_J_PER_MOL_K * T));
        verify_derivative("salt-chemical-potential", x[1] * ((a.gradient[3] - b.gradient[3]) + (a.gradient[4] - b.gradient[4])),
            energy_difference, {2}, x[1] * (std::abs(a.gradient[3]) + std::abs(a.gradient[4])));
        verify_derivative("active-diameter", d0 * (a.gradient[5] - b.gradient[5]), energy_difference, {3},
            d0 * std::max(std::abs(a.gradient[5]), std::abs(b.gradient[5])));
        verify_derivative("mixed-hessian", V * x[1] * (a.hessian[9] + a.hessian[10] - b.hessian[9] - b.hessian[10]),
            energy_difference, {1, 2}, V * x[1] * (std::abs(a.hessian[9]) + std::abs(a.hessian[10])));
        double floor = 0;
        const auto mt = properties::density_volume_taylor(m, T, rho, x, floor);
        const auto dt = properties::density_volume_taylor(d, T, rho, x, floor);
        verify_derivative("density-taylor", 2 * V * V * (mt.c[2] - dt.c[2]), energy_difference, {1, 1},
            2 * V * V * std::max(std::abs(mt.c[2]), std::abs(dt.c[2])));
        const std::vector<double> v{0, V, 0, 0, 0, 0}, size{0, 0, 0, 0, 0, d0};
        const auto third = m.residual_energy_action(T, V, x, {v, v, size});
        const auto third_dh = d.residual_energy_action(T, V, x, {v, v, size});
        const auto reverse = m.repeated_direction_energy_gradient(T, V, x, v, {5});
        const auto reverse_dh = d.repeated_direction_energy_gradient(T, V, x, v, {5});
        CHECK(third.status == PropertyActionStatus::Available && third_dh.status == PropertyActionStatus::Available
            && reverse.status == PropertyActionStatus::Available && reverse_dh.status == PropertyActionStatus::Available);
        verify_derivative("mixed-forward-third", third.action.value() - third_dh.action.value(), energy_difference, {1, 1, 3},
            std::abs(third.action.value()) + std::abs(third_dh.action.value()));
        verify_derivative("mixed-reverse-third", d0 * (reverse.gradient[0] - reverse_dh.gradient[0]),
            energy_difference, {1, 1, 3}, d0 * (std::abs(reverse.gradient[0]) + std::abs(reverse_dh.gradient[0])));
        const auto q0 = physical({1, 1, 1, 1});
        const auto exact = primitive_reference(q0.first, q0.second);
        const auto eps = [](const State& s) { return s.bulk_relative_permittivity; };
        CHECK(eps(ma) == eps(db));
        // At fixed dielectric and size, vary volume and the amounts along the neutral ion ratio.
        const auto rho_number = rho * NUMBER_DENSITY_PER_ANGSTROM3;
        const double analytic_mean = static_cast<double>((exact.energy + exact.pressure) / (2 * q0.second[0].n));
        const double ell = static_cast<double>(q0.first);
        const double n = static_cast<double>(q0.second[0].n);
        std::vector<Taylor<1>> fixed{ell, n, d0, 1, n, d1, -1};
        fixed[1].c[1] = fixed[4].c[1] = n;
        const auto pressure_energy = msa_helmholtz_density(msa_gamma(fixed), fixed);
        compare("derivative/primitive-pressure", "beta P / total number density", "1",
            (pressure_energy.c[1] - pressure_energy.c[0]) / rho_number, static_cast<double>(exact.pressure / rho_number),
            0, roundoff, "Blum-Hoye-3.7");
        compare("derivative/primitive-mean", "beta mu mean", "1", pressure_energy.c[1] / (2 * n),
            analytic_mean, 0, roundoff, "Blum-Hoye-3.9");
        verify_derivative("Gamma-density", msa_gamma(fixed).c[1], [&](const Point& q) {
            const Wide gamma = primitive_reference(ell, {{n * q[2], d0, 1}, {n * q[2], d1, -1}}).gamma;
            return Difference{gamma, 32 * std::numeric_limits<Wide>::epsilon() * gamma};
        }, {2}, 0, "1/angstrom");
        using Inner = Taylor<1, Wide>;
        using Outer = Taylor<1, Inner>;
        std::vector<Outer> nested{Inner(ell), Inner(n), Inner(d0), Inner(1), Inner(n), Inner(d1), Inner(-1)};
        nested[0].c[0].c[1] = ell;
        nested[2].c[1].c[0] = d0;
        verify_derivative("Gamma-nested", static_cast<double>(msa_gamma(nested).c[1].c[1]), [&](const Point& q) {
            const Wide gamma = primitive_reference(ell * q[0], {{n, d0 * q[3], 1}, {n, d1, -1}}).gamma;
            return Difference{gamma, 32 * std::numeric_limits<Wide>::epsilon() * gamma};
        }, {0, 3}, 0, "1/angstrom");
    }
    const std::vector<double> zero{1, 0, 0};
    const double T = 410, rho = 800;
    const auto ma = m.state(T, rho, zero), db = d.state(T, rho, zero);
    const auto a = m.evaluate(T, 1 / rho, zero, DerivativeOrder::Hessian);
    const auto b = d.evaluate(T, 1 / rho, zero, DerivativeOrder::Hessian);
    const auto same = [](double lhs, double rhs) { return std::memcmp(&lhs, &rhs, sizeof(double)) == 0; };
    const auto columns = [&](const auto& lhs, const auto& rhs) {
        return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin(), same);
    };
    const bool identical = same(ma.msa, 0.0) && same(ma.residual_helmholtz, db.residual_helmholtz)
        && same(ma.pressure, db.pressure)
        && same(ma.residual_chemical_potential_over_rt.front().value(), db.residual_chemical_potential_over_rt.front().value())
        && ma.component_property_status == db.component_property_status && columns(a.gradient, b.gradient)
        && columns(a.hessian, b.hessian);
    compare("neutral/bit-identity", "different values", "1", !identical, 0, 0, 0, "binary64-bit-comparison");
}

}  // namespace

void primitive_msa() {
    primitive_values();
    dilute_limit();
    derivatives_and_neutral_support();
}

}  // namespace properties_checks
