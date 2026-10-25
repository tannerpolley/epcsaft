// D1/D6/D7 of #658: manufactured aqueous liquid, exact actions, faces and dilute orders.
#include "properties_checks_shared.hpp"
#include "regression/fit.hpp"

namespace properties_checks {
void cation_mixing() {
    using namespace epcsaft;
    auto base = properties_support::parameters(4);
    base.segment_count = {1.2047, 1, 1, 1};
    base.sigma_constant = {2.7927, 2.8, 3, 2.756};
    base.dispersion_energy_over_k = {353.95, 230, 230, 170};
    base.molar_mass_kg_per_mol = {0.01801528, 0.00694, 0.024305, 0.035453};
    base.charge_number = {0, 1, 2, -1};
    base.packing_diameter_angstrom = {0, 2.5, 2.5, 2.5};
    base.debye_huckel_diameter_angstrom = {0, 3, 3, 3};
    base.relative_permittivity = {78.09, 0, 0, 0};
    base.electrolyte_form = ElectrolyteForm::DebyeHuckel;
    base.ionic_dispersion = IonicDispersion::ExcludeSameSign;
    base.permittivity_form = PermittivityForm::IonSuppression;
    auto p = base;
    p.cation_mixing = {{0, 1, 2, 3}};
    p.cation_pair_mixing = {-0.2};
    p.cation_triple_mixing = {0.02};
    const std::vector<ActiveParameter> active = {{ParameterFamily::CationPairMixing, {0}}, {ParameterFamily::CationTripleMixing, {0}}};
    const double T = 298.15, mass = p.molar_mass_kg_per_mol[0];
    const HelmholtzModel plain(base), mixed(p, active);
    std::printf("a,b,quantity,step,value\n");
    // #781: an extra neutral uses normalization, never extra solvent mass.
    auto neutral = properties_support::parameters(5);
    neutral.segment_count = {1.2047, 1, 1, 1, 5.306};
    neutral.sigma_constant = {2.7927, 2.8, 3, 2.756, 3.8959};
    neutral.dispersion_energy_over_k = {353.95, 230, 230, 170, 249.21};
    neutral.molar_mass_kg_per_mol = {mass, 0.00694, 0.024305, 0.035453, 0.17034};
    neutral.charge_number = {0, 1, 2, -1, 0};
    neutral.relative_permittivity = {78.09, 0, 0, 0, 2};
    neutral.packing_diameter_angstrom = {0, 2.5, 2.5, 2.5, 0};
    neutral.debye_huckel_diameter_angstrom = {0, 3, 3, 3, 0};
    neutral.electrolyte_form = p.electrolyte_form;
    neutral.ionic_dispersion = p.ionic_dispersion;
    neutral.permittivity_form = p.permittivity_form;
    neutral.cation_mixing = p.cation_mixing;
    neutral.cation_pair_mixing = p.cation_pair_mixing;
    neutral.cation_triple_mixing = p.cation_triple_mixing;
    const std::vector<double> aqueous = {1 / mass, 0.2, 0.4, 1.0};
    auto extended = aqueous;
    extended.push_back(0);
    const double volume = std::accumulate(aqueous.begin(), aqueous.end(), 0.0) / 30000;
    const auto four = HelmholtzModel(p).evaluate(T, volume, aqueous, DerivativeOrder::Hessian);
    const auto five = HelmholtzModel(neutral).evaluate(T, volume, extended, DerivativeOrder::Hessian);
    CHECK(four.value == five.value);
    for (std::size_t row = 0; row < 6; ++row) {
        CHECK(four.gradient[row] == five.gradient[row]);
        for (std::size_t col = 0; col < 6; ++col) CHECK(four.hessian[row*6+col] == five.hessian[row*7+col]);
    }
    // #819: fixed dry supports have zero correction; dense transverse entries are storage only.
    auto dry_plain = neutral;
    dry_plain.cation_mixing.clear(); dry_plain.cation_pair_mixing.clear(); dry_plain.cation_triple_mixing.clear();
    const HelmholtzModel dry_model(neutral, active), dry_off(dry_plain);
    const std::vector<std::vector<double>> supports = {{0, 0, 0, 0, 1}, {0, .1, 0, .1, .8},
        {0, 0, .1, .2, .7}, {.9, .02, .02, .06, 0}};
    // Exercise both cached positions and eviction by a third support, in both transition orders.
    for (const std::size_t index : {3U, 0U, 3U, 0U, 1U, 2U, 3U, 2U, 0U}) {
        const auto& x = supports[index];
        const auto cached = dry_model.evaluate(T, 1/1000., x, DerivativeOrder::Hessian);
        const auto fresh = HelmholtzModel(neutral, active).evaluate(T, 1/1000., x, DerivativeOrder::Hessian);
        CHECK(cached.value == fresh.value && cached.gradient == fresh.gradient && cached.hessian == fresh.hessian);
        if (index == 3) continue;
        const auto old = dry_off.evaluate(T, 1/1000., x, DerivativeOrder::Hessian);
        CHECK(cached.value == old.value);
        for (std::size_t row = 0; row < 9; ++row) {
            CHECK(cached.gradient[row] == (row < 7 ? old.gradient[row] : 0));
            for (std::size_t col = 0; col < 9; ++col)
                CHECK(cached.hessian[row*9+col] == (row < 7 && col < 7 ? old.hessian[row*7+col] : 0));
        }
        std::vector<detail::AD> variables(x.begin(), x.end());
        CppAD::Independent(variables);
        CppAD::ADFun<double> tape(variables, {cation_mixing_helmholtz(variables, neutral)});
        CHECK(tape.Forward(0, x)[0] == 0);
        for (const auto value : tape.Jacobian(x)) CHECK(value == 0);
        for (const auto value : tape.Hessian(x, 0)) CHECK(value == 0);
        const auto state = dry_model.state_at_pressure(T, 101325, x, DensityBranch::Liquid);
        CHECK(state.cation_mixing == 0);
        const auto missing = index == 1 ? 2U : 1U;
        CHECK(!state.residual_chemical_potential_over_rt[missing]);
        CHECK(state.component_property_status[missing] == PropertyActionStatus::IonicInsertion);
        CHECK(state.log_fugacity_coefficient && !(*state.log_fugacity_coefficient)[missing]);
        if (index != 0) CHECK(state.residual_chemical_potential_over_rt[3-missing]);
        const auto actions = dry_model.pressure_actions(T, 101325, x, DensityBranch::Liquid,
            {PropertyObservable::LogFugacityCoefficient}, {{1, 0, {0, 0, 0, 0, 0}, {}}});
        CHECK(actions.results[missing].status == PropertyActionStatus::IonicInsertion);
        CHECK(!actions.results[missing].value && !actions.results[missing].action);
        std::vector<double> insertion(9, 0); insertion[2+missing] = 1;
        const auto refused = dry_model.residual_energy_action(T, 1/1000., x, {insertion, insertion, insertion});
        CHECK(!refused.action && refused.status == PropertyActionStatus::IonicInsertion);
        std::vector<std::vector<double>> directions(3, std::vector<double>(9, 0));
        for (auto& direction : directions) direction[7] = 1;
        const auto zero = dry_model.residual_energy_action(T, 1/1000., x, directions);
        std::printf("dry-action,%zu,status,%d,value,%.17g\n", index, static_cast<int>(zero.status), zero.action.value_or(-999));
        CHECK(zero.action && *zero.action == 0);
        std::printf("dry,%zu,correction,%.17g,dense-difference,0,insertion-status,%d\n",
            index, state.cation_mixing, static_cast<int>(state.component_property_status[missing]));
    }
    for (const auto& x : {std::vector<double>{0, .1, .1, .3, .5}, {0, 1e-200, 1e-200, 3e-200, 1}}) {
        for (const bool zero_coefficients : {false, true}) {
            auto rejected = neutral;
            if (zero_coefficients) rejected.cation_pair_mixing[0] = rejected.cation_triple_mixing[0] = 0;
            bool refused = false;
            try { (void)cation_mixing_helmholtz(x, rejected); } catch (const std::invalid_argument&) { refused = true; }
            CHECK(refused);
        }
        bool refused = false;
        try { (void)dry_model.state(T, 1000, x); } catch (const std::invalid_argument&) { refused = true; }
        CHECK(refused);
    }
    // Wet faces retain their analytic insertion and mixed derivative even at very small positive solvent.
    for (const double solvent : {1.0, 1e-20}) for (const std::size_t missing : {1U, 2U}) {
        std::vector<double> x = {solvent, .1, .2, .3, .4}; x[missing] = 0;
        std::vector<detail::AD> variables(x.begin(), x.end());
        CppAD::Independent(variables);
        CppAD::ADFun<double> tape(variables, {cation_mixing_helmholtz(variables, neutral)});
        const double W = mass*solvent, cross = -.4/W + .02*x[3]/(W*W);
        CHECK(close(tape.Jacobian(x)[missing], x[3-missing]*cross, 1e-12));
        CHECK(close(tape.Hessian(x, 0)[1*5+2], cross, 1e-12));
        CHECK(tape.Jacobian(x)[missing] != 0);
        std::printf("wet-face,%.17g,%zu,insertion,%.17g,mixed,%.17g\n", solvent, missing,
            tape.Jacobian(x)[missing], tape.Hessian(x, 0)[1*5+2]);
    }
    neutral.cation_pair_mixing = {-0.2205};
    neutral.cation_triple_mixing = {0.024};
    extended[3] = 1.35;
    std::vector<double> reference_mu;
    for (const double amount : {1.0, 2.0}) for (const double neutral_mass : {0.17034, 0.34068}) {
        extended[4] = amount;
        neutral.molar_mass_kg_per_mol[4] = neutral_mass;
        const double F = cation_mixing_helmholtz(extended, neutral);
        CHECK(std::abs((F - (-0.032688)) / F) <= 1e-15);
        std::vector<detail::AD> amounts(extended.begin(), extended.end());
        CppAD::Independent(amounts);
        CppAD::ADFun<double> tape(amounts, {cation_mixing_helmholtz(amounts, neutral)});
        const auto mu = tape.Jacobian(extended);
        if (reference_mu.empty()) reference_mu = mu;
        CHECK(mu == reference_mu && mu[4] == 0);
        const auto hessian = tape.Hessian(extended, 0);
        for (std::size_t row = 0; row < 5; ++row) CHECK(hessian[row*5+4] == 0 && hessian[4*5+row] == 0);
        auto x = extended;
        const double total = std::accumulate(x.begin(), x.end(), 0.0);
        for (auto& value : x) value /= total;
        const auto state = HelmholtzModel(neutral).state(T, 30000, x);
        CHECK(std::abs((state.cation_mixing - F/total) / (F/total)) <= 1e-15);
        std::printf("neutral,%.17g,mass,%.17g,F,%.17g,molar,%.17g,mu,%.17g,%.17g,%.17g,%.17g\n",
            amount, neutral_mass, F, state.cation_mixing, mu[0], mu[1], mu[2], mu[3]);
    }
    // #814: hand-derived Li-Mg-Cl-SO4 reference, 1 kg water; undeclared sulfate psi is zero.
    auto sulfate = neutral;
    const double M = 0.0180153;
    sulfate.molar_mass_kg_per_mol = {M, 0.00694, 0.024305, 0.035453, 0.09606};
    sulfate.charge_number[4] = -2;
    sulfate.relative_permittivity[4] = 0;
    sulfate.packing_diameter_angstrom[4] = 2.5;
    sulfate.debye_huckel_diameter_angstrom[4] = 3;
    auto without = sulfate;
    without.cation_mixing.clear(); without.cation_pair_mixing.clear(); without.cation_triple_mixing.clear();
    const HelmholtzModel sulfate_model(sulfate), sulfate_plain(without);
    const std::vector<double> n = {1/M, 0.2, 0.4, 0.6, 0.2}, scales = {1/M, 1, 1, 1, 1};
    const double N = std::accumulate(n.begin(), n.end(), 0.0), V = N/30000, F = -0.034128;
    const std::vector<double> g = {0.032976*M, -0.17064, -0.08532, 0.00192, 0};
    const std::vector<double> H = {-0.063648*M*M, 0.16488*M, 0.08244*M, -0.00384*M, 0,
        0.16488*M, 0, -0.4266, 0.0096, 0, 0.08244*M, -0.4266, 0, 0.0048, 0,
        -0.00384*M, 0.0096, 0.0048, 0, 0, 0, 0, 0, 0, 0};
    const auto verify = [](const char* quantity, std::size_t row, std::size_t col,
                           double value, double reference, double scale, double tolerance) {
        const double error = (value-reference)*scale;
        std::printf("sulfate,%s,%zu,%zu,%.17g,%.17g,%.17g\n", quantity, row, col, value, reference, error);
        CHECK(std::abs(error) <= tolerance*std::max(1.0, std::abs(reference*scale)));
    };
    std::vector<double> inputs = {T, V}; inputs.insert(inputs.end(), n.begin(), n.end());
    std::vector<detail::AD> variables(inputs.begin(), inputs.end());
    CppAD::Independent(variables);
    CppAD::ADFun<double> isolated(variables, {cation_mixing_helmholtz(
        std::vector<detail::AD>(variables.begin()+2, variables.end()), sulfate)});
    const auto gradient = isolated.Jacobian(inputs), hessian = isolated.Hessian(inputs, 0);
    verify("F", 0, 0, isolated.Forward(0, inputs)[0], F, 1/std::abs(F), 1e-12);
    const auto on = sulfate_model.evaluate(T, V, n, DerivativeOrder::Hessian);
    const auto off = sulfate_plain.evaluate(T, V, n, DerivativeOrder::Hessian);
    verify("full-F", 0, 0, on.value-off.value, F, 1, 1e-10);
    for (std::size_t row = 0; row < 7; ++row) {
        const double reference = row < 2 ? 0 : g[row-2], scale = row < 2 ? 1 : scales[row-2];
        verify("gradient", row, 0, gradient[row], reference, scale, 1e-12);
        verify("full-gradient", row, 0, on.gradient[row]-off.gradient[row], reference, scale, 1e-10);
        if (row < 2 || row == 6) CHECK(gradient[row] == 0);
        for (std::size_t col = 0; col < 7; ++col) {
            const double expected = row < 2 || col < 2 ? 0 : H[(row-2)*5+col-2];
            const double factor = scale*(col < 2 ? 1 : scales[col-2]);
            verify("hessian", row, col, hessian[row*7+col], expected, factor, 1e-12);
            verify("full-hessian", row, col, on.hessian[row*7+col]-off.hessian[row*7+col], expected, factor, 1e-10);
            if (row < 2 || col < 2 || row == 6 || col == 6) CHECK(hessian[row*7+col] == 0);
        }
    }
    std::vector<double> composition = {T, 30000};
    for (const double amount : n) composition.push_back(amount/N);
    variables.assign(composition.begin(), composition.end());
    CppAD::Independent(variables);
    CppAD::ADFun<double> molar(variables, {cation_mixing_helmholtz(
        std::vector<detail::AD>(variables.begin()+2, variables.end()), sulfate)});
    const auto gx = molar.Jacobian(composition), hx = molar.Hessian(composition, 0);
    verify("molar", 0, 0, sulfate_model.state(T, 30000, std::vector<double>(composition.begin()+2, composition.end())).cation_mixing,
        F/N, N/std::abs(F), 1e-12);
    const std::vector<double> u = {0, 0, 0, 1, 0, -3, 2};
    double first = 0, second = 0;
    for (std::size_t row = 0; row < 7; ++row) {
        first += gx[row]*u[row];
        for (std::size_t col = 0; col < 7; ++col) {
            second += u[row]*hx[row*7+col]*u[col];
            if (row < 2 || col < 2) CHECK(hx[row*7+col] == 0);
        }
        if (row < 2) CHECK(gx[row] == 0);
    }
    verify("composition-first", 0, 0, first, -0.1764, 1, 1e-12);
    verify("composition-second", 0, 0, second, -0.0576*N, 1, 1e-12);
    for (std::size_t col = 0; col < 7; ++col) {
        double previous = std::numeric_limits<double>::infinity();
        for (const double step : {1e-3, 5e-4, 2.5e-4}) {
            auto plus = inputs, minus = inputs;
            const double scale = col < 2 ? inputs[col] : scales[col-2], h = step*scale;
            plus[col] += h; minus[col] -= h;
            const double fd = (isolated.Forward(0, plus)[0]-isolated.Forward(0, minus)[0])/(2*h);
            const double error = std::abs((fd-gradient[col])*scale);
            std::printf("sulfate,fd-gradient,%zu,%.17g,%.17g\n", col, step, error);
            if (col == 2) CHECK(error <= previous*0.4);
            if (step == 2.5e-4) CHECK(error <= 1e-7);
            previous = error;
            const auto gp = isolated.Jacobian(plus), gm = isolated.Jacobian(minus);
            for (std::size_t row = 0; row < 7; ++row) {
                const double e = std::abs(((gp[row]-gm[row])/(2*h)-hessian[row*7+col])*scale*(row < 2 ? inputs[row] : scales[row-2]));
                std::printf("sulfate,fd-hessian,%zu,%zu,%.17g,%.17g\n", row, col, step, e);
                if (step == 2.5e-4) CHECK(e <= 1e-7);
            }
        }
    }
    auto spectator = inputs; spectator[6] += 0.7;
    CHECK(isolated.Forward(0, spectator)[0] == isolated.Forward(0, inputs)[0]);
    CHECK(isolated.Jacobian(spectator) == gradient && isolated.Hessian(spectator, 0) == hessian);
    auto spectator_x = n; spectator_x[4] += 0.7;
    for (auto& x : spectator_x) x /= N+0.7;
    verify("spectator-molar", 0, 0, sulfate_model.state(T, 30000, spectator_x).cation_mixing, F/(N+0.7), 1, 1e-12);
    auto chloride_zero = inputs; chloride_zero[5] = 0; chloride_zero[6] = 0.5;
    verify("chloride-zero-F", 0, 0, isolated.Forward(0, chloride_zero)[0], -0.03528, 1/0.03528, 1e-12);
    verify("chloride-insertion", 0, 0, isolated.Jacobian(chloride_zero)[5], 0.00192, 1, 1e-12);
    std::vector<regression::Observation> observations;
    for (const auto& [a, b] : {std::pair{0.2, 0.4}, std::pair{1.2, 2.0}}) {
        const double c = a + 2 * b;
        const std::vector<double> n = {1 / mass, a, b, c};
        const double N = std::accumulate(n.begin(), n.end(), 0.0), V = N / 30000;
        const auto old = plain.evaluate(T, V, n, DerivativeOrder::Hessian);
        const auto now = mixed.evaluate(T, V, n, DerivativeOrder::Hessian);
        const double F = -0.4 * a * b + 0.02 * a * b * c;
        CHECK(close(now.value - old.value, F, 1e-12));
        std::printf("%.17g,%.17g,extensive-helmholtz-error,0,%.17g\n", a, b, now.value - old.value - F);
        const std::vector<double> mu = {-mass * (-0.4 * a * b + 0.04 * a * b * c),
            -0.4 * b + 0.02 * b * c, -0.4 * a + 0.02 * a * c, 0.02 * a * b};
        double euler = 0;
        for (std::size_t r = 0; r < 4; ++r) {
            CHECK(close(now.gradient[2 + r] - old.gradient[2 + r], mu[r], 1e-10));
            euler += n[r] * mu[r];
            std::printf("%.17g,%.17g,mu-%zu-error,0,%.17g\n", a, b, r, now.gradient[2+r]-old.gradient[2+r]-mu[r]);
        }
        CHECK(close(euler, F, 1e-10));
        std::printf("%.17g,%.17g,Euler-error,0,%.17g\n", a, b, euler-F);
        std::printf("%.17g,%.17g,pressure-increment,0,%.17g\n", a, b, now.gradient[1]-old.gradient[1]);
        CHECK(now.gradient[0] == old.gradient[0] && now.gradient[1] == old.gradient[1]);
        CHECK(close(now.gradient[6], 2 * a * b, 1e-12));
        CHECK(close(now.gradient[7], a * b * c, 1e-12));
        auto scaled = n;
        for (auto& v : scaled) v *= 3;
        const auto larger = mixed.evaluate(T, 3 * V, scaled, DerivativeOrder::Gradient);
        CHECK(close(larger.value, 3 * now.value, 1e-10));
        std::printf("%.17g,%.17g,homogeneity-error,0,%.17g\n", a, b, larger.value-3*now.value);
        CHECK(close(larger.gradient[0], 3 * now.gradient[0], 1e-10));
        for (std::size_t r = 1; r < 6; ++r) CHECK(close(larger.gradient[r], now.gradient[r], 1e-10));
        for (std::size_t col = 0; col < 8; ++col) {
            double previous = std::numeric_limits<double>::infinity();
            for (const double step : {1e-3, 5e-4, 2.5e-4}) {
                auto np = n, nm = n;
                auto pp = p, pm = p;
                double tp = T, tm = T, vp = V, vm = V;
                const double h = step * (col == 0 ? T : col == 1 ? V : col < 6 ? n[col - 2] : 1.0);
                if (col == 0) tp += h, tm -= h;
                else if (col == 1) vp += h, vm -= h;
                else if (col < 6) np[col - 2] += h, nm[col - 2] -= h;
                else {
                    auto family = col == 6 ? ParameterFamily::CationPairMixing : ParameterFamily::CationTripleMixing;
                    detail::family_values(pp, family).first[0] += h;
                    detail::family_values(pm, family).first[0] -= h;
                }
                const auto plus = HelmholtzModel(pp, active).evaluate(tp, vp, np, DerivativeOrder::Gradient);
                const auto minus = HelmholtzModel(pm, active).evaluate(tm, vm, nm, DerivativeOrder::Gradient);
                const double fd = (plus.value - minus.value) / (2 * h);
                const double error = std::abs(fd - now.gradient[col]) / std::max(1.0, std::abs(now.gradient[col]));
                CHECK(error <= 1e-6 || step != 2.5e-4);
                CHECK(error <= previous * 0.4 || error < 1e-8);
                previous = error;
                std::printf("%.17g,%.17g,gradient-%zu,%.17g,%.17g\n", a, b, col, step, error);
                for (std::size_t row = 0; row < 6; ++row) {
                    const double exact = now.hessian[row * 8 + col];
                    const double fd2 = (plus.gradient[row] - minus.gradient[row]) / (2 * h);
                    const double error2 = std::abs(fd2 - exact) / std::max(1.0, std::abs(exact));
                    std::printf("%.17g,%.17g,hessian-%zu-%zu,%.17g,%.17g\n", a, b, row, col, step, error2);
                    if (step == 2.5e-4) CHECK(error2 <= 1e-6);
                }
            }
        }
        for (const auto& indices : {std::vector<std::size_t>{3, 4, 6}, {3, 4, 5, 7}, {2, 3, 4, 5, 7}}) {
            std::vector<std::vector<double>> directions;
            for (const auto index : indices) {
                directions.emplace_back(8, 0.0);
                directions.back()[index] = 1;
            }
            const auto higher = mixed.residual_energy_action(T, V, n, directions);
            const double expected = indices.size() == 3 ? 2.0 : indices.size() == 4 ? 1.0 : -2 * mass;
            CHECK(higher.action && close(*higher.action, expected, 1e-10));
            std::printf("%.17g,%.17g,E%zu,0,%.17g\n", a, b, indices.size(), *higher.action - expected);
            const auto variable = indices.size() == 5 ? indices.front() : indices.back();
            auto lower = directions;
            if (indices.size() == 5) lower.erase(lower.begin()); else lower.pop_back();
            double previous_error = std::numeric_limits<double>::infinity();
            for (const double step : {1e-3, 5e-4, 2.5e-4}) {
                auto pp = p, pm = p; auto np = n, nm = n;
                const double h = step * (variable == 2 ? n[0] : 1.0);
                if (variable == 2) np[0] += h, nm[0] -= h;
                else {
                    const auto family = variable == 6 ? ParameterFamily::CationPairMixing : ParameterFamily::CationTripleMixing;
                    detail::family_values(pp, family).first[0] += h;
                    detail::family_values(pm, family).first[0] -= h;
                }
                const HelmholtzModel plus(pp, active), minus(pm, active);
                const double up = lower.size() == 2 ? plus.evaluate(T,V,np,DerivativeOrder::Hessian).hessian[3*8+4]
                    : plus.residual_energy_action(T,V,np,lower).action.value();
                const double down = lower.size() == 2 ? minus.evaluate(T,V,nm,DerivativeOrder::Hessian).hessian[3*8+4]
                    : minus.residual_energy_action(T,V,nm,lower).action.value();
                const double error = std::abs((up-down)/(2*h)-*higher.action) / std::max(1.0,std::abs(*higher.action));
                std::printf("%.17g,%.17g,E%zu-refinement,%.17g,%.17g\n",a,b,indices.size(),step,error);
                CHECK(error <= previous_error*0.4 || error < 1e-8);
                if (step == 2.5e-4) CHECK(error <= 1e-6);
                previous_error = error;
            }
        }
        std::vector<double> x = n;
        for (auto& v : x) v /= N;
        const auto s0 = plain.state_at_pressure(T, 1e5, x, DensityBranch::Liquid);
        const auto s1 = mixed.state_at_pressure(T, 1e5, x, DensityBranch::Liquid);
        CHECK(s0.molar_density == s1.molar_density && s0.pressure == s1.pressure);
        CHECK(close(s1.cation_mixing, F / N, 1e-12));
        regression::Observation ob;
        ob.kind = regression::ObservationKind::MixtureOsmoticCoefficient;
        ob.temperature = T; ob.pressure = 1e5; ob.composition = x; ob.branch = DensityBranch::Liquid;
        const auto y0 = regression::predict(plain, ob), y1 = regression::predict(mixed, ob);
        std::printf("prediction-causes,%s,%s\n", y0.message.c_str(), y1.message.c_str());
        CHECK(y0.value && y1.value);
        CHECK(close(*y1.value - *y0.value, 2 * a * b * (-0.2 + 0.02 * c) / (a + b + c), 1e-12));
        std::printf("%.17g,%.17g,Phi-error,0,%.17g\n", a, b, *y1.value-*y0.value-2*a*b*(-0.2+0.02*c)/(a+b+c));
        ob.observed = *y1.value;
        observations.push_back(ob);
        for (const std::size_t absent : {1U, 2U}) {
            auto face = n; face[absent] = 0; face[3] = face[1] + 2*face[2];
            const auto f0 = plain.evaluate(T, V, face, DerivativeOrder::Hessian);
            const auto f1 = mixed.evaluate(T, V, face, DerivativeOrder::Hessian);
            CHECK(f0.value == f1.value);
            auto xf = face;
            const double nf = std::accumulate(xf.begin(), xf.end(), 0.0);
            for (auto& v : xf) v /= nf;
            const auto sf0 = plain.state_at_pressure(T, 1e5, xf, DensityBranch::Liquid);
            const auto sf1 = mixed.state_at_pressure(T, 1e5, xf, DensityBranch::Liquid);
            CHECK(sf0.molar_density == sf1.molar_density && sf0.pressure == sf1.pressure);
            auto face_ob = ob; face_ob.composition = xf;
            const auto pf0 = regression::predict(plain, face_ob), pf1 = regression::predict(mixed, face_ob);
            CHECK(pf0.value && pf1.value && *pf0.value == *pf1.value);
            std::printf("%.17g,%.17g,face-%zu-Phi-difference,0,%.17g\n", a, b, absent, *pf1.value-*pf0.value);
            for (std::size_t r = 0; r < 6; ++r) if (r != absent + 2) CHECK(f0.gradient[r] == f1.gradient[r]);
            for (std::size_t row = 0; row < 6; ++row) for (std::size_t col = 0; col < 6; ++col)
                if (row != absent + 2 && col != absent + 2)
                    CHECK(f0.hessian[row * 6 + col] == f1.hessian[row * 8 + col]);
            CHECK(close(f1.gradient[absent + 2] - f0.gradient[absent + 2], face[3-absent]*(-0.4 + 0.02*face[3]), 1e-10));
        }
        auto zero = p; zero.cation_pair_mixing[0] = zero.cation_triple_mixing[0] = 0;
        const auto z = HelmholtzModel(zero, active).evaluate(T, V, n, DerivativeOrder::Hessian);
        CHECK(z.value == old.value && z.gradient[6] != 0 && z.gradient[7] != 0);
        for (std::size_t row = 0; row < 6; ++row) {
            CHECK(z.gradient[row] == old.gradient[row]);
            for (std::size_t col = 0; col < 6; ++col) CHECK(z.hessian[row*8+col] == old.hessian[row*6+col]);
        }
        std::printf("%.17g,%.17g,zero-start-theta-column,0,%.17g\n", a, b, z.gradient[6]);
        std::printf("%.17g,%.17g,zero-start-psi-column,0,%.17g\n", a, b, z.gradient[7]);
    }
    for (const bool triple : {false, true}) {
        auto isolated = p;
        (triple ? isolated.cation_pair_mixing : isolated.cation_triple_mixing)[0] = 0;
        double previous = 0, previous_phi = 0, previous_solvent = 0, previous_ion = 0;
        for (const double scale : {1e-2, 1e-3, 1e-4}) {
            const std::vector<double> n = {1 / mass, 0.2 * scale, 0.4 * scale, scale};
            const double F = cation_mixing_helmholtz(n, isolated);
            std::vector<detail::AD> amounts(n.begin(), n.end());
            CppAD::Independent(amounts);
            std::vector<detail::AD> output = {cation_mixing_helmholtz(amounts, isolated)};
            CppAD::ADFun<double> tape(amounts, output);
            const auto mu = tape.Jacobian(n);
            const double phi = -mu[0] / (mass * 1.6 * scale);
            if (previous != 0) CHECK(close(F / previous, triple ? 1e-3 : 1e-2, 1e-12));
            if (previous_phi != 0) {
                CHECK(close(phi / previous_phi, triple ? 1e-2 : 1e-1, 1e-12));
                CHECK(close(mu[0] / previous_solvent, triple ? 1e-3 : 1e-2, 1e-12));
                CHECK(close(mu[1] / previous_ion, triple ? 1e-2 : 1e-1, 1e-12));
            }
            std::printf("0,0,%s-dilute-F,%.17g,%.17g\n", triple ? "triple" : "pair", scale, F);
            std::printf("0,0,%s-dilute-Phi,%.17g,%.17g\n", triple ? "triple" : "pair", scale, phi);
            std::printf("0,0,%s-dilute-solvent-mu,%.17g,%.17g\n", triple ? "triple" : "pair", scale, mu[0]);
            std::printf("0,0,%s-dilute-ion-mu,%.17g,%.17g\n", triple ? "triple" : "pair", scale, mu[1]);
            previous = F;
            previous_phi = phi; previous_solvent = mu[0]; previous_ion = mu[1];
        }
    }
    const std::vector<double> pure = {1,0,0,0};
    const auto pure_reference = solvent_reference(mixed,T,1e5,pure);
    const auto pure_activity = solvent_activity(mixed,T,1e5,pure,pure_reference);
    CHECK(pure_activity.solvent_activity && *pure_activity.solvent_activity == 1.0
        && pure_activity.osmotic_cause == "OsmoticUndefinedPureSolvent");
    std::printf("0,0,pure-water-activity,0,%.17g\n",*pure_activity.solvent_activity);
    regression::FitRequest request;
    request.parameters = p;
    request.parameters.cation_pair_mixing[0] = request.parameters.cation_triple_mixing[0] = 0;
    for (const auto& target : active) {
        regression::Coordinate coordinate;
        coordinate.targets = {{target.family, target.indices, 1.0}};
        coordinate.origin = 0; coordinate.scale = 1;
        request.coordinates.coordinates.push_back(coordinate);
    }
    request.observations = observations;
    request.controls.maximum_iterations = 30; request.controls.maximum_elapsed_time_seconds = 30;
    request.controls.workers = 1;
    const auto recovered = regression::fit(request);
    CHECK(recovered.usable);
    CHECK(close(recovered.parameters.cation_pair_mixing[0], -0.2, 1e-8));
    CHECK(close(recovered.parameters.cation_triple_mixing[0], 0.02, 1e-8));
    std::printf("recovered,theta,%.17g,psi,%.17g\n", recovered.parameters.cation_pair_mixing[0], recovered.parameters.cation_triple_mixing[0]);
}
}  // namespace properties_checks
