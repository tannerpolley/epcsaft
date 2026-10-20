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
