// P2: association actions.
#include "properties_checks_shared.hpp"

namespace properties_checks {

enum class P2AssociationDirection { Temperature, Volume, AssociationEnergy };

struct P2AssociationAction final {
    double residual = 0.0;
    std::size_t rank = 0;
    double pivot_ratio = 1.0;
};

P2AssociationAction p2_association_action(
    double temperature, double volume, const std::vector<double>& composition,
    const epcsaft::Parameters& parameters, const std::vector<double>& sites,
    P2AssociationDirection direction
) {
    const auto conditioning = p2_association_conditioning(
        temperature, 1.0 / volume, composition, parameters, sites
    );
    using AD = CppAD::AD<double>;
    const double coordinate = direction == P2AssociationDirection::Temperature
        ? temperature
        : direction == P2AssociationDirection::Volume
            ? volume : parameters.association_energy_over_k[1];
    std::vector<AD> independent = {AD(coordinate)};
    CppAD::Independent(independent);
    const AD temperature_ad = direction == P2AssociationDirection::Temperature
        ? independent.front() : AD(temperature);
    const AD volume_ad = direction == P2AssociationDirection::Volume
        ? independent.front() : AD(volume);
    std::vector<AD> composition_ad(composition.begin(), composition.end());
    std::vector<AD> sites_ad(sites.begin(), sites.end());
    std::vector<epcsaft::detail::ActiveParameter> active;
    if (direction == P2AssociationDirection::AssociationEnergy) {
        active.push_back({
            epcsaft::detail::ParameterFamily::AssociationEnergyOverK, {1, 4}
        });
    }
    const auto tape_parameters = epcsaft::detail::make_tape_parameters(
        parameters, active, independent
    );
    const auto sigma = epcsaft::segment_sigma(temperature_ad, tape_parameters);
    const auto diameter = epcsaft::segment_diameter(
        temperature_ad, tape_parameters, sigma
    );
    const auto density_ad = AD(1.0) / volume_ad;
    const auto delta = epcsaft::association_delta_values(
        temperature_ad, density_ad, composition_ad, tape_parameters, sigma, diameter,
        tape_parameters.association_energy_over_k, tape_parameters.association_volume
    );
    const auto residuals = epcsaft::association_residuals(
        density_ad, composition_ad, tape_parameters, sites_ad, delta
    );
    CppAD::ADFun<double> function(independent, residuals);
    const auto jacobian = function.Jacobian(std::vector<double>{coordinate});
    Eigen::VectorXd f_direction(static_cast<Eigen::Index>(sites.size()));
    for (std::size_t site = 0; site < sites.size(); ++site) {
        f_direction(static_cast<Eigen::Index>(site)) =
            jacobian[site];
    }
    Eigen::FullPivLU<Eigen::MatrixXd> factorization(conditioning.jacobian);
    const Eigen::VectorXd x_direction = -factorization.solve(f_direction);
    const double residual = (conditioning.jacobian * x_direction + f_direction)
        .lpNorm<Eigen::Infinity>();
    return {residual, conditioning.rank, conditioning.pivot_ratio};
}

double p2_pressure_slope_floor_reference(
    double temperature, double volume, const std::vector<double>& composition,
    const epcsaft::Parameters& parameters, const std::vector<double>& sites
) {
    const double density = 1.0 / volume;
    const auto conditioning = p2_association_conditioning(
        temperature, density, composition, parameters, sites
    );
    if (sites.empty()) return 0.0;
    const double step = volume * 1.0e-6;
    const auto plus = p2_association_residuals(
        temperature, volume + step, composition, parameters, sites
    );
    const auto minus = p2_association_residuals(
        temperature, volume - step, composition, parameters, sites
    );
    Eigen::VectorXd potential_cross(static_cast<Eigen::Index>(sites.size()));
    for (std::size_t site = 0; site < sites.size(); ++site) {
        const auto component = static_cast<std::size_t>(parameters.site_component[site]);
        const double mass = composition[component] * parameters.site_multiplicity[site];
        const double f_volume = (plus[site] - minus[site]) / (2.0 * step);
        potential_cross(static_cast<Eigen::Index>(site)) =
            -mass / sites[site] * f_volume;
    }
    const Eigen::VectorXd ones =
        Eigen::VectorXd::Ones(static_cast<Eigen::Index>(sites.size()));
    Eigen::FullPivLU<Eigen::MatrixXd> factorization(conditioning.jacobian);
    // The floor propagates the achieved site residual, bounded below by eps per site.
    double achieved = std::numeric_limits<double>::epsilon() * static_cast<double>(sites.size());
    for (const double value : p2_association_residuals(temperature, volume, composition, parameters, sites))
        achieved = std::max(achieved, std::abs(value));
    return std::abs(potential_cross.dot(factorization.solve(ones))) * achieved;
}

struct P2PressureChain final {
    double measured = 0.0;
    double predicted = 0.0;
    double pressure_density_derivative = 0.0;
    double reduced_residual = 0.0;
    std::size_t stable_root_count = 0;
    epcsaft::DensityBranch branch = epcsaft::DensityBranch::Unspecified;
    std::size_t association_rank = 0;
    double association_pivot_ratio = 1.0;
};

P2PressureChain p2_fixed_pressure_chain(
    const epcsaft::HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, epcsaft::DensityBranch branch,
    double step, bool enthalpy
) {
    const auto central = model.state_at_pressure(
        temperature, pressure, composition, branch
    );
    const auto plus = model.state_at_pressure(
        temperature + step, pressure, composition, branch
    );
    const auto minus = model.state_at_pressure(
        temperature - step, pressure, composition, branch
    );
    const auto fixed_plus = model.state(
        temperature + step, central.molar_density, composition
    );
    const auto fixed_minus = model.state(
        temperature - step, central.molar_density, composition
    );
    const double pressure_temperature =
        (fixed_plus.pressure - fixed_minus.pressure) / (2.0 * step);
    const double density_temperature =
        -pressure_temperature / central.pressure_density_derivative;
    const double density_step = central.molar_density * 1.0e-5;
    const auto density_plus = model.state(
        temperature, central.molar_density + density_step, composition
    );
    const auto density_minus = model.state(
        temperature, central.molar_density - density_step, composition
    );
    const auto quantity = [enthalpy](const epcsaft::State& state) {
        return enthalpy ? state.residual_enthalpy : state.molar_density;
    };
    const double quantity_temperature = quantity(fixed_plus) - quantity(fixed_minus);
    const double quantity_density = quantity(density_plus) - quantity(density_minus);
    const double predicted = quantity_temperature / (2.0 * step)
        + quantity_density / (2.0 * density_step) * density_temperature;
    const double measured = (quantity(plus) - quantity(minus)) / (2.0 * step);
    const auto conditioning = p2_association_conditioning(
        temperature, central.molar_density, composition,
        model.parameters(), central.site_fractions
    );
    return {
        measured, predicted, central.pressure_density_derivative,
        central.density_convergence.has_value()
            ? central.density_convergence->residual : std::numeric_limits<double>::infinity(),
        central.stable_root_count.value_or(0), central.density_branch,
        conditioning.rank, conditioning.pivot_ratio
    };
}

// ALG-DENSITY: independent tape derivatives and verification-only log-volume stencils.
void p2_density_series() {
    struct Case { const char* name; epcsaft::Parameters parameters; double temperature, density; std::vector<double> composition; };
    auto density_temperature = properties_support::aqueous_electrolyte();
    density_temperature.permittivity_form = epcsaft::PermittivityForm::DensityTemperature;
    density_temperature.dielectric_volume_m3_per_mol = {3.777e-4, 0.0, 0.0};
    density_temperature.dielectric_temperature_k = {1403.0, 0.0, 0.0};
    auto dipolar = properties_support::neutral_binary();
    dipolar.dipole_moment_debye = {1.5, 0.0};
    dipolar.dipole_count = {1.0, 0.0};
    const std::vector<Case> cases = {
        {"pure-neutral", properties_support::ethane(), 280.0, 800.0, {1.0}},
        {"neutral", properties_support::neutral_binary(), 250.0, 800.0, {0.5, 0.5}},
        {"association", properties_support::associating_binary(), 300.0, 35000.0, {0.1, 0.9}},
        {"absent-species", properties_support::associating_binary(), 300.0, 35000.0, {0.0, 1.0}},
        {"electrolyte", properties_support::aqueous_electrolyte(), 298.15, 55000.0, {0.98, 0.01, 0.01}},
        {"trace-ion", properties_support::aqueous_electrolyte(), 298.15, 55000.0, {1.0-2e-10, 1e-10, 1e-10}},
        {"density-permittivity", density_temperature, 298.15, 55000.0, {0.98, 0.01, 0.01}},
        {"dipolar", dipolar, 250.0, 800.0, {0.5, 0.5}}};
    for (const auto& item : cases) {
        epcsaft::HelmholtzModel model(item.parameters);
        const double volume = 1.0/item.density;
        double floor = 0.0;
        const auto series = epcsaft::properties::density_volume_taylor(model, item.temperature, item.density, item.composition, floor);
        const auto tape = model.evaluate(item.temperature, volume, item.composition, epcsaft::DerivativeOrder::Hessian);
        const auto smaller = model.evaluate(item.temperature, volume, item.composition, epcsaft::DerivativeOrder::PressureSlope);
        std::vector<double> direction(item.composition.size()+2, 0.0); direction[1] = 1.0;
        const auto third = model.residual_energy_action(item.temperature, volume, item.composition, {direction,direction,direction});
        CHECK(third.action.has_value());
        const std::array<double,3> exact = {volume*tape.gradient[1], volume*volume*tape.hessian[tape.gradient.size()+1],
            volume*volume*volume*third.action.value_or(std::numeric_limits<double>::quiet_NaN())};
        const std::array<double,3> calculated = {volume*series.c[1], 2*volume*volume*series.c[2], 6*volume*volume*volume*series.c[3]};
        const auto inputs = recording_support::fields({model_inputs(item.parameters,item.temperature,volume,item.composition),
            recording_support::field("site_fractions",tape.site_fractions), recording_support::field("series_gradient_floor",floor),
            recording_support::field("tape_gradient_floor",tape.gradient_floor)});
        for (std::size_t k=0;k<3;++k) {
            recording_support::compare((std::string("volume-series/")+item.name+"/"+std::to_string(k+1)).c_str(),
                "scaled_volume_derivative","1",calculated[k],exact[k],inputs,"Identity","independent-tape",
                "normalized Taylor coefficients agree with implicit tape differentiation", "reference_or_one",0.0,1e-10,
                close(calculated[k],exact[k],1e-10),"p2-volume-series",item.name,"V^k Phi_V...V","T,V,n","series=tape","order");
            double best = std::numeric_limits<double>::infinity(), best_step = 0.0, previous = best;
            bool decreasing = false;
            std::vector<double> steps, errors;
            for (double h : {1e-2,3e-3,1e-3,3e-4,1e-4,3e-5}) {
                const auto plus=model.evaluate(item.temperature,volume*std::exp(h),item.composition,epcsaft::DerivativeOrder::Hessian);
                const auto minus=model.evaluate(item.temperature,volume*std::exp(-h),item.composition,epcsaft::DerivativeOrder::Hessian);
                const auto lower = [&](const epcsaft::HelmholtzDerivatives& d) {
                    return k==0 ? d.value : k==1 ? d.gradient[1] : d.hessian[d.gradient.size()+1]; };
                const double measured=std::pow(volume,static_cast<int>(k))*(lower(plus)-lower(minus))/(2*h);
                const double error=std::abs(measured-exact[k])/std::max(1.0,std::abs(exact[k]));
                decreasing = decreasing || (std::isfinite(previous) && std::isfinite(error) && error < previous); previous=error;
                steps.push_back(h); errors.push_back(error);
                if(error<best) { best=error; best_step=h; }
            }
            const double limit=k==2 ? 1e-5 : 1e-6;
            recording_support::compare((std::string("volume-stencil/")+item.name+"/"+std::to_string(k+1)).c_str(),
                "scaled_error","1",best,0.0,recording_support::fields({inputs,recording_support::field("log_volume_step",best_step),
                    recording_support::field("log_volume_steps",steps), recording_support::field("scaled_errors",errors)}),
                "Relationship","symmetric-log-volume-stencil","successive lower-derivative differences verify the exact volume derivative",
                "absolute",limit,0.0,decreasing&&best<=limit,"p2-volume-stencil",item.name,"scaled error","T,V,n","d/dlnV=V d/dV","order");
        }
        const std::array<double,3> series_pressure={1-calculated[0],1+calculated[1],-calculated[2]-2*calculated[1]};
        const std::array<double,3> tape_pressure={1-volume*smaller.gradient[1],1+volume*volume*smaller.volume_second_derivative,-exact[2]-2*exact[1]};
        if (item.composition.size()==1) {
            const auto pressure=model.neutral_property(item.temperature,item.density,epcsaft::NeutralProperty::Pressure,epcsaft::DerivativeOrder::Hessian);
            const double curvature=item.density*pressure.hessian[3]/(epcsaft::GAS_CONSTANT_J_PER_MOL_K*item.temperature);
            recording_support::compare("pressure-hessian/pure-neutral", "scaled_pressure_curvature", "1",
                series_pressure[2], curvature, inputs, "Identity", "neutral-pressure-Hessian",
                "the exact pressure Hessian agrees with the volume-series curvature", "reference_or_one", 0.0, 1e-10,
                close(series_pressure[2],curvature,1e-10), "p2-pressure-hessian", item.name,
                "rho P_rhorho/(RT)", "T,rho", "series=pressure Hessian", "case");
        }
        for(std::size_t k=0;k<3;++k)
            recording_support::compare((std::string("pressure-series/")+item.name+"/"+std::to_string(k)).c_str(),
                "scaled_pressure_derivative","1",series_pressure[k],tape_pressure[k],inputs,"Identity","pressure-chain",
                "pressure, slope and curvature follow the volume derivative chain", "reference_or_one",0.0,1e-10,
                close(series_pressure[k],tape_pressure[k],1e-10),"p2-pressure-series",item.name,"scaled P derivative","T,V,n","pressure-chain","order");
    }
}

// The finer native diagnostic resolves two positive-slope roots on this isotherm.
void p2_decane_fallback() {
    auto parameters=properties_support::parameters(1);
    parameters.segment_count={4.6627}; parameters.sigma_constant={3.8384};
    parameters.dispersion_energy_over_k={243.87}; parameters.molar_mass_kg_per_mol={0.142286};
    epcsaft::HelmholtzModel model(parameters);
    const std::vector<double> composition{1.0};
    constexpr double reference=0.9036845021232383, dense_reference=6174.3912350690825;
    CHECK(model.state(135.0,dense_reference,composition).pressure_density_derivative>0.0);
    for(const auto branch : {epcsaft::DensityBranch::Vapor,epcsaft::DensityBranch::Liquid}) {
        for (const bool request_sibling : {false, true}) {
            epcsaft::core::Capture capture(true);
            epcsaft::DensityRoot sibling;
            const auto root=epcsaft::solve_density(model,135.0,1000.0,composition,branch,std::nullopt,request_sibling ? &sibling : nullptr);
            const auto trace=capture.take();
            bool fallback=false;
            for(const auto& e:trace.events) fallback=fallback||(std::string(e.name)=="structured_broken"&&e.reason=="curvature or slope order");
            const double expected = branch == epcsaft::DensityBranch::Vapor ? reference : dense_reference;
            const double other = branch == epcsaft::DensityBranch::Vapor ? dense_reference : reference;
            const bool ok=close(root.molar_density,expected,1e-10)&&root.stable_root_count==2
                &&root.branch==branch && root.pressure_density_derivative > 0.0
                &&(!request_sibling || (sibling.branch!=root.branch && sibling.stable_root_count==2
                    && sibling.pressure_density_derivative > 0.0 && close(sibling.molar_density,other,1e-10)))&&fallback;
            recording_support::compare(branch==epcsaft::DensityBranch::Vapor
                    ? (request_sibling ? "decane/vapor-sibling" : "decane/vapor")
                    : (request_sibling ? "decane/liquid-sibling" : "decane/liquid"),
                "molar_density","mol/m3",root.molar_density,expected,recording_support::fields({
                    model_inputs(parameters,135.0,1.0/root.molar_density,composition),recording_support::field("P_star",1000.0),
                    recording_support::field("dense_diagnostic_root",dense_reference),recording_support::field("stable_root_count",root.stable_root_count.value_or(0)),
                    recording_support::field("fallback",fallback), recording_support::field("request_sibling",request_sibling)}),"Relationship","current-main-scan",
                "both branches and the actual sibling recover the diagnostic roots", "reference_or_one",0.0,1e-10,ok,
                "p2-decane",branch==epcsaft::DensityBranch::Vapor ? "vapor" : "liquid","rho","T,P,x,branch","resolved positive-slope roots","case");
        }
    }
    std::string cause;
    try { epcsaft::solve_density(model,135.0,1000.0,composition,epcsaft::DensityBranch::Unspecified,std::nullopt); }
    catch (const epcsaft::DensityStructureError& error) { cause = error.what(); }
    CHECK(cause.find("ambiguous") != std::string::npos);
    const auto batch = model.pressure_actions(135.0,1000.0,composition,epcsaft::DensityBranch::Unspecified, {},
        {epcsaft::PressureActionDirection{0.0,1.0,{0.0},{}}});
    CHECK(batch.density_structure_refused && batch.density_status == epcsaft::PropertyActionStatus::DensityUnavailable
        && !batch.molar_density && batch.density_message == cause);

    // #513 R1: the refined crossing can publish a nonzero-residual interval endpoint.
    constexpr double endpoint_pressure = 118697451.59547403, endpoint_reference = 6744.38451784218;
    for (const auto branch : {epcsaft::DensityBranch::Vapor, epcsaft::DensityBranch::Liquid}) {
        const auto root = epcsaft::solve_density(model, 135.0, endpoint_pressure, composition, branch, std::nullopt);
        recording_support::compare(branch == epcsaft::DensityBranch::Vapor ? "decane/endpoint-vapor" : "decane/endpoint-liquid",
            "molar_density", "mol/m3", root.molar_density, endpoint_reference,
            recording_support::fields({model_inputs(parameters, 135.0, 1.0 / root.molar_density, composition),
                recording_support::field("P_star", endpoint_pressure), recording_support::field("stable_root_count", root.stable_root_count.value_or(0))}),
            "Relationship", "native-4096-and-16384-cell-diagnostics", "the resolved crossing survives the local and ordinary root merge",
            "reference_or_one", 0.0, 1e-10, close(root.molar_density, endpoint_reference, 1e-10)
                && root.stable_root_count == 1 && root.branch == epcsaft::DensityBranch::Unspecified
                && root.pressure_density_derivative > 0.0,
            "p2-decane", "endpoint", "rho", "T,P,x,branch", "one positive-slope root", "case");
    }
}

void properties_p2() {
    p2_density_series();
    p2_decane_fallback();

    const auto associating_parameters = properties_support::associating_binary();
    epcsaft::HelmholtzModel associating(associating_parameters);
    {
        const double temperature = 330.0, volume = 6e-5;
        const std::vector<double> composition{0.3, 0.7};
        const auto at = associating.evaluate(temperature, volume, composition, epcsaft::DerivativeOrder::Hessian);
        const auto errors = associating.residual_coefficient_errors(
            {{temperature, volume, 0.3, 0.7}, {0.0, 1.0, 0.0, 0.0}}, at.site_fractions);
        const bool covered = errors[1] >= at.gradient_floor;
        recording_support::compare_custom("coefficient-volume-floor", "volume_derivative_error", "mol/m^3",
            errors[1], at.gradient_floor,
            recording_support::fields({model_inputs(associating_parameters, temperature, volume, composition)}),
            "Relationship", "achieved-site-residual", "coefficient errors retain the achieved site-residual contribution",
            "absolute", covered, "coefficient-volume-floor", "volume", "error", "T,V,n",
            "error>=achieved-site-residual-floor", "case");
        CHECK(covered);
    }
    const auto record_zero_support = [&](const char* point,
                                         const epcsaft::HelmholtzModel& model,
                                         double temperature_value, double density_value,
                                         const std::vector<double>& composition) {
        const auto state = model.state(temperature_value, density_value, composition);
        const auto residuals = p2_association_residuals(
            temperature_value, 1.0 / density_value, composition,
            model.parameters(), state.site_fractions
        );
        double residual = 0.0;
        bool finite = finite_state_without_cp(state);
        for (const double value : state.site_fractions) finite = finite && std::isfinite(value);
        for (const double value : residuals) residual = std::max(residual, std::abs(value));
        recording_support::compare_custom(
            (std::string("support-residual/") + point).c_str(),
            "association_residual_infinity", "1", residual, 0.0,
            recording_support::fields({
                model_inputs(model.parameters(), temperature_value, 1.0 / density_value, composition),
                recording_support::field("support", point),
                recording_support::field("site_fractions", state.site_fractions),
            }), "Relationship", "mass-action",
            "exact neutral/zero-ion support remains finite and converged", "absolute",
            residual <= 2e-13 && finite, "p2-support-residual", point,
            "||F_a||_inf", "T,rho,x", "||F_a||_inf<=2e-13", "support");
        CHECK(finite && residual <= 2e-13);
    };
    const std::array<std::vector<double>, 2> supports = {{{1.0, 0.0}, {0.0, 1.0}}};
    for (std::size_t support = 0; support < supports.size(); ++support) {
        const auto& composition = supports[support];
        constexpr double temperature = 300.0;
        constexpr double density = 5000.0;
        const auto state = associating.state(temperature, density, composition);
        const auto residuals = p2_association_residuals(
            temperature, 1.0 / density, composition,
            associating.parameters(), state.site_fractions
        );
        double residual = 0.0;
        bool finite = finite_state_without_cp(state);
        for (const double value : state.site_fractions) finite = finite && std::isfinite(value);
        for (const double value : residuals) residual = std::max(residual, std::abs(value));
        const auto pressure_slope = associating.evaluate(
            temperature, 1.0 / density, composition, epcsaft::DerivativeOrder::PressureSlope
        );
        const double expected_floor = p2_pressure_slope_floor_reference(
            temperature, 1.0 / density, composition,
            associating.parameters(), state.site_fractions
        );
        // The PressureSlope adjoint gradient and the Gradient-order sensitivity gradient are
        // algebraically distinct formulas for the same reduced gradient.
        const auto gradient_order = associating.evaluate(
            temperature, 1.0 / density, composition, epcsaft::DerivativeOrder::Gradient);
        double gradient_error = 0.0;
        for (std::size_t i = 0; i < gradient_order.gradient.size(); ++i)
            gradient_error = std::max(gradient_error, std::abs(pressure_slope.gradient[i] - gradient_order.gradient[i])
                / std::max(1.0, std::abs(gradient_order.gradient[i])));
        CHECK(pressure_slope.gradient.size() == gradient_order.gradient.size() && gradient_error <= 1e-12);
        const auto inputs = recording_support::fields({
            model_inputs(associating.parameters(), temperature, 1.0 / density, composition),
            recording_support::field("support", support),
            recording_support::field("pressure_slope_gradient_identity_error", gradient_error),
            recording_support::field("site_fractions", state.site_fractions),
            recording_support::field("association_floor", pressure_slope.gradient_floor),
            recording_support::field("expected_floor", expected_floor),
        });
        const auto point = std::to_string(support);
        recording_support::compare_custom(
            (std::string("support-residual/") + point).c_str(),
            "association_residual_infinity", "1", residual, 0.0, inputs,
            "Relationship", "mass-action", "exact support remains finite and converged",
            "absolute", residual <= 2e-13 && finite,
            "p2-support-residual", point.c_str(), "||F_a||_inf", "T,rho,x",
            "||F_a||_inf<=2e-13", "support");
        CHECK(finite && residual <= 2e-13);
        const double floor_scale = std::max(expected_floor, 1.0e-300);
        const double floor_error = std::abs(pressure_slope.gradient_floor - expected_floor);
        const bool floor_ok = floor_error <= 1.0e-6 * floor_scale;
        recording_support::compare_custom(
            (std::string("support-floor/") + point).c_str(), "gradient_floor", "1",
            pressure_slope.gradient_floor, expected_floor, inputs, "Relationship",
            "mass-action-residual-bound", "PressureSlope floor uses U^{-1} one",
            "relative", floor_ok, "p2-support-floor", point.c_str(), "floor",
            "T,rho,x", "floor=bound(U^{-1}1)", "support");
    }
    record_zero_support(
        "aqueous-zero-ion", epcsaft::HelmholtzModel(properties_support::aqueous_electrolyte()),
        298.15, 5000.0, {1.0, 0.0, 0.0}
    );
    record_zero_support(
        "nine-zero-ion", epcsaft::HelmholtzModel(properties_support::nine_species()),
        313.15, 5000.0, {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0}
    );

    const std::array<double, 5> association_trace = {
        1e-2, 1e-3, 1e-4, 1e-5, 0.0
    };
    std::array<double, 5> trace_to_zero_errors{};
    const auto zero_derivatives = associating.evaluate(
        300.0, 1.0 / 5000.0, {1.0, 0.0}, epcsaft::DerivativeOrder::Hessian
    );
    for (std::size_t trace_index = 0; trace_index < association_trace.size(); ++trace_index) {
        const double trace = association_trace[trace_index];
        const std::vector<double> composition = {1.0 - trace, trace};
        const auto state = associating.state(300.0, 5000.0, composition);
        const auto derivatives = associating.evaluate(
            300.0, 1.0 / 5000.0, composition, epcsaft::DerivativeOrder::Hessian
        );
        const auto reference = properties_support::reference_tensor(
            associating.parameters(), 300.0, 1.0 / 5000.0, composition,
            state.site_fractions
        );
        double hessian_error = 0.0;
        double trace_to_zero_error = 0.0;
        for (std::size_t index = 0; index < derivatives.hessian.size(); ++index) {
            hessian_error = std::max(
                hessian_error,
                std::abs(derivatives.hessian[index] - reference.hessian[index])
                    / std::max(1.0, std::abs(reference.hessian[index]))
            );
            trace_to_zero_error = std::max(
                trace_to_zero_error,
                std::abs(derivatives.hessian[index] - zero_derivatives.hessian[index])
                    / std::max(1.0, std::abs(zero_derivatives.hessian[index]))
            );
        }
        trace_to_zero_errors[trace_index] = trace_to_zero_error;
        const auto point = std::to_string(trace);
        recording_support::compare(
            (std::string("support-hessian/") + point).c_str(), "relative_hessian_error", "1",
            hessian_error, 0.0,
            recording_support::fields({
                model_inputs(associating.parameters(), 300.0, 1.0 / 5000.0, composition),
                recording_support::field("trace_amount", trace),
                recording_support::field("trace_to_zero_error", trace_to_zero_error),
                recording_support::field("scaled_trace_to_zero_error",
                    trace == 0.0 ? 0.0 : trace_to_zero_error / trace),
            }), "Relationship", "AD-reference",
            "association implicit Hessian remains support-reduced", "absolute", 1e-12, 0.0,
            hessian_error <= 1e-12, "p2-support-hessian", point.c_str(), "Hres",
            "T,rho,x", "Himplicit=Hreference", "trace");
    }
    const bool trace_refines = trace_to_zero_errors[0] >= trace_to_zero_errors[1]
        && trace_to_zero_errors[1] >= trace_to_zero_errors[2]
        && trace_to_zero_errors[2] >= trace_to_zero_errors[3]
        && trace_to_zero_errors[3] >= trace_to_zero_errors[4];
    recording_support::compare_custom(
        "support-hessian/trace-to-zero", "refinement", "bool",
        trace_refines ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("trace_errors", std::vector<double>(
                trace_to_zero_errors.begin(), trace_to_zero_errors.end())),
            recording_support::field("terminal_error", trace_to_zero_errors.back()),
            recording_support::field("scaled_terminal_error",
                trace_to_zero_errors[3] / association_trace[3]),
        }), "Relationship", "scaled-trace-ladder",
        "trace-to-zero Hessian error refines to the exact-zero row", "absolute", trace_refines,
        "p2-support-hessian", "trace-to-zero", "refinement", "trace",
        "error(trace)<=error(previous);error(0)=0", "trace");

    const std::vector<double> action_composition = {0.4, 0.6};
    const std::array<P2AssociationDirection, 3> directions = {
        P2AssociationDirection::Temperature,
        P2AssociationDirection::Volume,
        P2AssociationDirection::AssociationEnergy
    };
    const std::array<const char*, 3> direction_names = {
        "temperature", "volume", "association-energy"
    };
    const std::array<double, 3> direction_steps = {
        1e-3, 1.0 / 5000.0 * 1e-3, 1e-3
    };
    const auto action_state = associating.state(300.0, 5000.0, action_composition);
    for (std::size_t index = 0; index < directions.size(); ++index) {
        std::vector<epcsaft::ActiveParameter> active;
        if (directions[index] == P2AssociationDirection::AssociationEnergy) {
            active.push_back({epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}});
        }
        const epcsaft::HelmholtzModel production_model(
            associating_parameters, active
        );
        const auto production = production_model.evaluate(
            300.0, 1.0 / 5000.0, action_composition, epcsaft::DerivativeOrder::Hessian
        );
        const auto reference = properties_support::reference_tensor(
            associating_parameters, 300.0, 1.0 / 5000.0, action_composition,
            action_state.site_fractions, active
        );
        const std::size_t coordinate = directions[index] == P2AssociationDirection::Temperature
            ? 0 : directions[index] == P2AssociationDirection::Volume
                ? 1 : 2 + action_composition.size();
        const double production_error = std::abs(
            production.hessian[coordinate * production.gradient.size() + coordinate]
            - reference.hessian[coordinate * reference.gradient.size() + coordinate]
        ) / std::max(1.0, std::abs(reference.hessian[coordinate * reference.gradient.size() + coordinate]));
        const auto action = p2_association_action(
            300.0, 1.0 / 5000.0, action_composition,
            associating.parameters(), action_state.site_fractions,
            directions[index]
        );
        const bool action_ok = production_error <= 1e-10
            && action.rank == action_state.site_fractions.size()
            && std::isfinite(action.pivot_ratio)
            && action.pivot_ratio > epcsaft::policy::ASSOCIATION_JACOBIAN_CONDITION.value
            && std::isfinite(action.residual) && action.residual <= 1e-10;
        const auto inputs = recording_support::fields({
            model_inputs(associating.parameters(), 300.0, 1.0 / 5000.0, action_composition),
            recording_support::field("direction", direction_names[index]),
            recording_support::field("production_hessian_coordinate", coordinate),
            recording_support::field("production_reference_error", production_error),
            recording_support::field("step", direction_steps[index]),
            recording_support::field("rank", action.rank),
            recording_support::field("pivot_ratio", action.pivot_ratio),
        });
        recording_support::compare_custom(
            (std::string("association-action/") + direction_names[index]).c_str(),
            "production_action_error", "1", production_error, 0.0, inputs, "Relationship",
            "production-vs-raw-residual", "production Hessian matches independent raw-residual elimination",
            "absolute", action_ok, "p2-association-action", direction_names[index],
            "H_q,||UX_q+F_q||inf", "T,rho,V,theta", "H_production=H_raw;||UX_q+F_q||inf<=1e-10", "direction");
    }

    const auto dense_state = associating.state_at_pressure(
        353.15, 1e5, {0.25, 0.75}, epcsaft::DensityBranch::Liquid);
    const double dense_pressure_closure = std::abs(dense_state.pressure - 1e5) / 1e5;
    recording_support::compare(
        "root-residual/associating-dense", "relative_pressure_residual", "1",
        dense_pressure_closure, 0.0, recording_support::fields({
            model_inputs(associating.parameters(), 353.15, 1.0 / dense_state.molar_density, {0.25, 0.75}),
            recording_support::field("P_star", 1e5),
            recording_support::field("P_published", dense_state.pressure),
        }), "Relationship", "density-publication",
        "dense-liquid publication closes the imposed pressure", "absolute", 1e-10, 0.0,
        dense_pressure_closure <= 1e-10, "p2-root-residual", "associating-dense",
        "|P(rho)-P_star|/P_star", "T,P,x,liquid", "|P(rho)-P_star|/P_star<=1e-10", "branch");

    epcsaft::HelmholtzModel ethane(properties_support::ethane());
    constexpr double temperature = 200.0;
    constexpr double pressure = 217800.69494692178;
    const std::vector<double> composition = {1.0};
    const auto vapor_root = epcsaft::solve_density(
        ethane, temperature, pressure, composition, epcsaft::DensityBranch::Vapor, std::nullopt
    );
    const auto liquid_root = epcsaft::solve_density(
        ethane, temperature, pressure, composition, epcsaft::DensityBranch::Liquid, std::nullopt
    );
    epcsaft::DensityRoot sibling;
    const auto paired_vapor = epcsaft::solve_density(
        ethane, temperature, pressure, composition, epcsaft::DensityBranch::Vapor, std::nullopt, &sibling
    );
    const auto sibling_state = ethane.state(temperature, sibling.molar_density, composition);
    const double sibling_residual = std::abs(sibling_state.pressure - pressure)
        / (sibling.molar_density * epcsaft::GAS_CONSTANT_J_PER_MOL_K * temperature);
    const bool sibling_refined = close(paired_vapor.molar_density, vapor_root.molar_density, 1e-10)
        && close(sibling.molar_density, liquid_root.molar_density, 1e-10)
        && sibling.branch == epcsaft::DensityBranch::Liquid && sibling.stable_root_count == 2
        && paired_vapor.branch == epcsaft::DensityBranch::Vapor && paired_vapor.stable_root_count == 2
        && sibling_state.pressure_density_derivative > 0.0 && sibling_residual <= 1e-12;
    recording_support::compare("root-sibling/liquid", "molar_density", "mol/m3",
        sibling.molar_density, liquid_root.molar_density, recording_support::fields({
            model_inputs(ethane.parameters(), temperature, 1.0 / sibling.molar_density, composition),
            recording_support::field("P_star", pressure), recording_support::field("reduced_residual", sibling_residual)}),
        "Relationship", "pressure-closure", "an actual sibling is the refined opposite stable root",
        "reference_or_one", 0.0, 1e-10, sibling_refined, "p2-root-sibling", "liquid", "rho", "T,P,x,branch",
        "rho_sibling=rho_liquid", "case");
    const auto vapor_state = ethane.state_at_pressure(
        temperature, pressure, composition, epcsaft::DensityBranch::Vapor
    );
    const auto liquid_state = ethane.state_at_pressure(
        temperature, pressure, composition, epcsaft::DensityBranch::Liquid
    );
    const std::array<const epcsaft::DensityRoot*, 2> roots = {&vapor_root, &liquid_root};
    const std::array<const epcsaft::State*, 2> states = {&vapor_state, &liquid_state};
    const std::array<const char*, 2> names = {"vapor", "liquid"};
    const std::array<epcsaft::DensityBranch, 2> branches = {
        epcsaft::DensityBranch::Vapor, epcsaft::DensityBranch::Liquid
    };
    for (std::size_t index = 0; index < roots.size(); ++index) {
        const auto& root = *roots[index];
        const auto& state = *states[index];
        const auto root_state = ethane.state(temperature, root.molar_density, composition);
        const double root_pressure = root_state.pressure - pressure;
        const double pressure_error = std::abs(root.pressure_residual - root_pressure);
        const auto inputs = recording_support::fields({
            model_inputs(ethane.parameters(), temperature, 1.0 / root.molar_density, composition),
            recording_support::field("P_star", pressure),
            recording_support::field("requested_branch", names[index]),
            recording_support::field("stable_root_count", root.stable_root_count.value_or(0)),
            recording_support::field("density", root.molar_density),
            recording_support::field("root_pressure_residual_Pa", root.pressure_residual),
            recording_support::field("reduced_convergence_residual", root.convergence.residual),
            recording_support::field("pressure_density_derivative", root_state.pressure_density_derivative),
            recording_support::field("D", root_state.pressure_density_derivative
                / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * temperature)),
            recording_support::field("stop", epcsaft::core::stop_name(root.convergence.stop)),
            recording_support::field("criterion", root.convergence.criterion->name),
            recording_support::field("association_rank", std::size_t{0}),
            recording_support::field("association_pivot_ratio", 1.0),
        });
        recording_support::compare(
            (std::string("root-residual/") + names[index]).c_str(), "pressure_residual", "Pa",
            root.pressure_residual, root_pressure, inputs, "Identity", "state-replay",
            "DensityRoot pressure residual is P(root)-P", "absolute", 1e-7, 0.0,
            pressure_error <= 1e-7, "p2-root-residual", names[index], "P_{root}-P",
            "T,P,x,branch", "rP=P(root)-P", "branch");
        recording_support::compare(
            (std::string("root-convergence/") + names[index]).c_str(), "reduced_residual", "1",
            root.convergence.residual, 0.0, inputs, "Availability", "density-policy",
            "density convergence uses the reduced pressure residual", "absolute", 1e-12, 0.0,
            root.convergence.residual <= 1e-12, "p2-root-convergence", names[index],
            "r_reduced", "T,P,x,branch", "r_reduced<=1e-12", "branch");
        const double reuse_error = std::abs(state.molar_density - root.molar_density);
        const bool metadata = root.branch == branches[index]
            && root.stable_root_count == std::optional<std::size_t>{2}
            && state.density_branch == branches[index]
            && state.stable_root_count == std::optional<std::size_t>{2}
            && reuse_error <= 1e-12 * std::max(1.0, root.molar_density);
        const auto metadata_inputs = recording_support::fields({
            inputs,
            recording_support::field("state_branch", names[index]),
            recording_support::field("state_count", state.stable_root_count.value_or(0)),
            recording_support::field("state_density", state.molar_density),
            recording_support::field("density_reuse_error", reuse_error),
        });
        recording_support::compare(
            (std::string("root-metadata/") + names[index]).c_str(), "molar_density", "mol/m3",
            state.molar_density, root.molar_density, metadata_inputs, "Identity",
            "state-propagation", "selected root and State retain density, branch and count",
            "reference_or_one", 0.0, 1e-12, metadata, "p2-root-metadata", names[index],
            "rho", "T,P,x,branch", "rho_state=rho_root;branch=count", "branch");
    }

    bool ambiguous = false;
    try {
        (void)epcsaft::solve_density(
            ethane, temperature, pressure, composition, epcsaft::DensityBranch::Unspecified,
            std::nullopt
        );
    } catch (const std::runtime_error& error) {
        ambiguous = std::string(error.what()).find("ambiguous") != std::string::npos;
    }
    recording_support::compare_custom(
        "density-selection/ambiguous", "selected", "bool", ambiguous ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("T", temperature),
            recording_support::field("P_star", pressure),
            recording_support::field("requested_branch", "unspecified"),
            recording_support::field("stable_root_count", 2U),
            recording_support::field("outcome", "ambiguous"),
        }), "Availability", "bounded-scan", "unspecified multi-root request fails closed",
        "absolute", ambiguous, "p2-density-selection", "ambiguous", "selected",
        "T,P,x,branch,count", "no selected density", "case");

    constexpr double coalescing_temperature = 320.0;
    const double coalescing_upper = ethane.molar_density_upper_bound(
        coalescing_temperature, composition
    ) * epcsaft::policy::DENSITY_UPPER_FRACTION.value;
    double coalescing_pressure = 1e6;
    double coalescing_density = 0.0;
    for (std::size_t iteration = 0; iteration < 40; ++iteration) {
        const double ideal_density = coalescing_pressure
            / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * coalescing_temperature);
        const double lower_density = ideal_density
            * epcsaft::policy::DENSITY_LOWER_FRACTION.value;
        coalescing_density = std::exp(
            std::log(lower_density) + 0.5
            * (std::log(coalescing_upper) - std::log(lower_density))
        );
        coalescing_pressure = ethane.state(
            coalescing_temperature, coalescing_density, composition
        ).pressure;
    }
    const auto coalescing_root = epcsaft::solve_density(
        ethane, coalescing_temperature, coalescing_pressure, composition,
        epcsaft::DensityBranch::Vapor,
        std::nullopt
    );
    const double coalescing_alignment = std::abs(
        coalescing_root.molar_density - coalescing_density
    ) / std::max(1.0, coalescing_density);
    const bool coalesced = close(
        coalescing_root.molar_density, coalescing_density, 1e-12
    ) && coalescing_root.stable_root_count == std::optional<std::size_t>{1};
    recording_support::compare_custom(
        "root-coalescing/grid", "grid_alignment", "1", coalescing_alignment, 0.0,
        recording_support::fields({
            recording_support::field("T", coalescing_temperature),
            recording_support::field("P_star", coalescing_pressure),
            recording_support::field("grid_index", std::size_t{20}),
            recording_support::field("grid_density", coalescing_density),
            recording_support::field("root_density", coalescing_root.molar_density),
            recording_support::field("stable_root_count",
                coalescing_root.stable_root_count.value_or(0)),
        }), "Relationship", "bounded-scan", "grid-point root is coalesced once",
        "absolute", coalesced, "p2-root-coalescing", "grid", "rho",
        "T,P,x,grid_index", "rho_root=rho_grid", "case");

    const auto supercritical_vapor = ethane.state_at_pressure(
        320.0, 1e6, composition, epcsaft::DensityBranch::Vapor
    );
    const auto supercritical_liquid = ethane.state_at_pressure(
        320.0, 1e6, composition, epcsaft::DensityBranch::Liquid
    );
    const auto supercritical_root = epcsaft::solve_density(
        ethane, 320.0, 1e6, composition, epcsaft::DensityBranch::Vapor, std::nullopt
    );
    const auto supercritical_inputs = recording_support::fields({
        model_inputs(ethane.parameters(), 320.0, 1.0 / supercritical_vapor.molar_density, composition),
        recording_support::field("P", 1e6),
        recording_support::field("vapor_density", supercritical_vapor.molar_density),
        recording_support::field("liquid_density", supercritical_liquid.molar_density),
        recording_support::field("root_density", supercritical_root.molar_density),
        recording_support::field("root_branch", "unspecified"),
        recording_support::field("root_count", supercritical_root.stable_root_count.value_or(0)),
        recording_support::field("state_branch", "unspecified"),
        recording_support::field("state_count", supercritical_vapor.stable_root_count.value_or(0)),
    });
    recording_support::compare(
        "single-metadata/root", "molar_density", "mol/m3", supercritical_root.molar_density,
        supercritical_vapor.molar_density, supercritical_inputs, "Relationship", "single-root",
        "one stable root is shared by DensityRoot and both State requests", "reference_or_one", 0.0, 1e-12,
        supercritical_root.branch == epcsaft::DensityBranch::Unspecified
            && supercritical_root.stable_root_count == std::optional<std::size_t>{1}
            && supercritical_vapor.density_branch == epcsaft::DensityBranch::Unspecified
            && supercritical_vapor.stable_root_count == std::optional<std::size_t>{1}
            && supercritical_liquid.density_branch == epcsaft::DensityBranch::Unspecified
            && supercritical_liquid.stable_root_count == std::optional<std::size_t>{1}
            && supercritical_vapor.requested_density_branch == epcsaft::DensityBranch::Vapor
            && supercritical_liquid.requested_density_branch == epcsaft::DensityBranch::Liquid
            && supercritical_vapor.returned_density_branch() == epcsaft::DensityBranch::Unspecified
            && supercritical_liquid.returned_density_branch() == epcsaft::DensityBranch::Unspecified
            && !supercritical_vapor.density_branch_mismatch()
            && !supercritical_liquid.density_branch_mismatch()
            && close(supercritical_root.molar_density, supercritical_vapor.molar_density, 1e-12),
        "p2-single-metadata", "root", "rho", "T,P,x", "rho_root=rho_state;branch=count", "state");

    const auto cold_vapor = ethane.state_at_pressure(
        280.0, 1e6, composition, epcsaft::DensityBranch::Vapor
    );
    const auto cold_liquid = ethane.state_at_pressure(
        280.0, 1e6, composition, epcsaft::DensityBranch::Liquid
    );
    const auto anchored_vapor = ethane.state_at_pressure(
        280.0, 1e6, composition, epcsaft::DensityBranch::Vapor,
        cold_vapor.molar_density
    );
    const auto anchored_liquid = ethane.state_at_pressure(
        280.0, 1e6, composition, epcsaft::DensityBranch::Liquid,
        cold_liquid.molar_density
    );
    const auto anchor_inputs = recording_support::fields({
        recording_support::field("T", 280.0), recording_support::field("P", 1e6),
        recording_support::field("requested_branch", "vapor/liquid"),
        recording_support::field("returned_branch", "vapor/liquid"),
        recording_support::field("stable_root_count", "absent"),
        recording_support::field("envelope_log_half_width",
            epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value),
        recording_support::field("scan_intervals",
            static_cast<std::size_t>(epcsaft::policy::DENSITY_ANCHOR_SCAN_INTERVALS.value)),
        recording_support::field("cold_vapor_density", cold_vapor.molar_density),
        recording_support::field("cold_liquid_density", cold_liquid.molar_density),
        recording_support::field("vapor_density", anchored_vapor.molar_density),
        recording_support::field("liquid_density", anchored_liquid.molar_density),
        recording_support::field("vapor_reduced_residual",
            anchored_vapor.density_convergence->residual),
        recording_support::field("liquid_reduced_residual",
            anchored_liquid.density_convergence->residual),
        recording_support::field("association_rank", std::size_t{0}),
        recording_support::field("association_pivot_ratio", 1.0),
    });
    recording_support::compare(
        "anchor/vapor", "molar_density", "mol/m3", anchored_vapor.molar_density,
        cold_vapor.molar_density, anchor_inputs, "Relationship", "local-continuation",
        "vapor anchor replays the local vapor root", "reference_or_one", 0.0, 1e-12,
        close(anchored_vapor.molar_density, cold_vapor.molar_density, 1e-12),
        "p2-anchor", "vapor", "rho_v", "T,P,x",
        "rho_anchor_v=rho_v", "case");
    recording_support::compare(
        "anchor/liquid", "molar_density", "mol/m3", anchored_liquid.molar_density,
        cold_liquid.molar_density, anchor_inputs, "Relationship", "local-continuation",
        "liquid anchor replays the local liquid root", "reference_or_one", 0.0, 1e-12,
        close(anchored_liquid.molar_density, cold_liquid.molar_density, 1e-12),
        "p2-anchor", "liquid", "rho_l", "T,P,x",
        "rho_anchor_l=rho_l", "case");

    // The certificate covers only the offered local window, including its final bracket.
    for (const double offset : {-1e-11, 1e-13}) {
        const double offered=cold_vapor.molar_density*std::exp(-epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value+offset);
        const double edge=offered*std::exp(epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value);
        bool returned=false, confined=false;
        double density=0.0, bracket_lower=0.0, bracket_upper=0.0, residual=0.0;
        epcsaft::core::Capture capture(true);
        try {
            const auto root=epcsaft::solve_density(ethane,280.0,1e6,composition,epcsaft::DensityBranch::Vapor,offered);
            returned=true; density=root.molar_density; residual=root.convergence.residual;
            confined=density<=edge&&close(density,cold_vapor.molar_density,1e-10);
            for(const auto& event:capture.trace().events)
                if(std::string(event.name)=="density_root_bracket") {
                    bracket_lower=event.values[0]; bracket_upper=event.values[1];
                    confined=confined&&event.values[1]<=edge;
                }
        } catch(const std::runtime_error&) {}
        const bool inside=offset>0.0, ok=inside ? returned&&confined : !returned;
        recording_support::compare_custom((std::string("anchor/window-edge/")+(offset>0.0?"inside":"outside")).c_str(),
            "admitted","bool",returned?1.0:0.0,inside?1.0:0.0,recording_support::fields({anchor_inputs,
                recording_support::field("offered_density",offered),recording_support::field("upper_edge",edge),
                recording_support::field("offset",offset),recording_support::field("returned_density",density),
                recording_support::field("bracket_lower",bracket_lower),recording_support::field("bracket_upper",bracket_upper),
                recording_support::field("residual",residual)}),
            "Availability","local-window","the selected crossing and its refine bracket stay inside the local window",
            "absolute",ok,"p2-anchor-window",inside?"inside":"outside","admitted","T,P,x,anchor","root in window","case");
    }

    // An above-critical isotherm cannot use positive slope tangents to certify no hidden loop.
    for(const double target:{1e6*(1-1e-6),1e6,1e6*(1+1e-6)}) {
        double lo=1.0,hi=2000.0;
        CHECK(ethane.state(320.0,lo,composition).pressure<target);
        CHECK(ethane.state(320.0,hi,composition).pressure>target);
        for(std::size_t i=0;i<100&&std::nextafter(lo,hi)<hi;++i) {
            const double mid=0.5*(lo+hi);
            (ethane.state(320.0,mid,composition).pressure<target?lo:hi)=mid;
        }
        const double reference=0.5*(lo+hi);
        for(const auto branch:{epcsaft::DensityBranch::Vapor,epcsaft::DensityBranch::Liquid}) {
            epcsaft::core::Capture capture(true);
            const auto root=epcsaft::solve_density(ethane,320.0,target,composition,branch,std::nullopt);
            bool fallback=false;
            for(const auto& event:capture.trace().events)
                fallback=fallback||(std::string(event.name)=="structured_broken"&&event.reason=="loop absence unresolved");
            const bool ok=fallback&&root.branch==epcsaft::DensityBranch::Unspecified&&root.stable_root_count==1
                &&root.pressure_density_derivative>0.0&&close(root.molar_density,reference,1e-10);
            recording_support::compare((std::string("sole-root/fallback/")+std::to_string(target)+
                (branch==epcsaft::DensityBranch::Vapor?"/vapor":"/liquid")).c_str(),"molar_density","mol/m3",root.molar_density,
                reference,recording_support::fields({recording_support::field("T",320.0),recording_support::field("P_star",target),
                    recording_support::field("fallback",fallback),recording_support::field("stable_root_count",root.stable_root_count.value_or(0)),
                    recording_support::field("slope",root.pressure_density_derivative)}),"Relationship","fixed-density-bisection",
                "unresolved loop absence scans and preserves the sole unlabelled root", "reference_or_one",0.0,1e-10,ok,
                "p2-sole-root-fallback",branch==epcsaft::DensityBranch::Vapor?"vapor":"liquid","rho","T,P,x,branch","rho=rho_ref","case");
        }
    }

    const auto association_anchor = associating.state_at_pressure(
        376.88, 101330.0, {0.1385, 0.8615}, epcsaft::DensityBranch::Liquid
    );
    bool branch_unreachable = false;
    bool typed_failure = false;
    double failure_envelope = 0.0;
    std::size_t failure_intervals = 0;
    epcsaft::core::Stop failure_stop = epcsaft::core::Stop::stalled;
    double failure_residual = 0.0;
    try {
        (void)associating.state_at_pressure(
            376.88, 101330.0, {0.99, 0.01}, epcsaft::DensityBranch::Liquid,
            association_anchor.molar_density
        );
    } catch (const epcsaft::core::ConvergenceFailure& error) {
        typed_failure = true;
        branch_unreachable = error.record.stop == epcsaft::core::Stop::left_domain;
        failure_envelope = epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value;
        failure_intervals = static_cast<std::size_t>(
            epcsaft::policy::DENSITY_ANCHOR_SCAN_INTERVALS.value);
        failure_stop = error.record.stop;
        failure_residual = error.record.residual;
    } catch (const std::runtime_error&) {
        branch_unreachable = false;
    }
    const auto unreachable_inputs = recording_support::fields({
        recording_support::field("T", 376.88), recording_support::field("P", 101330.0),
        recording_support::field("anchor_x", std::vector<double>{0.1385, 0.8615}),
        recording_support::field("target_x", std::vector<double>{0.99, 0.01}),
        recording_support::field("anchor_density", association_anchor.molar_density),
        recording_support::field("typed_failure", typed_failure ? 1U : 0U),
        recording_support::field("stop", epcsaft::core::stop_name(failure_stop)),
        recording_support::field("residual", failure_residual),
        recording_support::field("envelope_log_half_width", failure_envelope),
        recording_support::field("scan_intervals", failure_intervals),
    });
    recording_support::compare_custom(
        "anchor/unreachable", "branch_reachable", "bool", branch_unreachable ? 1.0 : 0.0, 1.0,
        unreachable_inputs, "Availability", "local-continuation", "unreachable anchored branch fails closed",
        "absolute", branch_unreachable && typed_failure
            && failure_envelope == epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value
            && failure_intervals == static_cast<std::size_t>(
                epcsaft::policy::DENSITY_ANCHOR_SCAN_INTERVALS.value),
        "p2-anchor", "unreachable", "reachable",
        "T,P,x_anchor,x_target", "no cold fallback", "case");
    CHECK(branch_unreachable && typed_failure);

    bool invalid_anchor = false;
    try {
        (void)epcsaft::solve_density(
            ethane, temperature, pressure, composition, epcsaft::DensityBranch::Vapor, -1.0
        );
    } catch (const std::invalid_argument&) {
        invalid_anchor = true;
    }
    recording_support::compare_custom(
        "anchor-invalid/negative", "accepted", "bool", invalid_anchor ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("T", temperature),
            recording_support::field("P_star", pressure),
            recording_support::field("anchor_density", -1.0),
            recording_support::field("outcome", "invalid-anchor"),
        }), "Availability", "trust-boundary", "nonpositive anchor is rejected",
        "absolute", invalid_anchor, "p2-anchor-invalid", "negative", "accepted",
        "T,P,x,anchor", "anchor rejected", "case");

    struct P2ChainCase final {
        const epcsaft::HelmholtzModel* model;
        double temperature;
        double pressure;
        std::vector<double> composition;
        epcsaft::DensityBranch branch;
        const char* branch_name;
        const char* point;
        const char* quantity;
        const char* unit;
        bool enthalpy;
    };
    const std::array<double, 3> delta_temperatures = {1e-2, 3e-3, 1e-3};
    const std::array<P2ChainCase, 3> chain_cases = {{
        {&ethane, 320.0, 1e6, composition, epcsaft::DensityBranch::Vapor,
            "vapor", "one-root-density", "density_temperature_derivative", "mol/m3/K", false},
        {&ethane, 200.0, pressure, composition, epcsaft::DensityBranch::Liquid,
            "liquid", "multi-root-enthalpy", "residual_enthalpy_temperature_derivative", "J/mol/K", true},
        {&associating, 313.15, 1e5, {0.1385, 0.8615}, epcsaft::DensityBranch::Liquid,
            "liquid", "associating-enthalpy", "residual_enthalpy_temperature_derivative", "J/mol/K", true},
    }};
    const auto record_chain = [&](const char* point, const char* quantity,
                                  const char* unit, double temperature_value,
                                  double pressure_value,
                                  const char* requested_branch,
                                  double delta_temperature,
                                  const P2PressureChain& chain) {
        const bool ok = close(chain.measured, chain.predicted, 5e-6)
            && chain.pressure_density_derivative > 0.0
            && std::isfinite(chain.reduced_residual);
        recording_support::compare(
            (std::string("fixed-pressure/") + point).c_str(), quantity, unit,
            chain.measured, chain.predicted,
            recording_support::fields({
                recording_support::field("T", temperature_value),
                recording_support::field("P_star", pressure_value),
                recording_support::field("requested_branch", requested_branch),
                recording_support::field("returned_branch",
                    chain.branch == epcsaft::DensityBranch::Unspecified
                        ? "unspecified" : requested_branch),
                recording_support::field("stable_root_count", chain.stable_root_count),
                recording_support::field("P_rho", chain.pressure_density_derivative),
                recording_support::field("D", chain.pressure_density_derivative
                    / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * temperature_value)),
                recording_support::field("reduced_residual", chain.reduced_residual),
                recording_support::field("association_rank", chain.association_rank),
                recording_support::field("association_pivot_ratio",
                    chain.association_pivot_ratio),
                recording_support::field("delta_T", delta_temperature),
            }), "Relationship", "implicit-chain",
            "fixed-pressure property chain includes density response",
            "reference_or_one", 0.0, 5e-6, ok, "p2-fixed-pressure", point,
            quantity, "T,P,x,branch,delta_T", "Q_P=Q_T+Q_rho rho_P", "case");
    };
    for (const auto& item : chain_cases) {
        for (const double delta_temperature : delta_temperatures) {
            const auto chain = p2_fixed_pressure_chain(
                *item.model, item.temperature, item.pressure, item.composition,
                item.branch, delta_temperature, item.enthalpy
            );
            const std::string point = std::string(item.point) + "/"
                + std::to_string(delta_temperature);
            record_chain(point.c_str(), item.quantity, item.unit, item.temperature,
                item.pressure, item.branch_name, delta_temperature, chain);
        }
    }

    // #89: a van der Waals loop narrower than one 40-cell scan interval. The reference is
    // independent of solve_density: sign changes of P(T,rho)-P* on a 4000-point log grid at
    // fixed density, bisected to machine precision, classified by the analytic slope.
    // The 300 K state sits within 0.05 % of the liquid spinodal pressure: the liquid root lies
    // 1.1 % above the unstable one and the residual dip between them is far narrower than a cell.
    // The eighth state is 1e-4 K below the critical point: the negative-slope dip is about 0.002 wide
    // in log density, narrower than any fixed search width, and the pressure extrema differ by 0.021 Pa.
    // The ninth lies 7.4 Pa above it, near the vapor-side spinodal pressure: the vapor and unstable
    // roots are 3 mol/m3 apart (hence the denser reference grid) and the vapor root's residual
    // tolerance band is 0.05 mol/m3 wide, where a shared cell endpoint inside the band was once
    // reported as a third stable root. Near-tangent roots are compared at the solver's first-order
    // accuracy tol |d rho/d r| when that exceeds 1e-7 relative, and that accuracy must stay far
    // below the branch separation.
    // The final three inputs put the liquid root exactly on the former scan grid and one
    // pressure ULP on either side. A grid coincidence must not erase or duplicate a branch.
    struct LoopCase final {
        double temperature, pressure; std::size_t reference_points; const char* suffix = "";
        double reference_lower = 1000.0, reference_upper = 0.0;
    };
    const std::array<LoopCase, 13> loop_cases = {{
        {300.0, 3916073.568, 4000}, {305.0, 4758116.154, 4000}, {306.0, 4862786.389, 4000},
        {307.0, 4965901.440, 4000}, {308.0, 5067464.963, 4000}, {308.8, 5147500.0, 4000},
        {308.9, 5157549.110, 4000}, {308.9554032139963, 5163053.740696311, 4000},
        {308.95547769995335, 5163061.1390785752, 16000},
        {308.0, 5057915.461281232, 4000, "/grid"},
        {308.0, std::nextafter(5057915.461281232, 0.0), 4000, "/grid-below"},
        {308.0, std::nextafter(5057915.461281232, 1e7), 4000, "/grid-above"},
        // Residual-only stopping overshot its local sensitivity radius on the liquid branch.
        {308.9554732140194, 5163060.690898801, 1000, "/radius", 6380.0, 6400.0}}};
    for (const auto& item : loop_cases) {
        const double reference_upper = item.reference_upper > 0.0 ? item.reference_upper
            : ethane.molar_density_upper_bound(item.temperature, composition)
                * epcsaft::policy::DENSITY_UPPER_FRACTION.value;
        const auto residual_at = [&](double density) {
            return ethane.state(item.temperature, density, composition).pressure - item.pressure;
        };
        std::vector<std::pair<double, double>> references;
        double previous_density = item.reference_lower, previous_residual = residual_at(previous_density);
        for (std::size_t index = 1; index <= item.reference_points; ++index) {
            const double density = item.reference_lower * std::pow(reference_upper / item.reference_lower,
                                                     static_cast<double>(index) / static_cast<double>(item.reference_points));
            const double residual = residual_at(density);
            if (previous_residual * residual < 0.0) {
                double lower = previous_density, upper = density, lower_residual = previous_residual;
                for (std::size_t bisection = 0; bisection < 80; ++bisection) {
                    const double middle = 0.5 * (lower + upper), value = residual_at(middle);
                    if (lower_residual * value <= 0.0) upper = middle; else { lower = middle; lower_residual = value; }
                }
                const double root = 0.5 * (lower + upper);
                references.emplace_back(root, ethane.state(item.temperature, root, composition).pressure_density_derivative);
            }
            previous_density = density;
            previous_residual = residual;
        }
        CHECK(references.size() == 3 && references[0].second > 0.0 && references[1].second < 0.0
              && references[2].second > 0.0);
        const std::string label = std::to_string(item.temperature) + item.suffix;
        const auto loop_inputs = recording_support::fields({
            recording_support::field("T", item.temperature), recording_support::field("P_star", item.pressure),
            recording_support::field("reference_vapor", references[0].first),
            recording_support::field("reference_unstable", references[1].first),
            recording_support::field("reference_liquid", references[2].first),
            recording_support::field("loop_log_width", std::log(references[2].first / references[0].first)),
        });
        const std::array<std::pair<epcsaft::DensityBranch, const char*>, 2> branches = {{
            {epcsaft::DensityBranch::Vapor, "loop-vapor/"}, {epcsaft::DensityBranch::Liquid, "loop-liquid/"}}};
        for (std::size_t side = 0; side < 2; ++side) {
            const auto root = epcsaft::solve_density(
                ethane, item.temperature, item.pressure, composition, branches[side].first, std::nullopt);
            const double reference = references[side == 0 ? 0 : 2].first;
            const double accuracy = std::max(1e-7 * reference, epcsaft::policy::DENSITY_ROOT_RESIDUAL.value
                * root.convergence.amplification.value_or(0.0));
            CHECK(accuracy < 1e-2 * (references[2].first - references[0].first));
            const bool ok = std::abs(root.molar_density - reference) <= accuracy && root.branch == branches[side].first
                && root.stable_root_count == std::optional<std::size_t>{2} && root.pressure_density_derivative > 0.0
                && root.convergence.stop == epcsaft::core::Stop::converged;
            recording_support::compare((branches[side].second + label).c_str(), "molar_density", "mol/m3",
                root.molar_density, reference, recording_support::fields({loop_inputs,
                    recording_support::field("stable_root_count", root.stable_root_count.value_or(0)),
                    recording_support::field("pressure_density_derivative", root.pressure_density_derivative),
                    recording_support::field("first_order_accuracy", epcsaft::policy::DENSITY_ROOT_RESIDUAL.value
                        * root.convergence.amplification.value_or(0.0))}),
                "Relationship", "independent-bisection", "declared branch selects the extreme stable root",
                "absolute", accuracy, 0.0, ok, side == 0 ? "p2-loop-vapor" : "p2-loop-liquid",
                label.c_str(), "rho", "T,P,x,branch", "rho_branch=rho_ref", "case");
        }
        bool loop_ambiguous = false;
        try {
            (void)epcsaft::solve_density(ethane, item.temperature, item.pressure, composition,
                                         epcsaft::DensityBranch::Unspecified, std::nullopt);
        } catch (const std::runtime_error& error) {
            loop_ambiguous = std::string(error.what()).find("ambiguous") != std::string::npos;
        }
        recording_support::compare_custom(("loop-ambiguous/" + label).c_str(), "selected", "bool",
            loop_ambiguous ? 1.0 : 0.0, 1.0, loop_inputs, "Availability", "bounded-scan",
            "unspecified request on two stable roots fails closed", "absolute", loop_ambiguous,
            "p2-loop-ambiguous", label.c_str(), "selected", "T,P,x", "no selected density", "case");
        if (item.temperature != 308.8) continue;
        for (std::size_t side = 0; side < 2; ++side) {
            const auto batch = ethane.pressure_actions(item.temperature, item.pressure, composition,
                branches[side].first, {}, {epcsaft::PressureActionDirection{0.0, 1.0, {0.0}, {}}});
            const double reference = 1.0 / references[side == 0 ? 0 : 2].second;
            const bool ok = batch.density_status == epcsaft::PropertyActionStatus::Available
                && batch.density_action.size() == 1 && close(batch.density_action[0], reference, 1e-6);
            recording_support::compare((std::string("loop-action/") + (side == 0 ? "vapor" : "liquid")).c_str(),
                "density_action", "mol/(m3 Pa)", batch.density_action.empty() ? 0.0 : batch.density_action[0],
                reference, loop_inputs, "Identity", "reciprocal-slope",
                "d rho/dP on the selected branch is 1/P_rho at the independent root",
                "reference_or_one", 0.0, 1e-6, ok, "p2-loop-action", side == 0 ? "vapor" : "liquid",
                "d rho/dP", "T,P,x,branch", "rho_P=1/P_rho", "branch");
        }
    }
    // Independently locate the vapor-side fold at 308 K from its analytic fixed-density slope.
    // A pressure indistinguishable from this maximum cannot certify two nonsingular branches.
    double fold_lower = 5000.0, fold_upper = 6400.0;
    CHECK(ethane.state(308.0, fold_lower, composition).pressure_density_derivative > 0.0);
    CHECK(ethane.state(308.0, fold_upper, composition).pressure_density_derivative < 0.0);
    for (std::size_t i = 0; i < 80; ++i) {
        const double middle = 0.5 * (fold_lower + fold_upper);
        if (ethane.state(308.0, middle, composition).pressure_density_derivative > 0.0)
            fold_lower = middle;
        else fold_upper = middle;
    }
    const auto fold = ethane.state(308.0, 0.5 * (fold_lower + fold_upper), composition);
    CHECK(std::abs(fold.pressure_density_derivative) < 1e-8);
    bool unresolved = false;
    try {
        (void)epcsaft::solve_density(ethane, 308.0, fold.pressure, composition,
                                    epcsaft::DensityBranch::Vapor, std::nullopt);
    } catch (const std::runtime_error& error) {
        unresolved = std::string(error.what()).find("unresolved at spinodal") != std::string::npos;
    }
    recording_support::compare_custom("fold-resolution/vapor", "unresolved", "bool",
        unresolved ? 1.0 : 0.0, 1.0, recording_support::fields({
            recording_support::field("T", 308.0), recording_support::field("P_star", fold.pressure),
            recording_support::field("fold_density", fold.molar_density),
            recording_support::field("fold_bracket_lower", fold_lower),
            recording_support::field("fold_bracket_upper", fold_upper),
            recording_support::field("pressure_density_derivative", fold.pressure_density_derivative)}),
        "Availability", "independent-slope-bisection", "a mechanical fold is not a resolved stable branch",
        "absolute", unresolved, "p2-fold-resolution", "vapor", "unresolved", "T,P,x", "P_rho=0", "case");

    // #391: Gross & Sadowski (2001) propylene at its coexistence pressure. The last cold-scan cell,
    // 9394-27094 mol/m3, has negative slopes at both ends (spinodals 14133 and 24347 mol/m3) but its
    // pressure rises through the target; the liquid root inside it was once discarded. References
    // are independent bisections of P(rho)-P* inside each positive-slope interval.
    auto propylene_parameters = properties_support::parameters(1);
    propylene_parameters.segment_count[0] = 1.9597;
    propylene_parameters.sigma_constant[0] = 3.5356;
    propylene_parameters.dispersion_energy_over_k[0] = 207.19;
    propylene_parameters.molar_mass_kg_per_mol[0] = 0.042081;
    epcsaft::HelmholtzModel propylene(propylene_parameters);
    constexpr double cold_temperature = 87.48019017432647, cold_pressure = 7.837213701088098e-4;
    const auto bisect = [&](double lower, double upper) {
        for (std::size_t i = 0; i < 200 && std::nextafter(lower, upper) < upper; ++i) {
            const double middle = 0.5 * (lower + upper);
            (propylene.state(cold_temperature, middle, composition).pressure < cold_pressure ? lower : upper) = middle;
        }
        return lower;
    };
    const std::array<double, 2> cold_references = {bisect(1e-7, 1e-5), bisect(15000.0, 24000.0)};
    const auto cold_inputs = recording_support::fields({recording_support::field("T", cold_temperature),
        recording_support::field("P_star", cold_pressure), recording_support::field("reference_vapor", cold_references[0]),
        recording_support::field("reference_liquid", cold_references[1])});
    for (std::size_t side = 0; side < 2; ++side) {
        const double reference = cold_references[side];
        bool ok = propylene.state(cold_temperature, reference, composition).pressure_density_derivative > 0.0;
        double density = 0.0;
        for (const auto anchor : {std::optional<double>{}, std::optional<double>{reference}}) {
            const auto root = epcsaft::solve_density(propylene, cold_temperature, cold_pressure, composition,
                                                     side ? epcsaft::DensityBranch::Liquid : epcsaft::DensityBranch::Vapor, anchor);
            ok = ok && std::abs(root.molar_density - reference) <= std::max(1e-7 * reference,
                    epcsaft::policy::DENSITY_ROOT_RESIDUAL.value * root.convergence.amplification.value_or(0.0))
                && root.branch == (side ? epcsaft::DensityBranch::Liquid : epcsaft::DensityBranch::Vapor)
                && root.pressure_density_derivative > 0.0
                && (anchor || root.stable_root_count == std::optional<std::size_t>{2});
            if (!anchor) density = root.molar_density;
        }
        recording_support::compare((std::string("hidden-crossing/") + (side ? "liquid" : "vapor")).c_str(),
            "molar_density", "mol/m3", density, reference, cold_inputs, "Relationship", "independent-bisection",
            "cold and anchored branch requests return the extreme stable root", "absolute", 1e-7 * reference, 0.0, ok,
            "p2-hidden-crossing", side ? "liquid" : "vapor", "rho", "T,P,x,branch", "rho_branch=rho_ref", "case");
    }
    bool cold_ambiguous = false;
    try {
        (void)epcsaft::solve_density(propylene, cold_temperature, cold_pressure, composition,
                                     epcsaft::DensityBranch::Unspecified, std::nullopt);
    } catch (const std::runtime_error& error) {
        cold_ambiguous = std::string(error.what()).find("ambiguous") != std::string::npos;
    }
    recording_support::compare_custom("hidden-crossing/ambiguous", "selected", "bool", cold_ambiguous ? 1.0 : 0.0,
        1.0, cold_inputs, "Availability", "bounded-scan", "unspecified request on two stable roots fails closed",
        "absolute", cold_ambiguous, "p2-hidden-crossing", "ambiguous", "selected", "T,P,x", "no selected density", "case");
}

}  // namespace properties_checks
