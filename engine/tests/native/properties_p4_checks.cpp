// P4, electrolyte coordinates and Issue10: property actions on parameter and state charts.
#include "properties_checks_shared.hpp"

namespace properties_checks {

double p4_reference(
    const epcsaft::Parameters& parameters,
    const std::vector<epcsaft::ActiveParameter>& active,
    const std::optional<epcsaft::ThermochemistryRecord>& thermochemistry,
    double temperature, double molar_density, const std::vector<double>& composition,
    const std::vector<epcsaft::DensityActionDirection>& directions,
    epcsaft::PropertyObservable observable, double step,
    std::vector<double>* admission = nullptr
) {
    epcsaft::core::Capture capture;
    const epcsaft::HelmholtzModel central_model(parameters, {}, thermochemistry);
    (void)central_model.state(temperature, molar_density, composition);
    const auto central_signature = p4_take_signature(capture, thermochemistry, temperature);
    const std::size_t combinations = std::size_t{1} << directions.size();
    double result = 0.0;
    bool support_unchanged = true, predicate_unchanged = true;
    for (std::size_t mask = 0; mask < combinations; ++mask) {
        double signed_temperature = temperature;
        double signed_density = molar_density;
        std::vector<double> signed_composition = composition;
        std::vector<double> parameter_shift(active.size(), 0.0);
        int sign_product = 1;
        for (std::size_t index = 0; index < directions.size(); ++index) {
            const double sign = (mask & (std::size_t{1} << index)) != 0 ? 1.0 : -1.0;
            sign_product *= sign > 0.0 ? 1 : -1;
            signed_temperature += step * sign * directions[index].temperature;
            signed_density += step * sign * directions[index].molar_density;
            for (std::size_t component = 0; component < composition.size(); ++component)
                signed_composition[component] += step * sign * directions[index].composition[component];
            for (const auto& parameter : directions[index].parameter) {
                for (std::size_t coordinate = 0; coordinate < active.size(); ++coordinate) {
                    if (active[coordinate].family == parameter.family
                        && active[coordinate].indices == parameter.indices) {
                        parameter_shift[coordinate] += step * sign * parameter.value;
                    }
                }
            }
        }
        const auto shifted = p4_shifted_parameters(parameters, active, parameter_shift);
        const epcsaft::HelmholtzModel model(shifted, {}, thermochemistry);
        const auto shifted_state = model.state(
            signed_temperature, signed_density, signed_composition);
        const bool predicate_match =
            p4_take_signature(capture, thermochemistry, signed_temperature)
                == central_signature;
        predicate_unchanged = predicate_unchanged && predicate_match;
        CHECK(predicate_match);
        for (std::size_t i = 0; i < composition.size(); ++i) {
            const bool support_match =
                (composition[i] > 0.0) == (signed_composition[i] > 0.0);
            support_unchanged = support_unchanged && support_match;
            CHECK(support_match);
        }
        result += sign_product * p4_state_value(shifted_state, observable);
    }
    if (admission != nullptr) admission->insert(admission->end(), {
        step, static_cast<double>(support_unchanged),
        static_cast<double>(predicate_unchanged)});
    return result / std::pow(2.0 * step, static_cast<int>(directions.size()));
}

double p4_pressure_reference(
    const epcsaft::Parameters& parameters,
    const std::vector<epcsaft::ActiveParameter>& active,
    const std::optional<epcsaft::ThermochemistryRecord>& thermochemistry,
    double temperature, double pressure, const std::vector<double>& composition,
    const std::vector<epcsaft::PressureActionDirection>& directions,
    epcsaft::DensityBranch branch, epcsaft::PropertyObservable observable, double step,
    std::vector<double>* admission = nullptr
) {
    epcsaft::core::Capture capture;
    const epcsaft::HelmholtzModel central_model(parameters, {}, thermochemistry);
    const auto central = central_model.state_at_pressure(temperature, pressure, composition, branch);
    const auto central_signature = p4_take_signature(capture, thermochemistry, temperature);
    double minimum_slope = std::numeric_limits<double>::infinity();
    double maximum_closure = 0.0, maximum_site_residual = 0.0, minimum_pivot_ratio = 1.0;
    double result = 0.0;
    const std::size_t combinations = std::size_t{1} << directions.size();
    for (std::size_t mask = 0; mask < combinations; ++mask) {
        double signed_temperature = temperature;
        double signed_pressure = pressure;
        std::vector<double> signed_composition = composition;
        std::vector<double> parameter_shift(active.size(), 0.0);
        const int sign_product = [&] {
            int product = 1;
            for (std::size_t index = 0; index < directions.size(); ++index) {
                const double sign = (mask & (std::size_t{1} << index)) != 0 ? 1.0 : -1.0;
                product *= sign > 0.0 ? 1 : -1;
                signed_temperature += step * sign * directions[index].temperature;
                signed_pressure += step * sign * directions[index].pressure;
                for (std::size_t component = 0; component < composition.size(); ++component)
                    signed_composition[component] += step * sign * directions[index].composition[component];
                for (const auto& parameter : directions[index].parameter) {
                    for (std::size_t coordinate = 0; coordinate < active.size(); ++coordinate)
                        if (active[coordinate].family == parameter.family
                            && active[coordinate].indices == parameter.indices)
                            parameter_shift[coordinate] += step * sign * parameter.value;
                }
            }
            return product;
        }();
        const auto shifted_parameters = p4_shifted_parameters(parameters, active, parameter_shift);
        const epcsaft::HelmholtzModel shifted(shifted_parameters, {}, thermochemistry);
        const auto shifted_state = shifted.state_at_pressure(
            signed_temperature, signed_pressure, signed_composition, branch);
        const auto shifted_signature = p4_take_signature(
            capture, thermochemistry, signed_temperature);
        const double closure = std::abs(shifted_state.pressure - signed_pressure) / signed_pressure;
        minimum_slope = std::min(minimum_slope, shifted_state.pressure_density_derivative);
        maximum_closure = std::max(maximum_closure, closure);
        // Density closure is owned in reduced units; relative pressure error is amplified by 1/Z.
        const double reduced_closure = std::abs(shifted_state.pressure - signed_pressure)
            / (shifted_state.molar_density * epcsaft::GAS_CONSTANT_J_PER_MOL_K * signed_temperature);
        CHECK(std::isfinite(shifted_state.pressure_density_derivative)
            && shifted_state.pressure_density_derivative > 0.0
            && reduced_closure <= epcsaft::policy::DENSITY_ROOT_RESIDUAL.value);
        CHECK(shifted_state.density_branch == central.density_branch
            && shifted_state.stable_root_count == central.stable_root_count
            && shifted_signature == central_signature);
        for (std::size_t i = 0; i < composition.size(); ++i)
            CHECK((composition[i] > 0.0) == (signed_composition[i] > 0.0));
        const auto residuals = p2_association_residuals(signed_temperature,
            1.0 / shifted_state.molar_density, signed_composition,
            shifted_parameters, shifted_state.site_fractions);
        for (const auto residual : residuals) {
            CHECK(std::isfinite(residual) && std::abs(residual) <= 2.0e-13);
            maximum_site_residual = std::max(maximum_site_residual, std::abs(residual));
        }
        const auto conditioning = p2_association_conditioning(signed_temperature,
            shifted_state.molar_density, signed_composition,
            shifted_parameters, shifted_state.site_fractions);
        CHECK(conditioning.rank == shifted_state.site_fractions.size()
            && conditioning.pivot_ratio > epcsaft::policy::ASSOCIATION_JACOBIAN_CONDITION.value);
        minimum_pivot_ratio = std::min(minimum_pivot_ratio, conditioning.pivot_ratio);
        result += sign_product * p4_state_value(shifted_state, observable);
    }
    if (admission != nullptr) {
        admission->insert(admission->end(), {step, minimum_slope, maximum_closure,
            maximum_site_residual, minimum_pivot_ratio,
            static_cast<double>(central.density_branch),
            central.stable_root_count.has_value()
                ? static_cast<double>(*central.stable_root_count) : -1.0});
    }
    return result / std::pow(2.0 * step, static_cast<int>(directions.size()));
}

// EOS sections 6--7: charged-coordinate higher actions on a regular positive-ionic-strength state.
void electrolyte_parameter_actions() {
    using namespace epcsaft;
    auto parameters = properties_support::aqueous_electrolyte();
    parameters.born_c_shell = parameters.born_c_dielectric = 1.0;
    const std::vector<ActiveParameter> active = {
        {ParameterFamily::BornDiameter, {1}}, {ParameterFamily::DebyeHuckelDiameter, {1}},
        {ParameterFamily::PackingDiameter, {1}}, {ParameterFamily::IonSuppressionCoefficient, {1}}};
    const std::vector<double> x = {0.98, 0.01, 0.01};
    constexpr double T = 298.15, rho = 55000.0;
    const HelmholtzModel model(parameters, active);
    const auto observables = std::vector<PropertyObservable>{PropertyObservable::ResidualEnthalpy};
    for (std::size_t coordinate = 0; coordinate < active.size(); ++coordinate) {
        const DensityActionDirection parameter{0.0, 0.0, {0.0, 0.0, 0.0},
            {{active[coordinate].family, active[coordinate].indices, 1.0}}};
        for (std::size_t order = 2; order <= 3; ++order) {
            std::vector<DensityActionDirection> directions(order, parameter);
            if (coordinate != 0) directions[1] = {1.0, 0.0, {0.0, 0.0, 0.0}, {}};
            const auto result = model.density_actions(T, rho, x, observables, directions).results.front();
            CHECK(result.status == PropertyActionStatus::Available && result.action);
            const std::string id = "electrolyte-coordinate/" + std::to_string(coordinate) + "/A" + std::to_string(order);
            double reference, uncertainty = 0.0, allowance;
            std::string inputs;
            if (coordinate == 0) {
                // Independently evaluated corrected radial charging formula, h_B = -K/d_B.
                const double epsilon_bulk = parameters.relative_permittivity[0]
                    / (1.0 + x[1] * parameters.ion_suppression_coefficient[1]
                        + x[2] * parameters.ion_suppression_coefficient[2]);
                double f_mix = 0.0;
                for (std::size_t i = 0; i < x.size(); ++i) f_mix += x[i] * parameters.solvation_factor[i];
                const double shell_ratio = 1.0 + (f_mix - 1.0) / std::abs(parameters.charge_number[1]);
                const double outer = (1.0 - 1.0 / epsilon_bulk) / shell_ratio;
                const double inner = (1.0 - 1.0 / parameters.ionic_region_relative_permittivity) * (1.0 - 1.0 / shell_ratio);
                const double prefactor = GAS_CONSTANT_J_PER_MOL_K * ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C
                    * x[1] * parameters.charge_number[1] * parameters.charge_number[1]
                    / (4.0 * PI * VACUUM_PERMITTIVITY_F_PER_M * BOLTZMANN_J_PER_K * ANGSTROM_M);
                const double d = parameters.born_diameter_angstrom[1];
                const double derivative_factor = order == 2 ? -2.0 / std::pow(d, 3) : 6.0 / std::pow(d, 4);
                reference = derivative_factor * prefactor * (inner + outer);
                const double source_scale = std::abs(derivative_factor * prefactor) * (std::abs(inner) + std::abs(outer));
                allowance = 1.0e-10 * std::abs(reference) + 32.0 * std::numeric_limits<double>::epsilon() * source_scale;
                inputs = recording_support::fields({recording_support::field("K_J_A_per_mol", prefactor * (inner + outer)),
                    recording_support::field("d_B_A", d), recording_support::field("source_scale", source_scale)});
            } else {
                auto lower = directions;
                lower.pop_back();
                const auto stencil = [&](double h) {
                    const auto value = [&](double shift) {
                        std::vector<double> delta(active.size(), 0.0);
                        delta[coordinate] = order == 3 ? shift : 0.0;
                        const HelmholtzModel rebuilt(p4_shifted_parameters(parameters, active, delta), active);
                        const auto action = rebuilt.density_actions(order == 2 ? T + shift : T, rho, x, observables, lower).results.front();
                        CHECK(action.status == PropertyActionStatus::Available && action.action);
                        return *action.action;
                    };
                    return (value(h) - value(-h)) / (2.0 * h);
                };
                // A 0.01 suppression step loses the small T/parameter action
                // to cancellation in the preceding action; resolve it at larger nodes.
                const double step = coordinate == 3 && order == 3 ? 4.0 : 0.4;
                const double coarse = stencil(step), middle = stencil(step / 2.0), fine = stencil(step / 4.0);
                const double rich_coarse = (4.0 * middle - coarse) / 3.0;
                reference = (4.0 * fine - middle) / 3.0;
                uncertainty = std::abs(reference - rich_coarse);
                allowance = (order == 2 ? 5.0e-5 : 5.0e-4) * std::abs(reference);
                inputs = recording_support::fields({recording_support::field("coarse", coarse),
                    recording_support::field("middle", middle), recording_support::field("fine", fine),
                    recording_support::field("reference_uncertainty", uncertainty)});
            }
            const bool ok = std::abs(*result.action - reference) <= allowance
                && uncertainty <= allowance && std::abs(reference) > allowance;
            recording_support::compare_custom(id.c_str(), "enthalpy_parameter_action", "native", *result.action, reference,
                recording_support::fields({inputs, recording_support::field("T_K", T), recording_support::field("rho_mol_per_m3", rho),
                    recording_support::field("x", x), recording_support::field("allowance", allowance)}),
                "Relationship", coordinate == 0 ? "Born-charging-identity" : "preceding-action-Richardson-stencil",
                "nonzero second and third actions in charged EOS coordinates", "magnitude_scaled_with_resolved_reference", ok);
        }
    }
}

epcsaft::Parameters p4_association_stress(std::size_t sites, double q) {
    constexpr double temperature = 300.0;
    constexpr double density = 1000.0;
    auto parameters = properties_support::ethane();
    parameters.site_component.assign(sites, 0);
    parameters.site_multiplicity.assign(sites, 1);
    parameters.association_active.assign(sites * sites, 0);
    parameters.association_energy_over_k.assign(sites * sites, 0.0);
    parameters.association_volume.assign(sites * sites, 0.0);
    parameters.association_active[0] = 1;
    parameters.association_volume[0] = 0.1;
    parameters.association_energy_over_k[0] = temperature * std::log(2.0);
    const auto sigma = epcsaft::segment_sigma(temperature, parameters);
    const auto diameter = epcsaft::segment_diameter(temperature, parameters, sigma);
    const double unit_delta = epcsaft::association_delta(
        temperature, density, {1.0}, parameters, sigma, diameter)[0];
    const double number_density = density * epcsaft::NUMBER_DENSITY_PER_ANGSTROM3;
    parameters.association_energy_over_k[0] = temperature
        * std::log1p(q / (number_density * unit_delta));
    return parameters;
}

double p4_residual_enthalpy_reference(const epcsaft::Parameters& parameters) {
    constexpr double temperature = 300.0;
    constexpr double volume = 1.0e-3;
    constexpr double temperature_step = 1.0e-3;
    constexpr double volume_step = 1.0e-8;
    epcsaft::HelmholtzModel model(parameters);
    const auto potential = [&](double t, double v) {
        return model.evaluate(t, v, {1.0}, epcsaft::DerivativeOrder::Value).value;
    };
    const double temperature_derivative = (
        potential(temperature + temperature_step, volume)
        - potential(temperature - temperature_step, volume)) / (2.0 * temperature_step);
    const double volume_derivative = (
        potential(temperature, volume + volume_step)
        - potential(temperature, volume - volume_step)) / (2.0 * volume_step);
    return -epcsaft::GAS_CONSTANT_J_PER_MOL_K * temperature
        * (temperature * temperature_derivative + volume * volume_derivative);
}

void properties_p4() {
    const std::vector<epcsaft::ActiveParameter> active = {
        {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}},
        {epcsaft::ParameterFamily::AssociationVolume, {1, 4}},
    };
    epcsaft::HelmholtzModel model(
        properties_support::associating_binary(), active, p3_binary_record());
    const std::vector<double> composition = {0.25, 0.75};
    const auto observables = std::vector<epcsaft::PropertyObservable>{
        epcsaft::PropertyObservable::IdealEnthalpy,
        epcsaft::PropertyObservable::IdealInternalEnergy,
        epcsaft::PropertyObservable::IdealEntropy,
        epcsaft::PropertyObservable::IdealIsobaricHeatCapacity,
        epcsaft::PropertyObservable::IdealIsochoricHeatCapacity,
        epcsaft::PropertyObservable::ResidualEnthalpy,
        epcsaft::PropertyObservable::ResidualInternalEnergy,
        epcsaft::PropertyObservable::ResidualEntropy,
        epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity,
        epcsaft::PropertyObservable::ResidualIsochoricHeatCapacity,
        epcsaft::PropertyObservable::TotalEnthalpy,
        epcsaft::PropertyObservable::TotalInternalEnergy,
        epcsaft::PropertyObservable::TotalEntropy,
        epcsaft::PropertyObservable::TotalIsobaricHeatCapacity,
        epcsaft::PropertyObservable::TotalIsochoricHeatCapacity,
    };
    const std::array<epcsaft::DensityActionDirection, 3> directions = {{
        {0.2, 5.0, {-0.01, 0.01}, {
            {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}, 0.03},
            {epcsaft::ParameterFamily::AssociationVolume, {1, 4}, 0.001},
        }},
        {-0.1, -3.0, {0.005, -0.005}, {
            {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}, -0.02},
            {epcsaft::ParameterFamily::AssociationVolume, {1, 4}, -0.0005},
        }},
        {0.05, 2.0, {0.002, -0.002}, {
            {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}, 0.01},
            {epcsaft::ParameterFamily::AssociationVolume, {1, 4}, 0.0002},
        }},
    }};
    const auto state = model.state(313.15, 42441.09427799721, composition);
    const auto pressure_root = model.state_at_pressure(
        313.15, state.pressure, composition, epcsaft::DensityBranch::Unspecified);
    const std::array<epcsaft::PressureActionDirection, 3> pressure_directions_all = {{
        {directions[0].temperature, 15000.0, directions[0].composition, directions[0].parameter},
        {directions[1].temperature, -10000.0, directions[1].composition, directions[1].parameter},
        {directions[2].temperature, 6000.0, directions[2].composition, directions[2].parameter},
    }};
    for (std::size_t order = 1; order <= directions.size(); ++order) {
        const std::vector<epcsaft::DensityActionDirection> request(
            directions.begin(), directions.begin() + static_cast<std::ptrdiff_t>(order));
        const auto batch = model.density_actions(
            313.15, 42441.09427799721, composition, observables, request);
        CHECK(batch.results.size() == observables.size());
        const auto reference = property_action_reference(model.parameters(), active, p3_binary_record(),
            313.15, 42441.09427799721, composition, observables, request);
        for (std::size_t index = 0; index < observables.size(); ++index) {
            const auto& result = batch.results[index];
            const bool available = result.status == epcsaft::PropertyActionStatus::Available
                && result.value.has_value() && result.action.has_value()
                && std::isfinite(*result.value) && std::isfinite(*result.action);
            const double fine = reference.value[index], middle = reference.middle[index], coarse = reference.coarse[index];
            const double tolerance = std::array<double, 3>{5.0e-6, 5.0e-5, 5.0e-4}[order - 1];
            const bool action_ok = available && std::abs(*result.action - fine)
                <= tolerance * reference.scale[index];
            const bool refinement_ok = reference.uncertainty[index] <= tolerance * reference.scale[index];
            recording_support::compare_custom(
                (std::string("direct/action/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "action", "native",
                available ? *result.action : std::numeric_limits<double>::quiet_NaN(), fine,
                recording_support::fields({recording_support::field("observable", index),
                    recording_support::field("order", order),
                    recording_support::field("central", p4_state_value(state, observables[index])),
                    recording_support::field("energy_order", result.energy_order),
                    recording_support::field("reference_fine", fine),
                    recording_support::field("reference_middle", middle),
                    recording_support::field("reference_coarse", coarse),
                    recording_support::field("tolerance", tolerance),
                    recording_support::field("term_scale", reference.scale[index]),
                    recording_support::field("reference_uncertainty", reference.uncertainty[index]),
                    recording_support::field("allowance", tolerance * reference.scale[index])}),
                "Relationship", "central-difference-refinement", "direct action matches independent refinement",
                "tolerance_times_chart_term_scale", action_ok, "p4-direct",
                (std::string("action/") + std::to_string(index) + "/" + std::to_string(order)).c_str(),
                "action", "T,rho,x,theta", "available", "observable,order");
            recording_support::compare_custom(
                (std::string("direct/refinement/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "reference_refinement", "native",
                fine, middle,
                recording_support::fields({recording_support::field("observable", index),
                    recording_support::field("order", order),
                    recording_support::field("reference_fine", fine),
                    recording_support::field("reference_middle", middle),
                    recording_support::field("reference_coarse", coarse),
                    recording_support::field("tolerance", tolerance),
                    recording_support::field("term_scale", reference.scale[index]),
                    recording_support::field("reference_uncertainty", reference.uncertainty[index]),
                    recording_support::field("allowance", tolerance * reference.scale[index])}),
                "Relationship", "central-difference-ladder", "direct reference refines",
                "tolerance_times_chart_term_scale", refinement_ok, "p4-direct",
                (std::string("refinement/") + std::to_string(index) + "/" + std::to_string(order)).c_str(),
                "reference", "T,rho,x,theta,h", "refined", "observable,order");
        }
    }
    // Pressure rows: the 15 calorics, ln(phi_i) per component (#73), then the density action d^k rho.
    auto pressure_observables = observables;
    pressure_observables.push_back(epcsaft::PropertyObservable::LogFugacityCoefficient);
    for (std::size_t order = 1; order <= directions.size(); ++order) {
        const std::vector<epcsaft::PressureActionDirection> request(
            pressure_directions_all.begin(),
            pressure_directions_all.begin() + static_cast<std::ptrdiff_t>(order));
        const auto batch = model.pressure_actions(
            313.15, state.pressure, composition, pressure_root.density_branch,
            pressure_observables, request);
        const auto reference = property_action_reference(model.parameters(), active, p3_binary_record(),
            313.15, state.pressure, composition, pressure_observables, request, pressure_root.density_branch);
        CHECK(reference.value.size() == observables.size() + composition.size() + 1);
        for (std::size_t index = 0; index < reference.value.size(); ++index) {
            const bool density = index == batch.results.size();
            const auto& result = batch.results[std::min(index, batch.results.size() - 1)];
            const bool component = !density && result.component == (index < observables.size()
                ? std::nullopt : std::optional<std::size_t>(index - observables.size()));
            const auto action = density ? (order == 1 ? std::optional<double>(batch.density_action.at(0))
                : batch.density_mixed_action) : result.action;
            const double fine = reference.value[index], middle = reference.middle[index], coarse = reference.coarse[index];
            const double tolerance = std::array<double, 3>{5.0e-6, 5.0e-5, 5.0e-4}[order - 1];
            const bool action_ok = (density ? batch.density_status : result.status) == epcsaft::PropertyActionStatus::Available
                && (density || component) && action.has_value() && std::abs(*action - fine)
                    <= tolerance * reference.scale[index];
            const bool refinement_ok = reference.uncertainty[index] <= tolerance * reference.scale[index];
            recording_support::compare_custom(
                (std::string("pressure/action/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "action", "native",
                action.value_or(std::numeric_limits<double>::quiet_NaN()), fine,
                recording_support::fields({recording_support::field("observable", index),
                    recording_support::field("order", order), recording_support::field("central",
                        index < observables.size() ? p4_state_value(pressure_root, observables[index])
                        : density ? pressure_root.molar_density
                        : (*pressure_root.log_fugacity_coefficient)[index - observables.size()].value()),
                    recording_support::field("energy_order", density ? order : result.energy_order),
                    recording_support::field("reference_fine", fine),
                    recording_support::field("reference_middle", middle),
                    recording_support::field("reference_coarse", coarse),
                    recording_support::field("tolerance", tolerance),
                    recording_support::field("term_scale", reference.scale[index]),
                    recording_support::field("reference_uncertainty", reference.uncertainty[index]),
                    recording_support::field("allowance", tolerance * reference.scale[index])}),
                "Relationship", "central-difference-refinement",
                "fixed-pressure action matches independent refinement",
                "tolerance_times_chart_term_scale", action_ok,
                "p4-pressure", (std::string("action/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "action", "T,P,x,theta", "available",
                "observable,order");
            recording_support::compare_custom(
                (std::string("pressure/refinement/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "reference_refinement", "native", fine,
                middle, recording_support::fields({recording_support::field("observable", index),
                    recording_support::field("order", order),
                    recording_support::field("reference_fine", fine),
                    recording_support::field("reference_middle", middle),
                    recording_support::field("reference_coarse", coarse),
                    recording_support::field("tolerance", tolerance),
                    recording_support::field("term_scale", reference.scale[index]),
                    recording_support::field("reference_uncertainty", reference.uncertainty[index]),
                    recording_support::field("allowance", tolerance * reference.scale[index])}),
                "Relationship", "central-difference-ladder", "pressure reference refines",
                "tolerance_times_chart_term_scale", refinement_ok, "p4-pressure",
                (std::string("refinement/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "reference", "T,P,x,theta,h", "refined",
                "observable,order");
        }
    }
    const std::array<epcsaft::PropertyObservable, 5> thermal_observables = {
        epcsaft::PropertyObservable::IdealEnthalpy,
        epcsaft::PropertyObservable::IdealInternalEnergy,
        epcsaft::PropertyObservable::IdealEntropy,
        epcsaft::PropertyObservable::IdealIsobaricHeatCapacity,
        epcsaft::PropertyObservable::IdealIsochoricHeatCapacity};
    const std::vector<epcsaft::PropertyObservable> thermal_request_observables(
        thermal_observables.begin(), thermal_observables.end());
    const double gas_constant = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const std::array<std::pair<const char*, epcsaft::IdealCorrelation>, 5> thermal_forms = {{
        {"polynomial", epcsaft::IdealCorrelation{
            epcsaft::IdealPolynomial{{30.0, 0.02, 1.0e-4, 1.0e-4}, 300.0, 1000.0, 50.0},
            1.0e5}},
        {"shomate", epcsaft::IdealCorrelation{
            epcsaft::IdealShomate{{30.0, 20.0, 10.0, 5.0, 2.0, 0.0, 0.0, 0.0}, 1000.0},
            1.0e5}},
        {"nasa7", epcsaft::IdealCorrelation{
            epcsaft::IdealNASA7{{30.0 / gas_constant, 0.02 / gas_constant,
                2.0e-5 / gas_constant, 1.0e-7 / gas_constant, 1.0e-9 / gas_constant,
                0.0, 0.0}}, 1.0e5}},
        {"nasa9", epcsaft::IdealCorrelation{
            epcsaft::IdealNASA9{{0.0, 0.0, 30.0 / gas_constant, 0.02 / gas_constant,
                2.0e-5 / gas_constant, 1.0e-7 / gas_constant, 1.0e-9 / gas_constant,
                0.0, 0.0}}, 1.0e5}},
        {"dippr107", epcsaft::IdealCorrelation{
            epcsaft::IdealDIPPR107{{0.0, 30.0, 250.0, 10.0, 600.0},
                300.0, 1000.0, 50.0}, 1.0e5}}
    }};
    const std::array<epcsaft::DensityActionDirection, 3> thermal_directions = {{
        {20.0, 0.0, {0.0}, {}}, {-10.0, 0.0, {0.0}, {}}, {5.0, 0.0, {0.0}, {}}}};
    for (const auto& [name, correlation] : thermal_forms) {
        const auto thermal_record = p3_record(correlation);
        epcsaft::HelmholtzModel thermal_model(
            properties_support::ethane(), {}, thermal_record);
        for (std::size_t order = 1; order <= thermal_directions.size(); ++order) {
            const std::vector<epcsaft::DensityActionDirection> request(
                thermal_directions.begin(), thermal_directions.begin() + static_cast<std::ptrdiff_t>(order));
            const auto batch = thermal_model.density_actions(
                320.0, 1000.0, {1.0}, thermal_request_observables, request);
            for (std::size_t index = 0; index < thermal_observables.size(); ++index) {
                const double fine = p4_reference(
                    thermal_model.parameters(), {}, thermal_record, 320.0, 1000.0,
                    {1.0}, request, thermal_observables[index],
                    std::array<double, 3>{1.0e-3, 3.0e-3, 1.0e-2}[order - 1]);
                const bool nonconstant_cp = order != 3 || index != 3 || std::abs(fine) > 5.0e-4;
                const bool ok = batch.results[index].status == epcsaft::PropertyActionStatus::Available
                    && batch.results[index].action.has_value()
                    && close(*batch.results[index].action, fine, 5.0e-4)
                    && nonconstant_cp;
                recording_support::compare_custom(
                    (std::string("supplemental/thermal-form/") + name + "/"
                        + std::to_string(index) + "/" + std::to_string(order)).c_str(),
                    "action", "native",
                    batch.results[index].action.value_or(std::numeric_limits<double>::quiet_NaN()), fine,
                    recording_support::fields({recording_support::field("form", name),
                        recording_support::field("observable", index), recording_support::field("order", order),
                        recording_support::field("reference", fine),
                        recording_support::field("nonconstant_cp", nonconstant_cp ? 1 : 0)}),
                    "Relationship", "central-difference-refinement", "thermal form action",
                    "relative", ok, "p4-supplemental", (std::string("thermal-form/") + name
                        + "/" + std::to_string(index) + "/" + std::to_string(order)).c_str(),
                    "action", "T,rho,x", "available", "form,observable,order");
            }
        }
    }
    const auto multiroot_directions = std::array<epcsaft::PressureActionDirection, 2>{
        epcsaft::PressureActionDirection{0.1, 2000.0, {0.0}, {}},
        epcsaft::PressureActionDirection{-0.05, -1000.0, {0.0}, {}}};
    epcsaft::HelmholtzModel ethane_model(properties_support::ethane());
    for (const auto& [site_count, q, expected_status, name]
         : std::array<std::tuple<std::size_t, double, epcsaft::PropertyActionStatus,
             const char*>, 2>{{
             {2, 1.0e26, epcsaft::PropertyActionStatus::AssociationIllConditioned,
              "association-ill-conditioned"},
             {128, 1.0e27, epcsaft::PropertyActionStatus::AssociationSingular,
              "association-singular"}}}) {
        const auto parameters = p4_association_stress(site_count, q);
        const auto sigma = epcsaft::segment_sigma(300.0, parameters);
        const auto diameter = epcsaft::segment_diameter(300.0, parameters, sigma);
        const auto sites = epcsaft::solve_site_fractions(
            300.0, 1000.0, {1.0}, parameters, sigma, diameter);
        const double exact_site = 2.0 / (1.0 + std::sqrt(1.0 + 4.0 * q));
        const auto residuals = p2_association_residuals(
            300.0, 1.0e-3, {1.0}, parameters, sites);
        const auto conditioning = p2_association_conditioning(
            300.0, 1000.0, {1.0}, parameters, sites);
        epcsaft::HelmholtzModel probe(
            parameters, {}, p3_record(p3_forms()[0].second));
        const auto batch = probe.density_actions(300.0, 1000.0, {1.0},
            {epcsaft::PropertyObservable::ResidualEnthalpy,
             epcsaft::PropertyObservable::TotalEnthalpy,
             epcsaft::PropertyObservable::IdealEnthalpy},
            {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
        const double expected_residual_enthalpy = p4_residual_enthalpy_reference(parameters);
        const double maximum_residual = std::accumulate(residuals.begin(), residuals.end(), 0.0,
            [](double value, double residual) { return std::max(value, std::abs(residual)); });
        const bool inactive_exact = std::all_of(sites.begin() + 1, sites.end(),
            [](double value) { return value == 1.0; });
        const bool ok = batch.results.size() == 3 && sites.size() == site_count
            && close(sites[0], exact_site, 2.0e-12) && inactive_exact
            && maximum_residual <= epcsaft::policy::ASSOCIATION_RESIDUAL.value
            && batch.results[0].status == expected_status
            && batch.results[0].value.has_value() && !batch.results[0].action.has_value()
            && close(*batch.results[0].value, expected_residual_enthalpy, 2.0e-6)
            && batch.results[1].status == expected_status
            && batch.results[1].value.has_value() && !batch.results[1].action.has_value()
            && batch.results[2].status == epcsaft::PropertyActionStatus::Available
            && batch.results[2].value.has_value() && batch.results[2].action.has_value()
            && close(*batch.results[1].value,
                *batch.results[0].value + *batch.results[2].value, 2.0e-12)
            && batch.diagnostics.site_rank == conditioning.rank;
        recording_support::compare_custom(
            (std::string("supplemental/") + name).c_str(), "status_and_preservation", "bool",
            ok ? 1.0 : 0.0, 1.0,
            recording_support::fields({recording_support::field("site_count", site_count),
                recording_support::field("q", q),
                recording_support::field("site_fraction", sites[0]),
                recording_support::field("exact_site_fraction", exact_site),
                recording_support::field("maximum_primal_residual", maximum_residual),
                recording_support::field("rank", conditioning.rank),
                recording_support::field("pivot_ratio", conditioning.pivot_ratio),
                recording_support::field("residual_enthalpy", batch.results[0].value.value_or(
                    std::numeric_limits<double>::quiet_NaN())),
                recording_support::field("reference_residual_enthalpy", expected_residual_enthalpy),
                recording_support::field("central_ideal_enthalpy",
                    batch.results[2].value.value_or(std::numeric_limits<double>::quiet_NaN()))}),
            "Availability", "analytic-association-stress",
            "typed association failure preserves central residual, total and ideal values", "absolute", ok,
            "p4-supplemental", name, "status", "T,rho,x,q", "typed failure");
    }
    const auto missing_thermochemistry = ethane_model.density_actions(
        300.0, 1000.0, {1.0},
        {epcsaft::PropertyObservable::IdealEnthalpy},
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool missing_ok = missing_thermochemistry.results[0].status
        == epcsaft::PropertyActionStatus::MissingThermochemistry;
    recording_support::compare_custom("supplemental/missing-thermochemistry", "status", "bool",
        missing_ok ? 1.0 : 0.0, 1.0, "ethane without source", "Availability", "thermo",
        "ideal source absence is typed", "absolute", missing_ok,
        "p4-supplemental", "missing-thermochemistry", "status", "thermo", "missing_thermochemistry");
    const auto partial = ethane_model.density_actions(
        300.0, 1000.0, {1.0},
        {epcsaft::PropertyObservable::IdealEnthalpy, epcsaft::PropertyObservable::ResidualEnthalpy},
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool partial_ok = partial.results[0].status == epcsaft::PropertyActionStatus::MissingThermochemistry
        && partial.results[1].status == epcsaft::PropertyActionStatus::Available
        && partial.results[1].value.has_value();
    recording_support::compare_custom("supplemental/partial-availability", "sibling_preserved", "bool",
        partial_ok ? 1.0 : 0.0, 1.0, "ideal+residual without source", "Availability", "partial",
        "residual sibling remains available", "absolute", partial_ok,
        "p4-supplemental", "partial-availability", "sibling_preserved", "sibling", "partial");
    const auto negative_slope = ethane_model.density_actions(
        280.0, 4000.0, {1.0},
        {epcsaft::PropertyObservable::ResidualEnthalpy,
         epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity,
         epcsaft::PropertyObservable::ResidualIsochoricHeatCapacity},
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool negative_ok = negative_slope.results[0].status == epcsaft::PropertyActionStatus::Available
        && negative_slope.results[1].status == epcsaft::PropertyActionStatus::PressureSingular
        && negative_slope.results[0].value.has_value()
        && negative_slope.results[2].status == epcsaft::PropertyActionStatus::Available
        && negative_slope.results[2].action.has_value();
    recording_support::compare_custom("supplemental/negative-slope", "status", "bool",
        negative_ok ? 1.0 : 0.0, 1.0, "ethane T=280,rho=4000", "Availability", "pressure",
        "negative pressure slope suppresses residual cp only", "absolute", negative_ok,
        "p4-supplemental", "negative-slope", "status", "pressure", "pressure_singular");
    const auto density_unavailable = ethane_model.pressure_actions(
        200.0, 217800.69494692178, {1.0}, epcsaft::DensityBranch::Unspecified,
        {epcsaft::PropertyObservable::ResidualEnthalpy},
        {epcsaft::PressureActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool density_unavailable_ok = density_unavailable.results.size() == 1
        && density_unavailable.results[0].status
            == epcsaft::PropertyActionStatus::DensityUnavailable;
    recording_support::compare_custom(
        "supplemental/density-unavailable", "status", "bool",
        density_unavailable_ok ? 1.0 : 0.0, 1.0,
        "ethane two-root state with unspecified branch", "Availability", "density-selection",
        "ambiguous fixed-pressure density remains unavailable", "absolute", density_unavailable_ok,
        "p4-supplemental", "density-unavailable", "status", "T,P,x,branch",
        "density_unavailable");
    const std::vector<epcsaft::ActiveParameter> entropy_active = {
        {epcsaft::ParameterFamily::AssociationVolume, {1, 4}}
    };
    const epcsaft::PressureActionDirection entropy_direction{
        0.0, 0.0, {0.0, 0.0}, {
            {epcsaft::ParameterFamily::AssociationVolume, {1, 4}, 1.0}
        }};
    const auto entropy_batch = model.pressure_actions(
        313.15, state.pressure, composition, pressure_root.density_branch,
        {epcsaft::PropertyObservable::IdealEntropy}, {entropy_direction});
    const epcsaft::HelmholtzModel entropy_model(
        model.parameters(), entropy_active, p3_binary_record());
    const auto entropy_derivatives = entropy_model.evaluate(
        313.15, 1.0 / state.molar_density, composition, epcsaft::DerivativeOrder::Hessian);
    const std::size_t entropy_parameter_column = 2 + composition.size();
    const std::size_t entropy_point_width = entropy_parameter_column + entropy_active.size();
    const double pressure_parameter = -epcsaft::GAS_CONSTANT_J_PER_MOL_K * 313.15
        * entropy_derivatives.hessian[entropy_point_width + entropy_parameter_column];
    const double entropy_chain = epcsaft::GAS_CONSTANT_J_PER_MOL_K * pressure_parameter
        / (state.molar_density * state.pressure_density_derivative);
    const double entropy_fine = p4_pressure_reference(
        model.parameters(), entropy_active, p3_binary_record(), 313.15, state.pressure,
        composition, {entropy_direction}, pressure_root.density_branch,
        epcsaft::PropertyObservable::IdealEntropy, 1.0e-3);
    const double entropy_middle = p4_pressure_reference(
        model.parameters(), entropy_active, p3_binary_record(), 313.15, state.pressure,
        composition, {entropy_direction}, pressure_root.density_branch,
        epcsaft::PropertyObservable::IdealEntropy, 3.0e-3);
    const double entropy_coarse = p4_pressure_reference(
        model.parameters(), entropy_active, p3_binary_record(), 313.15, state.pressure,
        composition, {entropy_direction}, pressure_root.density_branch,
        epcsaft::PropertyObservable::IdealEntropy, 1.0e-2);
    const bool entropy_ok = entropy_batch.results.size() == 1
        && entropy_batch.results[0].status == epcsaft::PropertyActionStatus::Available
        && entropy_batch.results[0].action.has_value()
        && close(*entropy_batch.results[0].action, entropy_chain, 5.0e-6)
        && close(*entropy_batch.results[0].action, entropy_fine, 5.0e-6)
        && close(entropy_fine, entropy_middle, 1.0e-5);
    recording_support::compare_custom(
        "supplemental/parameter-ideal-entropy", "action", "J/(mol K)",
        entropy_batch.results.empty() ? std::numeric_limits<double>::quiet_NaN()
            : entropy_batch.results[0].action.value_or(std::numeric_limits<double>::quiet_NaN()),
        entropy_fine,
        recording_support::fields({recording_support::field("T", 313.15),
            recording_support::field("P", state.pressure),
            recording_support::field("rho", state.molar_density),
            recording_support::field("P_theta", pressure_parameter),
            recording_support::field("chain", entropy_chain),
            recording_support::field("reference_fine", entropy_fine),
            recording_support::field("reference_middle", entropy_middle),
            recording_support::field("reference_coarse", entropy_coarse)}),
        "Relationship", "implicit-pressure-chain", "parameter Ideal entropy retains P_theta/(rho P_rho)",
        "relative", entropy_ok, "p4-supplemental", "parameter-ideal-entropy", "action",
        "T,P,x,theta", "available", "reference_fine,reference_middle,reference_coarse");
    for (const auto branch : {epcsaft::DensityBranch::Vapor, epcsaft::DensityBranch::Liquid}) {
        const auto root = ethane_model.state_at_pressure(
            200.0, 217800.69494692178, {1.0}, branch);
        const std::vector<epcsaft::PropertyObservable> rows = {
            epcsaft::PropertyObservable::ResidualEnthalpy,
            epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity};
        for (std::size_t order = 1; order <= 2; ++order) {
            const std::vector<epcsaft::PressureActionDirection> request(
                multiroot_directions.begin(), multiroot_directions.begin() + static_cast<std::ptrdiff_t>(order));
            const auto branch_batch = ethane_model.pressure_actions(
                200.0, 217800.69494692178, {1.0}, branch, rows, request);
            for (std::size_t index = 0; index < rows.size(); ++index) {
                const double fine = p4_pressure_reference(
                    properties_support::ethane(), {}, std::nullopt, 200.0,
                    217800.69494692178, {1.0}, request, branch, rows[index],
                    std::array<double, 2>{3.0e-3, 1.0e-2}[order - 1]);
                const double coarse = p4_pressure_reference(
                    properties_support::ethane(), {}, std::nullopt, 200.0,
                    217800.69494692178, {1.0}, request, branch, rows[index],
                    std::array<double, 2>{1.0e-2, 3.0e-2}[order - 1]);
                const double tolerance = std::array<double, 2>{5.0e-6, 5.0e-5}[order - 1];
                const bool action_ok = branch_batch.results[index].status == epcsaft::PropertyActionStatus::Available
                    && branch_batch.results[index].action.has_value()
                    && close(*branch_batch.results[index].action, fine, tolerance);
                const bool refinement_ok = close(fine, coarse, 2.0 * tolerance);
                const std::string point = std::string(branch == epcsaft::DensityBranch::Vapor ? "vapor" : "liquid")
                    + "/" + std::to_string(index) + "/" + std::to_string(order);
                recording_support::compare_custom(("supplemental/multiroot/action/" + point).c_str(),
                    "action", "native",
                    branch_batch.results[index].action.value_or(std::numeric_limits<double>::quiet_NaN()), fine,
                    recording_support::fields({recording_support::field("branch", point), recording_support::field("order", order)}),
                    "Relationship", "central-difference-refinement", "multiroot pressure action",
                    "relative", action_ok, "p4-supplemental", ("multiroot/action/" + point).c_str(),
                    "action", "T,P,x", "available", "branch,observable,order");
                recording_support::compare_custom(("supplemental/multiroot/refinement/" + point).c_str(),
                    "reference_refinement", "native", fine, coarse,
                    recording_support::fields({recording_support::field("branch", point), recording_support::field("order", order)}),
                    "Relationship", "central-difference-ladder", "multiroot reference refines",
                    "relative", refinement_ok, "p4-supplemental", ("multiroot/refinement/" + point).c_str(),
                    "reference", "T,P,x,h", "refined", "branch,observable,order");
            }
        }
        CHECK(root.density_branch == branch);
    }
    const std::array<epcsaft::DensityActionDirection, 3> electrolyte_directions = {{
        {0.1, 1.0, {-0.02, 0.01, 0.01}, {}},
        {-0.05, -0.5, {0.01, -0.005, -0.005}, {}},
        {0.03, 0.2, {0.004, -0.002, -0.002}, {}}}};
    const std::vector<epcsaft::PropertyObservable> electrolyte_observables = {
        epcsaft::PropertyObservable::ResidualEnthalpy,
        epcsaft::PropertyObservable::ResidualEntropy,
        epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity};
    epcsaft::HelmholtzModel electrolyte_model(properties_support::aqueous_electrolyte());
    for (std::size_t order = 1; order <= electrolyte_directions.size(); ++order) {
        const std::vector<epcsaft::DensityActionDirection> request(
            electrolyte_directions.begin(), electrolyte_directions.begin() + static_cast<std::ptrdiff_t>(order));
        const auto batch = electrolyte_model.density_actions(
            300.0, 1000.0, {0.98, 0.01, 0.01}, electrolyte_observables, request);
        for (std::size_t index = 0; index < electrolyte_observables.size(); ++index) {
            const double fine = p4_reference(
                electrolyte_model.parameters(), {}, std::nullopt, 300.0, 1000.0,
                {0.98, 0.01, 0.01}, request, electrolyte_observables[index],
                std::array<double, 3>{1.0e-3, 3.0e-3, 1.0e-2}[order - 1]);
            const bool ok = batch.results[index].status == epcsaft::PropertyActionStatus::Available
                && batch.results[index].action.has_value()
                && close(*batch.results[index].action, fine,
                    std::array<double, 3>{5.0e-6, 5.0e-5, 5.0e-4}[order - 1]);
            const std::string point = std::to_string(index) + "/" + std::to_string(order);
            recording_support::compare_custom(("supplemental/electrolyte/" + point).c_str(),
                "action", "native",
                batch.results[index].action.value_or(std::numeric_limits<double>::quiet_NaN()), fine,
                recording_support::fields({recording_support::field("observable", index), recording_support::field("order", order)}),
                "Relationship", "central-difference-refinement", "charge-balanced electrolyte action",
                "relative", ok, "p4-supplemental", ("electrolyte/" + point).c_str(),
                "action", "T,rho,x", "available", "observable,order");
        }
    }
    epcsaft::ThermochemistryRecord electrolyte_record;
    electrolyte_record.standard_pressure = 1.0e5;
    for (std::size_t component = 0; component < 3; ++component) {
        electrolyte_record.component_intervals.push_back({
            epcsaft::IdealInterval{250.0, 350.0, true, true, p3_forms()[0].second}
        });
    }
    const std::vector<double> electrolyte_composition = {0.98, 0.01, 0.01};
    const epcsaft::DensityActionDirection electrolyte_total_direction{
        0.1, 1.0, {-0.02, 0.01, 0.01}, {}};
    epcsaft::HelmholtzModel electrolyte_with_ideal(
        properties_support::aqueous_electrolyte(), {}, electrolyte_record);
    const auto electrolyte_total = electrolyte_with_ideal.density_actions(
        300.0, 1000.0, electrolyte_composition,
        {epcsaft::PropertyObservable::TotalEnthalpy}, {electrolyte_total_direction});
    const std::optional<epcsaft::ThermochemistryRecord> electrolyte_record_optional = electrolyte_record;
    const double electrolyte_total_fine = p4_reference(
        electrolyte_with_ideal.parameters(), {}, electrolyte_record_optional, 300.0, 1000.0,
        electrolyte_composition, {electrolyte_total_direction},
        epcsaft::PropertyObservable::TotalEnthalpy, 1.0e-3);
    const double electrolyte_total_middle = p4_reference(
        electrolyte_with_ideal.parameters(), {}, electrolyte_record_optional, 300.0, 1000.0,
        electrolyte_composition, {electrolyte_total_direction},
        epcsaft::PropertyObservable::TotalEnthalpy, 3.0e-3);
    const double electrolyte_total_coarse = p4_reference(
        electrolyte_with_ideal.parameters(), {}, electrolyte_record_optional, 300.0, 1000.0,
        electrolyte_composition, {electrolyte_total_direction},
        epcsaft::PropertyObservable::TotalEnthalpy, 1.0e-2);
    const bool electrolyte_total_ok = electrolyte_total.results.size() == 1
        && electrolyte_total.results[0].status == epcsaft::PropertyActionStatus::Available
        && electrolyte_total.results[0].value.has_value()
        && electrolyte_total.results[0].action.has_value()
        && close(*electrolyte_total.results[0].action, electrolyte_total_fine, 5.0e-6)
        && close(electrolyte_total_fine, electrolyte_total_middle, 1.0e-5);
    recording_support::compare_custom(
        "supplemental/electrolyte-total", "action", "J/mol", electrolyte_total.results.empty()
            ? std::numeric_limits<double>::quiet_NaN()
            : electrolyte_total.results[0].action.value_or(std::numeric_limits<double>::quiet_NaN()), electrolyte_total_fine,
        recording_support::fields({recording_support::field("T", 300.0),
            recording_support::field("rho", 1000.0), recording_support::field("x", electrolyte_composition),
            recording_support::field("reference_fine", electrolyte_total_fine),
            recording_support::field("reference_middle", electrolyte_total_middle),
            recording_support::field("reference_coarse", electrolyte_total_coarse)}),
        "Relationship", "central-difference-refinement", "sourced Ideal plus residual total action",
        "relative", electrolyte_total_ok, "p4-supplemental", "electrolyte-total", "action",
        "T,rho,x", "available", "reference_fine,reference_middle,reference_coarse");
    const auto full_one = model.density_actions(
        313.15, 42441.09427799721, composition, observables, {directions[0]});
    const auto selected = model.density_actions(
        313.15, 42441.09427799721, composition,
        {observables[5], observables[13]}, {directions[0]});
    const auto result_matches = [](const epcsaft::PropertyActionResult& left,
                                   const epcsaft::PropertyActionResult& right) {
        if (left.observable != right.observable || left.status != right.status
            || left.energy_order != right.energy_order
            || left.value.has_value() != right.value.has_value()
            || left.action.has_value() != right.action.has_value()) return false;
        return (!left.value.has_value() || close(*left.value, *right.value, 2.0e-12))
            && (!left.action.has_value() || close(*left.action, *right.action, 2.0e-12));
    };
    const bool selected_h = selected.results.size() == 2 && full_one.results.size() == observables.size()
        && result_matches(selected.results[0], full_one.results[5]);
    const bool selected_cp = selected.results.size() == 2 && full_one.results.size() == observables.size()
        && result_matches(selected.results[1], full_one.results[13]);
    CHECK(selected_h && selected_cp);
    recording_support::compare_custom("supplemental/selected-h", "selected_full", "bool", selected_h ? 1.0 : 0.0,
        1.0, "residual_h", "Identity", "selected/full", "selected equals full", "absolute", selected_h,
        "p4-supplemental", "selected-h", "selected_full", "subset", "selected/full");
    recording_support::compare_custom("supplemental/selected-total-cp", "selected_full", "bool", selected_cp ? 1.0 : 0.0,
        1.0, "total_cp", "Identity", "selected/full", "selected equals full", "absolute", selected_cp,
        "p4-supplemental", "selected-total-cp", "selected_full", "subset", "selected/full");
    const auto full_missing = ethane_model.density_actions(
        300.0, 1000.0, {1.0}, observables,
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool failure_match = partial.results.size() == 2 && full_missing.results.size() == observables.size()
        && result_matches(partial.results[0], full_missing.results[0])
        && result_matches(partial.results[1], full_missing.results[5]);
    recording_support::compare_custom(
        "supplemental/selected-failure", "selected_full", "bool",
        failure_match ? 1.0 : 0.0, 1.0, "ideal+residual without source",
        "Identity", "selected/full", "selected preserves typed failure and valid sibling",
        "absolute", failure_match, "p4-supplemental", "selected-failure", "selected_full",
        "subset", "selected/full");
    const auto single_h = model.density_actions(
        313.15, 42441.09427799721, composition, {observables[5]}, {directions[0]});
    const auto single_cp = model.density_actions(
        313.15, 42441.09427799721, composition, {observables[13]}, {directions[0]});
    const bool output_selection_return = single_h.results.size() == 1
        && single_cp.results.size() == 1
        && result_matches(single_h.results[0], full_one.results[5])
        && result_matches(single_cp.results[0], full_one.results[13]);
    recording_support::compare_custom(
        "supplemental/output-selection-return", "selected_full", "bool",
        output_selection_return ? 1.0 : 0.0, 1.0,
        "residual-only then ideal+residual selection at the same state",
        "Identity", "selected/full", "returning to ideal output selection preserves values and actions",
        "absolute", output_selection_return, "p4-supplemental", "output-selection-return",
        "selected_full", "subset", "selected/full");
    epcsaft::HelmholtzModel cold_model(
        model.parameters(), active, p3_binary_record());
    const auto cold_batch = cold_model.density_actions(
        313.15, 42441.09427799721, composition,
        {epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity}, {directions[0]});
    const auto warm_batch = cold_model.density_actions(
        313.15, 42441.09427799721, composition,
        {epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity}, {directions[0]});
    const bool work_reuse = cold_batch.results.size() == 1 && warm_batch.results.size() == 1
        && cold_batch.results[0].status == epcsaft::PropertyActionStatus::Available
        && result_matches(warm_batch.results[0], cold_batch.results[0]);
    recording_support::compare_custom(
        "supplemental/work-reuse", "warm_equals_cold", "bool", work_reuse ? 1.0 : 0.0, 1.0,
        "q3 residual cp, repeated request", "Identity", "shared-batch-work",
        "a repeated request on the recorded tape reproduces the cold batch",
        "absolute", work_reuse, "p4-supplemental", "work-reuse", "warm_equals_cold", "case", "warm=cold");
    const auto state_values_match = [](
        const epcsaft::State& left, const epcsaft::State& right,
        const std::vector<epcsaft::PropertyObservable>& rows, double tolerance) {
        if (!close(left.molar_density, right.molar_density, tolerance)
            || !close(left.pressure, right.pressure, tolerance)
            || !close(left.pressure_density_derivative, right.pressure_density_derivative, tolerance))
            return false;
        for (const auto observable : rows) {
            const double left_value = p4_state_value(left, observable);
            const double right_value = p4_state_value(right, observable);
            if (std::isnan(left_value) && std::isnan(right_value)) continue;
            if (!std::isfinite(left_value) || !std::isfinite(right_value)
                || !close(left_value, right_value, tolerance)) return false;
        }
        return true;
    };
    const std::vector<double> pure_first = {1.0, 0.0};
    const std::vector<epcsaft::PropertyObservable> freshness_rows = {
        epcsaft::PropertyObservable::ResidualEnthalpy,
        epcsaft::PropertyObservable::ResidualIsochoricHeatCapacity,
        epcsaft::PropertyObservable::IdealEntropy};
    for (std::size_t order = 1; order <= 3; ++order) {
        const std::vector<epcsaft::DensityActionDirection> request(order, {0.2, 5.0, {0.0, 0.0}, {}});
        for (std::size_t point = 0; point < 3; ++point) {
            const bool mixed = point == 1;
            const auto& x = mixed ? composition : pure_first;
            const double t = point == 1 ? 320.0 : point == 2 ? 307.15 : 313.15;
            const double rho = point == 1 ? 6000.0 : point == 2 ? 4500.0 : 5000.0;
            const auto retained = model.density_actions(t, rho, x, freshness_rows, request);
            const auto fresh = epcsaft::HelmholtzModel(model.parameters(), active, p3_binary_record())
                .density_actions(t, rho, x, freshness_rows, request);
            for (std::size_t row = 0; row < freshness_rows.size(); ++row) {
                const auto& actual = retained.results[row];
                const auto& expected = fresh.results[row];
                const bool agrees = actual.status == epcsaft::PropertyActionStatus::Available
                    && expected.status == actual.status && actual.value && expected.value
                    && actual.action && expected.action
                    && close(*actual.value, *expected.value, 2.0e-12)
                    && close(*actual.action, *expected.action, 2.0e-12);
                const auto id = "supplemental/action-freshness/" + std::to_string(order)
                    + "/" + std::to_string(point) + "/" + std::to_string(row);
                recording_support::compare_custom(id.c_str(), "action", "native",
                    actual.action.value_or(std::numeric_limits<double>::quiet_NaN()),
                    expected.action.value_or(std::numeric_limits<double>::quiet_NaN()),
                    recording_support::fields({recording_support::field("T", t),
                        recording_support::field("rho", rho), recording_support::field("x", x),
                        recording_support::field("order", order), recording_support::field("row", row),
                        recording_support::field("value", actual.value.value_or(std::numeric_limits<double>::quiet_NaN())),
                        recording_support::field("fresh_value", expected.value.value_or(std::numeric_limits<double>::quiet_NaN()))}),
                    "Identity", "fresh-model-action", "retained action refreshes all state dependencies",
                    "relative", agrees, "p4-supplemental", id.c_str(), "action", "T,rho,x", "fresh-model");
            }
        }
    }
    const auto retained_pure = model.state(313.15, 5000.0, pure_first);
    const auto retained_mixed = model.state(320.0, 6000.0, composition);
    const auto retained_pure_again = model.state(313.15, 5000.0, pure_first);
    const auto fresh_pure = epcsaft::HelmholtzModel(
        model.parameters(), active, p3_binary_record()).state(
            313.15, 5000.0, pure_first);
    const auto fresh_mixed = epcsaft::HelmholtzModel(
        model.parameters(), active, p3_binary_record()).state(320.0, 6000.0, composition);
    const auto fresh_pure_again = epcsaft::HelmholtzModel(
        model.parameters(), active, p3_binary_record()).state(
            313.15, 5000.0, pure_first);
    const bool retained_freshness = state_values_match(retained_pure, fresh_pure, observables, 2.0e-12)
        && state_values_match(retained_mixed, fresh_mixed, observables, 2.0e-12)
        && state_values_match(retained_pure_again, fresh_pure_again, observables, 2.0e-12)
        && close(p4_state_value(retained_pure, epcsaft::PropertyObservable::IdealEntropy),
                 p4_state_value(retained_pure_again, epcsaft::PropertyObservable::IdealEntropy), 2.0e-12);
    recording_support::compare_custom(
        "supplemental/retained-freshness", "state", "bool", retained_freshness ? 1.0 : 0.0,
        1.0, recording_support::fields({recording_support::field("pure_x", pure_first),
            recording_support::field("mixed_x", composition), recording_support::field("T_pure", 313.15),
            recording_support::field("T_mixed", 320.0),
            recording_support::field("rho_pure", 5000.0),
            recording_support::field("rho_mixed", 6000.0),
            recording_support::field("pure_entropy", p4_state_value(retained_pure,
                epcsaft::PropertyObservable::IdealEntropy)),
            recording_support::field("mixed_entropy", p4_state_value(retained_mixed,
                epcsaft::PropertyObservable::IdealEntropy)),
            recording_support::field("repeat_entropy", p4_state_value(retained_pure_again,
                epcsaft::PropertyObservable::IdealEntropy))}),
        "Identity", "retained-state-replay", "changing T and x refreshes ideal support and state values",
        "absolute", retained_freshness, "p4-supplemental", "retained-freshness", "state",
        "T,rho,x", "fresh-state", "ideal_entropy,all_properties");
    const std::vector<epcsaft::PropertyObservable> residual_rows = {
        epcsaft::PropertyObservable::ResidualEnthalpy,
        epcsaft::PropertyObservable::ResidualInternalEnergy,
        epcsaft::PropertyObservable::ResidualEntropy,
        epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity,
        epcsaft::PropertyObservable::ResidualIsochoricHeatCapacity};
    bool branch_route_ok = true;
    std::vector<double> branch_densities;
    std::vector<double> branch_pressures;
    for (const auto branch : {epcsaft::DensityBranch::Vapor, epcsaft::DensityBranch::Liquid}) {
        const auto pressure_state = ethane_model.state_at_pressure(
            200.0, 217800.69494692178, {1.0}, branch);
        const auto direct_state = ethane_model.state(
            200.0, pressure_state.molar_density, {1.0});
        const auto fresh_state = epcsaft::HelmholtzModel(properties_support::ethane()).state_at_pressure(
            200.0, 217800.69494692178, {1.0}, branch);
        branch_densities.push_back(pressure_state.molar_density);
        branch_pressures.push_back(pressure_state.pressure);
        branch_route_ok = branch_route_ok
            && pressure_state.density_branch == branch
            && fresh_state.density_branch == branch
            && pressure_state.stable_root_count == fresh_state.stable_root_count
            && state_values_match(pressure_state, fresh_state, residual_rows, 2.0e-12)
            && state_values_match(pressure_state, direct_state, residual_rows, 1.0e-9);
    }
    recording_support::compare_custom(
        "supplemental/branch-route", "state", "bool", branch_route_ok ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", 200.0),
            recording_support::field("P", 217800.69494692178),
            recording_support::field("rho", branch_densities),
            recording_support::field("achieved_P", branch_pressures)}),
        "Identity", "resolved-state-route", "same branch pressure and density routes agree",
        "absolute", branch_route_ok, "p4-supplemental", "branch-route", "state",
        "T,P,x,branch", "residual_properties,branch,root_count");
    const auto fixed = model.density_actions(
        313.15, 42441.09427799721, composition, {observables[5]},
        {epcsaft::DensityActionDirection{0.0, 0.0, {0.0, 0.0},
            {{epcsaft::ParameterFamily::AssociationVolume, {0}, 1.0}}}});
    const bool fixed_zero = fixed.results[0].status == epcsaft::PropertyActionStatus::Available
        && fixed.results[0].action.has_value() && std::abs(*fixed.results[0].action) <= 1.0e-12;
    recording_support::compare_custom("supplemental/fixed-slot", "fixed_slot_zero", "bool",
        fixed_zero ? 1.0 : 0.0, 1.0, "association_volume[0]", "Availability", "fixed-slot",
        "known fixed slot is available structural zero", "absolute", fixed_zero,
        "p4-supplemental", "fixed-slot", "fixed_slot_zero", "slot", "fixed/zero");
    bool unknown_rejected = false;
    try {
        (void)model.density_actions(313.15, 42441.09427799721, composition, {observables[5]},
            {epcsaft::DensityActionDirection{0.0, 0.0, {0.0, 0.0},
                {{epcsaft::ParameterFamily::AssociationVolume, {999}, 1.0}}}});
    } catch (const std::invalid_argument&) { unknown_rejected = true; }
    recording_support::compare_custom("supplemental/unknown-slot", "unknown_rejected", "bool",
        unknown_rejected ? 1.0 : 0.0, 1.0, "association_volume[999]", "Availability", "unknown-slot",
        "unknown native slot rejects", "absolute", unknown_rejected,
        "p4-supplemental", "unknown-slot", "unknown_rejected", "slot", "invalid_argument");
    const auto rejects = [](auto&& call) {
        try { call(); }
        catch (const std::invalid_argument&) { return true; }
        return false;
    };
    const auto& valid_direction = directions[0];
    std::size_t malformed_rejections = 0;
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {}, {valid_direction}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0], observables[0]},
        {valid_direction}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]}, {}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]},
        std::vector<epcsaft::DensityActionDirection>(4, valid_direction)); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]},
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]},
        {epcsaft::DensityActionDirection{std::numeric_limits<double>::infinity(),
            0.0, {0.0, 0.0}, {}}}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]},
        {epcsaft::DensityActionDirection{0.0, 0.0, {0.1, 0.0}, {}}}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]},
        {epcsaft::DensityActionDirection{0.0, 0.0, {0.0, 0.0},
            {{epcsaft::ParameterFamily::Kij, {2, 1}, 1.0}}}}); });
    epcsaft::HelmholtzModel ionic_validation(properties_support::aqueous_electrolyte());
    malformed_rejections += rejects([&] { ionic_validation.density_actions(
        300.0, 1000.0, {0.98, 0.01, 0.01},
        {epcsaft::PropertyObservable::ResidualEnthalpy},
        {epcsaft::DensityActionDirection{0.0, 0.0, {-0.01, 0.01, 0.0}, {}}}); });
    const bool malformed_ok = malformed_rejections == 9;
    recording_support::compare_custom(
        "supplemental/malformed-requests", "invalid_argument", "count",
        static_cast<double>(malformed_rejections), 9.0,
        "empty,duplicate,order,dimension,nonfinite,balance,slot,charge",
        "Availability", "request-validation", "malformed requests fail at the trust boundary",
        "absolute", malformed_ok, "p4-supplemental", "malformed-requests",
        "invalid_argument", "request", "nine rejected cases");
    // A neutral insertion is its exact infinite-dilution limit (#73): f(0) = 2 f(eps) - f(2 eps) + O(eps^2).
    epcsaft::HelmholtzModel neutral(properties_support::neutral_binary());
    const auto inserted = [&](double eps) {
        return neutral.density_actions(300.0, 1000.0, {1.0 - eps, eps}, {observables[5]},
            {epcsaft::DensityActionDirection{0.0, 0.0, {-0.01, 0.01}, {}}}).results[0];
    };
    const auto support = inserted(0.0);
    const double support_limit = 2.0 * inserted(1.0e-7).action.value() - inserted(2.0e-7).action.value();
    const bool support_change = support.status == epcsaft::PropertyActionStatus::Available && support.action
        && std::abs(*support.action - support_limit) <= std::max(1.0e-12, 1.0e-9 * std::abs(support_limit));
    recording_support::compare_custom("supplemental/support-change", "action", "J/mol",
        support.action.value_or(std::numeric_limits<double>::quiet_NaN()), support_limit,
        "neutral x=(1,0), dx=(-.01,.01), eps=1e-7", "Relationship", "infinite-dilution-extrapolation",
        "a neutral insertion action is its infinite-dilution limit", "symmetric_close_1e-9_1e-12", support_change,
        "p4-supplemental", "support-change", "action", "support", "support_limit");
    epcsaft::HelmholtzModel aqueous(properties_support::aqueous_electrolyte());
    const auto ionic = aqueous.density_actions(300.0, 1000.0, {1.0, 0.0, 0.0}, {observables[5]},
        {epcsaft::DensityActionDirection{0.0, 0.0, {-0.02, 0.01, 0.01}, {}}});
    const bool ionic_insertion = ionic.results[0].status == epcsaft::PropertyActionStatus::IonicInsertion;
    recording_support::compare_custom("supplemental/ionic-insertion", "status", "bool",
        ionic_insertion ? 1.0 : 0.0, 1.0, "aqueous x=(1,0,0), dx=(-.02,.01,.01)", "Availability", "ionic",
        "ionic insertion unavailable", "absolute", ionic_insertion,
        "p4-supplemental", "ionic-insertion", "status", "ionic", "ionic_insertion");
    epcsaft::HelmholtzModel endpoint_model(
        properties_support::ethane(), {}, p3_record(p3_forms()[0].second));
    const auto endpoint_action = endpoint_model.density_actions(250.0, 1000.0, {1.0},
        {epcsaft::PropertyObservable::IdealEnthalpy},
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool thermal_boundary = endpoint_action.results[0].status == epcsaft::PropertyActionStatus::ThermalBoundary;
    recording_support::compare_custom("supplemental/thermal-boundary", "status", "bool",
        thermal_boundary ? 1.0 : 0.0, 1.0, "T=250 endpoint, dT=.1", "Availability", "thermal",
        "thermal boundary action unavailable", "absolute", thermal_boundary,
        "p4-supplemental", "thermal-boundary", "status", "thermal", "thermal_boundary");
    const auto counter_batch = model.density_actions(313.15, 42441.09427799721, composition,
        {epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity},
        {directions[0], directions[1], directions[2]});
    const auto pressure_counter_batch = model.pressure_actions(
        313.15, state.pressure, composition, pressure_root.density_branch,
        {epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity},
        {pressure_directions_all[0], pressure_directions_all[1], pressure_directions_all[2]});
    const bool counter_ok = counter_batch.results.size() == 1
        && pressure_counter_batch.results.size() == 1
        && counter_batch.results[0].status == epcsaft::PropertyActionStatus::Available
        && pressure_counter_batch.results[0].status == epcsaft::PropertyActionStatus::Available
        && pressure_counter_batch.diagnostics.pressure_density_derivative > 0.0
        && pressure_counter_batch.diagnostics.stable_root_count == pressure_root.stable_root_count
        && pressure_counter_batch.diagnostics.branch == pressure_root.density_branch;
    recording_support::compare_custom("supplemental/counters", "counter_contract", "bool",
        counter_ok ? 1.0 : 0.0, 1.0, "q3 residual cp", "Structure", "counters",
        "pressure-chart actions carry the state's root branch, root count and positive slope", "absolute", counter_ok,
        "p4-supplemental", "counters", "counter_contract", "counters", "E5/no-Hessian");
}

void neutral_property_pressure_hessian() {
    auto parameters = properties_support::ethane();
    parameters.segment_count[0] = 2.002;
    parameters.sigma_constant[0] = 3.6184;
    parameters.dispersion_energy_over_k[0] = 208.11;
    parameters.molar_mass_kg_per_mol[0] = 0.044096;
    epcsaft::HelmholtzModel model(parameters);
    constexpr double temperature = 400.0;
    constexpr double density = 60.0;
    const auto inputs = model_inputs(parameters, temperature, 1.0 / density, {1.0});
    const auto actual = model.neutral_property(
        temperature, density, epcsaft::NeutralProperty::Pressure,
        epcsaft::DerivativeOrder::Hessian
    );
    const auto value_after_hessian = model.neutral_property(
        temperature, density, epcsaft::NeutralProperty::Pressure,
        epcsaft::DerivativeOrder::Value
    );
    const auto state = model.state(temperature, density, {1.0});
    recording_support::compare(
        "issue10-pressure/value", "pressure", "Pa", actual.value, state.pressure, inputs,
        "Relationship", "state", "neutral pressure equals state pressure",
        "reference_or_one", 0.0, 1.0e-12,
        close(actual.value, state.pressure, 1.0e-12), "issue10-pressure", "value",
        "P", "T,rho", "P_neutral=P_state", "T,rho"
    );
    recording_support::compare(
        "issue10-pressure/value-after-hessian", "pressure", "Pa",
        value_after_hessian.value, actual.value, inputs,
        "Relationship", "order-transition", "value remains valid after Hessian",
        "reference_or_one", 0.0, 1.0e-12,
        close(value_after_hessian.value, actual.value, 1.0e-12),
        "issue10-pressure", "value-after-hessian", "P", "T,rho",
        "P_value=P_hessian", "T,rho"
    );
    const double delta = 1.0e-4;
    const auto plus = model.neutral_property(
        temperature + 400.0 * delta, density + 100.0 * delta,
        epcsaft::NeutralProperty::Pressure, epcsaft::DerivativeOrder::Gradient
    );
    const auto minus = model.neutral_property(
        temperature - 400.0 * delta, density - 100.0 * delta,
        epcsaft::NeutralProperty::Pressure, epcsaft::DerivativeOrder::Gradient
    );
    constexpr const char* units[] = {"Pa/K", "Pa m3/mol"};
    for (std::size_t row = 0; row < 2; ++row) {
        const double centered = (plus.gradient[row] - minus.gradient[row]) / (2.0 * delta);
        const double expected = actual.hessian[row] * 400.0
            + actual.hessian[row * 2 + 1] * 100.0;
        recording_support::compare(
            (std::string("issue10-hessian/") + std::to_string(row)).c_str(),
            "pressure_hessian_action", units[row], centered, expected,
            recording_support::fields({inputs, recording_support::field("delta", delta),
                                       recording_support::field("direction", row)}),
            "Relationship", "centered_gradient", "pressure Hessian action",
            "reference_or_one", 0.0, 1.0e-6,
            close(centered, expected, 1.0e-6), "issue10-hessian",
            std::to_string(row).c_str(), "H*d", "T,rho", "H*d=finite_difference", "T,rho"
        );
    }
}


// #566: State's Joule-Thomson coefficient against two independent pressure-chart action routes:
// (T v_T - v)/c_p with v_T from the volume action, and -(dh/dP)_T/(dh/dT)_P from total-enthalpy actions.
void joule_thomson_identities() {
    using namespace epcsaft;
    const HelmholtzModel model(properties_support::neutral_binary(), {}, properties_support::p3_binary_record());
    const double T = 300.0;
    const std::vector<double> x{0.4, 0.6};
    using Point = std::tuple<const char*, double, DensityBranch>;
    for (const auto& [name, P, branch] : std::array<Point, 2>{{
             {"vapor", 4.87e5, DensityBranch::Vapor}, {"dense", 9.52e6, DensityBranch::Liquid}}}) {
        const auto state = model.state_at_pressure(T, P, x, branch, std::nullopt, true);
        const auto actions = [&, P = P, branch = branch](const PressureActionDirection& direction) {
            return model.pressure_actions(T, P, x, branch,
                {PropertyObservable::MolarVolumeTemperatureDerivative, PropertyObservable::TotalEnthalpy}, {direction});
        };
        const auto t = actions({1.0, 0.0, {0.0, 0.0}, {}}), p = actions({0.0, 1.0, {0.0, 0.0}, {}});
        const double mu = state.joule_thomson_coefficient.value();
        const std::array<std::pair<const char*, double>, 2> routes{{
            {"volume", (T * t.results[0].value.value() - 1.0 / state.molar_density)
                / state.total_isobaric_heat_capacity.value()},
            {"enthalpy", -p.results[1].action.value() / t.results[1].action.value()}}};
        const auto inputs = recording_support::fields({recording_support::field("T_K", T),
            recording_support::field("P_Pa", P), recording_support::field("x", x),
            recording_support::field("rho", state.molar_density)});
        for (const auto& [route, reference] : routes)
            recording_support::compare(("joule-thomson/" + std::string(name) + "/" + route).c_str(),
                "joule_thomson_coefficient", "K/Pa", mu, reference, inputs, "Relationship", "exact-identity",
                "State mu_JT equals the isenthalpic pressure-action identity", "symmetric_close", 0.0, 1.0e-12,
                recording_support::scalar_verdict(mu, reference, "symmetric_close", 0.0, 1.0e-12), "joule-thomson");
    }
}

}  // namespace properties_checks
