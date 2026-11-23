#include "properties/helmholtz_model.hpp"
#include "domain_reports.hpp"

#include "equilibrium_record.hpp"
#include "fields.hpp"
#include "solver_options.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace nb = nanobind;
using epcsaft::PROPERTY_SCOPE;

namespace {

template <class Getter>
nb::list state_column(const epcsaft::StateBatch& batch, Getter get) {
    nb::list result;
    for (const auto& row : batch.rows) result.append(row ? nb::cast(get(*row)) : nb::none());
    return result;
}
template <class Value>
void state_field(nb::class_<epcsaft::State>& scalar, nb::class_<epcsaft::StateBatch>& batch,
    const char* name, Value epcsaft::State::*member) {
    scalar.def_ro(name, member);
    batch.def_prop_ro(name, [member](const epcsaft::StateBatch& value) {
        return state_column(value, [member](const epcsaft::State& row) { return row.*member; });
    });
}
#define STATE_FIELD(T, name) state_field(state, batch, #name, &T::name);

// Physical input records: every field is required and transported unchanged; the native
// model copies them at construction and native admission owns their consistency.
void register_thermal_inputs(nb::module_& module) {
    using namespace epcsaft;
    module.attr("GAS_CONSTANT_J_PER_MOL_K") = GAS_CONSTANT_J_PER_MOL_K;
    nb::class_<IdealPolynomial>(module, "IdealPolynomial")
        .def(nb::init<std::vector<double>, double, double, double>(), nb::arg("coefficients"), nb::arg("reference_temperature"),
            nb::arg("reference_enthalpy"), nb::arg("reference_entropy"))
        EPCSAFT_READWRITE(IdealPolynomial, coefficients, reference_temperature, reference_enthalpy, reference_entropy);
    nb::class_<IdealShomate>(module, "IdealShomate")
        .def(nb::init<std::array<double, 8>, double>(), nb::arg("coefficients"), nb::arg("formation_enthalpy"))
        EPCSAFT_READWRITE(IdealShomate, coefficients, formation_enthalpy);
    nb::class_<IdealNASA7>(module, "IdealNASA7")
        .def(nb::init<std::array<double, 7>>(), nb::arg("coefficients"))
        .def_rw("coefficients", &IdealNASA7::coefficients);
    nb::class_<IdealNASA9>(module, "IdealNASA9")
        .def(nb::init<std::array<double, 9>>(),
             nb::arg("coefficients"))
        .def_rw("coefficients", &IdealNASA9::coefficients);
    nb::class_<IdealDIPPR107>(module, "IdealDIPPR107")
        .def(nb::init<std::array<double, 5>, double, double, double>(), nb::arg("coefficients"), nb::arg("reference_temperature"),
            nb::arg("reference_enthalpy"), nb::arg("reference_entropy"))
        EPCSAFT_READWRITE(IdealDIPPR107, coefficients, reference_temperature, reference_enthalpy, reference_entropy);
    nb::class_<IdealCorrelation>(module, "IdealCorrelation")
        .def(nb::init<decltype(IdealCorrelation::form), double>(), nb::arg("form"), nb::arg("source_pressure"))
        EPCSAFT_READWRITE(IdealCorrelation, form, source_pressure);
    nb::class_<IdealInterval>(module, "IdealInterval")
        .def(nb::init<double, double, bool, bool, IdealCorrelation>(), nb::arg("lower"), nb::arg("upper"),
            nb::arg("lower_inclusive"), nb::arg("upper_inclusive"), nb::arg("correlation"))
        EPCSAFT_READWRITE(IdealInterval, lower, upper, lower_inclusive, upper_inclusive, correlation);
    nb::class_<ThermochemistryRecord>(module, "ThermochemistryRecord")
        .def(nb::init<double, std::vector<std::vector<IdealInterval>>>(), nb::arg("standard_pressure"), nb::arg("component_intervals"))
        EPCSAFT_READWRITE(ThermochemistryRecord, standard_pressure) EPCSAFT_COPIED(ThermochemistryRecord, component_intervals);
    nb::class_<ReactionLogPolynomial>(module, "ReactionLogPolynomial")
        .def(nb::init<double, double, double, double, double, double, double, std::string, double, double>(), nb::arg("a"), nb::arg("b"), nb::arg("c"), nb::arg("d"),
            nb::arg("reference_temperature"), nb::arg("temperature_min"),
            nb::arg("temperature_max"), nb::arg("standard_state_id"), nb::arg("e") = 0.0, nb::arg("f") = 0.0)
        EPCSAFT_READWRITE(ReactionLogPolynomial, a, b, c, d, reference_temperature, temperature_min, temperature_max,
            standard_state_id, e, f);
    nb::enum_<ParameterFamily>(module, "ParameterFamily")
        EPCSAFT_VALUES(ParameterFamily, BornCDielectric, IonicRegionRelativePermittivity, BornDiameter,
            DebyeHuckelDiameter, PackingDiameter, IonSuppressionCoefficient, SegmentCount, SigmaConstant,
            DispersionEnergyOverK, Kij, Lij, AssociationEnergyOverK, AssociationVolume, KijReciprocalTemperatureSlope,
            DielectricVolume, DielectricTemperature, TemperatureLawCoefficient, CationPairMixing, CationTripleMixing, PathwayWeight);
    // Native slot indices of one coordinate; a reciprocal pair lists both of its slots.
    nb::class_<ActiveParameter>(module, "ActiveParameter")
        .def(nb::init<ParameterFamily, std::vector<std::size_t>>(), nb::arg("family"), nb::arg("indices"))
        EPCSAFT_READWRITE(ActiveParameter, family, indices);
}

}  // namespace

void register_eos(nb::module_& module) {
    using namespace epcsaft;
    using bindings::domain_point;
    using bindings::domain_result;
    register_thermal_inputs(module);
    nb::enum_<DomainQuantity>(module, "DomainQuantity")
        EPCSAFT_VALUES(DomainQuantity, SoluteMolality, IonicStrength, IonMoleFraction, Temperature, Pressure);
    nb::enum_<DomainLimit>(module, "DomainLimit") EPCSAFT_VALUES(DomainLimit, Minimum, Maximum);
    nb::class_<DomainBound>(module, "DomainBound")
        .def(nb::init<std::string, std::string, std::string, std::string, DomainQuantity, DomainLimit, std::size_t, std::vector<std::size_t>, double>());
    nb::class_<DomainBounds>(module, "DomainBounds");
    module.def("admit_domains", &admit_domains, nb::arg("parameters"), nb::arg("bounds"));
    nb::class_<OutsideDomain>(module, "OutsideDomain")
        EPCSAFT_READONLY(OutsideDomain, packet, domain_id, variable, side, unit, solute, value, bound, row, phase);
    module.def("outside_domains", [](const Parameters& p, const DomainBounds* bounds, const std::vector<double>& x,
        std::optional<double> temperature, std::optional<double> pressure) {
        std::vector<OutsideDomain> rows;
        outside_domains(rows, p, bounds, x, temperature, pressure, 0, 0, true);
        return rows;
    }, nb::arg("parameters"), nb::arg("bounds").none(), nb::arg("x"),
        nb::arg("temperature") = std::nullopt, nb::arg("pressure") = std::nullopt);
    nb::class_<epcsaft::ParameterEntry>(module, "ParameterEntry")
        .def(nb::init<>())
        EPCSAFT_READWRITE(epcsaft::ParameterEntry, identity, field, index, magnitude, unit);
    nb::class_<epcsaft::ParameterDocument>(module, "ParameterDocument")
        .def(nb::init<>())
        EPCSAFT_READWRITE(epcsaft::ParameterDocument, component_ids, site_ids, site_component, site_multiplicity, base_choice,
            association_choice, electrolyte_choice, ionic_dispersion_choice, permittivity_choice, polar_choice)
        EPCSAFT_COPIED(epcsaft::ParameterDocument, entries);
    module.def("build_parameters", &epcsaft::build_parameters, nb::arg("document"));
    module.def("prepare_thermochemistry", &epcsaft::prepare_thermochemistry, nb::arg("record"), nb::arg("source_gas_constant"));
    module.def("fitted_document", &epcsaft::fitted_document, nb::arg("start"), nb::arg("fitted"));

    nb::class_<epcsaft::Parameters>(module, "Parameters")
        .def_prop_ro("component_count", &epcsaft::Parameters::component_count)
        .def_prop_ro("site_count", &epcsaft::Parameters::site_count)
        .def_ro("ideal_gas", &epcsaft::Parameters::ideal_gas)
        .def_ro("molar_masses", &epcsaft::Parameters::molar_mass_kg_per_mol, "kg/mol, in component_ids order")
        .def_ro("charge_number", &epcsaft::Parameters::charge_number, "elementary charges, in component_ids order")
        .def_ro("temperature_law_coefficient_identities", &epcsaft::Parameters::temperature_law_coefficient_identity,
            "source identity of each TemperatureLawCoefficient slot");
    nb::enum_<epcsaft::DerivativeOrder>(module, "DerivativeOrder")
        EPCSAFT_VALUES(epcsaft::DerivativeOrder, Value, Gradient, PressureSlope, Hessian);
    nb::enum_<epcsaft::NeutralProperty>(module, "NeutralProperty")
        EPCSAFT_VALUES(epcsaft::NeutralProperty, ResidualHelmholtz, Pressure);
    nb::class_<epcsaft::NeutralPropertyDerivatives>(module, "NeutralPropertyDerivatives")
        EPCSAFT_READONLY(epcsaft::NeutralPropertyDerivatives, value, gradient, hessian);
    nb::enum_<epcsaft::PropertyObservable>(module, "PropertyObservable")
        EPCSAFT_VALUES(epcsaft::PropertyObservable, IdealEnthalpy, IdealInternalEnergy, IdealEntropy,
            IdealIsobaricHeatCapacity, IdealIsochoricHeatCapacity, ResidualEnthalpy, ResidualInternalEnergy,
            ResidualEntropy, ResidualIsobaricHeatCapacity, ResidualIsochoricHeatCapacity, TotalEnthalpy,
            TotalInternalEnergy, TotalEntropy, TotalIsobaricHeatCapacity, TotalIsochoricHeatCapacity,
            LogFugacityCoefficient);
    nb::enum_<epcsaft::PropertyActionStatus>(module, "PropertyActionStatus")
        EPCSAFT_VALUES(epcsaft::PropertyActionStatus, Available, MissingThermochemistry, OutsideInterval,
            ThermalBoundary, SupportChange, IonicInsertion, AssociationSingular, AssociationIllConditioned,
            AssociationResidual, DensityUnavailable, PressureSingular, CentralPropertyUnavailable, Nonfinite);
    // A native family/slot coordinate increment in its ALG-PARAMETERS unit.
    nb::class_<epcsaft::ParameterActionDirection>(module, "ParameterActionDirection")
        .def(nb::init<epcsaft::ParameterFamily, std::vector<std::size_t>, double>(), nb::arg("family"), nb::arg("indices"), nb::arg("value"))
        EPCSAFT_READONLY(epcsaft::ParameterActionDirection, family, indices, value);
    nb::class_<epcsaft::DensityActionDirection>(module, "DensityActionDirection")
        .def(nb::init<double, double, std::vector<double>, std::vector<epcsaft::ParameterActionDirection>>(), nb::arg("temperature") = 0.0, nb::arg("molar_density") = 0.0,
            nb::arg("composition") = std::vector<double>{},
            nb::arg("parameter") = std::vector<epcsaft::ParameterActionDirection>{})
        EPCSAFT_READONLY(epcsaft::DensityActionDirection, temperature, molar_density, composition, parameter);
    nb::class_<epcsaft::PressureActionDirection>(module, "PressureActionDirection")
        .def(nb::init<double, double, std::vector<double>, std::vector<epcsaft::ParameterActionDirection>>(), nb::arg("temperature") = 0.0, nb::arg("pressure") = 0.0,
            nb::arg("composition") = std::vector<double>{},
            nb::arg("parameter") = std::vector<epcsaft::ParameterActionDirection>{})
        EPCSAFT_READONLY(epcsaft::PressureActionDirection, temperature, pressure, composition, parameter);
    nb::class_<epcsaft::PropertyActionResult>(module, "PropertyActionResult")
        EPCSAFT_READONLY(epcsaft::PropertyActionResult, observable, component, value, action, status, energy_order);
    nb::enum_<epcsaft::DensityBranch>(module, "DensityBranch")
        EPCSAFT_VALUES(epcsaft::DensityBranch, Unspecified, Vapor, Liquid)
        .def("__int__", [](epcsaft::DensityBranch value) { return static_cast<int>(value); });
    using Diagnostics = epcsaft::PropertyActionDiagnostics;
    nb::class_<Diagnostics>(module, "PropertyActionDiagnostics")
        EPCSAFT_READONLY(Diagnostics, site_conditioning, site_rank, site_factorizations, site_solves, site_pivot_ratio,
            association_residual, pressure_density_derivative, pressure_residual, stability_determinant,
            pressure_factorizations, pressure_solves, thermal_evaluations, recordings, directional_paths,
            forward_sweeps, branch, returned_branch, stable_root_count, first_numerical_failure);
    nb::class_<epcsaft::ResidualEnergyAction>(module, "ResidualEnergyAction")
        EPCSAFT_READONLY(epcsaft::ResidualEnergyAction, action, status, diagnostics);
    nb::class_<epcsaft::PropertyActionBatch>(module, "PropertyActionBatch", nb::dynamic_attr())
        EPCSAFT_DOMAIN_VIOLATIONS
        EPCSAFT_SOLVER_OPTIONS(epcsaft::PropertyActionBatch, PROPERTY_SCOPE)
        EPCSAFT_READONLY(epcsaft::PropertyActionBatch, results, maximum_energy_order, diagnostics, molar_density,
            density_action, density_mixed_action, density_message, density_convergence, density_status);
    nb::enum_<epcsaft::ThermochemistryStatus>(module, "ThermochemistryStatus")
        EPCSAFT_VALUES(epcsaft::ThermochemistryStatus, MissingRecord, Available, OutsideInterval, Nonfinite);
    nb::class_<epcsaft::ReactionThermochemistry>(module, "ReactionThermochemistry")
        EPCSAFT_READONLY(epcsaft::ReactionThermochemistry, gibbs, enthalpy, heat_capacity, pressure_gibbs, ln_k_pressure, ln_k_density,
            k_pressure, k_density, gibbs_derivatives, status, cause);
    module.def("reaction_thermochemistry_density", &epcsaft::reaction_thermochemistry_density,
        nb::arg("record"), nb::arg("stoichiometry"), nb::arg("temperature"));
    nb::class_<epcsaft::VirialCoefficients>(module, "VirialCoefficients", nb::dynamic_attr())
        EPCSAFT_DOMAIN_VIOLATIONS
        EPCSAFT_READONLY(epcsaft::VirialCoefficients, second, third);
    nb::class_<epcsaft::LogFugacityTerms>(module, "LogFugacityTerms", nb::dynamic_attr())
        EPCSAFT_DOMAIN_VIOLATIONS
        EPCSAFT_READONLY(epcsaft::LogFugacityTerms, hard_chain, dispersion, association, debye_huckel, msa, born, dipolar,
            cation_mixing, log_compressibility_factor);
    nb::class_<epcsaft::VirialCoefficientActions>(module, "VirialCoefficientActions", nb::dynamic_attr())
        EPCSAFT_DOMAIN_VIOLATIONS
        EPCSAFT_READONLY(epcsaft::VirialCoefficientActions, second, third, second_action, third_action);
    nb::enum_<epcsaft::StateStatus>(module, "StateStatus")
        EPCSAFT_VALUES(epcsaft::StateStatus, Available, ParameterError, DeclarationError, ConvergenceError, DensityStructureError, EosError, DomainError);
    nb::class_<epcsaft::StateBatch> batch(module, "StateBatch", nb::dynamic_attr());
    batch EPCSAFT_DOMAIN_VIOLATIONS EPCSAFT_READONLY(epcsaft::StateBatch, rows, status, cause, failure_record)
        .def_prop_ro("first_numerical_failure", [](const epcsaft::StateBatch& value) {
            std::vector<std::optional<epcsaft::core::Convergence>> result;
            for (std::size_t k = 0; k < value.rows.size(); ++k)
                result.push_back(value.rows[k] ? value.rows[k]->first_numerical_failure : value.failure_first[k]);
            return result;
        });
    nb::class_<epcsaft::State> state(module, "State", "One homogeneous state in SI molar units (K, Pa, mol/m^3, J/mol); component lists in component_ids order.", nb::dynamic_attr());
    state EPCSAFT_DOMAIN_VIOLATIONS.def_prop_ro("returned_density_branch", &epcsaft::State::returned_density_branch)
        .def_prop_ro("density_branch_mismatch", &epcsaft::State::density_branch_mismatch)
        .def_ro("first_numerical_failure", &epcsaft::State::first_numerical_failure)
        EPCSAFT_SOLVER_OPTIONS(epcsaft::State, PROPERTY_SCOPE)
        EPCSAFT_REPR(State, "pressure", "molar_density", "compressibility_factor", "density_branch", "residual_helmholtz");
    state_field(state, batch, "requested_density_branch", &epcsaft::State::requested_density_branch);
    batch.def_prop_ro("returned_density_branch", [](const epcsaft::StateBatch& v) {
        return state_column(v, [](const epcsaft::State& r) { return r.returned_density_branch(); });
    }).def_prop_ro("density_branch_mismatch", [](const epcsaft::StateBatch& v) {
        return state_column(v, [](const epcsaft::State& r) { return r.density_branch_mismatch(); });
    }).def_prop_ro("solver_options", [](const epcsaft::StateBatch& v) {
        nb::list result;
        for (const auto& row : v.rows) result.append(row ? epcsaft::bindings::solver_options_view(row->solver_options, PROPERTY_SCOPE) : nb::none());
        return result;
    });
        EPCSAFT_EACH(STATE_FIELD, epcsaft::State, hard_chain, dispersion, association, debye_huckel, msa, born, dipolar, cation_mixing,
            bulk_relative_permittivity, residual_helmholtz, compressibility_factor, pressure,
            // None for a charged component on an ion-free support; see component_property_status (IonicInsertion).
            residual_chemical_potential_over_rt, component_property_status,
            log_fugacity_coefficient, log_fugacity_parameter_derivative, residual_entropy, residual_internal_energy,
            residual_enthalpy, residual_isochoric_heat_capacity, residual_isobaric_heat_capacity, ideal_enthalpy,
            ideal_internal_energy, ideal_entropy, ideal_isobaric_heat_capacity, ideal_isochoric_heat_capacity,
            total_enthalpy, total_internal_energy, total_entropy, total_isobaric_heat_capacity,
            total_isochoric_heat_capacity, molar_density, pressure_density_derivative, site_fractions, pathway_weight,
            density_convergence, thermochemistry_status,
            // Unspecified and absent in the T-rho chart: no pressure-root search was performed.
            density_branch, stable_root_count)
        EPCSAFT_EACH(STATE_FIELD, epcsaft::State, unavailable_properties, rank_records, isothermal_compressibility, isobaric_expansivity, heat_capacity_ratio,
            isentropic_compressibility, speed_of_sound, joule_thomson_coefficient, partial_molar_volume,
            partial_molar_enthalpy, partial_molar_entropy)
    nb::class_<epcsaft::DensityDiagnosticPoint>(module, "DensityDiagnosticPoint")
        EPCSAFT_READONLY(epcsaft::DensityDiagnosticPoint, molar_density, packing_fraction, pressure,
            pressure_density_derivative, lower, upper, residual, mechanical_sign, linked, declared, cause, stop);
    nb::class_<epcsaft::DensityDiagnosticBracket>(module, "DensityDiagnosticBracket")
        EPCSAFT_READONLY(epcsaft::DensityDiagnosticBracket, lower, upper, lower_residual, upper_residual);
    nb::class_<epcsaft::DensityDiagnosticRequest>(module, "DensityDiagnosticRequest")
        .def(nb::init<>())
        EPCSAFT_READWRITE(epcsaft::DensityDiagnosticRequest, nodes, density_tolerance, log_tolerance, finite_upper,
            max_evaluations, max_seconds, finite_anchor, refine, adaptive, evaluation_order);
    nb::class_<epcsaft::DensityDiagnostics>(module, "DensityDiagnostics")
        EPCSAFT_SOLVER_OPTIONS(epcsaft::DensityDiagnostics, PROPERTY_SCOPE)
        EPCSAFT_READONLY(epcsaft::DensityDiagnostics, samples, pressure_roots, spinodals, declared_brackets, unresolved,
            evaluation_order, finite_domain, geometric_cap, declared_crossings, evaluations, stop, first_numerical_failure);
    nb::class_<epcsaft::SweepCounters>(module, "SweepCounters")
        EPCSAFT_READONLY(epcsaft::SweepCounters, forward, reverse, hessian_calls, recordings, directional_paths, maximum_energy_order);
    nb::class_<epcsaft::HelmholtzModel>(module, "HelmholtzModel")
        .def("__init__", [](epcsaft::HelmholtzModel* self, epcsaft::Parameters p, std::vector<epcsaft::ActiveParameter> a,
                         std::optional<epcsaft::ThermochemistryRecord> t, nb::object options, std::shared_ptr<const DomainBounds> domains) {
            new (self) epcsaft::HelmholtzModel(std::move(p), std::move(a), std::move(t), epcsaft::bindings::solver_options(options, PROPERTY_SCOPE));
            self->domains = std::move(domains);
        },
             nb::arg("parameters"), nb::arg("active_parameters") = std::vector<epcsaft::ActiveParameter>{},
             nb::arg("thermochemistry") = std::nullopt, nb::arg("solver_options") = nb::none(), nb::arg("domains") = nb::none())
        .def_prop_ro("solver_options", [](const epcsaft::HelmholtzModel& v) { return epcsaft::bindings::solver_options_view(v.solver_options(), PROPERTY_SCOPE); })
        .def_prop_ro("active_parameter_values", &epcsaft::HelmholtzModel::active_parameter_values)
        .def("virial_coefficients", [](const HelmholtzModel& model, double t, const std::vector<double>& x) {
            return domain_point(model.virial_coefficients(t, x), model, x, t);
        }, nb::arg("temperature"), nb::arg("composition"))
        .def("virial_coefficient_actions", [](const HelmholtzModel& model, double t, const std::vector<double>& x,
            double dt, const std::vector<double>& dx) {
            return domain_point(model.virial_coefficient_actions(t, x, dt, dx), model, x, t);
        }, nb::arg("temperature"), nb::arg("composition"), nb::arg("temperature_direction"), nb::arg("composition_direction"))
        .def("state", [](const HelmholtzModel& model, double t, double rho, const std::vector<double>& x,
                         bool derived) {
            model.validate_electroneutrality(x);
            return domain_result(model.state(t, rho, x, derived), &model,
                [&](auto& rows, const State& r) { outside_domains(rows, model.parameters(), model.domains.get(), x, t, r.pressure); });
        }, nb::arg("temperature"), nb::arg("molar_density"), nb::arg("mole_fractions"),
             nb::arg("derived_properties") = false)
        .def("log_fugacity_terms", [](const HelmholtzModel& model, double t, double rho, const std::vector<double>& x) {
            return domain_point(model.log_fugacity_terms(t, rho, x), model, x, t, std::nullopt, rho);
        }, nb::arg("temperature"), nb::arg("molar_density"), nb::arg("mole_fractions"))
        .def("states", [](const HelmholtzModel& model, const std::vector<double>& t, const std::vector<double>& input,
            const std::vector<std::vector<double>>& x, bool at_pressure, DensityBranch branch, bool derived) {
            return domain_result(model.states(t, input, x, at_pressure, branch, derived), &model, [&](auto& rows, const StateBatch& r) {
                for (std::size_t i = 0; i < r.rows.size(); ++i) if (r.rows[i]) {
                    const double pressure = at_pressure ? input[input.size() == 1 ? 0 : i] : r.rows[i]->pressure;
                    outside_domains(rows, model.parameters(), model.domains.get(), x[x.size() == 1 ? 0 : i], t[i], pressure, i);
                }
            });
        }, nb::arg("temperature"), nb::arg("input"), nb::arg("composition"), nb::arg("at_pressure"), nb::arg("branch"),
             nb::arg("derived_properties"))
        .def("neutral_property", &epcsaft::HelmholtzModel::neutral_property, nb::arg("temperature"), nb::arg("molar_density"),
             nb::arg("property"), nb::arg("order"))
        .def("density_actions", [](const epcsaft::HelmholtzModel& model,
                                    double temperature, double molar_density,
                                    const std::vector<double>& composition,
                                    const std::vector<epcsaft::PropertyObservable>& observables,
                                    const std::vector<epcsaft::DensityActionDirection>& directions) {
            epcsaft::core::Scope scope("density_actions");
            model.validate_electroneutrality(composition);
            auto result = model.density_actions(temperature, molar_density, composition, observables, directions);
            if (!model.domains || std::none_of(result.results.begin(), result.results.end(), [](const auto& row) { return row.value.has_value(); }))
                return nb::cast(std::move(result));
            return domain_point(std::move(result), model, composition, temperature, std::nullopt, molar_density);
        }, nb::arg("temperature"), nb::arg("molar_density"), nb::arg("mole_fractions"), nb::arg("observables"),
             nb::arg("directions"))
        .def("pressure_actions", [](const epcsaft::HelmholtzModel& model,
                                     double temperature, double pressure,
                                     const std::vector<double>& composition,
                                     epcsaft::DensityBranch branch,
                                     const std::vector<epcsaft::PropertyObservable>& observables,
                                     const std::vector<epcsaft::PressureActionDirection>& directions,
                                     std::optional<double> anchor_density,
                                     const std::optional<std::vector<double>>& active_parameter_values) {
            epcsaft::core::Scope scope("pressure_actions");
            const auto shifted = active_parameter_values && *active_parameter_values != model.active_parameter_values()
                ? epcsaft::equilibrium::model_at(model, *active_parameter_values) : nullptr;
            const auto& current = shifted ? *shifted : model;
            current.validate_electroneutrality(composition);
            return domain_result(current.pressure_actions(temperature, pressure, composition, branch, observables, directions, anchor_density),
                &model, [&](auto& rows, const PropertyActionBatch& r) {
                    if (r.molar_density) outside_domains(rows, current.parameters(), model.domains.get(), composition, temperature, pressure);
                });
        },
             nb::arg("temperature"), nb::arg("pressure"), nb::arg("mole_fractions"),
             nb::arg("branch"), nb::arg("observables"), nb::arg("directions"),
             nb::arg("anchor_density") = std::nullopt, nb::kw_only(), nb::arg("active_parameter_values") = std::nullopt)
        .def("pressure_action_difference", &epcsaft::HelmholtzModel::pressure_action_difference, nb::arg("temperature"),
             nb::arg("pressure"), nb::arg("mole_fractions"), nb::arg("branch"), nb::arg("observables"), nb::arg("direction"),
             nb::arg("step"), nb::arg("tau"), nb::arg("kappa"))
        .def("packing_fraction_coefficient", [](const epcsaft::HelmholtzModel& model,
                                               double temperature, const std::vector<double>& composition) {
            epcsaft::core::Scope scope("packing_fraction_coefficient");
            return model.packing_fraction_coefficient(temperature, composition);
        }, nb::arg("temperature"), nb::arg("composition"))
        .def("state_at_pressure", [](const epcsaft::HelmholtzModel& model, double temperature, double pressure,
                                     const std::vector<double>& composition, epcsaft::DensityBranch branch,
                                     std::optional<double> anchor_density, bool derived_properties) {
            model.validate_electroneutrality(composition);
            return domain_point(model.state_at_pressure(temperature, pressure, composition, branch, anchor_density, derived_properties),
                model, composition, temperature, pressure);
        },
             nb::arg("temperature"), nb::arg("pressure"), nb::arg("mole_fractions"),
             nb::arg("branch"), nb::arg("anchor_density") = std::nullopt, nb::arg("derived_properties") = false)
        .def("molar_density_upper_bound", &epcsaft::HelmholtzModel::molar_density_upper_bound, nb::arg("temperature"),
             nb::arg("composition"))
        .def("density_diagnostics", &epcsaft::HelmholtzModel::density_diagnostics, nb::arg("temperature"),
             nb::arg("pressure"), nb::arg("composition"), nb::arg("request"))
        .def("counters", &epcsaft::HelmholtzModel::counters)
        .def("reset_counters", &epcsaft::HelmholtzModel::reset_counters);
    nb::class_<epcsaft::InfiniteDilutionReference>(module, "InfiniteDilutionReference")
        EPCSAFT_SOLVER_OPTIONS(epcsaft::InfiniteDilutionReference, PROPERTY_SCOPE)
        EPCSAFT_READONLY(epcsaft::InfiniteDilutionReference, log_fugacity, molality, convergence_error,
            parameter_derivative, pure_solvent_log_fugacity, pure_solvent_state, solvent_state, diagnostics,
            solvent_status, solvent_cause, solvent_convergence);
    nb::class_<epcsaft::SolventActivity>(module, "SolventActivity", nb::dynamic_attr())
        EPCSAFT_DOMAIN_VIOLATIONS
        EPCSAFT_READONLY(epcsaft::SolventActivity, solvent_activity, osmotic_coefficient,
            solvent_activity_parameter_derivative, osmotic_parameter_derivative, parameter_status, parameter_cause,
            status, cause, osmotic_cause, state);
    nb::class_<epcsaft::HenryConstant>(module, "HenryConstant", nb::dynamic_attr())
        EPCSAFT_DOMAIN_VIOLATIONS
        EPCSAFT_READONLY(epcsaft::HenryConstant, component_index, raw_log_reference, henry_constant_px, status, message);
    module.def("henry_constants", [](const HelmholtzModel& model, double t, double p, const std::vector<double>& x,
        const std::vector<std::size_t>& solutes) {
        return domain_point(henry_constants(model, t, p, x, solutes), model, x, t, p);
    }, nb::arg("model"), nb::arg("temperature"), nb::arg("pressure"), nb::arg("solvent_fractions"), nb::arg("solute_indices"));
    module.def("solvent_reference", [](const epcsaft::HelmholtzModel& model, double T, double P, const std::vector<double>& x) {
        return epcsaft::solvent_reference(model, T, P, x);
    }, nb::arg("model"), nb::arg("temperature"), nb::arg("pressure"), nb::arg("solvent_fractions"));
    module.def("solvent_activity", [](const epcsaft::HelmholtzModel& model, double T, double P,
        const std::vector<double>& x, const epcsaft::InfiniteDilutionReference& reference) {
        auto result = epcsaft::solvent_activity(model, T, P, x, reference);
        if (!result.solvent_activity) return nb::cast(std::move(result));
        return domain_point(std::move(result), model, x, T, P);
    }, nb::arg("model"), nb::arg("temperature"), nb::arg("pressure"), nb::arg("composition"), nb::arg("reference"));
    nb::class_<epcsaft::ElectrolyteActivity>(module, "ElectrolyteActivity", "Selected salt coefficients on the declared solvent infinite-dilution reference; mean_ionic_molality_scale is the molality-scale gamma+-.", nb::dynamic_attr())
        EPCSAFT_DOMAIN_VIOLATIONS
        EPCSAFT_SOLVER_OPTIONS(epcsaft::ElectrolyteActivity, PROPERTY_SCOPE)
        EPCSAFT_READONLY(epcsaft::ElectrolyteActivity, single_ion, mean_ionic_mole_fraction_scale,
            mean_ionic_log_molality_scale, mean_ionic_molality_scale, osmotic_coefficient, mean_ionic_parameter_derivative,
            reference_solvent_fractions)
        EPCSAFT_REPR(ElectrolyteActivity, "mean_ionic_molality_scale", "osmotic_coefficient", "single_ion");
    module.def("infinite_dilution_reference", [](const epcsaft::HelmholtzModel& model, double temperature, double pressure,
        const std::vector<double>& solvent, const std::array<std::size_t, 2>& ions, const std::array<int, 2>& stoichiometry) {
        return epcsaft::infinite_dilution_reference(model, temperature, pressure, solvent, ions, stoichiometry);
    }, nb::arg("model"), nb::arg("temperature"), nb::arg("pressure"), nb::arg("solvent_fractions"), nb::arg("ion_indices"),
       nb::arg("stoichiometry"));
    module.def("electrolyte_activity", [](const epcsaft::HelmholtzModel& model, double temperature, double pressure,
        const std::vector<double>& solvent, const std::array<std::size_t, 2>& ions, const std::array<int, 2>& stoichiometry,
        double molality, const epcsaft::InfiniteDilutionReference& reference,
        const std::optional<std::vector<double>>& solute_molalities) {
        const auto evaluate = [&](std::vector<double>* composition) {
            return epcsaft::electrolyte_activity(model, temperature, pressure, solvent, ions, stoichiometry, molality, reference,
                std::nullopt, solute_molalities, composition);
        };
        if (!model.domains) return nb::cast(evaluate(nullptr));
        std::vector<double> x;
        auto result = evaluate(&x);
        return domain_point(std::move(result), model, x, temperature, pressure);
    }, nb::arg("model"), nb::arg("temperature"), nb::arg("pressure"), nb::arg("solvent_fractions"), nb::arg("ion_indices"),
       nb::arg("stoichiometry"), nb::arg("molality"), nb::arg("reference"), nb::arg("solute_molalities") = nb::none());
}
