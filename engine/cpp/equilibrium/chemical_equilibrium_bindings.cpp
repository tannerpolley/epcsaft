#include "chemical_equilibrium.hpp"
#include "grepe.hpp"
#include "grepe_pybind.hpp"
#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace epcsaft_engine::equilibrium {
namespace {

DenseMatrix dense_matrix(const py::handle& value, const char* field) {
    const std::vector<std::vector<double>> rows = py::cast<std::vector<std::vector<double>>>(
        value
    );
    if (rows.empty()) {
        return {};
    }
    const std::size_t columns = rows.front().size();
    if (columns == 0) {
        throw py::value_error(std::string(field) + " rows must not be empty");
    }
    DenseMatrix result{rows.size(), columns, {}};
    result.values.reserve(rows.size() * columns);
    for (const std::vector<double>& row : rows) {
        if (row.size() != columns) {
            throw py::value_error(std::string(field) + " must be rectangular");
        }
        result.values.insert(result.values.end(), row.begin(), row.end());
    }
    return result;
}

DenseMatrix phase_direction_matrix(
    const std::vector<std::vector<double>>& rows,
    std::size_t expected_columns,
    const char* field
) {
    if (rows.empty()) {
        return {0, expected_columns, {}};
    }
    const std::size_t columns = rows.front().size();
    DenseMatrix result{rows.size(), columns, {}};
    result.values.reserve(rows.size() * columns);
    for (const std::vector<double>& row : rows) {
        if (row.size() != columns) {
            throw py::value_error(std::string(field) + " must be rectangular");
        }
        result.values.insert(result.values.end(), row.begin(), row.end());
    }
    if (columns != expected_columns) {
        throw py::value_error(
            std::string(field) + " columns must match the phase-species dimension"
        );
    }
    return result;
}

std::vector<std::vector<double>> matrix_rows(const DenseMatrix& matrix) {
    std::vector<std::vector<double>> result(
        matrix.rows,
        std::vector<double>(matrix.columns)
    );
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t column = 0; column < matrix.columns; ++column) {
            result[row][column] = matrix.values[row * matrix.columns + column];
        }
    }
    return result;
}

py::dict compile_phase_direction_system_binding(
    const std::vector<std::string>& phase_ids,
    const std::vector<std::string>& species_ids,
    const std::vector<double>& charges,
    const std::vector<std::vector<double>>& balance_rows,
    const std::vector<std::vector<double>>& phase_reference_rows,
    const std::vector<std::vector<double>>& direction_rows,
    const std::vector<double>& driving_constants,
    double rank_tolerance,
    double conservation_tolerance
) {
    const std::size_t dimension = phase_ids.size() * species_ids.size();
    const PhaseDirectionCompilation compilation =
        compile_phase_direction_system(
            phase_ids,
            species_ids,
            charges,
            phase_direction_matrix(
                balance_rows,
                species_ids.size(),
                "balance_matrix"
            ),
            phase_direction_matrix(
                phase_reference_rows,
                species_ids.size(),
                "phase_reference_amounts"
            ),
            phase_direction_matrix(
                direction_rows,
                dimension,
                "direction_rows"
            ),
            driving_constants,
            rank_tolerance,
            conservation_tolerance
        );
    py::dict result;
    result["phase_ids"] = compilation.phase_ids;
    result["species_ids"] = compilation.species_ids;
    result["charges"] = compilation.charges;
    result["coordinate_count"] = dimension;
    py::list coordinate_order;
    for (const std::string& phase_id : compilation.phase_ids) {
        for (const std::string& species_id : compilation.species_ids) {
            coordinate_order.append(py::make_tuple(phase_id, species_id));
        }
    }
    result["coordinate_order"] = std::move(coordinate_order);
    result["direction_matrix"] = matrix_rows(compilation.direction_matrix);
    result["directions"] = matrix_rows(compilation.direction_matrix);
    result["invariant_rows"] = matrix_rows(compilation.invariant_matrix);
    result["invariant_matrix"] = matrix_rows(compilation.invariant_matrix);
    result["phase_reference_amounts"] = compilation.phase_reference_amounts;
    result["driving_constants"] = compilation.driving_constants;
    result["invariant_totals"] = compilation.invariant_totals;
    result["reference_potentials"] = compilation.reference_potentials;
    result["g"] = compilation.reference_potentials;
    result["global_balance_residuals"] = compilation.global_balance_residuals;
    result["phase_charge_residuals"] = compilation.phase_charge_residuals;
    result["reference_phase_charge_residuals"] =
        compilation.reference_phase_charge_residuals;
    result["direction_rank"] = compilation.direction_rank;
    result["rank"] = compilation.direction_rank;
    result["invariant_rank"] = compilation.invariant_rank;
    result["direction_count"] = compilation.direction_matrix.rows;
    result["invariant_count"] = compilation.invariant_matrix.rows;
    result["direction_scale"] = compilation.direction_scale;
    result["rank_tolerance"] = compilation.rank_tolerance;
    result["conservation_tolerance"] = compilation.conservation_tolerance;
    result["invariant_orthonormal_residual_inf_norm"] =
        compilation.invariant_orthonormal_residual_inf;
    result["invariant_direction_residual_inf_norm"] =
        compilation.invariant_direction_residual_inf;
    result["potential_equation_residual_inf_norm"] =
        compilation.potential_equation_residual_inf;
    result["minimum_norm_potential_norm"] =
        compilation.minimum_norm_potential_norm;
    result["global_balance_residual_inf_norm"] =
        compilation.global_balance_residual_inf;
    result["phase_charge_residual_inf_norm"] =
        compilation.phase_charge_residual_inf;
    result["reference_neutrality_residual_inf_norm"] =
        compilation.reference_neutrality_residual_inf;
    return result;
}

double quantity_magnitude(const py::handle& value, const char* unit) {
    return py::cast<double>(value.attr("to")(unit).attr("magnitude"));
}

GrepeIntensive grepe_intensive(
    const py::handle& value,
    const char* unit
) {
    GrepeIntensive result;
    result.solved = py::hasattr(value, "initial");
    if (!result.solved) {
        result.value = quantity_magnitude(value.attr("value"), unit);
        return result;
    }
    result.value = quantity_magnitude(value.attr("initial"), unit);
    const py::tuple bounds = py::cast<py::tuple>(value.attr("bounds"));
    result.bounds = {
        quantity_magnitude(bounds[0], unit),
        quantity_magnitude(bounds[1], unit),
    };
    for (const py::handle start : py::cast<py::tuple>(value.attr("starts"))) {
        result.starts.push_back(quantity_magnitude(start, unit));
    }
    return result;
}

epcsaft::NativeModel grepe_active_model(
    const epcsaft::NativeModel& base,
    const py::object& active_parameters
) {
    if (active_parameters.is_none()) {
        return base;
    }
    std::vector<epcsaft::NativeActiveCoordinate> coordinates;
    for (const py::handle coordinate : py::cast<py::tuple>(
             active_parameters.attr("coordinates")
         )) {
        const py::object applicability = coordinate.attr("applicability");
        if (applicability.is_none()) {
            throw py::value_error(
                "active GREPE parameters require resolved physical applicability"
            );
        }
        coordinates.push_back({
            py::cast<std::string>(coordinate.attr("identity")),
            py::cast<std::string>(coordinate.attr("family")),
            py::cast<std::string>(applicability.attr("application")),
            py::cast<std::vector<std::size_t>>(
                applicability.attr("component_indices")
            ),
            py::cast<std::vector<std::array<std::size_t, 2>>>(
                applicability.attr("association_edges")
            ),
            py::cast<std::vector<double>>(
                applicability.attr("association_energy_scales")
            ),
            py::cast<double>(coordinate.attr("value")),
            py::cast<double>(applicability.attr("reference_temperature_k")),
        });
    }
    return base.with_active_coordinates(coordinates);
}

std::vector<std::size_t> grepe_support_indices(
    const py::handle& support,
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges,
    std::string& support_name
) {
    const std::string kind = py::cast<std::string>(
        py::type::of(support).attr("__name__")
    );
    std::vector<std::size_t> result;
    if (kind == "AllComponents") {
        support_name = "all_components";
        result.resize(component_ids.size());
        std::iota(result.begin(), result.end(), 0);
    } else if (kind == "NeutralComponents") {
        support_name = "neutral_components";
        for (std::size_t index = 0; index < charges.size(); ++index) {
            if (charges[index] == 0) {
                result.push_back(index);
            }
        }
    } else {
        support_name = kind == "DeclaredComponents"
            ? "declared_components"
            : "declared_neutral_components";
        for (const std::string& identity : py::cast<std::vector<std::string>>(
                 support.attr("component_ids")
             )) {
            const auto found = std::find(
                component_ids.begin(), component_ids.end(), identity
            );
            if (found == component_ids.end()) {
                throw py::value_error("GREPE phase support names an unknown component");
            }
            result.push_back(static_cast<std::size_t>(found - component_ids.begin()));
        }
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
    }
    return result;
}

ReactionSystemInput grepe_reaction_input(
    const py::handle& reaction,
    double temperature_k,
    double pressure_pa
) {
    ReactionSystemInput input;
    input.species_ids = py::cast<std::vector<std::string>>(
        reaction.attr("species_ids")
    );
    input.charges = py::cast<std::vector<int>>(reaction.attr("charges"));
    input.molar_masses_kg_per_mol = py::cast<std::vector<double>>(
        reaction.attr("molar_masses_kg_per_mol")
    );
    input.balance_matrix = dense_matrix(reaction.attr("balance_matrix"), "balance matrix");
    input.conserved_totals = py::cast<std::vector<double>>(
        reaction.attr("conserved_totals")
    );
    input.reaction_matrix = dense_matrix(
        reaction.attr("reaction_matrix"), "reaction matrix"
    );
    if (input.reaction_matrix.rows == 0) {
        input.reaction_matrix.columns = input.species_ids.size();
    }
    input.feed_amounts = py::cast<std::vector<double>>(
        reaction.attr("feed_amounts_mol")
    );
    const py::object standard = reaction.attr("source_standard_state");
    const double record_pressure = standard.is_none()
        ? pressure_pa
        : py::cast<double>(standard.attr("reference_pressure_pa"));
    for (const py::handle constant : py::cast<py::tuple>(
             reaction.attr("equilibrium_constants")
         )) {
        const py::object reaction_id = constant.attr("reaction_id");
        const py::object correlation_kind = constant.attr("correlation_kind");
        input.ln_k.push_back(py::cast<double>(constant.attr("ln_value")));
        input.equilibrium_constant_records.push_back({
            py::cast<std::string>(constant.attr("source_id")),
            py::cast<std::string>(constant.attr("reference_id")),
            py::cast<std::string>(constant.attr("reaction_orientation")),
            py::cast<std::string>(constant.attr("conversion_id")),
            py::cast<bool>(constant.attr("dimensionless")),
            reaction_id.is_none() ? "" : py::cast<std::string>(reaction_id),
            correlation_kind.is_none()
                ? "" : py::cast<std::string>(correlation_kind),
            py::cast<std::vector<std::string>>(
                constant.attr("coefficient_identities")
            ),
            py::cast<std::vector<double>>(
                constant.attr("coefficient_values")
            ),
            temperature_k,
            record_pressure,
        });
    }
    input.temperature_k = temperature_k;
    input.pressure_pa = pressure_pa;
    return input;
}

std::optional<GrepeSourceStandardState> grepe_source(
    const py::handle& reaction,
    const std::vector<std::string>& component_ids,
    bool active_parameters_requested,
    const std::string& reference_fingerprint
) {
    const py::object standard = reaction.attr("source_standard_state");
    if (standard.is_none()) {
        return std::nullopt;
    }
    GrepeSourceStandardState result;
    result.id = py::cast<std::string>(standard.attr("id"));
    result.log_activity_scale_factors = py::cast<std::vector<double>>(
        standard.attr("log_activity_scale_factors")
    );
    result.reference_pressure_pa = py::cast<double>(
        standard.attr("reference_pressure_pa")
    );
    SourceReferenceRequest& request = result.reference;
    request.component_ids = py::cast<std::vector<std::string>>(
        standard.attr("source_reference_component_ids")
    );
    request.source_solvent_mole_fractions = py::cast<std::vector<double>>(
        standard.attr("source_reference_solvent_composition")
    );
    for (const py::handle pair : py::cast<py::tuple>(
             standard.attr("source_reference_ion_pairs")
         )) {
        const std::vector<std::string> identities = py::cast<std::vector<std::string>>(
            pair.attr("component_ids")
        );
        request.ion_pair_indices.push_back({
            static_cast<std::size_t>(std::find(
                component_ids.begin(), component_ids.end(), identities[0]
            ) - component_ids.begin()),
            static_cast<std::size_t>(std::find(
                component_ids.begin(), component_ids.end(), identities[1]
            ) - component_ids.begin()),
        });
        request.ion_pair_stoichiometry.push_back(
            py::cast<std::array<int, 2>>(pair.attr("stoichiometry"))
        );
    }
    for (const py::handle pure : py::cast<py::tuple>(
             standard.attr("source_reference_pure_components")
         )) {
        const std::string identity = py::cast<std::string>(pure.attr("component_id"));
        request.pure_component_indices.push_back(static_cast<std::size_t>(std::find(
            component_ids.begin(), component_ids.end(), identity
        ) - component_ids.begin()));
        request.pure_component_phases.push_back(
            py::cast<std::string>(pure.attr("phase"))
        );
    }
    request.phase = py::cast<std::string>(standard.attr("source_reference_phase"));
    request.reference_convention = py::cast<std::string>(
        standard.attr("source_reference_convention")
    );
    request.reference_state_id = result.id;
    request.activity_convention_id = py::cast<std::string>(
        standard.attr("source_reference_activity_convention_id")
    );
    request.reference_state_fingerprint = reference_fingerprint;
    request.standard_molality_mol_per_kg = py::cast<double>(
        standard.attr("source_reference_standard_molality_mol_per_kg")
    );
    request.source_reference_pressure_pa = result.reference_pressure_pa;
    request.required_derivatives = active_parameters_requested ? 10u : 2u;
    return result;
}

std::optional<GrepeDirectionSystem> grepe_directions(
    const py::object& phase_reactions,
    const ReactionSystemInput& reaction,
    const std::vector<GrepePhase>& phases
) {
    if (phase_reactions.is_none()) {
        return std::nullopt;
    }
    GrepeDirectionSystem result;
    result.phase_ids = py::cast<std::vector<std::string>>(
        phase_reactions.attr("phase_ids")
    );
    result.species_ids = py::cast<std::vector<std::string>>(
        phase_reactions.attr("species_ids")
    );
    result.charges.reserve(reaction.charges.size());
    for (int charge : reaction.charges) {
        result.charges.push_back(static_cast<double>(charge));
    }
    result.balance_matrix = reaction.balance_matrix;
    const std::size_t columns = phases.size() * reaction.species_ids.size();
    const py::tuple rows = py::cast<py::tuple>(phase_reactions.attr("reactions"));
    result.direction_matrix = {
        rows.size(), columns, std::vector<double>(rows.size() * columns, 0.0)
    };
    std::unordered_map<std::string, std::size_t> phase_index;
    std::unordered_map<std::string, std::size_t> species_index;
    for (std::size_t index = 0; index < phases.size(); ++index) {
        phase_index.emplace(phases[index].identity, index);
    }
    for (std::size_t index = 0; index < reaction.species_ids.size(); ++index) {
        species_index.emplace(reaction.species_ids[index], index);
    }
    for (std::size_t row = 0; row < rows.size(); ++row) {
        const py::handle declared = rows[row];
        for (const py::handle term : py::cast<py::tuple>(declared.attr("terms"))) {
            result.direction_matrix(
                row,
                phase_index.at(py::cast<std::string>(term.attr("phase_id")))
                    * reaction.species_ids.size()
                    + species_index.at(py::cast<std::string>(term.attr("species_id")))
            ) += py::cast<double>(term.attr("coefficient"));
        }
        const py::object ln_k = declared.attr("ln_equilibrium_constant");
        result.driving_constants.push_back(
            ln_k.is_none() ? 0.0 : py::cast<double>(ln_k)
        );
    }
    return result;
}

GrepeRequest grepe_request(
    const py::handle& problem,
    bool active_parameters_requested,
    const std::string& source_reference_fingerprint,
    const GrepeIntensive* temperature = nullptr,
    const GrepeIntensive* pressure = nullptr
) {
    GrepeRequest request;
    request.identity = py::cast<std::string>(problem.attr("identity"));
    request.temperature = temperature == nullptr
        ? grepe_intensive(problem.attr("temperature"), "kelvin")
        : *temperature;
    request.pressure = pressure == nullptr
        ? grepe_intensive(problem.attr("pressure"), "pascal")
        : *pressure;
    const py::handle reaction = problem.attr("reaction_system");
    request.reaction = grepe_reaction_input(
        reaction,
        request.temperature.value,
        request.pressure.value
    );
    request.trace_floor = py::cast<double>(
        reaction.attr("strict_interior_amount_floor_mol")
    );
    request.reaction_phase_ids = py::cast<std::vector<std::string>>(
        problem.attr("reaction_phase_ids")
    );
    std::unordered_map<std::string, std::size_t> phase_index;
    for (const py::handle phase : py::cast<py::tuple>(problem.attr("phases"))) {
        GrepePhase declared;
        declared.identity = py::cast<std::string>(phase.attr("identity"));
        phase_index.emplace(declared.identity, request.phases.size());
        declared.role = py::cast<std::string>(phase.attr("fluid_role")) == "liquid"
            ? PhaseRole::Liquid
            : PhaseRole::Vapor;
        declared.finite = py::cast<std::string>(phase.attr("amount_role")) == "finite";
        const py::handle phase_model = phase.attr("model");
        declared.ideal_gas = !py::hasattr(
            phase_model, "admissible_packing_fraction_interval"
        );
        declared.model_reference_id = py::cast<std::string>(
            phase_model.attr("reference_id")
        );
        if (!declared.ideal_gas) {
            declared.packing_fraction_interval = py::cast<std::array<double, 2>>(
                phase_model.attr("admissible_packing_fraction_interval")
            );
        }
        declared.support_indices = grepe_support_indices(
            phase.attr("support"),
            request.reaction.species_ids,
            request.reaction.charges,
            declared.support_name
        );
        const py::object start = phase.attr("start");
        if (!start.is_none()) {
            declared.start = ReactiveFinitePhaseStart{
                py::cast<std::vector<double>>(start.attr("amounts_mol")),
                py::cast<double>(start.attr("molar_volume_m3_per_mol")),
            };
        }
        request.phases.push_back(std::move(declared));
    }
    for (const py::handle output : py::cast<py::tuple>(problem.attr("outputs"))) {
        const py::object phase_identity = output.attr("phase_identity");
        std::vector<double> coefficients = py::cast<std::vector<double>>(
            output.attr("coefficients")
        );
        if (coefficients.empty()) {
            coefficients.assign(request.reaction.species_ids.size(), 0.0);
        }
        request.outputs.push_back({
            py::cast<std::string>(output.attr("selector")),
            phase_identity.is_none()
                ? 0u
                : phase_index.at(py::cast<std::string>(phase_identity)),
            std::move(coefficients),
            py::cast<std::vector<double>>(
                output.attr("solvent_mass_coefficients_kg_per_mol")
            ),
        });
    }
    request.directions = grepe_directions(
        problem.attr("phase_reactions"), request.reaction, request.phases
    );
    request.source = grepe_source(
        reaction,
        request.reaction.species_ids,
        active_parameters_requested,
        source_reference_fingerprint
    );
    const py::object feed = problem.attr("feed");
    if (!feed.is_none()) {
        GrepeFeed parsed;
        parsed.affine = py::hasattr(feed, "base_amounts_mol");
        parsed.amounts_or_base_mol = py::cast<std::vector<double>>(
            feed.attr(parsed.affine ? "base_amounts_mol" : "amounts_mol")
        );
        if (parsed.affine) {
            parsed.component_coordinate_matrix = dense_matrix(
                feed.attr("basis_matrix_mol_per_coordinate"), "affine feed basis"
            );
            for (const py::handle coordinate : py::cast<py::tuple>(
                     feed.attr("coordinates")
                 )) {
                const std::string identity = py::cast<std::string>(
                    coordinate.attr("identity")
                );
                parsed.coordinate_ids.push_back(identity);
                parsed.coordinates.push_back({
                    identity,
                    py::cast<std::array<double, 2>>(coordinate.attr("bounds")),
                });
            }
        } else {
            parsed.coordinate_ids = py::cast<std::vector<std::string>>(
                feed.attr("coordinate_ids")
            );
            if (!parsed.coordinate_ids.empty()) {
                parsed.component_coordinate_matrix = dense_matrix(
                    feed.attr("amount_jacobian_mol_per_coordinate"),
                    "resolved feed derivative"
                );
            }
        }
        request.feed = std::move(parsed);
    }
    const py::tuple boundaries = py::cast<py::tuple>(problem.attr("intensive_boundaries"));
    if (!boundaries.empty()) {
        const py::handle boundary = boundaries[0];
        const py::object identity = boundary.attr("phase_identity");
        request.boundary = GrepeIntensiveBoundary{
            py::cast<std::string>(boundary.attr("identity")),
            py::cast<std::string>(boundary.attr("coordinate_identity")),
            identity.is_none() ? "" : py::cast<std::string>(identity),
            py::cast<std::string>(boundary.attr("selector")),
            py::cast<std::vector<double>>(boundary.attr("coefficients")),
            py::cast<double>(boundary.attr("target")),
            py::cast<double>(boundary.attr("initial")),
        };
    }
    const py::object continuation = problem.attr("continuation_state");
    if (!continuation.is_none()) {
        const std::vector<std::string> identities = py::cast<std::vector<std::string>>(
            continuation.attr("phase_identities")
        );
        const py::tuple phases = py::cast<py::tuple>(continuation.attr("phases"));
        for (std::size_t index = 0; index < identities.size(); ++index) {
            request.continuation_phases.push_back({
                identities[index],
                py::cast<std::vector<double>>(phases[index].attr("mole_fractions")),
                py::cast<double>(phases[index].attr("molar_volume_m3_per_mol")),
            });
        }
        request.continuation_pressure_pa = py::cast<double>(
            continuation.attr("pressure_pa")
        );
    }
    return request;
}

GrepeResult solve_grepe_binding(
    const epcsaft::NativeModel& base_model,
    const epcsaft::NativeModel& incipient_model,
    const py::handle& problem,
    const py::object& active_parameters,
    const std::string& source_reference_fingerprint,
    bool perform_derivative_audits,
    bool declared_starts_only
) {
    epcsaft::NativeModel model = grepe_active_model(base_model, active_parameters);
    GrepeRequest request = grepe_request(
        problem,
        !active_parameters.is_none(),
        source_reference_fingerprint
    );
    request.perform_derivative_audits = perform_derivative_audits;
    request.declared_starts_only = declared_starts_only;
    py::gil_scoped_release release;
    return solve_grepe(model, incipient_model, request);
}

GrepeResult initialize_grepe_local_binding(
    const epcsaft::NativeModel& base_model,
    const py::handle& problem,
    const py::object& active_parameters,
    const std::string& source_reference_fingerprint
) {
    epcsaft::NativeModel model = grepe_active_model(base_model, active_parameters);
    GrepeRequest request = grepe_request(
        problem,
        !active_parameters.is_none(),
        source_reference_fingerprint
    );
    py::gil_scoped_release release;
    return initialize_grepe_local(model, request);
}

GrepeResult evaluate_grepe_local_binding(
    const epcsaft::NativeModel& base_model,
    const py::handle& problem,
    const py::object& active_parameters,
    const std::string& source_reference_fingerprint,
    const std::vector<double>& local_variables,
    const std::vector<double>& local_lower_bounds,
    const std::vector<double>& local_upper_bounds
) {
    epcsaft::NativeModel model = grepe_active_model(base_model, active_parameters);
    GrepeRequest request = grepe_request(
        problem,
        !active_parameters.is_none(),
        source_reference_fingerprint
    );
    py::gil_scoped_release release;
    return evaluate_grepe_local(
        model,
        request,
        local_variables,
        local_lower_bounds,
        local_upper_bounds
    );
}

}  // namespace

GrepeRequest parse_grepe_request(
    const py::handle& problem,
    bool active_parameters_requested,
    const std::string& source_reference_fingerprint,
    const GrepeIntensive* temperature,
    const GrepeIntensive* pressure
) {
    return grepe_request(
        problem,
        active_parameters_requested,
        source_reference_fingerprint,
        temperature,
        pressure
    );
}

void bind_chemical_equilibrium(py::module_& module) {
    py::class_<ReactiveFinitePhaseResult>(module, "_ReactiveFinitePhaseResult")
        .def_readonly("accepted", &ReactiveFinitePhaseResult::accepted)
        .def_readonly("solver_status", &ReactiveFinitePhaseResult::solver_status)
        .def_readonly("callback_error", &ReactiveFinitePhaseResult::callback_error)
        .def_readonly(
            "callback_failure_kind", &ReactiveFinitePhaseResult::callback_failure_kind
        )
        .def_readonly("failure_reason", &ReactiveFinitePhaseResult::failure_reason)
        .def_readonly("numerical_status", &ReactiveFinitePhaseResult::numerical_status)
        .def_readonly("physical_status", &ReactiveFinitePhaseResult::physical_status)
        .def_readonly(
            "provider_domain_status", &ReactiveFinitePhaseResult::provider_domain_status
        )
        .def_readonly("resolved_policy", &ReactiveFinitePhaseResult::resolved_policy)
        .def_readonly("phase_amounts_mol", &ReactiveFinitePhaseResult::phase_amounts_mol)
        .def_readonly("phase_volumes_m3", &ReactiveFinitePhaseResult::phase_volumes_m3)
        .def_readonly("phase_pressures_pa", &ReactiveFinitePhaseResult::phase_pressures_pa)
        .def_readonly(
            "phase_packing_fractions", &ReactiveFinitePhaseResult::phase_packing_fractions
        )
        .def_readonly(
            "phase_chemical_potentials_over_rt",
            &ReactiveFinitePhaseResult::phase_chemical_potentials_over_rt
        )
        .def_readonly(
            "phase_log_fugacities", &ReactiveFinitePhaseResult::phase_log_fugacities
        )
        .def_readonly("local_variables", &ReactiveFinitePhaseResult::local_variables)
        .def_readonly("local_lower_bounds", &ReactiveFinitePhaseResult::local_lower_bounds)
        .def_readonly("local_upper_bounds", &ReactiveFinitePhaseResult::local_upper_bounds)
        .def_readonly("input_order", &ReactiveFinitePhaseResult::input_order)
        .def_readonly(
            "kkt_input_jacobian", &ReactiveFinitePhaseResult::kkt_input_jacobian
        )
        .def_readonly(
            "input_amount_derivatives", &ReactiveFinitePhaseResult::input_amount_derivatives
        )
        .def_readonly(
            "input_volume_derivatives", &ReactiveFinitePhaseResult::input_volume_derivatives
        )
        .def_readonly(
            "input_pressure_derivatives", &ReactiveFinitePhaseResult::input_pressure_derivatives
        )
        .def_readonly(
            "input_chemical_potential_derivatives",
            &ReactiveFinitePhaseResult::input_chemical_potential_derivatives
        )
        .def_readonly(
            "input_log_fugacity_derivatives",
            &ReactiveFinitePhaseResult::input_log_fugacity_derivatives
        )
        .def_readonly("kkt_root_values", &ReactiveFinitePhaseResult::kkt_root_values)
        .def_readonly("kkt_root_jacobian", &ReactiveFinitePhaseResult::kkt_root_jacobian)
        .def_readonly("kkt_dimension", &ReactiveFinitePhaseResult::kkt_dimension)
        .def_readonly("kkt_rank", &ReactiveFinitePhaseResult::kkt_rank)
        .def_readonly("condition_number_inf", &ReactiveFinitePhaseResult::condition_number_inf)
        .def_readonly("output_values", &ReactiveFinitePhaseResult::output_values)
        .def_readonly(
            "output_input_derivatives", &ReactiveFinitePhaseResult::output_input_derivatives
        )
        .def_readonly(
            "output_local_derivatives", &ReactiveFinitePhaseResult::output_local_derivatives
        )
        .def_readonly(
            "output_parameter_partials", &ReactiveFinitePhaseResult::output_parameter_partials
        )
        .def_readonly("objective_over_rt", &ReactiveFinitePhaseResult::objective_over_rt)
        .def_readonly("balance_inf_norm", &ReactiveFinitePhaseResult::balance_inf_norm)
        .def_readonly(
            "pressure_relative_inf_norm", &ReactiveFinitePhaseResult::pressure_relative_inf_norm
        )
        .def_readonly(
            "reaction_affinity_inf_norm", &ReactiveFinitePhaseResult::reaction_affinity_inf_norm
        )
        .def_readonly(
            "transfer_chemical_potential_inf_norm",
            &ReactiveFinitePhaseResult::transfer_chemical_potential_inf_norm
        )
        .def_readonly(
            "kkt_stationarity_inf_norm", &ReactiveFinitePhaseResult::kkt_stationarity_inf_norm
        )
        .def_readonly("phase_distinctness", &ReactiveFinitePhaseResult::phase_distinctness)
        .def_readonly("minimum_amount_mol", &ReactiveFinitePhaseResult::minimum_amount_mol)
        .def_readonly("minimum_bound_margin", &ReactiveFinitePhaseResult::minimum_bound_margin)
        .def_readonly(
            "reduced_hessian_min_eigenvalue",
            &ReactiveFinitePhaseResult::reduced_hessian_min_eigenvalue
        )
        .def_readonly("reduced_hessian_scale", &ReactiveFinitePhaseResult::reduced_hessian_scale)
        .def_readonly(
            "reduced_hessian_threshold", &ReactiveFinitePhaseResult::reduced_hessian_threshold
        )
        .def_readonly(
            "declared_direction_rank", &ReactiveFinitePhaseResult::declared_direction_rank
        )
        .def_readonly("admissible_nullity", &ReactiveFinitePhaseResult::admissible_nullity)
        .def_readonly("optimizer_iterations", &ReactiveFinitePhaseResult::optimizer_iterations)
        .def_readonly(
            "provider_evaluation_count", &ReactiveFinitePhaseResult::provider_evaluation_count
        )
        .def_readonly(
            "derivative_audits_performed",
            &ReactiveFinitePhaseResult::derivative_audits_performed
        )
        .def_readonly(
            "derivative_audit_status", &ReactiveFinitePhaseResult::derivative_audit_status
        )
        .def_readonly(
            "objective_gradient_audit_max_relative_error",
            &ReactiveFinitePhaseResult::objective_gradient_audit_max_relative_error
        )
        .def_readonly(
            "constraint_jacobian_audit_max_relative_error",
            &ReactiveFinitePhaseResult::constraint_jacobian_audit_max_relative_error
        )
        .def_readonly(
            "lagrangian_hessian_audit_max_relative_error",
            &ReactiveFinitePhaseResult::lagrangian_hessian_audit_max_relative_error
        )
        .def_readonly("parameter_order", &ReactiveFinitePhaseResult::parameter_order)
        .def_readonly("sensitivity_status", &ReactiveFinitePhaseResult::sensitivity_status)
        .def_readonly(
            "sensitivity_failure_reason", &ReactiveFinitePhaseResult::sensitivity_failure_reason
        )
        .def_readonly("amount_derivatives", &ReactiveFinitePhaseResult::amount_derivatives)
        .def_readonly("volume_derivatives", &ReactiveFinitePhaseResult::volume_derivatives)
        .def_readonly("pressure_derivatives", &ReactiveFinitePhaseResult::pressure_derivatives)
        .def_readonly(
            "chemical_potential_derivatives",
            &ReactiveFinitePhaseResult::chemical_potential_derivatives
        )
        .def_readonly(
            "log_fugacity_derivatives", &ReactiveFinitePhaseResult::log_fugacity_derivatives
        );
    py::class_<GrepeResult>(module, "_GrepeResult")
        .def_readonly("formulation", &GrepeResult::formulation)
        .def_readonly("certificate_scope", &GrepeResult::certificate_scope)
        .def_readonly("equilibrium", &GrepeResult::equilibrium)
        .def_readonly("failure_code", &GrepeResult::failure_code)
        .def_readonly("temperature_k", &GrepeResult::temperature_k)
        .def_readonly("pressure_pa", &GrepeResult::pressure_pa)
        .def_readonly(
            "source_reference_representation_residual_inf_norm",
            &GrepeResult::source_reference_representation_residual_inf_norm
        )
        .def_readonly("evidence", &GrepeResult::evidence)
        .def_readonly("branch_kind", &GrepeResult::branch_kind)
        .def_readonly("branch_coordinates", &GrepeResult::branch_coordinates)
        .def_readonly(
            "phase_mechanical_classes", &GrepeResult::phase_mechanical_classes
        );
    module.def(
        "_solve_grepe",
        &solve_grepe_binding,
        py::arg("model"),
        py::arg("incipient_model"),
        py::arg("problem"),
        py::arg("active_parameters") = py::none(),
        py::arg("source_reference_fingerprint") = "",
        py::arg("perform_derivative_audits") = false,
        py::arg("declared_starts_only") = false
    );
    module.def(
        "_initialize_grepe_local",
        &initialize_grepe_local_binding,
        py::arg("model"),
        py::arg("problem"),
        py::arg("active_parameters"),
        py::arg("source_reference_fingerprint")
    );
    module.def(
        "_evaluate_grepe_local",
        &evaluate_grepe_local_binding,
        py::arg("model"),
        py::arg("problem"),
        py::arg("active_parameters"),
        py::arg("source_reference_fingerprint"),
        py::arg("local_variables"),
        py::arg("local_lower_bounds"),
        py::arg("local_upper_bounds")
    );
    module.def(
        "_compile_phase_direction_system",
        &compile_phase_direction_system_binding,
        py::arg("phase_ids"),
        py::arg("species_ids"),
        py::arg("charges"),
        py::arg("balance_matrix"),
        py::arg("phase_reference_amounts"),
        py::arg("direction_rows"),
        py::arg("driving_constants"),
        py::arg("rank_tolerance") = kPhaseDirectionRankTolerance,
        py::arg("conservation_tolerance") = kPhaseDirectionConservationTolerance
    );
    // Keep the private seam discoverable under the concise plural spelling as
    // well; both names call the same compiler and have no separate behavior.
    module.def(
        "_compile_phase_directions",
        &compile_phase_direction_system_binding,
        py::arg("phase_ids"),
        py::arg("species_ids"),
        py::arg("charges"),
        py::arg("balance_matrix"),
        py::arg("phase_reference_amounts"),
        py::arg("directions"),
        py::arg("driving_constants"),
        py::arg("rank_tolerance") = kPhaseDirectionRankTolerance,
        py::arg("conservation_tolerance") = kPhaseDirectionConservationTolerance
    );
}

}  // namespace epcsaft_engine::equilibrium
