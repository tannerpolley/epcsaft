#include "eos_regression_evaluator.hpp"

#include "../../../cpp/sha256.hpp"

#include "physical_constants.hpp"
#include "pure_saturation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace epcsaft {
namespace {

using epcsaft::detail::GAS_CONSTANT;
constexpr std::size_t kMaximumFailureCodeLength = 128;
constexpr std::size_t kMaximumFailureReasonLength = 1024;

using epcsaft_engine::interfaces::Evaluation;
using epcsaft_engine::interfaces::EvaluationFailure;
using epcsaft_engine::interfaces::EvaluationResult;

EvaluationFailure failure(std::string code, std::string reason) {
    if (code.empty()) {
        code = "native_residual_failure";
    }
    if (reason.empty()) {
        reason = "EOS owner evaluation failed";
    }
    code.resize(std::min(code.size(), kMaximumFailureCodeLength));
    reason.resize(std::min(reason.size(), kMaximumFailureReasonLength));
    return {std::move(code), std::move(reason)};
}

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) { return std::isfinite(value); });
}

void append_field(std::string& payload, const std::string& value) {
    payload += std::to_string(value.size());
    payload.push_back(':');
    payload += value;
    payload.push_back('|');
}

bool ends_with(const std::string& value, const char* suffix) {
    const std::size_t length = std::char_traits<char>::length(suffix);
    return value.size() >= length
        && value.compare(value.size() - length, length, suffix) == 0;
}

std::vector<std::string> canonical_parameter_units(
    const std::vector<NativeActiveCoordinate>& coordinates
) {
    std::vector<std::string> result;
    result.reserve(coordinates.size());
    for (const auto& coordinate : coordinates) {
        const auto& family = coordinate.family;
        if (family == "segment_diameter" || family == "born_diameter") {
            result.push_back("angstrom");
        } else if (family == "dispersion_energy_over_k"
                   || family == "association_energy_over_k"
                   || family == "k_ij_reciprocal_temperature_slope") {
            result.push_back("kelvin");
        } else if (family == "exponential_temperature_coefficient") {
            result.push_back("1 / kelvin");
        } else if (family.rfind("temperature_polynomial_coefficient_", 0) == 0) {
            result.push_back("1 / kelvin ** " + family.substr(35));
        } else if (family == "reaction_correlation") {
            const auto& identity = coordinate.identity;
            result.push_back(ends_with(identity, "_per_k") ? "1 / kelvin"
                : ends_with(identity, "_k") ? "kelvin"
                : "dimensionless");
        } else if (family == "segment_count" || family == "relative_permittivity"
                   || family == "solvation_factor" || family == "association_volume"
                   || family == "k_ij" || family == "l_ij" || family == "k_hb_ij"
                   || family == "ion_fraction_suppression_coefficient"
                   || family == "ionic_region_relative_permittivity"
                   || family == "ion_specific_suppression_coefficient"
                   || family == "correlation_integral_parameter") {
            result.push_back("dimensionless");
        } else {
            throw std::invalid_argument("unsupported EOS parameter unit family");
        }
    }
    return result;
}

std::vector<std::string> canonical_parameter_bases(
    const std::vector<NativeActiveCoordinate>& coordinates
) {
    return std::vector<std::string>(coordinates.size(), "physical");
}

std::vector<NativeActiveCoordinate> validated_active_parameters(
    const NativeModel& model,
    std::vector<NativeActiveCoordinate> coordinates
) {
    static_cast<void>(model.with_active_coordinates(coordinates));
    return coordinates;
}

std::string derived_parameter_contract_fingerprint(
    const std::string& model_fingerprint,
    const std::vector<NativeActiveCoordinate>& coordinates,
    const std::vector<std::string>& units,
    const std::vector<std::string>& bases
) {
    if (coordinates.size() != units.size() || coordinates.size() != bases.size()) {
        throw std::invalid_argument("EOS parameter contract metadata is inconsistent");
    }
    std::string payload = "epcsaft.eos-owner-parameter-contract.v1|";
    append_field(payload, model_fingerprint);
    for (std::size_t index = 0; index < coordinates.size(); ++index) {
        append_field(payload, coordinates[index].identity);
        std::ostringstream value;
        value << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10)
              << coordinates[index].value;
        append_field(payload, value.str());
        append_field(payload, units[index]);
        append_field(payload, bases[index]);
    }
    return "sha256:" + epcsaft::detail::sha256_digest(payload);
}

bool same_strings(
    const std::vector<std::string>& left,
    const std::vector<std::string>& right
) {
    return left.size() == right.size()
        && std::equal(left.begin(), left.end(), right.begin());
}

void require_finite(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

void require_digest(const std::string& value, const char* name) {
    if (!epcsaft::detail::sha256_fingerprint(value)) {
        throw std::invalid_argument(std::string(name) + " must be a sha256 fingerprint");
    }
}

const char* typed_failure_code(const std::string& reason) {
    if (reason == "association Jacobian is ill-conditioned") {
        return "ill_conditioned_sensitivity";
    }
    if (reason == "association Jacobian is singular") {
        return "singular_sensitivity";
    }
    if (reason == "association solve did not converge"
        || reason == "density root refinement did not converge"
        || reason == "density root pressure residual is too large"
        || reason == "maximum_iterations_exceeded"
        || reason == "evaluation_failed"
        || reason == "no finite in-domain saturation candidate"
        || reason == "coexistence feasibility or local physical checks failed") {
        return "numerical_convergence_failure";
    }
    if (reason == "association solve left physical domain"
        || reason == "density interval is empty"
        || reason == "packing fraction is outside the model domain") {
        return "eos_domain_rejection";
    }
    if (reason == "density residual is nonfinite"
        || reason == "mixture phase block contains nonfinite values"
        || reason == "fixed-topology phase block contains nonfinite values"
        || reason == "phase block contains nonfinite values"
        || reason == "active packing fraction block contains nonfinite values"
        || reason == "packing fraction block contains nonfinite values"
        || reason == "inverse packing geometry contains nonfinite values") {
        return "native_residual_failure";
    }
    return nullptr;
}

const char* exception_failure_code(const std::exception& error) {
    const auto code = typed_failure_code(error.what());
    return code == nullptr ? "native_residual_failure" : code;
}

const char* saturation_failure_code(const SaturationSolveResult& solve) {
    if (solve.rejection_status == "branch-loss") {
        return "branch_loss";
    }
    if (solve.rejection_status == "root-swap"
        || solve.rejection_status == "coalescence") {
        return "physical_rejection";
    }
    if (const auto code = typed_failure_code(solve.failure_reason);
        code != nullptr) {
        return code;
    }
    return "numerical_convergence_failure";
}

std::vector<EosRegressionOutputRow> require_rows(
    const EosRegressionSpec& specification,
    const std::vector<std::string>& component_ids
) {
    if (specification.family != EosRegressionFamily::FixedState
        && specification.family != EosRegressionFamily::PureSaturation) {
        throw std::invalid_argument("EOS evaluator observation family is unsupported");
    }
    const std::size_t rows = specification.row_identities.size();
    if (rows == 0) {
        throw std::invalid_argument("EOS evaluator rows are inconsistent");
    }
    std::vector<EosRegressionOutputRow> result;
    result.reserve(rows);
    for (std::size_t row = 0; row < rows; ++row) {
        const auto& identity = specification.row_identities[row];
        if (identity.empty()
            || std::count(specification.row_identities.begin(),
                          specification.row_identities.end(), identity) != 1) {
            throw std::invalid_argument("EOS evaluator row identities must be unique");
        }
        std::string unit;
        std::string basis;
        EosRegressionRowKind kind;
        if (identity == "helmholtz_over_rt_reference_amount") {
            unit = "1";
            basis = "reduced";
            kind = EosRegressionRowKind::FixedHelmholtz;
        } else if (identity == "pressure_pa") {
            unit = "pascal";
            basis = "absolute";
            kind = EosRegressionRowKind::FixedPressure;
        } else if (identity == "saturation_pressure_pa") {
            unit = "pascal";
            basis = "absolute";
            kind = EosRegressionRowKind::SaturationPressure;
        } else if (identity == "vapor_molar_density_mol_per_m3"
                   || identity == "liquid_molar_density_mol_per_m3") {
            unit = "mole / meter ** 3";
            basis = "molar";
            kind = identity == "vapor_molar_density_mol_per_m3"
                ? EosRegressionRowKind::VaporMolarDensity
                : EosRegressionRowKind::LiquidMolarDensity;
        } else if (identity == "vapor_mass_density_kg_per_m3"
                   || identity == "liquid_mass_density_kg_per_m3") {
            unit = "kilogram / meter ** 3";
            basis = "mass";
            kind = identity == "vapor_mass_density_kg_per_m3"
                ? EosRegressionRowKind::VaporMassDensity
                : EosRegressionRowKind::LiquidMassDensity;
        } else if (identity.rfind("chemical_potential_over_rt/", 0) == 0) {
            const auto component = identity.substr(27);
            if (std::find(component_ids.begin(), component_ids.end(), component)
                == component_ids.end()) {
                throw std::invalid_argument("EOS evaluator chemical row has an unknown component");
            }
            unit = "1";
            basis = "reduced";
            kind = EosRegressionRowKind::FixedChemicalPotential;
        } else {
            throw std::invalid_argument("EOS evaluator row identity is unsupported");
        }
        const bool fixed_kind = kind == EosRegressionRowKind::FixedHelmholtz
            || kind == EosRegressionRowKind::FixedPressure
            || kind == EosRegressionRowKind::FixedChemicalPotential;
        if (fixed_kind != (specification.family == EosRegressionFamily::FixedState)) {
            throw std::invalid_argument("EOS evaluator row is outside its observation family");
        }
        result.push_back({identity, std::move(unit), std::move(basis), kind});
    }
    return result;
}

std::vector<double> active_values(
    const std::vector<NativeActiveCoordinate>& coordinates
) {
    std::vector<double> values;
    values.reserve(coordinates.size());
    for (const auto& coordinate : coordinates) {
        values.push_back(coordinate.value);
    }
    return values;
}

}  // namespace

EosOwnerEvaluator::EosOwnerEvaluator(
    NativeModel model,
    EosRegressionSpec specification,
    std::vector<NativeActiveCoordinate> active_parameters
)
    : model_(std::move(model)),
      specification_(std::move(specification)),
      active_parameters_(validated_active_parameters(
          model_, std::move(active_parameters))),
      sidecar_{
          "eos",
          model_.parameter_fingerprint(),
          derived_parameter_contract_fingerprint(
              model_.parameter_fingerprint(), active_parameters_,
              canonical_parameter_units(active_parameters_),
              canonical_parameter_bases(active_parameters_)),
          specification_.source_spec_fingerprint,
          specification_.evaluator_fingerprint,
          active_parameter_identities(active_parameters_),
          active_values(active_parameters_),
          canonical_parameter_units(active_parameters_),
          canonical_parameter_bases(active_parameters_),
          require_rows(specification_, model_.component_ids()),
          "CppAD exact total physical parameter derivatives",
          specification_.family == EosRegressionFamily::FixedState
              ? "fixed-state EOS domain and finite exact Jacobian"
              : "local pure-saturation physical certificate and finite exact Jacobian",
          model_.eos_metadata().topology_fingerprint,
          model_.eos_metadata().component_order_fingerprint,
          model_.eos_metadata().domain_fingerprint,
          model_.eos_metadata().artifact_fingerprint,
          model_.eos_metadata().candidate_qualification,
          model_.eos_metadata().candidate_source_fingerprint,
          model_.eos_metadata().temperature_min_k,
          model_.eos_metadata().temperature_max_k,
          model_.eos_metadata().pressure_min_pa,
          model_.eos_metadata().pressure_max_pa,
          model_.eos_metadata().ion_mole_fraction_max,
      } {
    if (active_parameters_.empty()) {
        throw std::invalid_argument("EOS evaluator active parameter set must be nonempty");
    }
    if (specification_.active_parameter_units.size() != active_parameters_.size()
        || specification_.active_parameter_bases.size() != active_parameters_.size()
        || std::any_of(specification_.active_parameter_units.begin(),
                       specification_.active_parameter_units.end(),
                       [](const std::string& unit) { return unit.empty(); })) {
        throw std::invalid_argument("EOS evaluator active parameter units are inconsistent");
    }
    if (std::any_of(specification_.active_parameter_bases.begin(),
                    specification_.active_parameter_bases.end(),
                    [](const std::string& basis) { return basis.empty(); })) {
        throw std::invalid_argument("EOS evaluator active parameter bases are inconsistent");
    }
    if (!same_strings(specification_.active_parameter_units,
                      sidecar_.active_parameter_units)
        || !same_strings(specification_.active_parameter_bases,
                         sidecar_.active_parameter_bases)) {
        throw std::invalid_argument("EOS evaluator active parameter metadata is not canonical");
    }
    require_digest(specification_.source_spec_fingerprint,
                   "source specification fingerprint");
    require_digest(specification_.evaluator_fingerprint,
                   "evaluator fingerprint");
    const auto identities = active_parameter_identities(active_parameters_);
    if (!same_strings(identities, specification_.active_parameter_identities)) {
        throw std::invalid_argument("EOS evaluator active parameter order changed");
    }
    if (!same_strings(identities, model_.active_parameter_identities())
        && !model_.active_coordinates().empty()) {
        throw std::invalid_argument("EOS evaluator active parameter order changed");
    }
    if (std::any_of(identities.begin(), identities.end(),
                    [](const std::string& identity) { return identity.empty(); })) {
        throw std::invalid_argument("EOS evaluator active parameter identities are invalid");
    }
    for (const auto& coordinate : active_parameters_) {
        require_finite(coordinate.value, "active parameter value");
    }
    if (specification_.family == EosRegressionFamily::FixedState) {
        require_finite(specification_.temperature_k, "fixed-state temperature");
        require_finite(specification_.molar_density_mol_m3,
                       "fixed-state molar density");
        if (specification_.molar_density_mol_m3 <= 0.0) {
            throw std::invalid_argument("fixed-state molar density must be positive");
        }
        if (specification_.mole_fractions.size() != model_.component_count()
            || !finite(specification_.mole_fractions)) {
            throw std::invalid_argument("fixed-state composition is inconsistent");
        }
        const double composition_total = std::accumulate(
            specification_.mole_fractions.begin(),
            specification_.mole_fractions.end(), 0.0
        );
        if (std::abs(composition_total - 1.0) > 1.0e-12) {
            throw std::invalid_argument("fixed-state composition must be normalized");
        }
    } else if (specification_.family == EosRegressionFamily::PureSaturation) {
        if (model_.component_count() != 1) {
            throw std::invalid_argument("pure-saturation evaluator requires one component");
        }
        require_finite(specification_.temperature_k, "saturation temperature");
    } else {
        throw std::invalid_argument("EOS evaluator observation family is unsupported");
    }
}

std::size_t EosOwnerEvaluator::parameter_count() const noexcept {
    return active_parameters_.size();
}

std::size_t EosOwnerEvaluator::output_count() const noexcept {
    return specification_.row_identities.size();
}

const std::string& EosOwnerEvaluator::parameter_contract_fingerprint()
    const noexcept {
    return sidecar_.parameter_contract_fingerprint;
}

const std::string& EosOwnerEvaluator::evaluator_fingerprint()
    const noexcept {
    return sidecar_.evaluator_fingerprint;
}

const EosRegressionEvaluatorSidecar& EosOwnerEvaluator::sidecar()
    const noexcept {
    return sidecar_;
}

std::vector<NativeActiveCoordinate> EosOwnerEvaluator::coordinates(
    const std::vector<double>& physical_parameters
) const {
    std::vector<NativeActiveCoordinate> result = active_parameters_;
    for (std::size_t index = 0; index < result.size(); ++index) {
        result[index].value = physical_parameters[index];
    }
    return result;
}

EvaluationResult EosOwnerEvaluator::evaluate_fixed_state(
    const std::vector<double>& physical_parameters,
    bool exact_jacobian
) const {
    const auto current_coordinates = coordinates(physical_parameters);
    const auto current_model = model_.with_active_coordinates(current_coordinates);
    const auto block = current_model.evaluate_phase_derivatives(
        specification_.temperature_k,
        specification_.mole_fractions,
        1.0 / specification_.molar_density_mol_m3,
        current_coordinates,
        2
    );
    const std::size_t components = model_.component_count();
    const std::size_t parameters = active_parameters_.size();
    const std::size_t width = components + 1 + parameters;
    if (block.gradient.size() != width || block.hessian.size() != width * width) {
        return failure("owner_contract_failure", "EOS phase derivative tensor dimensions changed");
    }
    Evaluation result;
    result.predictions.reserve(output_count());
    if (exact_jacobian) {
        result.jacobian.reserve(output_count() * parameters);
    }
    auto append = [&](double value, const std::vector<double>& derivatives) {
        result.predictions.push_back(value);
        if (exact_jacobian) {
            result.jacobian.insert(result.jacobian.end(), derivatives.begin(), derivatives.end());
        }
    };
    for (const auto& row : sidecar_.output_rows) {
        std::vector<double> derivatives;
        if (exact_jacobian) {
            derivatives.resize(parameters);
        }
        double value = 0.0;
        switch (row.kind) {
        case EosRegressionRowKind::FixedHelmholtz:
            value = block.helmholtz_over_rt_reference_amount;
            if (exact_jacobian) {
                for (std::size_t parameter = 0; parameter < parameters; ++parameter) {
                    derivatives[parameter] = block.gradient[components + 1 + parameter];
                }
            }
            break;
        case EosRegressionRowKind::FixedPressure:
            value = block.pressure_pa;
            if (exact_jacobian) {
                for (std::size_t parameter = 0; parameter < parameters; ++parameter) {
                    derivatives[parameter] = -GAS_CONSTANT
                        * specification_.temperature_k
                        * block.hessian[components * width + components + 1 + parameter];
                }
            }
            break;
        case EosRegressionRowKind::FixedChemicalPotential: {
            constexpr std::size_t prefix_size = 27;
            const auto component_id = row.identity.substr(prefix_size);
            const auto component = std::find(model_.component_ids().begin(),
                                             model_.component_ids().end(),
                                             component_id);
            if (component == model_.component_ids().end()) {
                return failure("owner_contract_failure", "EOS chemical-potential row has an unknown component");
            }
            const std::size_t component_index = static_cast<std::size_t>(
                component - model_.component_ids().begin());
            value = block.gradient[component_index];
            if (exact_jacobian) {
                for (std::size_t parameter = 0; parameter < parameters; ++parameter) {
                    derivatives[parameter] = block.hessian[
                        component_index * width + components + 1 + parameter
                    ];
                }
            }
            break;
        }
        default:
            return failure("owner_contract_failure", "EOS fixed-state row kind changed");
        }
        append(value, derivatives);
    }
    if (!finite(result.predictions)) {
        return failure("native_residual_failure",
                       "EOS fixed-state prediction is nonfinite");
    }
    if (exact_jacobian && !finite(result.jacobian)) {
        return failure("native_residual_failure",
                       "EOS fixed-state derivative is nonfinite");
    }
    return result;
}

EvaluationResult EosOwnerEvaluator::evaluate_pure_saturation(
    const std::vector<double>& physical_parameters,
    bool exact_jacobian
) const {
    const auto current_coordinates = coordinates(physical_parameters);
    SaturationSensitivityResult solve;
    if (exact_jacobian) {
        solve = solve_local_saturation_active(model_, specification_.temperature_k,
                                              current_coordinates);
    } else {
        const auto current_model = model_.with_active_coordinates(current_coordinates);
        const double upper_density = current_model.molar_density_upper_bound(
            specification_.temperature_k, {1.0});
        solve = SaturationSensitivityResult{
            solve_local_saturation(current_model, specification_.temperature_k,
                                   upper_density),
            {},
            {},
        };
    }
    if (!solve.solve.accepted) {
        const auto reason = solve.solve.failure_reason.empty()
            ? solve.solve.rejection_status
            : solve.solve.failure_reason;
        return failure(saturation_failure_code(solve.solve), reason);
    }
    if (exact_jacobian && solve.sensitivity_status != "available") {
        return failure("sensitivity_unavailable",
                       solve.sensitivity_failure_reason.empty()
                           ? "saturation sensitivity is unavailable"
                           : solve.sensitivity_failure_reason);
    }
    const std::size_t parameters = active_parameters_.size();
    if (exact_jacobian
        && (solve.parameter_identities !=
                active_parameter_identities(active_parameters_)
            || solve.derivatives.size() != 3 * parameters)) {
        return failure("owner_contract_failure", "EOS saturation derivative dimensions changed");
    }
    const auto& result = solve.solve.evaluation;
    const double vapor_density = 1.0 / result.vapor.volume_m3;
    const double liquid_density = 1.0 / result.liquid.volume_m3;
    const auto masses = model_.molar_masses_kg_per_mol();
    const double molar_mass = masses.size() == 1 ? masses[0] : 0.0;
    Evaluation output;
    output.predictions.reserve(output_count());
    if (exact_jacobian) {
        output.jacobian.reserve(output_count() * parameters);
    }
    for (const auto& row : sidecar_.output_rows) {
        double value = 0.0;
        std::size_t source_row = 0;
        double scale = 1.0;
        switch (row.kind) {
        case EosRegressionRowKind::SaturationPressure:
            value = result.saturation_pressure_pa;
            source_row = 0;
            break;
        case EosRegressionRowKind::VaporMolarDensity:
            value = vapor_density;
            source_row = 1;
            break;
        case EosRegressionRowKind::LiquidMolarDensity:
            value = liquid_density;
            source_row = 2;
            break;
        case EosRegressionRowKind::VaporMassDensity:
            if (!(molar_mass > 0.0)) {
                return failure("unsupported_observable", "saturated mass density requires a positive molar mass");
            }
            value = vapor_density * molar_mass;
            source_row = 1;
            scale = molar_mass;
            break;
        case EosRegressionRowKind::LiquidMassDensity:
            if (!(molar_mass > 0.0)) {
                return failure("unsupported_observable", "saturated mass density requires a positive molar mass");
            }
            value = liquid_density * molar_mass;
            source_row = 2;
            scale = molar_mass;
            break;
        default:
            return failure("owner_contract_failure", "EOS saturation row kind changed");
        }
        output.predictions.push_back(value);
        if (exact_jacobian) {
            for (std::size_t parameter = 0; parameter < parameters; ++parameter) {
                output.jacobian.push_back(scale * solve.derivatives[
                    source_row * parameters + parameter
                ]);
            }
        }
    }
    if (!finite(output.predictions)) {
        return failure("native_residual_failure",
                       "EOS saturation prediction is nonfinite");
    }
    if (exact_jacobian && !finite(output.jacobian)) {
        return failure("native_residual_failure",
                       "EOS saturation derivative is nonfinite");
    }
    return output;
}

EvaluationResult EosOwnerEvaluator::evaluate(
    const std::vector<double>& physical_parameters,
    bool exact_jacobian
) const noexcept {
    try {
        if (physical_parameters.size() != parameter_count()) {
            return failure("parameter_topology_mismatch", "EOS evaluator physical parameter dimension changed");
        }
        if (!finite(physical_parameters)) {
            return failure("eos_domain_rejection", "EOS evaluator physical parameters must be finite");
        }
        return specification_.family == EosRegressionFamily::FixedState
            ? evaluate_fixed_state(physical_parameters, exact_jacobian)
            : evaluate_pure_saturation(physical_parameters, exact_jacobian);
    } catch (const EosDomainError& error) {
        return failure("eos_domain_rejection", error.what());
    } catch (const std::invalid_argument& error) {
        return failure("eos_domain_rejection", error.what());
    } catch (const std::domain_error& error) {
        return failure("eos_domain_rejection", error.what());
    } catch (const std::exception& error) {
        return failure(exception_failure_code(error), error.what());
    } catch (...) {
        return failure("native_residual_failure", "unknown EOS owner exception");
    }
}

CompiledEosRegressionEvaluator compile_eos_regression_evaluator(
    NativeModel model,
    EosRegressionSpec specification,
    std::vector<NativeActiveCoordinate> active_parameters
) {
    auto owner = std::make_shared<const EosOwnerEvaluator>(
        std::move(model), std::move(specification), std::move(active_parameters)
    );
    std::shared_ptr<const EosRegressionEvaluatorSidecar> sidecar(
        owner, &owner->sidecar()
    );
    return {std::move(owner), std::move(sidecar)};
}

}  // namespace epcsaft
