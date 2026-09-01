#pragma once

#include "eos.hpp"
#include "owner_evaluator.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace epcsaft {

enum class EosRegressionFamily {
    FixedState,
    PureSaturation,
};

struct EosRegressionSpec final {
    EosRegressionFamily family;
    double temperature_k = 0.0;
    double molar_density_mol_m3 = 0.0;
    std::vector<double> mole_fractions;
    std::vector<std::string> active_parameter_identities;
    std::vector<std::string> active_parameter_units;
    std::vector<std::string> active_parameter_bases;
    std::vector<std::string> row_identities;
    std::string parameter_contract_fingerprint;
    std::string source_spec_fingerprint;
    std::string evaluator_fingerprint;
};

enum class EosRegressionRowKind {
    FixedHelmholtz,
    FixedPressure,
    FixedChemicalPotential,
    SaturationPressure,
    VaporMolarDensity,
    LiquidMolarDensity,
    VaporMassDensity,
    LiquidMassDensity,
};

struct EosRegressionOutputRow final {
    std::string identity;
    std::string unit;
    std::string basis;
    EosRegressionRowKind kind = EosRegressionRowKind::FixedHelmholtz;
};

struct EosRegressionEvaluatorSidecar final {
    std::string owner_kind;
    std::string parameter_fingerprint;
    std::string parameter_contract_fingerprint;
    std::string source_spec_fingerprint;
    std::string evaluator_fingerprint;
    std::vector<std::string> active_parameter_identities;
    std::vector<double> active_parameter_values;
    std::vector<std::string> active_parameter_units;
    std::vector<std::string> active_parameter_bases;
    std::vector<EosRegressionOutputRow> output_rows;
    std::string derivative_provenance;
    std::string acceptance_policy;
    std::string topology_fingerprint;
    std::string component_order_fingerprint;
    std::string domain_fingerprint;
    std::string artifact_fingerprint;
    std::string candidate_qualification;
    std::string candidate_source_fingerprint;
    double temperature_min_k;
    double temperature_max_k;
    double pressure_min_pa;
    double pressure_max_pa;
    double ion_mole_fraction_max;
};

class EosOwnerEvaluator final
    : public epcsaft_engine::interfaces::OwnerEvaluator {
public:
    EosOwnerEvaluator(
        NativeModel model,
        EosRegressionSpec specification,
        std::vector<NativeActiveCoordinate> active_parameters
    );

    [[nodiscard]] std::size_t parameter_count() const noexcept override;
    [[nodiscard]] std::size_t output_count() const noexcept override;
    [[nodiscard]] const std::string& parameter_contract_fingerprint()
        const noexcept override;
    [[nodiscard]] const std::string& evaluator_fingerprint()
        const noexcept override;
    [[nodiscard]] epcsaft_engine::interfaces::EvaluationResult evaluate(
        const std::vector<double>& physical_parameters,
        bool exact_jacobian
    ) const noexcept override;
    [[nodiscard]] const EosRegressionEvaluatorSidecar& sidecar()
        const noexcept;

private:
    [[nodiscard]] epcsaft_engine::interfaces::EvaluationResult
    evaluate_fixed_state(
        const std::vector<double>& physical_parameters,
        bool exact_jacobian
    ) const;
    [[nodiscard]] epcsaft_engine::interfaces::EvaluationResult
    evaluate_pure_saturation(
        const std::vector<double>& physical_parameters,
        bool exact_jacobian
    ) const;

    [[nodiscard]] std::vector<NativeActiveCoordinate> coordinates(
        const std::vector<double>& physical_parameters
    ) const;

    const NativeModel model_;
    const EosRegressionSpec specification_;
    const std::vector<NativeActiveCoordinate> active_parameters_;
    const EosRegressionEvaluatorSidecar sidecar_;
};

struct CompiledEosRegressionEvaluator final {
    std::shared_ptr<const epcsaft_engine::interfaces::OwnerEvaluator> evaluator;
    std::shared_ptr<const EosRegressionEvaluatorSidecar> sidecar;
};

[[nodiscard]] CompiledEosRegressionEvaluator compile_eos_regression_evaluator(
    NativeModel model,
    EosRegressionSpec specification,
    std::vector<NativeActiveCoordinate> active_parameters
);

}  // namespace epcsaft
