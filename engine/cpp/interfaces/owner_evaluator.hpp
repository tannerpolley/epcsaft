#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <variant>
#include <vector>

namespace epcsaft_engine::interfaces {

struct Evaluation final {
  std::vector<double> predictions;
  std::vector<double> jacobian;
};

struct EvaluationFailure final {
  std::string code;
  std::string reason;
};

using EvaluationResult = std::variant<Evaluation, EvaluationFailure>;

class OwnerEvaluator {
public:
  virtual ~OwnerEvaluator() = default;

  [[nodiscard]] virtual std::size_t parameter_count() const noexcept = 0;
  [[nodiscard]] virtual std::size_t output_count() const noexcept = 0;
  [[nodiscard]] virtual const std::string &
  parameter_contract_fingerprint() const noexcept = 0;
  [[nodiscard]] virtual const std::string &
  evaluator_fingerprint() const noexcept = 0;
  [[nodiscard]] virtual EvaluationResult
  evaluate(const std::vector<double> &physical_parameters,
           bool exact_jacobian) const noexcept = 0;
};

using OwnerEvaluatorHandle = std::shared_ptr<const OwnerEvaluator>;

} // namespace epcsaft_engine::interfaces
