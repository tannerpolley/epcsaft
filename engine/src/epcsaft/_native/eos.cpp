#include "eos.hpp"

#include <cppad/cppad.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace epcsaft {
namespace {

bool positive_finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value) && value > 0.0;
    });
}

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

}  // namespace

NativeModel::NativeModel(
    std::vector<double> m,
    std::vector<double> sigma,
    std::vector<double> epsilon_over_k,
    std::vector<double> kij
)
    : m_(std::move(m)),
      sigma_(std::move(sigma)),
      epsilon_over_k_(std::move(epsilon_over_k)),
      kij_(std::move(kij)) {
    const auto size = m_.size();
    if (size == 0 || sigma_.size() != size || epsilon_over_k_.size() != size ||
        kij_.size() != size * size) {
        throw std::invalid_argument("inconsistent native model dimensions");
    }
    if (!positive_finite(m_) || !positive_finite(sigma_) ||
        !positive_finite(epsilon_over_k_) || !finite(kij_)) {
        throw std::invalid_argument("invalid native model parameters");
    }
}

}  // namespace epcsaft
