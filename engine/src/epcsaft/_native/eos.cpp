#include "eos.hpp"

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

NativeResult NativeModel::evaluate(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions
) const {
    if (!std::isfinite(temperature) || temperature <= 0.0 ||
        !std::isfinite(molar_density) || molar_density <= 0.0 ||
        mole_fractions.size() != m_.size() || !finite(mole_fractions)) {
        throw std::invalid_argument("invalid native state");
    }

    const auto values = detail::residual_helmholtz(
        temperature, molar_density, mole_fractions,
        m_, sigma_, epsilon_over_k_, kij_
    );

    using ADScalar = CppAD::AD<double>;
    std::vector<ADScalar> independent(1);
    independent[0] = molar_density;
    CppAD::Independent(independent);
    const auto taped = detail::residual_helmholtz(
        temperature, independent[0], mole_fractions,
        m_, sigma_, epsilon_over_k_, kij_
    );
    std::vector<ADScalar> dependent = {taped.hard_chain + taped.dispersion};
    CppAD::ADFun<double> function(independent, dependent);
    const std::vector<double> state = {molar_density};
    const double derivative = function.Jacobian(state)[0];

    const double residual = values.hard_chain + values.dispersion;
    const double compressibility = 1.0 + molar_density * derivative;
    return {
        values.hard_chain,
        values.dispersion,
        residual,
        compressibility,
        compressibility * molar_density * detail::GAS_CONSTANT * temperature,
    };
}

}  // namespace epcsaft
