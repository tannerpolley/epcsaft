#pragma once

#include "core/physical_constants.hpp"
#include "core/scalar.hpp"
#include "core/trace.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace epcsaft {

enum class MsaEvaluation { Initial, Residual, Slope, Energy };

template <typename Scalar>
Scalar msa_ion_volume(const Scalar& n, const Scalar& diameter) {
    return n * diameter * diameter * diameter;
}

template <typename Scalar>
struct MsaMoments final {
    Scalar denominator = Scalar(2.0 / PI), packing{}, coupling{}, q{}, q_gamma{}, h_moment{}, charge_square{}, energy{};
};

// The input sequence is the Bjerrum length followed by (number density, diameter, valence) for each ion.
// EqIDs primitive_msa_auxiliaries, primitive_msa_screening, primitive_msa_helmholtz
template <MsaEvaluation evaluation, typename Scalar>
MsaMoments<Scalar> msa_moments(const Scalar& gamma, const std::vector<Scalar>& inputs) {
    MsaMoments<Scalar> m;
    for (std::size_t i = 1; i < inputs.size(); i += 3) {
        const auto& n = inputs[i];
        const auto& diameter = inputs[i + 1];
        const auto& charge = inputs[i + 2];
        const Scalar h = Scalar(1.0) + gamma * diameter;
        const Scalar volume = msa_ion_volume(n, diameter);
        m.packing += volume;
        m.denominator += volume / h;
        m.coupling += n * charge * diameter / h;
        if constexpr (evaluation == MsaEvaluation::Initial) m.charge_square += n * charge * charge;
        if constexpr (evaluation == MsaEvaluation::Energy) m.energy += n * charge * charge / h;
    }
    m.denominator -= m.packing / Scalar(3.0);
    m.coupling /= m.denominator;
    if constexpr (evaluation == MsaEvaluation::Energy) {
        m.energy = gamma * m.energy + m.denominator * m.coupling * m.coupling;
    } else {
        Scalar moment{};
        for (std::size_t i = 1; i < inputs.size(); i += 3) {
            const auto& n = inputs[i];
            const auto& diameter = inputs[i + 1];
            const Scalar h = Scalar(1.0) + gamma * diameter;
            const Scalar effective = (inputs[i + 2] - m.coupling * diameter * diameter) / h;
            m.q += n * effective * effective;
            if constexpr (evaluation == MsaEvaluation::Slope) {
                m.h_moment += n * effective * diameter * diameter / h;
                moment += n * effective * effective * diameter / h;
            }
        }
        if constexpr (evaluation == MsaEvaluation::Slope)
            m.q_gamma = Scalar(2.0) * (m.h_moment * m.h_moment / m.denominator - moment);
    }
    return m;
}

template <typename Scalar>
// EqID: primitive_msa_screening
Scalar msa_residual(const Scalar& gamma, const std::vector<Scalar>& inputs) {
    const auto m = msa_moments<MsaEvaluation::Residual>(gamma, inputs);
    return gamma - scalar_sqrt(Scalar(PI) * inputs.front()) * scalar_sqrt(m.q);
}

template <typename Scalar>
Scalar solve_msa_gamma(const std::vector<Scalar>& inputs) {
    core::Scope scope("msa_screening");
    if (!(inputs.front() > 0.0) || !all_finite(inputs))
        throw std::domain_error("MSA requires finite inputs and a positive Bjerrum length");
    for (std::size_t i = 1; i < inputs.size(); i += 3)
        if (inputs[i] < 0.0 || !(inputs[i + 1] > 0.0))
            throw std::domain_error("MSA requires nonnegative ion densities and positive contact diameters");
    const auto initial = msa_moments<MsaEvaluation::Initial>(Scalar(0.0), inputs);
    if (!(initial.packing < Scalar(6.0 / PI)))
        throw std::domain_error("MSA ion packing requires positive Delta");
    const Scalar scale = scalar_sqrt(Scalar(PI) * inputs.front());
    Scalar lo = 0.0, hi = scale * scalar_sqrt(initial.q);
    if (!(hi > 0.0) || !wide_math::isfinite(hi))
        throw std::domain_error("MSA requires a finite positive screening bracket");
    Scalar gamma = std::min(hi, scale * scalar_sqrt(initial.charge_square));
    const auto maximum = static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value);
    const auto& target = policy::DENSITY_ROOT_RESIDUAL;
    double relative_residual = 0.0;
    std::size_t evaluations = 0;
    while (evaluations < maximum) {
        ++evaluations;
        const auto m = msa_moments<MsaEvaluation::Slope>(gamma, inputs);
        const Scalar root_q = scalar_sqrt(m.q);
        const Scalar slope = Scalar(1.0) - scale * m.q_gamma / (Scalar(2.0) * root_q);
        const Scalar residual = gamma - scale * root_q;
        relative_residual = static_cast<double>(wide_math::abs(residual) / gamma);
        if (relative_residual <= target.value || !std::isfinite(relative_residual) || evaluations == maximum) break;
        if (residual > 0.0) hi = gamma;
        else lo = gamma;
        const Scalar next = gamma - residual / slope;
        gamma = lo < next && next < hi ? next : lo + (hi - lo) / Scalar(2.0);
    }
    if (!scope.on() && relative_residual <= target.value) return gamma;
    core::Convergence record;
    record.stop = core::classify(relative_residual, target, std::nullopt, std::nullopt, evaluations == maximum);
    record.cause = "MSA screening";
    record.residual = relative_residual;
    record.iterations = evaluations - 1;
    record.evaluations = evaluations;
    record.criterion = std::make_shared<const policy::Tolerance>(target);
    scope.finish(record);
    if (record.stop != core::Stop::converged) throw core::ConvergenceFailure(record);
    return gamma;
}

template <std::size_t N, typename Base>
void msa_gamma_coefficients(const std::vector<Taylor<N, Base>>& inputs, Taylor<N, Base>& gamma,
                           std::size_t first, std::size_t last) {
    std::vector<Base> central;
    for (const auto& input : inputs) central.push_back(input.c[0]);
    const auto m = msa_moments<MsaEvaluation::Slope>(gamma.c[0], central);
    const Base slope = Base(1.0) - scalar_sqrt(Base(PI) * central.front()) * m.q_gamma
        / (Base(2.0) * scalar_sqrt(m.q));
    for (std::size_t k = first; k <= last; ++k)
        gamma.c[k] = -msa_residual(gamma, inputs).c[k] / slope;
}

template <typename Scalar>
Scalar msa_gamma(const std::vector<Scalar>& inputs) {
    if constexpr (is_taylor<Scalar>) {
        using Base = std::decay_t<decltype(Scalar{}.c[0])>;
        std::vector<Base> central;
        for (const auto& input : inputs) central.push_back(input.c[0]);
        Scalar gamma(msa_gamma(central));
        msa_gamma_coefficients(inputs, gamma, 1, gamma.c.size() - 1);
        return gamma;
    } else if constexpr (std::is_floating_point_v<Scalar> || std::is_same_v<Scalar, WiderScalar>) {
        return solve_msa_gamma(inputs);
    } else {
        throw std::logic_error("MSA screening requires the native Taylor or double CppAD scalar");
    }
}

template <typename Scalar>
// EqID: primitive_msa_helmholtz
Scalar msa_helmholtz_density(const Scalar& gamma, const std::vector<Scalar>& inputs) {
    return -inputs.front() * msa_moments<MsaEvaluation::Energy>(gamma, inputs).energy + gamma * gamma * gamma / Scalar(3.0 * PI);
}

}  // namespace epcsaft
