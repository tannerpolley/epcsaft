#pragma once

#include "core/physical_constants.hpp"
#include "core/scalar.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

namespace epcsaft {

// Constant heat capacity is the one-coefficient centered polynomial.
struct IdealPolynomial final {
    std::vector<double> coefficients;
    double reference_temperature, reference_enthalpy, reference_entropy;
};
struct IdealShomate final {
    std::array<double, 8> coefficients; // A, B, C, D, E, F, G, H; NIST units
    double formation_enthalpy; // J/mol; added to the unmodified fitted increment
};
struct IdealNASA7 final { std::array<double, 7> coefficients; };
struct IdealNASA9 final { std::array<double, 9> coefficients; };
struct IdealDIPPR107 final {
    std::array<double, 5> coefficients; // A, B, C, D, E; SI molar heat capacity
    double reference_temperature, reference_enthalpy, reference_entropy;
};
struct IdealCorrelation final {
    std::variant<IdealPolynomial, IdealShomate, IdealNASA7, IdealNASA9, IdealDIPPR107> form;
    double source_pressure;
};
struct IdealInterval final {
    double lower, upper;
    bool lower_inclusive, upper_inclusive;
    IdealCorrelation correlation;
};
template <typename Scalar>
struct IdealThermalValue final { Scalar cp, h, s; };

// EqExclude: explicit coefficient and domain admission
inline void validate_ideal_correlation(const IdealCorrelation& correlation) {
    if (!(correlation.source_pressure > 0.0) || !std::isfinite(correlation.source_pressure))
        throw std::invalid_argument("invalid ideal source pressure");
    std::visit([](const auto& c) {
        using C = std::decay_t<decltype(c)>;
        if (c.coefficients.empty()) throw std::invalid_argument("missing ideal coefficients");
        for (double value : c.coefficients)
            if (!std::isfinite(value)) throw std::invalid_argument("nonfinite ideal coefficient");
        if constexpr (std::is_same_v<C, IdealPolynomial> || std::is_same_v<C, IdealDIPPR107>) {
            if (!(c.reference_temperature > 0.0) || !std::isfinite(c.reference_temperature)
                || !std::isfinite(c.reference_enthalpy) || !std::isfinite(c.reference_entropy))
                throw std::invalid_argument("invalid ideal thermal anchor");
        } else if constexpr (std::is_same_v<C, IdealShomate>) {
            if (!std::isfinite(c.formation_enthalpy))
                throw std::invalid_argument("missing Shomate formation enthalpy");
        }
    }, correlation.form);
}

inline bool ideal_interval_contains(const IdealInterval& interval, double temperature) {
    return (temperature > interval.lower
            || (interval.lower_inclusive && temperature == interval.lower))
        && (temperature < interval.upper
            || (interval.upper_inclusive && temperature == interval.upper));
}

inline void validate_ideal_intervals(const std::vector<IdealInterval>& intervals) {
    if (intervals.empty()) throw std::invalid_argument("missing ideal intervals");
    std::vector<std::size_t> order(intervals.size());
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        order[i] = i;
        const auto& interval = intervals[i];
        if (!(interval.lower > 0.0) || !(interval.upper > interval.lower)
            || !std::isfinite(interval.upper))
            throw std::invalid_argument("invalid ideal interval");
        validate_ideal_correlation(interval.correlation);
        std::visit([&](const auto& c) {
            using C = std::decay_t<decltype(c)>;
            if constexpr (std::is_same_v<C, IdealPolynomial>
                          || std::is_same_v<C, IdealDIPPR107>) {
                if (!ideal_interval_contains(interval, c.reference_temperature))
                    throw std::invalid_argument("ideal anchor outside its correlation interval");
            }
        }, interval.correlation.form);
    }
    std::sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
        return intervals[left].lower < intervals[right].lower;
    });
    for (std::size_t i = 1; i < order.size(); ++i) {
        const auto& previous = intervals[order[i - 1]];
        const auto& current = intervals[order[i]];
        if (current.lower != previous.upper
            || current.lower_inclusive == previous.upper_inclusive)
            throw std::invalid_argument("overlapping or gapped ideal intervals");
    }
}

// EqID: cp_ideal_polynomial
// Integral of (T-Tref)^k/T; ALG-CALORICS bounds each representation.
// Multiply the coefficient before powers. Only immutable degree controls loops;
// clamped inactive branches remain finite during same-interval AD replay.
template <typename Scalar>
Scalar ideal_polynomial_entropy(const std::vector<double>& coefficients,
                                const Scalar& temperature, double reference) {
    const Scalar delta = temperature - Scalar(reference), u = delta / Scalar(reference);
    const Scalar magnitude = scalar_select_less(delta, Scalar(0.0), -delta, delta);
    const Scalar positive = scalar_select_less(u, Scalar(0.0), Scalar(0.0), u);
    const Scalar negative = scalar_select_less(u, Scalar(0.0), -u, Scalar(0.0));
    const Scalar near_u = scalar_select_less(
        u, Scalar(-0.5), Scalar(-0.5), scalar_select_less(u, Scalar(0.5), u, Scalar(0.5)));
    const Scalar log_ratio = scalar_select_less(
        magnitude, Scalar(reference * 0.5), scalar_log1p(near_u),
        scalar_log(temperature) - Scalar(std::log(reference)));
    Scalar result = Scalar(coefficients[0]) * log_ratio;
    for (std::size_t k = 1; k < coefficients.size(); ++k) {
        if (coefficients[k] == 0.0) continue;
        Scalar amplitude = Scalar(coefficients[k]);
        for (std::size_t p = 0; p < k; ++p) amplitude *= magnitude;
        const Scalar small_u = scalar_select_less(positive, Scalar(1.0), positive, Scalar(1.0));
        const Scalar ratio = small_u / (Scalar(1.0) + small_u);
        Scalar term = Scalar(1.0 / (k + 1)), positive_sum = term;
        for (int j = 0; j < std::numeric_limits<double>::digits + 2; ++j) {
            term *= ratio * Scalar(j + 1) / Scalar(k + j + 2);
            positive_sum += term;
        }
        const Scalar inverse_u = Scalar(reference) / scalar_select_less(
            delta, Scalar(reference), Scalar(reference), delta);
        Scalar large_positive = scalar_select_less(u, Scalar(1.0), Scalar(std::log(2.0)), log_ratio);
        for (std::size_t p = 1; p <= k; ++p)
            large_positive = Scalar(1.0 / p) - large_positive * inverse_u;
        const Scalar positive_moment = scalar_select_less(
            positive, Scalar(1.0), positive_sum * ratio, large_positive);

        // v^k >= exp(-1) bounds cancellation in the logarithmic tail.
        const double crossover = std::exp(-1.0 / k);
        const Scalar v = scalar_select_less(negative, Scalar(crossover), negative, Scalar(crossover));
        const std::size_t terms = (std::numeric_limits<double>::digits + 2
            + static_cast<std::size_t>(std::ceil(3.0 * std::log2(k + 1)))) * (k + 1);
        Scalar negative_sum = Scalar(0.0);
        for (std::size_t j = terms; j-- > 0;)
            negative_sum = Scalar(1.0 / (k + j + 1)) + v * negative_sum;
        const Scalar tail_v = scalar_select_less(
            negative, Scalar(crossover), Scalar(crossover), negative);
        Scalar tail = scalar_select_less(
            negative, Scalar(crossover), Scalar(-std::log1p(-crossover)), Scalar(-log_ratio));
        Scalar power = Scalar(1.0);
        for (std::size_t p = 1; p <= k; ++p) {
            power *= tail_v;
            tail -= power / Scalar(p);
        }
        const Scalar negative_moment = scalar_select_less(
            negative, Scalar(crossover), v * negative_sum, tail / power)
            * Scalar(k % 2 == 0 ? -1.0 : 1.0);
        result += amplitude * scalar_select_less(u, Scalar(0.0), negative_moment, positive_moment);
    }
    return result;
}

// EqID: cp_ideal_constant
// EqID: cp_ideal_polynomial
// EqID: cp_ideal_shomate
// EqID: cp_ideal_nasa7
// EqID: cp_ideal_nasa9
// EqID: cp_ideal_dippr_107
// Scalar recording is local to one selected source interval. Reselect the
// interval before reusing a tape; Properties owns that future integration.
template <typename Scalar>
IdealThermalValue<Scalar> evaluate_ideal_correlation(
    const IdealCorrelation& correlation, const Scalar& temperature, double pressure
) {
    validate_ideal_correlation(correlation);
    if (!(pressure > 0.0) || !std::isfinite(pressure))
        throw std::invalid_argument("invalid ideal standard pressure");
    if constexpr (std::is_same_v<Scalar, double>) {
        if (!(temperature > 0.0) || !std::isfinite(temperature))
            throw std::invalid_argument("invalid ideal temperature");
    }
    const double R = GAS_CONSTANT_J_PER_MOL_K;
    auto result = std::visit([&](const auto& c) -> IdealThermalValue<Scalar> {
        using C = std::decay_t<decltype(c)>;
        const auto& a = c.coefficients;
        const Scalar& T = temperature;
        if constexpr (std::is_same_v<C, IdealPolynomial>) {
            const Scalar dt = T - Scalar(c.reference_temperature);
            Scalar cp = Scalar(0.0), dh = Scalar(0.0);
            for (std::size_t k = a.size(); k-- > 0;) {
                cp = cp * dt + Scalar(a[k]);
                dh = dh * dt + Scalar(a[k] / static_cast<double>(k + 1));
            }
            return {cp, Scalar(c.reference_enthalpy) + dt * dh,
                    Scalar(c.reference_entropy) + ideal_polynomial_entropy(a, T, c.reference_temperature)};
        } else if constexpr (std::is_same_v<C, IdealShomate>) {
            const Scalar t = T / Scalar(1000.0), t2 = t*t, t3 = t2*t, t4 = t3*t;
            return {a[0] + a[1]*t + a[2]*t2 + a[3]*t3 + a[4]/t2,
                c.formation_enthalpy + 1000.0*(a[0]*t + a[1]*t2/2.0
                    + a[2]*t3/3.0 + a[3]*t4/4.0 - a[4]/t + a[5] - a[7]),
                a[0]*scalar_log(t) + a[1]*t + a[2]*t2/2.0 + a[3]*t3/3.0
                    - a[4]/(2.0*t2) + a[6]};
        } else if constexpr (std::is_same_v<C, IdealNASA7>) {
            const Scalar t2 = T*T, t3 = t2*T, t4 = t3*T;
            return {R*(a[0]+a[1]*T+a[2]*t2+a[3]*t3+a[4]*t4),
                R*(a[0]*T+a[1]*t2/2.0+a[2]*t3/3.0+a[3]*t4/4.0+a[4]*t4*T/5.0+a[5]),
                R*(a[0]*scalar_log(T)+a[1]*T+a[2]*t2/2.0+a[3]*t3/3.0+a[4]*t4/4.0+a[6])};
        } else if constexpr (std::is_same_v<C, IdealNASA9>) {
            const Scalar t2 = T*T, t3 = t2*T, t4 = t3*T;
            return {R*(a[0]/t2+a[1]/T+a[2]+a[3]*T+a[4]*t2+a[5]*t3+a[6]*t4),
                R*(-a[0]/T+a[1]*scalar_log(T)+a[2]*T+a[3]*t2/2.0
                    +a[4]*t3/3.0+a[5]*t4/4.0+a[6]*t4*T/5.0+a[7]),
                R*(-a[0]/(2.0*t2)-a[1]/T+a[2]*scalar_log(T)+a[3]*T
                    +a[4]*t2/2.0+a[5]*t3/3.0+a[6]*t4/4.0+a[8])};
        } else {
            // Even arguments; remove constant primitive offsets before taking
            // differences. exp(-2u) never overflows at large positive u.
            auto terms = [&](const Scalar& t) {
                IdealThermalValue<Scalar> value{Scalar(a[0]), a[0]*t, a[0]*scalar_log(t)};
                for (int term = 0; term != 2; ++term) {
                    const double energy = std::abs(a[2 + 2*term]), weight = a[1 + 2*term];
                    if (energy == 0.0) {
                        if (term == 0) { value.cp += weight; value.h += weight*t; value.s += weight*scalar_log(t); }
                        continue;
                    }
                    const Scalar u = Scalar(energy)/t, q = scalar_exp(Scalar(-2.0)*u);
                    const Scalar denominator = term == 0 ? Scalar(-scalar_expm1(Scalar(-2.0)*u)) : Scalar(1.0)+q;
                    const Scalar ratio = Scalar(2.0)*u*scalar_exp(Scalar(-u))/denominator;
                    const Scalar tail = Scalar(2.0)*u*q/denominator;
                    value.cp += weight*ratio*ratio;
                    value.h += weight*t*tail;
                    value.s += weight*(term == 0 ? Scalar(tail-scalar_log(denominator))
                                                : Scalar(tail+scalar_log1p(q)));
                }
                return value;
            };
            const auto value = terms(T), anchor = terms(Scalar(c.reference_temperature));
            return {value.cp, c.reference_enthalpy + value.h-anchor.h,
                    c.reference_entropy + value.s-anchor.s};
        }
    }, correlation.form);
    result.s -= R * (scalar_log(Scalar(pressure))
                     - scalar_log(Scalar(correlation.source_pressure)));
    if constexpr (std::is_same_v<Scalar, double>) {
        if (!std::isfinite(result.cp) || !std::isfinite(result.h) || !std::isfinite(result.s))
            throw std::domain_error("nonfinite ideal thermal result");
    }
    return result;
}

// EqExclude: source-declared interval selection over intervals admitted by validate_ideal_intervals
inline const IdealInterval* ideal_interval(const std::vector<IdealInterval>& intervals, double T) {
    const IdealInterval* selected = nullptr;
    for (const auto& interval : intervals)
        if (ideal_interval_contains(interval, T)) selected = &interval;
    return selected;
}

} // namespace epcsaft
