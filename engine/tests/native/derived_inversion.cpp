#include "core/numerical_policy.hpp"

inline constexpr epcsaft::policy::Tolerance INVERTED = epcsaft::policy::derived(
    epcsaft::policy::ASSOCIATION_RESIDUAL,
    0.5,
    "INVERTED",
    "must not compile"
);
