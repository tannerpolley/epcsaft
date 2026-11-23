#pragma once

#include <array>
#include <cstddef>

namespace epcsaft {

struct AssociationPlan {
    enum Strategy { General, ReciprocalPair, BalancedBinary };
    Strategy strategy = General;
    bool has_edges = false, one_partner_per_site = true;
    std::array<std::size_t, 4> sites{};
};

struct AssociationAction;

}  // namespace epcsaft
