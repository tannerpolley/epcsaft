#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

namespace epcsaft::equilibrium::stability_grid_detail {
using Integer = std::uint64_t;
using Node = std::vector<Integer>;
constexpr auto maximum = std::numeric_limits<Integer>::max();

// C(m+k-1,k-1), cancelling the divisor before each checked product.
inline Integer count(Integer m, Integer k) {
    if (!k) return 0;
    Integer total = 1;
    for (Integer j = 1; j < k; ++j) {
        if (m > maximum - j) throw std::overflow_error("raw grid count overflow");
        const auto divisor = std::gcd(total, j);
        const auto factor = (m + j) / (j / divisor);
        total /= divisor;
        if (total > maximum / factor) throw std::overflow_error("raw grid count overflow");
        total *= factor;
    }
    return total;
}

struct Vertex {
    std::size_t i, j;
    Integer left, denominator;
};

template<class Emit>
void enumerate(const std::vector<Vertex>& vertices, Integer m, Integer scale, std::size_t n, Emit emit) {
    Node weight(vertices.size());
    weight.front() = m;
    for (;;) {
        Node point(n);
        for (std::size_t v = 0; v < vertices.size(); ++v) {
            const auto [i, j, left, denominator] = vertices[v];
            const auto units = scale / denominator * weight[v];
            point[i] += units * left;
            point[j] += units * (denominator - left);
        }
        emit(std::move(point));
        if (weight.back() == m) return;
        auto cursor = weight.end() - 2;
        while (!*cursor) --cursor;
        const auto tail = weight.back();
        weight.back() = 0;
        --*cursor;
        cursor[1] = tail + 1;
    }
}

using Order = std::vector<std::size_t>::iterator;
inline void partition(const std::vector<Node>& nodes, Order first, Order last, std::size_t axis) {
    if (last - first < 2) return;
    const auto mid = first + (last - first) / 2;
    std::nth_element(first, mid, last, [&](auto a, auto b) { return nodes[a][axis] < nodes[b][axis]; });
    const auto next = (axis + 1) % nodes.front().size();
    partition(nodes, first, mid, next);
    partition(nodes, mid, last, next);
}

inline Integer difference(Integer a, Integer b) { return std::max(a, b) - std::min(a, b); }

// The median splits the current bounding box; child boxes inherit all other bounds.
inline void nearest(const std::vector<Node>& nodes, const Node& gap, Order first, Order last, std::size_t axis,
    Node& lower, Node& upper, Integer bound, Integer& best, std::size_t& index) {
    if (bound >= best) return;
    if (last - first == 1) {
        const auto distance = std::inner_product(gap.begin(), gap.end(), nodes[*first].begin(), Integer{0},
            std::plus<>{}, difference);
        if (distance < best) { best = distance; index = *first; }
        return;
    }
    const auto mid = first + (last - first) / 2;
    const auto cut = nodes[*mid][axis], x = gap[axis];
    const auto next = (axis + 1) % gap.size();
    const auto visit = [&](bool left) {
        auto& edge = left ? upper[axis] : lower[axis];
        const auto saved = edge;
        const auto prior = difference(x, std::clamp(x, lower[axis], upper[axis]));
        edge = cut;
        const auto distance = difference(x, std::clamp(x, lower[axis], upper[axis]));
        nearest(nodes, gap, left ? first : mid, left ? mid : last, next, lower, upper,
            bound - prior + distance, best, index);
        edge = saved;
    };
    visit(x <= cut);
    visit(x > cut);
}

struct Coverage { Integer distance = 0; std::size_t gap, admitted; };
inline std::optional<Coverage> coverage(const std::vector<Node>& nodes, const std::vector<bool>& admitted, Integer scale) {
    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < nodes.size(); ++i) if (admitted[i]) order.push_back(i);
    if (order.empty()) return {};
    partition(nodes, order.begin(), order.end(), 0);
    Node lower(nodes.front().size()), upper(lower.size(), scale);
    Coverage worst{0, order.front(), order.front()};
    for (std::size_t i = 0; i < nodes.size(); ++i) if (!admitted[i]) {
        Integer best = maximum;
        std::size_t index = 0;
        nearest(nodes, nodes[i], order.begin(), order.end(), 0, lower, upper, 0, best, index);
        if (best > worst.distance) worst = {best, i, index};
    }
    return worst;
}
}  // namespace epcsaft::equilibrium::stability_grid_detail
