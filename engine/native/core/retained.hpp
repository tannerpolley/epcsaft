#pragma once

#include <atomic>
#include <memory>
#include <utility>

namespace epcsaft::core {

template<class Key, class Value>
class Retained final {
    friend Key;

    struct Publication final {
        const Key key;
        const Value value;
        mutable std::atomic_bool live{true};

        Publication(Key k, Value v) : key(std::move(k)), value(std::move(v)) {}
    };
    std::shared_ptr<const Publication> publication_;

    static Retained publish(Key key, Value value) {
        Retained entry;
        entry.publication_ = std::make_shared<const Publication>(std::move(key), std::move(value));
        return entry;
    }

    template<class Request>
    std::shared_ptr<const Value> resolve(const Request& request) const {
        const auto& entry = publication_;
        if (!entry || !entry->live.load(std::memory_order_acquire) || !entry->key.matches(request)) return {};
        if (!Key::admits(request, entry->value)) {
            entry->live.store(false, std::memory_order_release);
            return {};
        }
        return {entry, &entry->value};
    }
};

}  // namespace epcsaft::core
