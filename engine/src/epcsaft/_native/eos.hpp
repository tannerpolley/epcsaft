#pragma once

#include <vector>

namespace epcsaft {

class NativeModel final {
public:
    NativeModel(
        std::vector<double> m,
        std::vector<double> sigma,
        std::vector<double> epsilon_over_k,
        std::vector<double> kij
    );

private:
    const std::vector<double> m_;
    const std::vector<double> sigma_;
    const std::vector<double> epsilon_over_k_;
    const std::vector<double> kij_;
};

}  // namespace epcsaft
