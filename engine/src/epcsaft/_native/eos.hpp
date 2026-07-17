#pragma once

#include "model.hpp"

#include <vector>

namespace epcsaft {

struct NativeResult final {
    double hard_chain;
    double dispersion;
    double association;
    double debye_huckel;
    double born_ssm_ds;
    double residual_helmholtz;
    double compressibility_factor;
    double pressure;
};

class NativeModel final {
public:
    explicit NativeModel(NativeParameters parameters);
    NativeResult evaluate(
        double temperature,
        double molar_density,
        const std::vector<double>& mole_fractions
    ) const;

private:
    const NativeParameters parameters_;
};

}  // namespace epcsaft
