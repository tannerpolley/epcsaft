#pragma once

#include "model.hpp"

#include <optional>
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
    std::vector<double> residual_chemical_potential_over_rt;
    std::optional<std::vector<double>> log_fugacity_coefficient;
    double residual_entropy;
    double residual_internal_energy;
    double residual_enthalpy;
    double residual_isochoric_heat_capacity;
    std::optional<double> residual_isobaric_heat_capacity;
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
