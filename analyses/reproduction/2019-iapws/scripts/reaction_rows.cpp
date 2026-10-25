#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
// Contract native Chemical rows; no activity quotient or equilibrium calculation.
int main() {
    std::array<double,3> chemical{}, x{}, n{};
    std::array<bool,3> seen{};
    double total, raw, retained;
    int count, kind, phase, index;
    std::cin >> total;
    for (int i=0;i<3;++i) { std::cin >> x[i]; n[i]=total*x[i]; }
    std::cin >> count;
    for (int i=0;i<count;++i) {
        std::cin >> kind >> phase >> index >> raw >> retained;
        if (!std::cin || !std::isfinite(raw) || raw!=retained) throw std::runtime_error("raw evaluation differs from retained residual");
        if (kind && phase==0) {
            if (index<0 || index>2 || seen[index]) throw std::runtime_error("invalid Chemical row identity");
            chemical[index]=raw; seen[index]=true;
        }
    }
    for (bool found:seen) if (!found) throw std::runtime_error("missing Chemical row");
    std::cout << std::setprecision(17) << "{\"species_amounts_mol\":[" << n[0] << ',' << n[1] << ',' << n[2]
              << "],\"reaction_raw\":" << -2*chemical[0]+chemical[1]+chemical[2]
              << ",\"oxygen_balance_mol\":" << n[0]+n[1]+n[2]-1
              << ",\"hydrogen_balance_mol\":" << 2*n[0]+3*n[1]+n[2]-2
              << ",\"charge_mol\":" << n[1]-n[2]
              << ",\"normalization\":" << x[0]+x[1]+x[2]-1 << '}';
}
