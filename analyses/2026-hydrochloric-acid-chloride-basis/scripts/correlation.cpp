#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    std::ifstream input(std::string(argv[1]) + "/fit-values.csv");
    std::ofstream output(std::string(argv[1]) + "/correlation.csv");
    output << std::setprecision(17) << "stage,start,row,column,correlation\n";
    std::map<std::pair<int, int>, std::vector<double>> matrices;
    std::string line; std::getline(input, line);
    while (std::getline(input, line)) {
        std::stringstream stream(line); std::string stage, start, field, index, value;
        std::getline(stream, stage, ','); std::getline(stream, start, ','); std::getline(stream, field, ',');
        std::getline(stream, index, ','); std::getline(stream, value);
        if (field == "C") matrices[{std::stoi(stage), std::stoi(start)}].push_back(std::stod(value));
    }
    for (const auto& [key, C] : matrices) {
        const auto [q, start] = key;
        if (C.size() != std::size_t(q*q)) return 1;
        for (int a = 0; a < q; ++a) for (int b = 0; b < q; ++b) {
            if (!(C[a*q+a] > 0 && C[b*q+b] > 0)) return 1;
            const double rho = C[a*q+b]/std::sqrt(C[a*q+a]*C[b*q+b]);
            if (!std::isfinite(rho) || (a == b && std::abs(rho-1) > 1e-12)) return 1;
            output << q << ',' << start << ',' << a << ',' << b << ',' << rho << '\n';
        }
    }
    return input.eof() && output.good() ? 0 : 1;
}
