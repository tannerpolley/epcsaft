#include "eos.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Arguments final {
    bool help = false;
    std::vector<std::pair<std::string, std::string>> values;
};

const std::string& required_value(
    const Arguments& arguments,
    const std::string& name
) {
    for (const auto& [key, value] : arguments.values) {
        if (key == name) {
            return value;
        }
    }
    throw std::invalid_argument("missing required option --" + name);
}

double parse_double(const std::string& text, const std::string& name) {
    std::size_t consumed = 0;
    double value = 0.0;
    try {
        value = std::stod(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument("invalid value for --" + name);
    }
    if (consumed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument("invalid value for --" + name);
    }
    return value;
}

int parse_int(const std::string& text, const std::string& name) {
    std::size_t consumed = 0;
    long long value = 0;
    try {
        value = std::stoll(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument("invalid value for --" + name);
    }
    if (consumed != text.size()
        || value < std::numeric_limits<int>::min()
        || value > std::numeric_limits<int>::max()) {
        throw std::invalid_argument("invalid value for --" + name);
    }
    return static_cast<int>(value);
}

template <typename Value, typename Parser>
std::vector<Value> parse_list(
    const std::string& text,
    const std::string& name,
    Parser parser
) {
    if (text.empty()) {
        return {};
    }
    std::vector<Value> values;
    std::size_t begin = 0;
    while (begin <= text.size()) {
        const std::size_t end = text.find(',', begin);
        const std::string token = text.substr(
            begin,
            end == std::string::npos ? std::string::npos : end - begin
        );
        if (token.empty()) {
            throw std::invalid_argument("invalid comma-separated value for --" + name);
        }
        values.push_back(parser(token, name));
        if (end == std::string::npos) {
            break;
        }
        begin = end + 1;
    }
    return values;
}

std::vector<double> doubles(const Arguments& arguments, const std::string& name) {
    return parse_list<double>(
        required_value(arguments, name),
        name,
        [](const std::string& text, const std::string& option) {
            return parse_double(text, option);
        }
    );
}

std::vector<int> integers(const Arguments& arguments, const std::string& name) {
    return parse_list<int>(
        required_value(arguments, name),
        name,
        [](const std::string& text, const std::string& option) {
            return parse_int(text, option);
        }
    );
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments arguments;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            arguments.help = true;
            continue;
        }
        if (argument.rfind("--", 0) != 0) {
            throw std::invalid_argument("options must start with --");
        }
        const std::size_t separator = argument.find('=');
        const std::string name = argument.substr(
            2,
            separator == std::string::npos ? std::string::npos : separator - 2
        );
        if (name.empty()) {
            throw std::invalid_argument("empty option name");
        }
        constexpr std::array known_options = {
            "temperature-k",
            "molar-density-mol-per-m3",
            "mole-fractions",
            "segment-count",
            "sigma-constant-angstrom",
            "sigma-term-component",
            "sigma-term-amplitude",
            "sigma-term-exponent",
            "dispersion-energy-over-k-k",
            "k-ij",
            "l-ij",
            "site-component",
            "site-multiplicity",
            "association-active",
            "association-energy-over-k-k",
            "association-volume",
            "molar-mass-kg-per-mol",
            "charge-number",
            "relative-permittivity",
            "dielectric-ion-suppression-coefficient",
            "born-diameter-angstrom",
            "solvation-factor",
            "ionic-region-relative-permittivity",
            "electrolyte-model",
        };
        bool known = false;
        for (const char* known_option : known_options) {
            if (name == known_option) {
                known = true;
                break;
            }
        }
        if (!known) {
            throw std::invalid_argument("unknown option --" + name);
        }
        std::string value;
        if (separator == std::string::npos) {
            if (index + 1 >= argc) {
                throw std::invalid_argument("missing value for --" + name);
            }
            value = argv[++index];
        } else {
            value = argument.substr(separator + 1);
        }
        for (const auto& [existing_name, unused] : arguments.values) {
            if (existing_name == name) {
                throw std::invalid_argument("duplicate option --" + name);
            }
        }
        arguments.values.emplace_back(name, std::move(value));
    }
    return arguments;
}

void print_help() {
    std::cout << R"(Usage: epcsaft_cpp_debug [options]

Evaluate one explicit-density state with a caller-supplied NativeParameters.
Scalar values use the units in their option names. Comma-separated lists are
row-major for pair and site matrices; an empty value is an explicit empty list.
All options below are required except --help.

State:
  --temperature-k=K
  --molar-density-mol-per-m3=mol/m^3
  --mole-fractions=x1,x2,...

NativeParameters:
  --segment-count=m1,m2,...
  --sigma-constant-angstrom=angstroms
  --sigma-term-component=i1,i2,...
  --sigma-term-amplitude=angstroms
  --sigma-term-exponent=1/K
  --dispersion-energy-over-k-k=K
  --k-ij=values
  --l-ij=values
  --site-component=i1,i2,...
  --site-multiplicity=n1,n2,...
  --association-active=0,1,...
  --association-energy-over-k-k=K
  --association-volume=dimensionless
  --molar-mass-kg-per-mol=kg/mol
  --charge-number=z1,z2,...
  --relative-permittivity=values
  --dielectric-ion-suppression-coefficient=value
  --born-diameter-angstrom=angstroms
  --solvation-factor=values
  --ionic-region-relative-permittivity=value
  --electrolyte-model=0, 1 (Figiel SSM+DS), or 2 (Ascani original Born)
)";
}

void print_vector(const std::vector<double>& values) {
    std::cout << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << ']';
}

void validate_state(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    std::size_t component_count
) {
    if (!(temperature > 0.0) || !(molar_density > 0.0)) {
        throw std::invalid_argument(
            "temperature and molar density must be positive"
        );
    }
    if (mole_fractions.size() != component_count) {
        throw std::invalid_argument(
            "mole-fraction count must equal the segment-count count"
        );
    }
    double total = 0.0;
    for (double fraction : mole_fractions) {
        if (fraction < 0.0) {
            throw std::invalid_argument("mole fractions must be nonnegative");
        }
        total += fraction;
    }
    if (std::abs(total - 1.0) > 1.0e-12) {
        throw std::invalid_argument("mole fractions must sum to one");
    }
}

epcsaft::NativeParameters parse_parameters(const Arguments& arguments) {
    return epcsaft::NativeParameters{
        doubles(arguments, "segment-count"),
        doubles(arguments, "sigma-constant-angstrom"),
        integers(arguments, "sigma-term-component"),
        doubles(arguments, "sigma-term-amplitude"),
        doubles(arguments, "sigma-term-exponent"),
        doubles(arguments, "dispersion-energy-over-k-k"),
        doubles(arguments, "k-ij"),
        doubles(arguments, "l-ij"),
        integers(arguments, "site-component"),
        integers(arguments, "site-multiplicity"),
        integers(arguments, "association-active"),
        doubles(arguments, "association-energy-over-k-k"),
        doubles(arguments, "association-volume"),
        doubles(arguments, "molar-mass-kg-per-mol"),
        integers(arguments, "charge-number"),
        doubles(arguments, "relative-permittivity"),
        parse_double(
            required_value(arguments, "dielectric-ion-suppression-coefficient"),
            "dielectric-ion-suppression-coefficient"
        ),
        doubles(arguments, "born-diameter-angstrom"),
        doubles(arguments, "solvation-factor"),
        parse_double(
            required_value(arguments, "ionic-region-relative-permittivity"),
            "ionic-region-relative-permittivity"
        ),
        parse_int(
            required_value(arguments, "electrolyte-model"),
            "electrolyte-model"
        ),
    };
}

int run(const Arguments& arguments) {
    const epcsaft::NativeParameters parameters = parse_parameters(arguments);
    const double temperature = parse_double(
        required_value(arguments, "temperature-k"),
        "temperature-k"
    );
    const double molar_density = parse_double(
        required_value(arguments, "molar-density-mol-per-m3"),
        "molar-density-mol-per-m3"
    );
    const std::vector<double> mole_fractions = doubles(arguments, "mole-fractions");
    validate_state(
        temperature,
        molar_density,
        mole_fractions,
        parameters.segment_count.size()
    );

    const epcsaft::NativeModel model(parameters);
    const epcsaft::NativeResult result = model.evaluate(
        temperature,
        molar_density,
        mole_fractions
    );

    std::cout << std::setprecision(17);
    std::cout << "temperature_k=" << temperature << '\n';
    std::cout << "molar_density_mol_per_m3=" << result.molar_density << '\n';
    std::cout << "mole_fractions=";
    print_vector(mole_fractions);
    std::cout << '\n';
    std::cout << "hard_chain=" << result.hard_chain << '\n';
    std::cout << "dispersion=" << result.dispersion << '\n';
    std::cout << "association=" << result.association << '\n';
    std::cout << "debye_huckel=" << result.debye_huckel << '\n';
    std::cout << "born_ssm_ds=" << result.born_ssm_ds << '\n';
    std::cout << "residual_helmholtz=" << result.residual_helmholtz << '\n';
    std::cout << "pressure_pa=" << result.pressure << '\n';
    std::cout << "compressibility_factor=" << result.compressibility_factor << '\n';
    std::cout << "residual_chemical_potential_over_rt=";
    print_vector(result.residual_chemical_potential_over_rt);
    std::cout << '\n';
    std::cout << "log_fugacity_coefficient=";
    if (result.log_fugacity_coefficient.has_value()) {
        print_vector(*result.log_fugacity_coefficient);
    } else {
        std::cout << "unavailable";
    }
    std::cout << '\n';
    std::cout << "density_diagnostics=unavailable (explicit-density state)\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Arguments arguments = parse_arguments(argc, argv);
        if (arguments.help) {
            print_help();
            return 0;
        }
        return run(arguments);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
