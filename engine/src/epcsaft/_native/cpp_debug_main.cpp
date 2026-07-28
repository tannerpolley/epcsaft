#include <epcsaft/native_model_v1.h>
#include <epcsaft/native_sdk_v1.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string option(int argc, char** argv, const std::string& name) {
    const std::string prefix = "--" + name + "=";
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument.rfind(prefix, 0) == 0) {
            return argument.substr(prefix.size());
        }
    }
    throw std::invalid_argument("missing required option --" + name);
}

double number(const std::string& text, const std::string& name) {
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument("invalid value for --" + name);
    }
    return value;
}

std::vector<double> numbers(const std::string& text) {
    std::vector<double> values;
    std::size_t begin = 0;
    while (begin < text.size()) {
        const std::size_t end = text.find(',', begin);
        values.push_back(number(
            text.substr(begin, end == std::string::npos ? end : end - begin),
            "mole-fractions"
        ));
        if (end == std::string::npos) {
            break;
        }
        begin = end + 1;
    }
    return values;
}

int run(int argc, char** argv) {
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE] = {};
    using Model = std::unique_ptr<
        epcsaft_native_model_handle_v1,
        decltype(&epcsaft_native_model_destroy_v1)
    >;
    Model model(
        epcsaft_native_model_load_v1(
            option(argc, argv, "model-config").c_str(), error, sizeof(error)
        ),
        &epcsaft_native_model_destroy_v1
    );
    if (!model) {
        throw std::runtime_error(error);
    }
    const auto* sdk = epcsaft_native_model_sdk_v1(model.get());
    const double temperature = number(
        option(argc, argv, "temperature-k"), "temperature-k"
    );
    const double density = number(
        option(argc, argv, "molar-density-mol-per-m3"),
        "molar-density-mol-per-m3"
    );
    const std::vector<double> amounts = numbers(
        option(argc, argv, "mole-fractions")
    );
    const std::size_t coordinates = amounts.size() + 1;
    std::vector<double> gradient(coordinates);
    std::vector<double> hessian(coordinates * coordinates);
    epcsaft_mixture_phase_block_result_v1 result{
        sizeof(result),
        EPCSAFT_NATIVE_STATUS_OK_V1,
        coordinates,
        gradient.size(),
        hessian.size(),
        gradient.data(),
        hessian.data(),
    };
    if (sdk == nullptr || sdk->evaluate_mixture_phase == nullptr
        || amounts.size() != sdk->component_count
        || sdk->evaluate_mixture_phase(
            sdk->model_context,
            temperature,
            amounts.data(),
            amounts.size(),
            1.0 / density,
            &result
        ) != EPCSAFT_NATIVE_STATUS_OK_V1) {
        throw std::runtime_error(result.error);
    }
    std::cout << std::setprecision(17)
              << "pressure_pa=" << result.pressure_pa << '\n'
              << "helmholtz_over_rt_reference_amount="
              << result.helmholtz_over_rt_reference_amount << '\n'
              << "parameter_fingerprint=" << result.parameter_fingerprint << '\n'
              << "capability_count=" << sdk->capability_count << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
